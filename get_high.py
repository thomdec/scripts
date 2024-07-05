import allel
import numpy as np
import pandas as pd
import gzip
import argparse
import sys
import subprocess
from uuid import uuid4
from warnings import filterwarnings
from os import path

#TO DO
#   - Read the VCF chunk by chunk, rather than reading in the entire genotype matrix at once

def main():

    filterwarnings('ignore') 

    help ="\nGET HIGH: Extract high impact alleles frequency and effect from a SNPEff annotated VCF"

    parser = argparse.ArgumentParser(description = get_high_image + help, prog='tool', formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog,max_help_position=40))

    parser._action_groups.pop()
    required = parser.add_argument_group('Required parameters')
    optional = parser.add_argument_group('Optional parameters')

    required.add_argument('-v', '--vcf', type=str, help='VCF file', required=True)
    required.add_argument('-o', '--output', type=str, help='Output name', required=True)
    optional.add_argument('-r', '--region', type=str, help='Region string as per bcftools: chr|chr:pos|chr:beg-end', required=False)
    optional.add_argument('-p', '--pop', type=str, help='Populations file.\nCsv file with id and population in each row: [id,pop]', required=False)

    args = parser.parse_args()

    # checks
    if path.exists(args.vcf) is not True:
        raise Exception('[YAWP] ERROR: The specified VCF ' + str(args.vcf) + ' does not exist') 

    # remove invariant sites
    print("\n[GET HIGH]: Removing invariant sites")
    tmp_vcf = "tmp.variant." + uuid4().hex + ".gz"

    if args.region is not None:
        run_command(f"""bcftools filter -r {args.region} -e 'N_ALT == 0' {args.vcf} -Oz -o {tmp_vcf}""")
    else:
        run_command(f"""bcftools filter -e 'N_ALT == 0' {args.vcf} -Oz -o {tmp_vcf}""")

    # create a temporary query of the VCF
    tmp_annotation_file = "tmp.annotation.file." + uuid4().hex + ".gz"
    run_command(f"""bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ANN\n" {tmp_vcf} | gzip > {tmp_annotation_file}""")

    # import genotype matrix
    callset = allel.read_vcf(tmp_vcf, fields=['calldata/GT'])
    g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))

    annotation_file = gzip.open(tmp_annotation_file, "rt")

    if args.output: output_file = gzip.open(args.output, "wt") if args.output.endswith(".gz") else open(args.output, "wt")

    if args.pop is not None:
        if path.exists(args.pop) is not True:
            raise Exception('[YAWP] ERROR: The specified populations file ' + str(args.pop) + ' does not exist')

        pop_df = pd.read_csv(args.pop, names=["id", "pop"])
        vcf_header = allel.read_vcf_headers(args.vcf)

        if not np.all(np.isin(vcf_header.samples, pop_df.id)):
            raise Exception('[YAWP] ERROR: Samples in the population file could not be found in the VCF')

        pops = {}
        unique_pops = np.unique(pop_df["pop"]).astype(str)
        for pop in unique_pops:
            ids = pop_df[pop_df['pop'] == pop]['id']
            pops[pop] = np.where(np.isin(vcf_header.samples, ids))[0]

        output_file.write('chr' + '\t' + 'pos' + '\t' + 'pop' + '\t' +  'n_alleles' + '\t' +  'n_high_alleles' + '\t' +  'effect' + '\n')

    else:
        output_file.write('chr' + '\t' + 'pos' + '\t' +  'n_alleles' + '\t' +  'n_high_alleles' + '\t' +  'effect' + '\n')

    print("[GET HIGH]: Extracting frequency and effect of high impact variants")

    for i, line in enumerate(annotation_file):
        chr, pos, ref, alt, annotations = line.split("\t")

        # get the alleles
        alleles = np.array(alt.split(","))
        alleles = np.insert(alleles, 0, str(ref))

        # get the effect impacts (MODIFIER, LOW, MODERATE, HIGH) and effects (stop_gained, frameshift etc)
        effects = [[] for y in range(len(alleles))]
        effect_impacts = np.zeros(len(alleles), dtype=int)

        for annotation_string in annotations.split(","):
            ann = annotation_string.split("|")

            if ann[2] == "HIGH":
                effect_impacts[int(np.where(alleles == ann[0])[0])] += 1
                effects[int(np.where(alleles == ann[0])[0])].append(ann[1])

        # go further if "HIGH" impact variants are observed
        if any(effect_impacts > 0):

            if args.pop is not None:

                # get the high impact effects
                high_effects = list(filter(None, effects))
                if len(high_effects) == 1: high_effects = ",".join(high_effects[0])

                for pop in pops:
                    # get the alleles count
                    ac = g[i:i+1][:, pops[pop]].count_alleles(max_allele = len(alleles) - 1)[0]

                    # get the number of "HIGH impact variants"
                    n_high = ",".join(map(str, ac[effect_impacts > 0]))

                    output_file.write(chr + '\t' + pos + '\t' + pop + '\t' +  str(np.sum(ac)) + '\t' +  str(n_high) + '\t' +  str(high_effects) + '\n')

            else:
                # get the alleles count
                ac = g[i:i+1].count_alleles(max_allele = len(alleles) - 1)[0]

                # get the number of "HIGH impact variants"
                n_high = ",".join(map(str, ac[effect_impacts > 0]))

                # get the high impact effects
                high_effects = list(filter(None, effects))
                if len(high_effects) == 1: high_effects = ",".join(high_effects[0])

                output_file.write(chr + '\t' + pos + '\t' +  str(np.sum(ac)) + '\t' +  str(n_high) + '\t' +  str(high_effects) + '\n')

    annotation_file.close()
    output_file.close()

    run_command(f"rm {tmp_annotation_file} && rm {tmp_vcf}")

    print("[GET HIGH]: Done. Got high")

def run_command(args):
    try:
        call = subprocess.run(args, text=True, capture_output=True, shell=True)
        call.check_returncode()
    except subprocess.CalledProcessError as cpe:
        if call.stdout:
            sys.exit('[X] The command following command failed with error code %r:\n[X] => %s\n[X] (STDOUT): %r\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stdout, cpe.stderr))
        sys.exit('[X] The command following command failed with error code %r:\n[X] => %s\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stderr.rstrip("\n")))

get_high_image = \
    " ░▒▓██████▓▒░░▒▓████████▓▒░▒▓████████▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░ \n" \
    "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░         ░▒▓█▓▒░          ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ \n" \
    "░▒▓█▓▒░      ░▒▓█▓▒░         ░▒▓█▓▒░          ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░ \n" \
    "░▒▓█▓▒▒▓███▓▒░▒▓██████▓▒░    ░▒▓█▓▒░          ░▒▓████████▓▒░▒▓█▓▒░▒▓█▓▒▒▓███▓▒░▒▓████████▓▒░ \n" \
    "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░         ░▒▓█▓▒░          ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ \n" \
    "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░         ░▒▓█▓▒░          ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ \n" \
    " ░▒▓██████▓▒░░▒▓████████▓▒░  ░▒▓█▓▒░          ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░ \n" \

if __name__ == "__main__":

   main()