import numpy as np
from cyvcf2 import VCF

from scipy.special import comb
from itertools import combinations

import os
import sys
import argparse
import subprocess

## TO DO
#   - Add a multiprocessing 
#   - Add region parameters! + only read in the samples supplied in the population file
#   - Add more windows options: mask, number of sites and info about the number of missing genotypes (% of genotypes covered?)
#   - Make more verbose: ADD TIME INFO

def main():
    
    sys.tracebacklimit = 0

    help = "\nYAYAWP: computes common popgen stats"

    parser = argparse.ArgumentParser(description=help, prog='tool', formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog,max_help_position=40))

    parser._action_groups.pop()
    required = parser.add_argument_group('Required parameters')
    optional = parser.add_argument_group('Optional parameters')

    #required
    required.add_argument('-v', '--vcf', type=str, help='VCF file', required=True)
    required.add_argument('-p', '--populations', type=str, help='Populations file.\nCsv file with id and population in each row: [id,pop]', required=True)
    required.add_argument('-o', '--output', type=str, help='Output prefix', required=True)

    # optional
    optional.add_argument('-w', '--window_size', type=int, help='Window size. If not supplied, will return global estimates', required=False)
    #optional.add_argument('-m', '--mask', type=str, help='Mask (e.g 4D sites). Must be in the BED format', required=False) #should I rename mask to filter?

    args = parser.parse_args()

    # check parameters
    print("[YAYAWP] Checking command line arguments")

    if os.path.exists(args.vcf) is not True:
        raise Exception(f"[YAYAWP] ERROR: The specified VCF {args.vcf} does not exist") 

    if not os.path.exists(args.vcf + ".csi") | os.path.exists(args.vcf + ".tbi"): 
        raise Exception('[YAYAWP] ERROR: The vcf is not indexed. The vcf can be indexed with "bcftools index" or "tabix -p vcf"') 

    if os.path.exists(args.populations) is not True:
        raise Exception(f"[YAYAWP] ERROR: The specified populations file {args.populations} does not exist")

    if os.access(os.path.dirname("./" + args.output), os.W_OK) is not True:
        raise Exception(f"[YAYAWP] ERROR: Cannot write output to {args.output}")

    pop_array = np.genfromtxt(args.populations, delimiter=',', dtype=str)

    _stdout = subprocess.check_output(f"bcftools query -l {args.vcf}", shell=True, text=True)
    samples = np.array([sample_id for sample_id in _stdout.split("\n") if sample_id], dtype=str)

    if not np.all(np.isin(pop_array[:, 0], samples)):
        raise Exception(f"[YAYAWP] ERROR: The following sample in the population file could not be found in the VCF:\n{', '.join(pop_array[:, 0][~np.isin(pop_array[:, 0], samples)])}")

    pops = {}
    for pop in np.unique(pop_array[:, 1]):
        ids = pop_array[pop_array[:, 1] == pop][:, 0]
        pops[pop] = np.where(np.isin(samples, ids))[0]


    if args.window_size is not None:

        print("[YAYAWP] YAYAWPing in windows")

        data = subprocess.check_output("bcftools index -s " + args.vcf, shell=True, text=True)
        chromosomes_info = [x.split('\t') for x in data[:-1].split('\n')]
        windows = list()
        for a in chromosomes_info:
            windows.extend(create_windows(args.window_size, a[1], a[0]))

        pi_file = open(f"{args.output}.pi.tsv", 'w')
        dxy_file = open(f'{args.output}.dxy.tsv', 'w')

        pi_file.write('chromosome' + '\t' + 'start' + '\t' +  'end' + '\t' +  'pop' + '\t' +  'pi' + '\t' +  'n_diff' + '\t' +  'n_comp' + '\n')
        dxy_file.write('chromosome' + '\t' + 'start' + '\t' +  'end' + '\t' +  'pop1' + '\t' +  'pop2' + '\t' +  'dxy' + '\t' + 'n_diff' + '\t' +  'n_comp' + '\n')

        for window in windows:

            pi_array, dxy_array = calc_pi_dxy(args.vcf, pops, window)

            for i, pop in enumerate(pops):
                pi_file.write(str(window[0]) + '\t' +  str(window[1]) + '\t' +  str(window[2]) + '\t' + pop + '\t' +  str(pi_array[i, 0] / pi_array[i, 1]) + '\t' +  str(pi_array[i, 0]) + '\t' +  str(pi_array[i, 1]) + '\n')

            for i, pop in enumerate(combinations(pops, 2)):
                dxy_file.write(str(window[0]) + '\t' +  str(window[1]) + '\t' +  str(window[2]) + '\t' + pop[0] + '\t' +  pop[1] + '\t' +  str(dxy_array[i, 0] / dxy_array[i, 1]) + '\t' +  str(dxy_array[i, 0]) + '\t' +  str(dxy_array[i, 1]) + '\n')

        pi_file.close()
        dxy_file.close()

    else:

        print("[YAYAWP] YAYAWPing globally")

        pi_file = open(f"{args.output}.pi.tsv", 'w')
        dxy_file = open(f'{args.output}.dxy.tsv', 'w')

        pi_file.write('pop' + '\t' +  'pi' + '\t' +  'n_diff' + '\t' +  'n_comp' + '\n')
        dxy_file.write('pop1' + '\t' +  'pop2' + '\t' +  'dxy' + '\t' + 'n_diff' + '\t' +  'n_comp' + '\n')

        pi_array, dxy_array = calc_pi_dxy(args.vcf, pops, None)

        for i, pop in enumerate(pops):
            pi_file.write(pop + '\t' +  str(pi_array[i, 0] / pi_array[i, 1]) + '\t' +  str(pi_array[i, 0]) + '\t' +  str(pi_array[i, 1]) + '\n')

        for i, pop in enumerate(combinations(pops, 2)):
            dxy_file.write(pop[0] + '\t' +  pop[1] + '\t' +  str(dxy_array[i, 0] / dxy_array[i, 1]) + '\t' +  str(dxy_array[i, 0]) + '\t' +  str(dxy_array[i, 1]) + '\n')

    print("[YAYAWP] Successfully ran YAYAWP")

def create_windows(window_size, chr_len, chr):
    """Create windows"""

    window_starts = np.array([*range(0, int(chr_len), window_size)])
    window_ends = np.array([*range(0 + window_size, int(chr_len) + window_size, window_size)])
    window_ends[-1] = int(chr_len)
    windows = list(zip([chr for i in range(len(window_starts))], window_starts, window_ends))
    return windows

def diff_comp_within(variant, pop):
    """Computes the number of nucleotide comparisons and differences for a single site within a population"""
    
    alleles_count = np.array([sum([variant.genotypes[ind][:2].count(allele) for ind in pop]) for allele in np.arange(len(variant.ALT) + 1)])
    n_allele = np.sum(alleles_count)
    n_comp = comb(n_allele, 2)
    n_same = np.sum(alleles_count * (alleles_count - 1) / 2)
    n_diff =  n_comp - n_same
    return n_diff, n_comp


def diff_comp_between(variant, pop1, pop2):
    """Computes the number of nucleotide comparisons and differences for a single site between two population"""
    
    alleles_count_pop1 = np.array([sum([variant.genotypes[ind][:2].count(allele) for ind in pop1]) for allele in np.arange(len(variant.ALT) + 1)])
    alleles_count_pop2 = np.array([sum([variant.genotypes[ind][:2].count(allele) for ind in pop2]) for allele in np.arange(len(variant.ALT) + 1)])
    n_allele_pop1 = np.sum(alleles_count_pop1)
    n_allele_pop2 = np.sum(alleles_count_pop2)
    n_comp = n_allele_pop1 * n_allele_pop2
    n_same = np.sum(alleles_count_pop1 * alleles_count_pop2)
    n_diff = n_comp - n_same
    return n_diff, n_comp

def calc_pi_dxy(vcf_file, pops, window):
    """Computes the number of diff and comps for all the populations in a defined window
       Returns two numpy array, one for pi and one for dxy, with the diffs and comps for all the populations"""
    
    tmp_pi = np.zeros((len(pops), 2), dtype=int)
    tmp_dxy = np.zeros((int(comb(len(pops), 2)), 2), dtype=int)
    
    window_string = f"{window[0]}:{window[1]}-{window[2]}" if window is not None else None
        
    vcf = VCF(vcf_file)
    
    for variant in vcf(window_string):

        for i, pop in enumerate(pops):
            
            n_diff, n_comp = diff_comp_within(variant, pops[pop])

            tmp_pi[i, 0] += n_diff
            tmp_pi[i, 1] += n_comp

        for i, pop in enumerate(combinations(pops, 2)):
            
            n_diff, n_comp = diff_comp_between(variant, pops[pop[0]], pops[pop[1]])
            
            tmp_dxy[i, 0] += n_diff
            tmp_dxy[i, 1] += n_comp

    return tmp_pi, tmp_dxy


if __name__ == "__main__":
            
   main()