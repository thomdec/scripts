import allel
import numpy as np
import pandas as pd
import argparse
from warnings import filterwarnings
import os

def main():

    filterwarnings('ignore') 

    parser = argparse.ArgumentParser(description = "PCA", prog='tool', formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog,max_help_position=50))

    parser._action_groups.pop()
    required = parser.add_argument_group('Required parameters')
    optional = parser.add_argument_group('Optional parameters')

    required.add_argument('-v', '--vcf', type=str, help='VCF file', required=True)
    required.add_argument('-o', '--output', type=str, help='Output name prefix', required=True)

    optional.add_argument('-n', '--n_components', type=int, help='Number of components [Default = 5]', default=5, required=False)
    optional.add_argument('-r', '--region', type=str, help='Region string: chr|chr:pos|chr:beg-end [Best with indexed VCF]', default=None, required=False)
    optional.add_argument('-s', '--samples', type=str, help='Samples file. One sample per row', default=None, required=False) #or populations

    args = parser.parse_args()

    # checks
    if os.path.exists(args.vcf) is not True:
        raise Exception('ERROR: The specified VCF ' + str(args.vcf) + ' does not exist') 

    if os.access(os.path.dirname("./" + args.output), os.W_OK) is not True:
        raise Exception(f"ERROR: Cannot write output to {args.output}.coords.csv")

    # parse optional arguments
    if args.samples is not None:
        if os.path.exists(args.samples) is not True:
            raise Exception('ERROR: The specified samples file ' + str(args.samples) + ' does not exist')

        sample_df = pd.read_csv(args.samples, names=["id"])
        vcf_header = allel.read_vcf_headers(args.vcf)

        if not np.all(np.isin(sample_df.id, vcf_header.samples)):
            raise Exception(f"ERROR: The following sample could not be found in the VCF:\n{', '.join(sample_df.id[~np.isin(sample_df.id, vcf_header.samples)])}")
        
        sample_string = sample_string = list(sample_df.id)

    else: 
        sample_string = None

    region_string = args.region if args.region is not None else None

    # import vcf 
    callset = allel.read_vcf(args.vcf,  fields=['calldata/GT', 'samples'], samples=sample_string, region=region_string)
    g = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT'])) # genotype array
    samples = callset['samples']

    # filter out singletons and multi-allelic sites
    ac = g.count_alleles()
    mask = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1) # change the multi-allelic command
    masked_g = g[mask]

    # create the PCA input (i.e. number of alternative alleles per individual at each site)
    gn = masked_g.to_n_alt()

    # run the PCA
    coords, model = allel.pca(gn, n_components=args.n_components, copy=True, scaler='patterson', ploidy=2)

    coords_df = pd.DataFrame(coords, columns=["PC" + str(i + 1) for i in range(coords.shape[1])])
    samples_df = pd.DataFrame(samples, columns=["id"])
    pca = pd.concat([samples_df, coords_df], axis=1)

    # output to file
    pca.to_csv(f"{args.output}.coords.csv", index=False, sep='\t')

    with open(f"{args.output}.explained_variance.txt", 'w') as f:
        f.write("component\tvariance_explained\n")
        for idx, line in enumerate(model.explained_variance_ratio_.tolist()):
            f.write(f"{idx + 1}\t{line}\n")

if __name__ == "__main__":

   main()