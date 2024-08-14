# Scripts for population genomics

## Yayawp

The script `yayawp.py` computes $\pi$ and $d_{xy}$ in windows or globally from an all-site VCF. It takes advantage of cyvcf2 to read a vcf site per site. As a result it is relatively slow, but takes up basicallly no RAM (it can be used to compute global estimates without causing memory issues, even on laptops).

Dependencies: python>=3.10, numpy and [cyvcf2](https://github.com/brentp/cyvcf2)

Arguments: 

    - Required:
        - An indexed all-site VCF
        - Population file (csv with id and population)
        - The output prefix
    - Optional:
        - Window size

Output:

    - a .pi.tsv with: [chromosome, start, end,] pop, pi, n_diff, n_comp
    - a .dxy.tsv with: [chromosome, start, end,] pop1, pop2, dxy, n_diff, n_comp


### Example command

```bash
python yayawp.py -v vcf.gz -p pops.csv -o global.tsv

python yayawp.py -v vcf.gz -p pops.csv -w 500_000 -o windows.tsv
```

## Get high 

The script `get_high.py` extracts the allele frequency and effect of so-called "HIGH impact variants" from a VCF annotated with SNPEff. 

Main dependencies: [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/) and [bcftools](https://samtools.github.io/bcftools/)

Arguments: 

    - Required:
        - A VCF annotated with [SNPEff](https://pcingola.github.io/SnpEff/)
        - The output name
    - Optional:
        - Region string as per bcftools: chr|chr:pos|chr:beg-end (vcfs need to be indexed with tabix or bcftools)
        - Population file (csv with id and population)

Output:

    - A tsv file with chromosome, position, [population,] total number of alleles, number of "high impact" alleles, effect of the "high impact" allele (eg. stop_gained, start_lost, splice_donor_variant...)


### Example command

```bash
python get_high.py -v annotated.vcf.gz -r chr1:100-1000 -p pops.csv -o output.tsv
```

## PCA

The script `pca.py` make use of [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/) to run a Principal Component Analysis.

Main dependency: [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/)

Arguments: 

    - Required:
        - A VCF
        - The output prefix
    - Optional:
        - Number of components [Default = 5]
        - Region string: chr|chr:pos|chr:beg-end (works best with vcfs indexed with tabix or bcftools)
        - Sample file (newline separated samples)

Output:

    - A tsv file with ids and PC coordinates
    - A tsv file with the proportion of variance explained by each components

### Example command

```bash
python pca.py -v vcf.gz -n 10 -r chr1:100-1000 -s samples.txt -o output_prefix
```

