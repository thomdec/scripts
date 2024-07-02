# Scripts for evolutionary genomics

## Get high 

The script `get_high.py` extracts the allele frequency and effect of so-called "HIGH impact variants" from a VCF annotated with SNPEff. 

Main dependencies: [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/) and [bcftools](https://samtools.github.io/bcftools/)

Arguments: 

    - Required:
        - A VCF annotated with SNPEff
        - The output name
    - Optional:
        - Region string as per bcftools: chr|chr:pos|chr:beg-end
        - Population file (Csv with id and population)

### Example command

```bash
python get_high.py -v annotated.vcf.gz -r chr1:100-1000 -p pops.csv -o output.tsv
```



