# mutect2 VCF 2 BSVI


## What does this do?

Script to normalise multiallelic variants in VCF output from mutect2 for
importing in to BSVI.

Required as BSVI is terrible and can't handle multiallelic variants with greater
than 2 entries in genotype field.


## What are typical use cases for this?

Takes output VCF of mutect2, splits multiallelic variants into individual
biallelic records using [bcftools norm][bcftools-url], then modifies genotype
entries to `0/1`.

This should correctly split all multi allelic records into biallelic ones, with the
alt specific records split out to the appropriate new record.

To add as part of DNAnexus workflow this should be used with the Swiss Army
Knife app (v.4.1.1). 


## What is required for this to run?

- Python 3.x
- pandas
- bcftools
- bgzip

To run with Swiss Army Knife app:
- Inputs:
    - python_packages.tar.gz (contains: pandas, numpy, & pytz)
    - mutect2_vcf_2_BSVI_v*_.py
    - swiss_army_cmd_line_v*_.sh
    - reference fasta file
    - mutect2 VCF

- CMD to run: `sh swiss_army_cmd_line_v*.sh`


## What does this output?

Modified VCF to import into BSVI.

[bcftools-url]: http://samtools.github.io/bcftools/bcftools.html#norm
