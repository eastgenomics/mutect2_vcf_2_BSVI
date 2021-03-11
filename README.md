# mutect2 VCF 2 BSVI


## What does this do?

Script to normalise multiallelic variants in VCF output from mutect2 for
importing in to BSVI.

Required as BSVI is terrible and can't handle multiallelic variants with greater
than 2 entries in genotype field.


## What are typical use cases for this?

Takes output VCF of mutect2, splits multiallelic variants into individual
biallelic lines using [bcftools norm][bcftools-url] then changes genotype entries to `0/1`.


## What is required for this to run?

- Python 3.x
- pandas
- bcftools


## What does this output?

Modified VCF to import into BSVI.

[bcftools-url]: http://samtools.github.io/bcftools/bcftools.html#norm