#!/bin/bash

# Bash command to call Python script that normalises multiallelic variants in VCF
# for importing in to BSVI

# Inputs python_packages.tar.gz, mutect2_vcf_2_BSVI.py script, reference fasta file
# and a VCF output from mutect2

# Output: modified VCF with normalised multiallelic variants

set -exo

echo "Starting app"

# install required python packages from local tar
tar xvf python_packages.tar.gz
pip3 install pytz-*.whl numpy-*.whl pandas-*.whl

# remove installation packages so they don't get uploaded
rm -f *.whl

ref=$(find . -name "*.fa.gz")
vcf=$(find . -name "*.vcf*")

# call Python script to normalise and generate new VCF
python3 mutect2_vcf_2_BSVI_v*.py $ref $vcf

echo "Done"
