"""
Script to modify VCFs from mutect2 for importing to BSVI.
Multiallelic sites need splitting to individual biallelic records, and
the genotype fields splitting from 0/0/1/0 -> 0/1 for BSVI to handle.

Requires bcftools be installed and on path.

Jethro Rainford
210311
"""
import io
from pathlib import Path
from shutil import which
import subprocess
import sys

import pandas as pd


def bcf_norm(input_vcf):
    """
    Normalise multiallelic records using bcftools norm

    Args:
        - input_vcf (file): vcf file passed at cmd line

    Returns:
        - vcf_header (list): header lines read in from VCF
        - vcf_df (df): df of variants
    """
    print("Calling bcftools")

    # call bcftools to normalise multiallelic sites
    process = subprocess.Popen(
        f'bcftools norm -m -both {input_vcf}', shell=True, stdout=subprocess.PIPE
    )

    vcf_data = io.StringIO()

    vcf_header = []

    for line in process.stdout:
        line = line.decode().strip('"\n') + '\n'
        vcf_data.write(line)

        if line.startswith('#'):
            # dump out header to list to write back
            vcf_header.append(line)

    vcf_data.seek(0)

    cols = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", "SAMPLE"
    ]

    # read normalised vcf into df
    vcf_df = pd.read_csv(vcf_data, sep="\t", comment='#', names=cols)

    print('Adjusting multiallelic genotypes')
    for index, row in vcf_df.iterrows():
        # loop over rows, change genotype if contains greater than 2 fields
        sample = row['SAMPLE'].split(':')

        # sense check that genotype field doesn't have anything funky,
        # if it does then it can be reviewed manually
        assert len(sample[0]) >= 3, \
            f'Genotype field has < 3 characters: {sample[0]}'

        if len(sample[0]) > 3:
            # >3 => not 0/1 => modify
            sample[0] = '0/1'
            sample = ':'.join(sample)

            # write new entry back to row
            vcf_df.at[index, 'SAMPLE'] = sample

    return vcf_header, vcf_df


def write_file(input_vcf, vcf_header, vcf_df):
    """
    Write modified vcf to file

    Args:
        - input_vcf (file): vcf file passed at cmd line
        - vcf_header (list): header lines read in from VCF
        - vcf_df (df): df of variants

    Outputs:
        - vcf file with modified multiallelic records
    """
    # set name for output vcf from input
    fname = str(Path(input_vcf).name).replace(
        '.vcf', '_multiallelic_split.vcf').rstrip('.gz')

    print(f'Writing to outfile: {fname}.gz')

    with open(fname, 'w') as f:
        for line in vcf_header:
            # write header to vcf
            f.write(line)

    # apend variants to vcf & compress
    with open(fname, 'a') as f:
        vcf_df.to_csv(f, sep='\t', header=False, index=False)

    subprocess.Popen(f'gzip {fname}', shell=True)


if __name__ == "__main__":

    assert which('bcftools')  # check bcftools is installed and on path

    input_vcf = sys.argv[1]
    vcf_header, vcf_df = bcf_norm(input_vcf)
    write_file(input_vcf, vcf_header, vcf_df)
