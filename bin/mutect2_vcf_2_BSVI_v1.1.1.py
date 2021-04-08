"""
Script to modify VCFs from mutect2 for importing to BSVI.
Multiallelic sites need splitting to individual biallelic records, and
the genotype fields splitting from 0/0/1/0 -> 0/1 for BSVI to handle.

Requires bgzip & bcftools be installed and on path.

Inputs:
    - reference fasta file - required by bcftools norm
    - input VCF to be modified

Jethro Rainford
210311
"""
import io
from pathlib import Path
from shutil import which
import subprocess
import sys

import pandas as pd


def bcf_norm(ref_fasta, input_vcf):
    """
    Normalise multiallelic records using bcftools norm

    Args:
        - ref_fasta (file): reference fasta file to use for bcftools norm
        - input_vcf (file): vcf file passed at cmd line

    Returns:
        - vcf_header (list): header lines read in from VCF
        - vcf_df (df): df of variants
    """
    print("Calling bcftools")

    # call bcftools to normalise multiallelic sites
    process = subprocess.Popen((
        f"zcat {input_vcf} | sed 's/AD,Number=./AD,Number=R/g' | "
        "sed 's/RPA,Number=./RPA,Number=R/g' | "
        f"bcftools norm -f {ref_fasta} -m -both"
    ), shell=True, stdout=subprocess.PIPE)

    vcf_data = io.StringIO()

    vcf_header = []

    for line in process.stdout:
        line = line.decode()
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
            f'Genotype field has < 3 characters: {sample[0]} for sample: {row}'

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
    """
    Output vcf name comments:
        - Output vcf name shortened to keep only meaningful characters
        - Removing ngr_ prefix and markdup_recalibrated 
        - Add _ms suffix to denote that vcf has multiallelic variants
          split to separate lines e.g. 
          Input vcf: ngr_1701389_S77_markdup_recalibrated_tnhaplotyper2.vcf.gz 
          Output vcf: 1701389_S77_tnhaplotyper2_ms.vcf.gz
    """
    
    fname = str(Path(input_vcf).name).replace('.vcf', '_ms.vcf').rstrip(
        '.gz').replace('ngr_', '').replace('_markdup_recalibrated', '')

    print(f'Writing to outfile: {fname}.gz')

    with open(fname, 'w') as f:
        for line in vcf_header:
            # write header to vcf
            f.write(line)

    # apend variants to vcf & compress
    with open(fname, 'a') as f:
        vcf_df.to_csv(f, sep='\t', header=False, index=False)

    subprocess.Popen(f'bgzip {fname}', shell=True)


if __name__ == "__main__":

    # check bcftools is installed and on path
    assert which('bcftools'), 'bcftools is not installed / on path'

    # check only reference fasta and one vcf passed
    assert len(sys.argv) == 3, 'Incorrect no. VCFs passed, requires one.'

    # check reference fasta and vcf are passed correctly
    assert '.fa' in sys.argv[1] and '.vcf' in sys.argv[2], (
        'Files passed incorrectly, it should be:'
        ' python3 mutect2_vcf_2_bsvi.py reference_fasta input_vcf'
    )

    ref_fasta = sys.argv[1]
    input_vcf = sys.argv[2]

    vcf_header, vcf_df = bcf_norm(ref_fasta, input_vcf)
    write_file(input_vcf, vcf_header, vcf_df)
