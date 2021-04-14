"""
Script to modify VCFs from mutect2 for importing to BSVI.
Multiallelic sites need splitting first with bcftools norm and then the
genotype fields are split from 0/0/1/0 -> 0/1 for BSVI to handle.

Outputs vcf with split multiallelics and a tsv with formatted INFO
column to individual columns.

Inputs:
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


def mod_genotype(input_vcf):
    """
    Overwrite GT with 0/1 when it's derived from a multiallelic
    e.g. 0/0/1/0 -> 0/1 and 0/0/0/1 -> 0/1

    Args:
        - input_vcf (file): vcf file passed at cmd line

    Returns:
        - vcf_header (list): header lines read in from VCF
        - vcf_df (df): df of variants
    """
    print("Calling bcftools")

    # stream the vcf
    process = subprocess.Popen(
        f"cat {input_vcf} ", shell=True, stdout=subprocess.PIPE
    )

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


def generate_tsv(vcf_df):
    """
    Generates tsv file from modified vcf with INFO column split out to
    individual columns. Expects min. 9 '|' separated fields in INFO
    column to split out.

    Args:
        - vcf_df (df): df of variants from vcf
    
    Returns:
        - tsv_df (df): df of variants with split info column
    """
    info_cols = [
        'GENE', 'VARIANT_CLASS', 'CONS', 'EXON', 'HGVSc', 'HGVSp',
        'gnomAD_AF', 'SIFT', 'POLYPHEN', 'DB'
    ]

    # sense check correct annotation has been added to all rows else it
    # gives an unhelpful pandas error on trying to split
    assert all(vcf_df.INFO.str.count('\|') > 8), \
        "Incorrectly formatted INFO field, some records have < 9 fields."

    # splits info column to cols defined in info_cols, any missing
    vcf_df[info_cols] = vcf_df['INFO'].str.split('|', 9, expand=True)

    # remove info id from gene
    vcf_df['GENE'] = vcf_df['GENE'].apply(lambda x: x.replace('CSQ=', ''))

    # split messy db annotation field out to clinvar, cosmic & dbsnp
    # cols have multiple fields and diff delimeters then join with ','
    # in case of having more than one entry
    vcf_df['COSMIC'] = vcf_df['DB'].str.split(r'\&|\||,').apply(
        lambda x: ','.join((y for y in x if y.startswith('COS')))
    )
    vcf_df['CLINVAR'] = vcf_df['DB'].str.split(r'\&|\||,').apply(
        lambda x: ','.join((y for y in x if y.startswith('CM')))
    )
    vcf_df['dbSNP'] = vcf_df['DB'].str.split(r'\&|\||,').apply(
        lambda x: ','.join((y for y in x if y.startswith('rs')))
    )

    # add interestingly formatted report text column
    vcf_df['Report_text'] = vcf_df[vcf_df.columns.tolist()].apply(
        lambda x: (
            f"{x['GENE']} {x['VARIANT_CLASS']} variant in {x['EXON']} \\n"
            f"{x['HGVSc']} \\n {x['HGVSp']} \\n COSMIC ID: {x['COSMIC']} \\n"
            f"Allele Frequency: {x['gnomAD_AF']}"
        ), axis=1
    )

    # drop unneeded columns
    vcf_df = vcf_df.drop(['INFO', 'DB'], axis=1)

    return vcf_df

def write_files(input_vcf, vcf_header, vcf_df, tsv_df):
    """
    Write modified vcf and tsv df with split info field to files

    Args:
        - input_vcf (file): vcf file passed at cmd line
        - vcf_header (list): header lines read in from VCF
        - vcf_df (df): df of variants
        - tsv_df (df): df of variants with split info column

    Outputs:
        - vcf file with modified multiallelic records
        - tsv file with modified multialleic records and split info field
    """
    vcf_fname = str(Path(input_vcf).name).replace('vepfilter', 'bsvi')
    tsv_fname = vcf_fname.replace('.vcf', '.tsv')

    print(f'Writing vcf to outfile: {vcf_fname}')
    print(f'Writing tsv to outfile: {tsv_fname}')

    with open(vcf_fname, 'w') as f:
        for line in vcf_header:
            # write header to vcf
            f.write(line)

    # apend variants to vcf
    with open(vcf_fname, 'a') as f:
        vcf_df.to_csv(f, sep='\t', header=False, index=False)
    
    # write tsv file
    with open(tsv_fname, 'w') as tsv:
        tsv_df.to_csv(tsv, sep='\t', header=True, index=False)


if __name__ == "__main__":

    # check only one vcf passed
    assert len(sys.argv) == 2, 'Incorrect no. VCFs passed, requires one.'

    input_vcf = sys.argv[1]

    vcf_header, vcf_df = mod_genotype(input_vcf)
    tsv_df = generate_tsv(vcf_df)
    write_files(input_vcf, vcf_header, vcf_df, tsv_df)
