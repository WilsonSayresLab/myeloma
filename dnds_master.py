#!/usr/bin/python
# -*- coding: utf-8 -*-

# IMPORTS
import argparse
import subprocess

"""This module uses an input VCF file and a reference genome
coding sequence FASTA file to create a new FASTA from a reference
FASTA and a VCF with GATK, creates a multisample FASTA file,
creates a codeml input file using fasta2paml.py, generates a codeml
input file, and runs PAML codeml to estimate dN/dS"""


def run_paml_codeml(paml_path, codeml_input, output_dir, output_file):
    """This function runs codeml with PAML to calculate dN/dS"""


def create_codeml_input(fasta2paml_path, output_dir, output_file):
    """This function creates a codeml input file from a multisample
    FASTA using fasta2paml.py"""


def concatenate_fastas(ref, mt_fasta, output_dir, output_file):
    """Concatenating the reference and mutation FASTA sequences into
    one multisample FASTA file"""

    concat_cmd = ('cat file1.fasta file2.fasta > combined.fasta').format(
        ref=ref,
        mt_fasta=mt_fasta,
        output_dir=output_dir,
        output_file=output_file
    )


def create_mt_fasta(gatk_path, vcf, ref, output_dir, output_file):
    """This function creates a new FASTA file from a reference FASTA
    and a VCF with GATK"""

   #  java -jar GenomeAnalysisTK.jar \
   # -T FastaAlternateReferenceMaker \
   # -R reference.fasta \
   # -o output.fasta \
   # -L input.intervals \
   # -V input.vcf \
   # [--snpmask mask.vcf]

    run_gatk_cmd = ('cd {output_dir}'
                    'java -jar {gatk_path} '
                    '-T FastaAlternateReferenceMaker '
                    '-R reference.fasta -o output.fasta '
                    '-L input.intervals -V input.vcf '
                    '[--snpmask mask.vcf];').format(
        output_dir=output_dir,
        gatk_path=gatk_path,
        )

    # Returning path to mt FASTA
    return mt_fasta


def parse_arguments():

    """Parsing command line arguments"""

    # Description for the help doc
    parser = argparse.ArgumentParser(
        description='This module uses an input VCF file and a reference '
                    'genome coding sequence FASTA file to create a new FASTA'
                    'from a reference FASTA and a VCF with GATK, creates a '
                    'multisample FASTA file, creates a codeml input file '
                    'using fasta2paml.py, generates a codeml input file, '
                    'and runs PAML codeml to estimate dN/dS')

    # Adding arguments
    parser.add_argument('--vcf', type=str, dest='vcf',
                        action='store',
                        help='Path to the input VCF file')
    parser.add_argument('--ref', type=str, dest='vcf',
                        action='store',
                        help='Path to the reference FASTA')
    parser.add_argument('--output_dir', type=str, dest='output_dir',
                        action='store',
                        help='Path to output directory')
    parser.add_argument('--output_file', type=str, dest='output_file',
                        action='store', help='Output file prefix')
    parser.add_argument('--gatk_path', action='store',
                        dest='gatk_path', default=None, nargs='?',
                        help='Path to GATK')
    parser.add_argument('--gtf', action='store',
                        dest='gtf', default=None, nargs='?',
                        help='Path to PAML')
    parser.add_argument('--fasta2paml_path', action='store',
                        dest='fasta2paml_path', default=None, nargs='?',
                        help='Path to fasta2paml.py')
    parser.add_argument('--paml_path', action='store',
                        dest='paml_path', default=None, nargs='?',
                        help='Path to PAML')

    # Convert argument strings to objects and assign them as attributes of the
    # namespace. Return the populated namespace.
    args = parser.parse_args()

    return (args.vcf, args.ref, args.output_dir, args.output_file,
            args.gatk_path, args.fasta2paml_path, args.paml_path)


def main():

    # Parsing arguments

    (vcf, ref, output_dir, output_file, gatk_path, fasta2paml_path,
     paml_path) = parse_arguments()


if __name__ == "__main__":
    main()