#!/usr/bin/python
# -*- coding: utf-8 -*-

# IMPORTS
import argparse
import os
import subprocess

import pandas as pd

"""This module uses an input multisample VCF, reference genome FASTA, and a
GTF file containing reference genome coding sequences and calculates dN/dS
using an R package PopGenome"""


def parse_dnds_result(chromosome, output_dir, output_file, gene_info_path,
                      nonsyn_taj_path, syn_taj_path):

    result_file = output_file + chromosome + "_dnds.txt"
    result_file_path = os.path.join(output_dir, result_file)

    # Opening result files as Pandas DataFrames
    gene_info_df = pd.read_table(gene_info_path, sep='\t', header=1,
                                 index_col=0, skip_blank_lines=True)
    gene_info_df.columns = ['Gene']

    nonsyn_taj_df = pd.read_table(nonsyn_taj_path, sep='\t', header=1,
                                  index_col=0, skip_blank_lines=True)
    nonsyn_taj_df.columns = ['nonsynTaj']

    syn_taj_df = pd.read_table(syn_taj_path, sep='\t', header=1, index_col=0,
                               skip_blank_lines=True)
    syn_taj_df.columns = ['synTaj']

    dfs = [gene_info_df, nonsyn_taj_df, syn_taj_df]
    result_df = pd.concat(dfs, axis=1)

    # Picking up gene names and IDs
    result_df['GeneID'] = (result_df.Gene.str.split(';').str.get(
        1)).str.split('=').str.get(1)
    result_df['GeneName'] = (result_df.Gene.str.split(';').str.get(
        3)).str.split('=').str.get(1)

    # Calculating dN/dS
    result_df['dN/dS'] = result_df.nonsynTaj / result_df.synTaj

    # Saving to gile
    result_df.to_csv(result_file_path, sep='\t')

    return result_df


def run_rscript(chromosome, rscript, output_dir, output_file, slurm_script):
    # Writing Rscript
    rscript_name = chromosome + output_file + "_Rscript.R"
    rscript_file = os.path.join(output_dir, rscript_name)
    with open(rscript_file, 'w') as ro:
        ro.write(rscript)

    # Writing SLURM file
    slurm_file_name = chromosome + output_file + "_SLURM.sh"
    slurm_file = os.path.join(output_dir, slurm_file_name)
    with open(slurm_file, 'w') as slurm:
        slurm.write(slurm_script)

    # Running Rscript
    cmd = "sbatch %s" % slurm_file
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    std, err = p.communicate()
    # std = std.decode("utf-8")
    # err = err.decode("utf-8")
    return rscript_file


def format_slurm(rscript_file_path, chromosome, output_file):

    template = """
    #!/bin/bash
    #SBATCH --job-name={chromosome}_{output_name} # Job name
    #SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
    #SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
    #SBATCH -n 1
    #SBATCH -t 96:00:00
    #SBATCH --qos=normal
    srun R CMD BATCH ./{rscript_file_path}
    """
    context = {
        rscript_file_path: rscript_file_path,
        chromosome: chromosome,
        output_file: output_file
    }

    return template.format(**context)


def format_popgenome_rscript(vcf, reference, gff, output_dir, output_file,
                             chromosome):

    template = """
    ##########################################################################
    # Requirements: PopGenome
    #!/usr/bin/env Rscript

    setwd("{output_dir}")

    library("PopGenome")

    # As long as a set.gff function is not implemented we have to read in the
    # data with the corresponding GFF file in order to verify syn & nonsyn
    # SNPs afterwards.
    # Reading in the data with the corresponding GFF file
    GENOME.class <- readVCF("{vcf}", 1000,"{chromosome}",1,60000000, include.unknown=TRUE, gffpath="{gff}")

    # Set syn & nonsyn SNPs: The results are stored
    # in the slot GENOME.class@region.data@synonymous
    # The input of the set.synnonsyn function is an object of
    # class GENOME and a reference chromosome in FASTA format.
    GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="{reference}")

    # number of synonymous changes
    sum(GENOME.class@region.data@synonymous[[1]]==1, na.rm=TRUE)

    # number of non-synonymous changes
    sum(GENOME.class@region.data@synonymous[[1]]==0, na.rm=TRUE)

    # Here, we have to define the parameter na.rm=TRUE because NaN values in this slot
    # indicate that the observed SNP is in a non-coding region

    # We now could split the data into gene regions again
    genePos  <- get_gff_info(gff.file="{gff}", chr="{chromosome}", feature="gene")
    genes <- splitting.data(GENOME.class, positions=genePos, type=2)
    genes <- splitting.data(GENOME.class, subsites="gene")

    # Now we perform The Tajima’s D statistic on the whole data set and
    # consider only nonsyn SNPs in each gene/region.
    genes <- neutrality.stats(genes, subsites="nonsyn", FAST=TRUE)
    nonsynTaj <- genes@Tajima.D
    write.table(nonsynTaj, "{chromosome}_nonsynTaj.txt", sep="\t")

    # The same now for synonymous SNPs
    genes <- neutrality.stats(genes, subsites="syn", FAST=TRUE)
    synTaj <- genes@Tajima.D
    write.table(synTaj, "{chromosome}_synTaj.txt", sep="\t")

    # To have a look at the differences of syn and nonsyn Tajima D values in each gene we
    # could do the following plot:
    plot(nonsynTaj, synTaj, main="Chr{chromosome}: Genes: Tajima’s D ")

    gene_info <- get_gff_info(gff.file="{gff}", chr="{chromosome}", feature="gene", extract.gene.names=TRUE)

    write.table(gene_info, "{chromosome}_gene_info.txt", sep="\t")

    ##########################################################################
    """
    context = {
        "input_dir": output_dir,
        "output_file": output_file,
        "chr": chromosome,
        "vcf": vcf,
        "reference": reference,
        "gff": gff
    }

    return template.format(**context)


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
    parser.add_argument('--reference, -r', type=str, dest='reference',
                        action='store',
                        help='Path to the reference FASTA')
    parser.add_argument('--gff', type=str, dest='gff',
                        action='store',
                        help='Path to a GFF file containing reference genome'
                             'coding sequences')
    parser.add_argument('--output_dir', type=str, dest='output_dir',
                        action='store',
                        help='Path to output directory')
    parser.add_argument('--output_file', type=str, dest='output_file',
                        action='store', help='Output file prefix')
    parser.add_argument('--chromosome', type=str, dest='chromosome',
                        action='store', help='Chromosome number')

    # Convert argument strings to objects and assign them as attributes of the
    # namespace. Return the populated namespace.
    args = parser.parse_args()

    return (args.vcf, args.reference, args.gff, args.output_dir,
            args.output_file, args.chromosome)


def main():

    # Parsing arguments
    (vcf, reference, gff, output_dir, output_file, chromosome
     ) = parse_arguments()

    # Formatting PopGenome Rscript
    rscript = format_popgenome_rscript(vcf, reference, gff, output_dir,
                                       output_file, chromosome)

    # Formatting SLURM script
    rscript_file_name = chromosome + output_file + "_Rscript.R"
    rscript_file_path = os.path.join(output_dir, rscript_file_name)
    slurm_script = format_slurm(rscript_file_path, chromosome, output_file)

    # Writing and intitiating Rscript
    run_rscript(chromosome, rscript, output_dir, output_file, slurm_script)

    # Formatting result output filenames
    gene_info = chromosome + "_gene_info.txt"
    gene_info_path = os.path.join(output_dir, output_file, gene_info)
    nonsyn_taj = chromosome + "_nonsynTaj.txt"
    nonsyn_taj_path = os.path.join(output_dir, output_file, nonsyn_taj)
    syn_taj = chromosome + "_synTaj.txt"
    syn_taj_path = os.path.join(output_dir, output_file, syn_taj)

    # Parsing PopGenome results
    result_df = parse_dnds_result(chromosome, output_dir, output_file,
                                  gene_info_path, nonsyn_taj_path,
                                  syn_taj_path)


if __name__ == "__main__":
    main()
