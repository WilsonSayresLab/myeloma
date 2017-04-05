#!/usr/bin/python
# -*- coding: utf-8 -*-

# IMPORTS
import argparse
import subprocess

import pandas as pd

"""This module uses an input multisample VCF, reference genome FASTA, and a
GTF file containing reference genome coding sequences and calculates dN/dS
using an R package PopGenome"""


def parse_dnds_result(chr, output_dir, output_file, gene_info_path,
                      nonsynTaj_path, synTaj_path):

    result_file = output_file + chr + "_dnds.txt"
    result_file_path = os.path.join(output_dir, result_file)

    # Opening result files as Pandas DataFrames
    gene_info_df = pd.read_table(gene_info_path, sep='\t', header=1,
                                 index_col=0, skip_blank_lines=True)
    gene_info_df.columns = ['Gene']

    nonsynTaj_df = pd.read_table(nonsynTaj_path, sep='\t', header=1,
                                 index_col=0, skip_blank_lines=True)
    nonsynTaj_df.columns = ['nonsynTaj']

    synTaj_df = pd.read_table(synTaj_path, sep='\t', header=1, index_col=0,
                              skip_blank_lines=True)
    synTaj_df.columns = ['synTaj']

    dfs = [gene_info_df, nonsynTaj_df, synTaj_df]
    result_df = pd.concat(dfs, axis=1)

    # Picking up gene names and IDs
    result_df['GeneID'] = (result_df.Gene.str.split(';').str.get(1)).str.split(
        '=').str.get(1)
    result_df['GeneName'] = (result_df.Gene.str.split(';').str.get(3)).str.split(
        '=').str.get(1)

    # Calculating dN/dS
    result_df['dN/dS'] = result_df.nonsynTaj / result_df.synTaj

    # Saving to gile
    result_df.to_csv(result_file_path, sep='\t')

    return result_df


def run_rscript(rscript):
    # Writing Rscript
    rscriptFile = os.path.join(output_dir, output_file)
    with open(rscriptFile, 'w') as ro:
        ro.write(rscript)
    # Running Rscript
    cmd = "Rscript %s" % rscriptFile
    p = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std, err = p.communicate()
    #std = std.decode("utf-8")
    #err = err.decode("utf-8")
    return rscriptFile


def format_popgenome_rscript(vcf, reference, gff, output_dir, output_file,
                             chr):

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
    GENOME.class <- readVCF("{vcf}", 1000,"{chr}",1,60000000, include.unknown=TRUE, gffpath="{gff}")

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
    genePos  <- get_gff_info(gff.file="{gff}", chr="{chr}", feature="gene")
    genes <- splitting.data(GENOME.class, positions=genePos, type=2)
    genes <- splitting.data(GENOME.class, subsites="gene")

    # Now we perform The Tajima’s D statistic on the whole data set and
    # consider only nonsyn SNPs in each gene/region.
    genes <- neutrality.stats(genes, subsites="nonsyn", FAST=TRUE)
    nonsynTaj <- genes@Tajima.D
    write.table(nonsynTaj, "{chr}_nonsynTaj.txt", sep="\t")

    # The same now for synonymous SNPs
    genes <- neutrality.stats(genes, subsites="syn", FAST=TRUE)
    synTaj <- genes@Tajima.D
    write.table(synTaj, "{chr}_synTaj.txt", sep="\t")

    # To have a look at the differences of syn and nonsyn Tajima D values in each gene we
    # could do the following plot:
    plot(nonsynTaj, synTaj, main="Chr{chr}: Genes: Tajima’s D ")

    gene_info <- get_gff_info(gff.file="{gff}", chr="{chr}", feature="gene", extract.gene.names=TRUE)

    write.table(gene_info, "{chr}_gene_info.txt", sep="\t")

    ##########################################################################
    """
    context = {
        "input_dir": output_dir,
        "output_file": output_file,
        "chr": chr,
        "vcf":vcf,
        "reference":reference,
        "gff":gff
    }

    return (template.format(**context))


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
    parser.add_argument('--chr', type=str, dest='chr',
                        action='store', help='Chromosome number')

    # Convert argument strings to objects and assign them as attributes of the
    # namespace. Return the populated namespace.
    args = parser.parse_args()

    return (args.vcf, args.ref, args.gff, args.output_dir, args.output_file,
            args.chr)


def main():

    # Parsing arguments
    (vcf, ref, gff, output_dir, output_file, chr) = parse_arguments()

    # Formatting PopGenome Rscript
    rscript = format_popgenome_rscript(vcf, reference, gff, output_dir, output_file, chr)

    # Writing and intitiating Rscript
    run_rscript(rscript)

    # Formatting result output filenames
    gene_info = chr + "_gene_info.txt"
    gene_info_path = os.path.join(output_dir, output_file, gene_info)
    nonsynTaj = chr + "_nonsynTaj.txt"
    nonsynTaj_path = os.path.join(output_dir, output_file, nonsynTaj)
    synTaj = chr + "_synTaj.txt"
    synTaj_path = os.path.join(output_dir, output_file, synTaj)

    # Parsing PopGenome results
    result_df = parse_dnds_result(gene_info_path, nonsynTaj_path, synTaj_path)


if __name__ == "__main__":
    main()
