#!/usr/bin/env python3

#########################################
# Author: Anan Ibrahim
# Date: May 2023
# Title: Filter the blast hits from BLASTn alignment
#########################################

# Required modules
import sys
import os
import pandas as pd
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning
    
# Define input arguments:
parser = argparse.ArgumentParser(prog = 'BLASTn HSPs Filter', formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=('''\
    ................................................................................
                                    *BLASTn HSPs Filter*
    ................................................................................
    This tool filters the BLASTn results given that per subject more than one HSPs
    is retrieved. Therfore, this tool retains only the HSPs with the longest
    alignment length. This is adjusted to only work on blastn format:

    #blastn -query $PROJECT_DIR/01-resources/mycin_nrps.fasta -db $NCBI/nt
    #-max_target_seqs 100
    #-num_threads 25
    #-task megablast 
    #-evalue 0.05
    #-outfmt "6 qseqid sseqid pident sacc mismatch sseq sstart send evalue sscinames nident mismatch length qcovs qcovhsp "
    #-out matches_nt_blastn_fna_qcov_local_5.txt

    To use you can run: 
    python3 blastn_filter_hsps.py 
    --blastn_file '/Net/Groups/ccdata/users/AIbrahim/Pseudomonas_spflanze/github/pseudomonas_Spflanze/04-analysis/blast/matches_nt_blastn_fna_qcov_local_5.txt' 
    --output_folder '/Net/Groups/ccdata/users/AIbrahim/Pseudomonas_spflanze/github/pseudomonas_Spflanze/04-analysis/blast/' 
    .................................................................................'''),
                                epilog='''Thank you for running GBK_extractor!''',
                                add_help=True)

parser.add_argument("--blastn_file", dest="file", nargs='?', help="Enter the path to the .txt blastn file. Eg. The txt file generated from blastn with the outfmt 6 described in the description. '/Net/Groups/ccdata/users/AIbrahim/Pseudomonas_spflanze/github/pseudomonas_Spflanze/04-analysis/blast/matches_nt_blastn_fna_qcov_local_5.txt'",
                    type=str)
parser.add_argument("--output_folder", dest="outputpath", nargs='?', help="Enter the path to where the new filtered blastn file. By default is called 'matches_nt_blastn_fna_qcov_local_filtered.txt'",
                    type=str)

# get command line arguments
args = parser.parse_args()

file = args.file
out_path = args.outputpath


#########################################
# Filter the blastn results 
#########################################
def filter(file, out_path):
    out_file='matches_nt_blastn_fna_qcov_local_filtered.txt'
    # open the input file
    input=pd.read_csv(file, sep='\t', header=None)
    # Note: that the qcov result is not to be trusted - i have double checked with the manual and the [5]blasthits and its correct nonetheless
    # Subset accession number by grabing the one with the higest alignment length
    grouped = input.groupby(3)[11].max().reset_index()
    # Then keep the subsets by merging the df
    input_filtered = pd.merge(input, grouped , how='right', on=[3,11])
    # Add headers to the file and write to a directory 
    input_filtered.columns = ["qseqid", 
                              "sseqid", 
                              "pident", 
                              "sacc",
                              "mismatch",
                              "sseq",
                              "sstart",
                              "send",
                              "evalue",
                              "sscinames",
                              "nident",
                              "length",
                              "qcovs",
                              "qcovhsp"]
    input_filtered.to_csv(out_path + out_file, sep='\t', index=False)

#########################################
# Main Function
#########################################

def main_workflow():
    filter(file, out_path)

if __name__ == "__main__":
    main_workflow()
