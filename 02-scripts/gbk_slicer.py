#!/usr/bin/env python3

# Required modules
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation

input = sys.argv[1]
ouput = sys.argv[2]
gene = sys.argv[3]
upstream = int(sys.argv[4])
downstream = int(sys.argv[5])

# change the following variables to match your file and gene of interest
#gb_file="/Net/Groups/ccdata/users/AIbrahim/Pseudomonas_spflanze/github/pseudomonas_Spflanze/04-analysis/blast/tblastn_corrected_sec/prokka/CP007410.1/CP007410.1.gbk"
#gene_id = "tycC_4"

gb_file=input
gene_id=gene
output_file=ouput

# open the GenBank file and read in the sequence record
with open(gb_file, "r") as input_handle:
    record = SeqIO.read(input_handle, "genbank")

# find the index of the CDS feature with the matching gene ID
matching_index = None
for i, feature in enumerate(record.features):
    if feature.type == "CDS" and feature.qualifiers.get("gene", [""])[0] == gene_id:
        matching_index = i
        break

# choose 20 CDS features before and after the matching one
features_to_copy = []
if matching_index is not None:
    start_index = max(0, matching_index - upstream) #41 means 20 !!
    end_index = min(len(record.features), matching_index + downstream)
    features_to_copy = record.features[start_index:end_index]

# create a new SeqRecord object for the extracted region and add all selected features to it
sub_seq = record.seq
output_record = SeqRecord(sub_seq, id=record.id, name=record.name, description="CDS features surrounding %s" % gene_id)
output_record.features = features_to_copy

# set the molecule type annotation to DNA
output_record.annotations["molecule_type"] = "DNA"

# print the output record for testing purposes
#print(output_record)

# write the output record to a new GenBank file
out_file = output_file
with open(out_file, "w") as output_handle:
    SeqIO.write(output_record, output_handle, "genbank")