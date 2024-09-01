#!/usr/bin/env python

from Bio import SeqIO
import sys
import os

# Input directory containing GenBank files
input_directory = "./MGCs/"

# Output file for combined FASTA sequences
output_file = "./MGCs/combined.fasta"

# Open output file for writing
with open(output_file, "w") as output_fh:
    # Iterate over GenBank files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".gbk"):
            filepath = os.path.join(input_directory, filename)
            # Parse GenBank file and write sequences in FASTA format to output file
            for record in SeqIO.parse(filepath, "genbank"):
                SeqIO.write(record, output_fh, "fasta")

print("Combined FASTA file created:", output_file)
