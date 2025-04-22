import csv
import sys
import os
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq

def csv_to_fasta(csv_file, fasta_file):
    with open(csv_file, 'r') as infile, open(fasta_file, 'w') as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            header = f">{row['Gene_ID']}${row['Annotation']}"
            # Convert nucleotide sequence to protein sequence
            nucleotide_sequence = Seq(row['Nucleotide_Sequence'])
            sequence = str(nucleotide_sequence.translate())
            outfile.write(f"{header}\n{sequence}\n")

def process_file(csv_file, input_dir, output_dir):
    relative_path = os.path.relpath(os.path.dirname(csv_file), input_dir)
    fasta_subdir = os.path.join(output_dir, relative_path)
    os.makedirs(fasta_subdir, exist_ok=True)
    fasta_file = os.path.join(fasta_subdir, f"{os.path.splitext(os.path.basename(csv_file))[0]}.fasta")
    csv_to_fasta(csv_file, fasta_file)

def process_directory(input_dir, output_dir):
    csv_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.csv'):
                csv_files.append(os.path.join(root, file))

    with ThreadPoolExecutor() as executor:
        for csv_file in csv_files:
            executor.submit(process_file, csv_file, input_dir, output_dir)

if __name__ == "__main__":
    
    input_dir = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations"
    output_dir = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_fasta"
    process_directory(input_dir, output_dir)