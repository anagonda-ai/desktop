import csv
import sys
import os
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq

def csv_to_fasta(csv_file, fasta_file):
    with open(csv_file, 'r') as infile, open(fasta_file, 'w') as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            first_part = row['Locus_Tag'] if row['Locus_Tag'] else row['Protein_ID'] if row['Protein_ID'] else row['Gene']
            header = f">{first_part} | {os.path.basename(csv_file).replace('.csv', '')} | {row['Start']} | {row['End']}"
            # Convert nucleotide sequence to protein sequence
            protein_sequence = row['Translation']
            outfile.write(f"{header}\n{protein_sequence}\n")

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
    
    input_dir = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/csv_files"
    output_dir = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files"
    process_directory(input_dir, output_dir)