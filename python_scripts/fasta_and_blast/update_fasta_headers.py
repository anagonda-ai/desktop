import os
import re

def extract_positions(header):
    # Extract the start and end positions from the FASTA header
    parts = header.split()
    if len(parts) > 2:
        chromosome_info = parts[2]
        chromosome_parts = chromosome_info.split(':')
        if len(chromosome_parts) > 4:
            start = chromosome_parts[3]
            end = chromosome_parts[4]
            print (start, end)
            return int(start), int(end)
    return None, None

def update_fasta_headers(fasta_file, output_fasta_file):
    with open(fasta_file, 'r') as f_in, open(output_fasta_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                # Extract the ID from the FASTA header
                fasta_id = line.split()[0][1:]
                start, end = extract_positions(line)
                if start and end:
                    # Update the FASTA header with start and end positions
                    updated_line = f">{fasta_id} | start={start} | end={end}\n"
                    f_out.write(updated_line)
                else:
                    f_out.write(line)
            else:
                f_out.write(line)

def process_fasta_files(base_dir):
    # Iterate through each .pep.all.fa file in the base directory and subdirectories
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".pep.all.fa"):
                fasta_path = os.path.join(root, file)
                output_fasta_file = os.path.join(root, f"updated_{file}")
                update_fasta_headers(fasta_path, output_fasta_file)
                print(f"Created updated FASTA file: {output_fasta_file}")

def main():
    base_dir = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl"
    process_fasta_files(base_dir)

if __name__ == "__main__":
    main()