import os
import sys

def merge_fasta_files(files, output_file):
    with open(output_file, 'w') as outfile:
        for file in files:
            with open(file, 'r') as infile:
                outfile.write(infile.read())

def find_fasta_files(base_dir):
    fasta_files = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))
                print("Found FASTA file:", os.path.join(root, file))
    return fasta_files

def main():
    if len(sys.argv) != 4:
        print("Usage: python merge_fasta_files.py <fasta_paths.txt> <output_dir> <file_name>")
        sys.exit(1)

    fasta_paths_file = sys.argv[1]
    output_dir = sys.argv[2]
    file_name = sys.argv[3]

    with open(fasta_paths_file, 'r') as f:
        files = [line.strip() for line in f if line.strip()]

    output_file = os.path.join(output_dir, f"{file_name}.fasta")
    
    # Merge the found FASTA files into one big file
    merge_fasta_files(files, output_file)
    print(f"Merged FASTA file written to: {output_file}")

if __name__ == '__main__':
    main()