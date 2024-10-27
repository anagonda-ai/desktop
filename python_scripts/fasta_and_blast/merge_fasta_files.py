import os

def merge_fasta_files(files, output_file):
    with open(output_file, 'w') as outfile:
        for file in files:
            with open(file, 'r') as infile:
                outfile.write(infile.read())

def find_fasta_files(base_dir, pattern):
    fasta_files = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.startswith(pattern) and file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))
                print("Found FASTA file:", os.path.join(root, file))
    return fasta_files

def main():
    base_dir = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza"
    pattern = "start_end_"
    output_file = os.path.join(base_dir, "plaza.fasta")
    
    # Find all start_end_*.fasta files in the base directory and subdirectories
    files = find_fasta_files(base_dir, pattern)
    
    # Merge the found FASTA files into one big file
    merge_fasta_files(files, output_file)
    print(f"Merged FASTA file written to: {output_file}")

if __name__ == '__main__':
    main()