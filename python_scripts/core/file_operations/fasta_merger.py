import os
import sys

def merge_fasta_files(files, output_file):
    with open(output_file, 'w') as outfile:
        for file in files:
            with open(file, 'r') as infile:
                for line in infile:
                    if line.startswith(">"):
                        filename = os.path.basename(file).replace(".fasta","")
                        outfile.write(f"{line.strip()}${filename}\n")
                    else:
                        outfile.write(line)

def find_fasta_files(base_dir):
    fasta_files = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))
                print("Found FASTA file:", os.path.join(root, file))
    return fasta_files

def main():


    fasta_paths_dir = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_test2_fasta"
    file_name = "merged_metabolic_pathways"

    files = find_fasta_files(fasta_paths_dir)

    output_file = os.path.join(fasta_paths_dir, f"{file_name}.fasta")
    
    # Merge the found FASTA files into one big file
    merge_fasta_files(files, output_file)
    print(f"Merged FASTA file written to: {output_file}")

if __name__ == '__main__':
    main()