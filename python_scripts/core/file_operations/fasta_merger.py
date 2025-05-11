import os

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

def find_fasta_files(input_file):
    fasta_files = []
    with open(input_file, 'r') as f:
        fasta_files = [file.strip() for file in f.readlines()]

    return fasta_files

def main():
    input_file = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_output/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/merged_list.txt"
    output_file_name = "merged_mgc_candidartes"

    files = find_fasta_files(input_file)

    output_dir = os.path.dirname(input_file)
    output_file = os.path.join(output_dir, f"{output_file_name}.fasta")
    
    # Merge the found FASTA files into one big file
    merge_fasta_files(files, output_file)
    print(f"Merged FASTA file written to: {output_file}")

if __name__ == '__main__':
    main()