import os
import glob

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

def main():
    search_dir = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files"
    # Recursively find all .fasta files in a directory and its subdirectories using os.walk
    fasta_files = []
    for root, dirs, files in os.walk(search_dir):
        for file in files:
            if file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))
                
                
                
    output_file_name = "merged_metabolic_pathways"

    output_file = os.path.join(search_dir, f"{output_file_name}.fasta")
    
    # Merge the found FASTA files into one big file
    merge_fasta_files(fasta_files, output_file)
    print(f"Merged FASTA file written to: {output_file}")

if __name__ == '__main__':
    main()