import os
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq

def convert_to_protein_fasta(input_file, output_file):
    """
    Converts a tab-delimited file with Gene_ID and Nucleotide_Sequence into FASTA format with protein sequences.
    
    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output FASTA file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip comments and empty lines
            if line.startswith("#") or not line.strip():
                continue
            
            # Split the line into Gene_ID and Nucleotide_Sequence
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue  # Skip malformed lines
            
            gene_id, nucleotide_sequence = parts
            
            # Convert nucleotide sequence to protein sequence
            try:
                protein_sequence = str(Seq(nucleotide_sequence).translate(to_stop=True))
            except Exception as e:
                print(f"Error translating sequence for {gene_id}: {e}")
                continue
            
            # Write in FASTA format
            outfile.write(f">{gene_id}\n")
            outfile.write(f"{protein_sequence}\n")


def process_file(input_file, output_file):
    """
    Wrapper function to process a single file.
    
    Args:
        input_file (str): Path to the input file.
        output_file (str): Path to the output FASTA file.
    """
    print(f"Processing {input_file}...")
    convert_to_protein_fasta(input_file, output_file)


def process_all_directories(base_input_dir, base_output_dir, max_workers=4):
    """
    Processes all subdirectories under the base input directory, converting all `.txt` files
    to FASTA format with protein sequences and saving them in corresponding subdirectories
    under the base output directory. Uses concurrency for faster processing.
    
    Args:
        base_input_dir (str): Path to the base input directory.
        base_output_dir (str): Path to the base output directory.
        max_workers (int): Maximum number of threads to use for concurrent processing.
    """
    tasks = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for subdir, _, files in os.walk(base_input_dir):
            # Determine the relative path of the subdirectory
            relative_path = os.path.relpath(subdir, base_input_dir)
            # Create the corresponding output subdirectory
            output_subdir = os.path.join(base_output_dir, relative_path)
            os.makedirs(output_subdir, exist_ok=True)
            
            for file in files:
                if file.endswith(".txt"):  # Process only .txt files
                    input_file = os.path.join(subdir, file)
                    output_file = os.path.join(output_subdir, file.replace(".txt", ".fasta"))
                    # Submit the task to the thread pool
                    tasks.append(executor.submit(process_file, input_file, output_file))
        
        # Wait for all tasks to complete
        for task in tasks:
            task.result()
    print("Conversion complete for all directories.")


# Define the base input and output directories
base_input_dir = "/groups/itay_mayrose/alongonda/datasets/KEGG"
base_output_dir = "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta"

# Process all directories with concurrency
process_all_directories(base_input_dir, base_output_dir, max_workers=8)