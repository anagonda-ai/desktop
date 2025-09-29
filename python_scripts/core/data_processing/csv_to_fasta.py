import os
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

def csv_to_fasta(csv_file, fasta_file):
    df = pd.read_csv(csv_file)
    with open(fasta_file, 'w') as f:
        for _, row in df.iterrows():
            gene_id = str(row['gene_id'])
            sequence = str(row['sequence'])
            f.write(f">{gene_id}\n{sequence}\n")


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
    
    input_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/random_mgc_candidates_csv_files"
    output_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/random_mgc_candidates_fasta_files"
    process_directory(input_dir, output_dir)