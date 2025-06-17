import os
import pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

# Input files
CSV_FILE = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_w10_g3/kegg_scanner_min_genes_based_metabolic/min_genes_3/potential_groups_w10_filtered.csv"
OUTPUT_DIR = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_w10_g3/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)
FINAL_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "final_fasta_files")
os.makedirs(FINAL_OUTPUT_DIR, exist_ok=True)

# Function to parse gene IDs from a comma-separated string
def parse_gene_ids(gene_id_field):
    """Parse gene IDs from a comma-separated string."""
    return gene_id_field.split(",")

# Function to read the CSV file and parse the pathway-gene map
def read_csv(file_path):
    """Read the CSV file and parse the pathway-gene map."""
    entries = []
    mgc_candidate_id = 1
    df = pd.read_csv(file_path)
    
    for _, row in df.iterrows():
        pathway = row["pathway"]
        source_file = row["source_file"]
        gene_ids = parse_gene_ids(row["metabolic_genes"])
        
        for gene_id in gene_ids:
            entries.append({
                "pathway": pathway,
                "metabolic_genes": gene_id,
                "source_file": source_file,
                "mgc_candidate_id": f"MGC_CANDIDATE_{mgc_candidate_id}",
            })
        mgc_candidate_id += 1
    return entries

# Function to load sequences from the source file into a dictionary
def load_sequences(source_file):
    """Load metabolic gene sequences into a dictionary from the source file."""
    df = pd.read_csv(source_file)
    return dict(zip(df["id"], df["sequence"]))

# Function to process each entry and return the FASTA content
def process_entry(entry):
    """Process an entry and return the FASTA content instead of writing directly."""
    gene_id = entry["metabolic_genes"]
    mgc_candidate_id = entry["mgc_candidate_id"]
    pathway = entry["pathway"]
    source_file = entry["source_file"]

    if not os.path.exists(source_file):
        return None, None  # Skip missing files

    sequences = load_sequences(source_file)
    if gene_id not in sequences:
        return None, None  # Skip missing genes

    sequence = sequences[gene_id]
    organism = os.path.basename(source_file).replace("_filtered.csv", "")
    new_header = f">{mgc_candidate_id} | {gene_id} | {pathway} | {organism}"
    print(mgc_candidate_id)
    
    return mgc_candidate_id, f"{new_header}\n{sequence}\n"

# Function to write FASTA files using multiprocessing
def write_fasta_files(entries, output_dir):
    """Write FASTA files using multiprocessing and batch writing."""
    fasta_batches = defaultdict(list)

    with ProcessPoolExecutor() as executor:
        results = executor.map(process_entry, entries)
        
        for mgc_candidate_id, fasta_content in results:
            if mgc_candidate_id and fasta_content:
                fasta_batches[mgc_candidate_id].append(fasta_content)

    # Batch write to minimize file I/O overhead
    for mgc_candidate_id, contents in fasta_batches.items():
        fasta_filename = os.path.join(output_dir, f"{mgc_candidate_id}.fasta")
        with open(fasta_filename, "a") as fasta_file:
            fasta_file.writelines(contents)
            
# Main function to read the CSV and write FASTA files
def main():
    entries = read_csv(CSV_FILE)
    write_fasta_files(entries, OUTPUT_DIR)
    print(f"FASTA files created in folder: {OUTPUT_DIR}")

# Entry point for the script
if __name__ == "__main__":
    main()
