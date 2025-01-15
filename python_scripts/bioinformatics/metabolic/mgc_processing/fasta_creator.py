import csv
import os
from Bio import SeqIO

# Input files
CSV_FILE = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/results/candidates.csv"
FASTA_FILE = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/pmn.fasta"
OUTPUT_DIR = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates_fasta_files_e2p2_filtered"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_gene_ids(gene_id_field):
    """Parse gene IDs from a semicolon-separated string."""
    return gene_id_field.split("; ")

def read_csv(file_path):
    """Read the CSV file and parse the pathway-gene map."""
    entries = []
    mgc_candidate_id = 1
    with open(file_path, "r") as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            pathway = row["Pathway (Occurrence)"].split(" (")[0]  # Remove occurrence info
            gene_ids = parse_gene_ids(row["Gene IDs"])
            for gene_id in gene_ids:
                entries.append({
                    "pathway": pathway,
                    "window_cluster_genes": gene_id,
                    "mgc_candidate_id": f"MGC_CANDIDATE_{mgc_candidate_id}",
                })
            mgc_candidate_id += 1
    return entries

def load_fasta_records(fasta_file):
    """Load all FASTA records into a dictionary for fast lookup."""
    return {record.id.lower(): record for record in SeqIO.parse(fasta_file, "fasta")}

def write_fasta_files(entries, fasta_records, output_dir):
    """Write individual FASTA files for each candidate."""
    for entry in entries:
        gene_id = entry["window_cluster_genes"]
        mgc_candidate_id = entry["mgc_candidate_id"]
        pathway = entry["pathway"]
        
        if gene_id.lower() in fasta_records:
            record = fasta_records[gene_id.lower()]
            organism = record.description.split(" | ")[1] if " | " in record.description else "unknown"
            new_header = f">{mgc_candidate_id} | {gene_id} | {pathway} | {organism}"
            fasta_filename = os.path.join(output_dir, f"{mgc_candidate_id}.fasta")
            
            with open(fasta_filename, "a") as fasta_file:
                fasta_file.write(f"{new_header}\n{str(record.seq)}\n")

def main():
    entries = read_csv(CSV_FILE)
    fasta_records = load_fasta_records(FASTA_FILE)
    write_fasta_files(entries, fasta_records, OUTPUT_DIR)
    print(f"FASTA files created in folder: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()