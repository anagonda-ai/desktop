import csv
import os
from Bio import SeqIO

# Input files
csv_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/unique_clusters_start_end.csv"
fasta_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/pmn.fasta"
output_dir = "/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates_fasta_files_without_e2p2_filtered"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Read the FASTA file and store records in a dictionary
fasta_records = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

# Read the CSV file and parse the pathway-gene map
entries = []
mgc_candidate_id = 1

with open(csv_file, "r") as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
        pathway = row["pathway"].split(" (")[0]  # Remove occurrence info
        gene_ids = row["window_cluster_genes"].split(",")
        for gene_id in gene_ids:
            entries.append({
                "pathway": pathway,
                "window_cluster_genes": gene_id,
                "mgc_candidate_id": f"MGC_CANDIDATE_{mgc_candidate_id}",
            })
        mgc_candidate_id += 1

# Process the entries and write to FASTA files
for entry in entries:
    gene_id = entry["window_cluster_genes"]
    mgc_candidate_id = entry["mgc_candidate_id"]
    pathway = entry["pathway"]
    
    if gene_id in fasta_records:
        record = fasta_records[gene_id]
        organism = record.description.split(" | ")[1] if " | " in record.description else "unknown"
        new_header = f">{mgc_candidate_id} | {gene_id} | {pathway} | {organism}"
        fasta_filename = os.path.join(output_dir, f"{mgc_candidate_id}.fasta")
        
        with open(fasta_filename, "a") as fasta_file:
            fasta_file.write(f"{new_header}\n{str(record.seq)}\n")