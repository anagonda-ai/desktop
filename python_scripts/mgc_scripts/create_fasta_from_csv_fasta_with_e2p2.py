import csv
import os
from Bio import SeqIO

# Input files
csv_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates.csv"
fasta_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/pmn.fasta"
output_dir = "/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates_fasta_files_e2p2_filtered"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Helper function to parse gene IDs from the CSV
def parse_gene_ids(gene_id_field):
    gene_ids = []
    for item in gene_id_field.split(" "):
        if ":" in item:
            gene_ids.append(item.split(":")[0])
    return gene_ids

# Read the CSV file and parse the pathway-gene map
entries = []
mgc_candidate_id = 1

with open(csv_file, "r") as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
        pathway = row["Pathway"].split(" (")[0]  # Remove occurrence info
        gene_ids = parse_gene_ids(row["Gene IDs and Predicted EC numbers"])
        for gene_id in gene_ids:
            entries.append({
                "mgc_candidate_id": f"MGC_CANDIDATE_{mgc_candidate_id}",
                "gene_id": gene_id,
                "pathway": pathway
            })
        mgc_candidate_id += 1

# Load all FASTA records into a dictionary for fast lookup
fasta_records = {
    record.id.split(" ")[0]: record for record in SeqIO.parse(fasta_file, "fasta")
}

# Write individual FASTA files for each candidate
for entry in entries:
    gene_id = entry["gene_id"]
    mgc_candidate_id = entry["mgc_candidate_id"]
    pathway = entry["pathway"]
    
    if gene_id in fasta_records:
        record = fasta_records[gene_id]
        organism = record.description.split(" | ")[1] if " | " in record.description else "unknown"
        new_header = f">{mgc_candidate_id} | {gene_id} | {pathway} | {organism}"
        fasta_filename = os.path.join(output_dir, f"{mgc_candidate_id}.fasta")
        
        with open(fasta_filename, "a") as fasta_file:
            fasta_file.write(f"{new_header}\n{str(record.seq)}\n")

print(f"FASTA files created in folder: {output_dir}")
