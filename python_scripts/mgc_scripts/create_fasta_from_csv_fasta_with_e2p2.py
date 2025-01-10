import csv
import os
from Bio import SeqIO

# Input files
csv_file = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/negative_training_set.csv"
fasta_file = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/pmn.fasta"
output_dir = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/negative_candidates_fasta_files"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Helper function to parse gene IDs from the CSV
def parse_gene_ids(gene_id_field):
    gene_ids = []
    for item in gene_id_field.split(","):
        gene_ids.append(item)
    return gene_ids

# Read the CSV file and parse the pathway-gene map
entries = []
negative_candidate_id = 1

with open(csv_file, "r") as csvfile:
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
        pathway = row["Pathway (Occurrence)"].split(" (")[0]  # Remove occurrence info
        gene_ids = parse_gene_ids(row["Window Genes"])
        for gene_id in gene_ids:
            entries.append({
                "negative_candidate_id": f"NEGATIVE_CANDIDATE_{negative_candidate_id}",
                "gene_id": gene_id,
                "pathway": pathway
            })
        negative_candidate_id += 1

# Load all FASTA records into a dictionary for fast lookup
fasta_records = {
    record.id.split(" ")[0]: record for record in SeqIO.parse(fasta_file, "fasta")
}

# Write individual FASTA files for each candidate
for entry in entries:
    gene_id = entry["gene_id"]
    negative_candidate_id = entry["negative_candidate_id"]
    pathway = entry["pathway"]
    
    if gene_id in fasta_records:
        record = fasta_records[gene_id]
        organism = record.description.split(" | ")[1] if " | " in record.description else "unknown"
        new_header = f">{negative_candidate_id} | {gene_id} | {pathway} | {organism}"
        fasta_filename = os.path.join(output_dir, f"{negative_candidate_id}.fasta")
        
        with open(fasta_filename, "a") as fasta_file:
            fasta_file.write(f"{new_header}\n{str(record.seq)}\n")

print(f"FASTA files created in folder: {output_dir}")
