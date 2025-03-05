import os
import subprocess
import csv
import pandas as pd
import concurrent.futures

# Directories containing CSV files
csv_dirs = [
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test_no_chloroplast_with_sequences",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences"
]

# Directory to store generated FASTA files and BLAST DBs
blast_db_dir = "/groups/itay_mayrose/alongonda/datasets/generated_blast_dbs"
os.makedirs(blast_db_dir, exist_ok=True)

# Query FASTA files
query_fastas = [
    "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/HAAAP_Transporters.fasta",
    "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/adcs.fasta",
    "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/cs.fasta",
    "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/pABA_transporter.fasta"
]

# Output directory for BLAST results
blast_results_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results"
os.makedirs(blast_results_dir, exist_ok=True)

# Function to convert CSV to FASTA, skipping empty files
def csv_to_fasta(csv_path):
    fasta_filename = os.path.splitext(os.path.basename(csv_path))[0] + ".fasta"
    fasta_path = os.path.join(blast_db_dir, fasta_filename)
    valid_entries = 0

    with open(csv_path, 'r') as csv_in, open(fasta_path, 'w') as fasta_out:
        reader = csv.DictReader(csv_in, delimiter=',')
        for row in reader:
            if all(k in row for k in ['id', 'sequence', 'chromosome', 'start', 'end']) and row['sequence'].strip():
                fasta_out.write(f">{row['id']}|{row['chromosome']}|{row['start']}|{row['end']}\n{row['sequence']}\n")
                valid_entries += 1

    if valid_entries == 0:
        os.remove(fasta_path)  # Remove empty FASTA
        print(f"⚠️ Skipped empty CSV: {csv_path}")
        return None

    print(f"✅ Converted {csv_path} → {fasta_path}")
    return fasta_path


# Convert CSV files in parallel
csv_to_fasta_map = {}
with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = {executor.submit(csv_to_fasta, os.path.join(csv_dir, csv_file)): csv_file for csv_dir in csv_dirs for csv_file in os.listdir(csv_dir) if csv_file.endswith(".csv")}
    for future in concurrent.futures.as_completed(futures):
        fasta_file = future.result()
        if fasta_file:
            csv_to_fasta_map[futures[future]] = fasta_file

# Function to create BLAST database
def create_blast_db(fasta_file):
    db_name = fasta_file + "_blastdb"
    cmd = f"makeblastdb -in {fasta_file} -dbtype prot -out {db_name}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        return db_name
    except subprocess.CalledProcessError:
        return None

# Create BLAST databases in parallel
blast_dbs = {}
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = {executor.submit(create_blast_db, fasta_file): csv_file for csv_file, fasta_file in csv_to_fasta_map.items()}
    for future in concurrent.futures.as_completed(futures):
        db_name = future.result()
        if db_name:
            blast_dbs[futures[future]] = db_name

# Function to run BLASTP
def run_blastp(query_fasta, csv_file, blast_db):
    output_file = os.path.join(blast_results_dir, f"{os.path.basename(query_fasta)}_{csv_file}_results.txt")
    cmd = f"blastp -query {query_fasta} -db {blast_db} -evalue 0.01 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 50 -gapopen 11 -gapextend 1 -out {output_file}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        return (query_fasta, csv_file, output_file)
    except subprocess.CalledProcessError:
        return None

# Run BLASTP in parallel
blast_results_map = {}
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = {executor.submit(run_blastp, query_fasta, csv_file, blast_db) for query_fasta in query_fastas for csv_file, blast_db in blast_dbs.items()}
    for future in concurrent.futures.as_completed(futures):
        result = future.result()
        if result:
            query_fasta, csv_file, output_file = result
            blast_results_map[(query_fasta, csv_file)] = output_file

# Function to parse BLAST results and extract best hit
def parse_blast_results(blast_file):
    best_hit = None
    best_score = float("-inf")

    with open(blast_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue  # Skip malformed lines
            
            query = cols[0]
            subject_info = cols[1].split("|")  # Extract first part before `;`
            
            # ✅ Extract chromosome, start, and end correctly
            subject_id = subject_info[0]
            chromosome = subject_info[1] if len(subject_info) > 1 else "Unknown"
            start = subject_info[2] if len(subject_info) > 2 else "Unknown"
            end = subject_info[3] if len(subject_info) > 3 else "Unknown"

            identity = float(cols[2])
            evalue = float(cols[10])  # Handles scientific notation properly
            bit_score = float(cols[11])  # Extracts bit score correctly
            
            if bit_score > best_score:
                best_score = bit_score
                best_hit = (subject_id, chromosome, start, end, identity, evalue, bit_score)  # ✅ Now includes Chromosome, Start, End

    return best_hit



# Extract best hits in parallel
best_hits_by_csv = {csv_file: [] for csv_file in blast_dbs.keys()}
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = {executor.submit(parse_blast_results, blast_file): (query_fasta, csv_file) for (query_fasta, csv_file), blast_file in blast_results_map.items()}
    for future in concurrent.futures.as_completed(futures):
        best_hit = future.result()
        if best_hit:
            query_fasta, csv_file = futures[future]
            best_hits_by_csv[csv_file].append([os.path.basename(query_fasta)] + list(best_hit))

# Save best hits summary in parallel
def save_best_hits(csv_file, best_hits):
    best_hits_dir = os.path.join(blast_results_dir, f"best_hits")
    os.makedirs(best_hits_dir, exist_ok=True)
    summary_file = os.path.join(best_hits_dir, f"best_hits_{csv_file}.csv")
    print(best_hits)
    with open(summary_file, "w") as f:
        f.write("query_gene,subject_gene,chromosome,start,end,identity,value,bit_score\n")
        for hit in best_hits:
            f.write(",".join(map(str, hit)) + "\n")
    print(f"✅ Best hits saved: {summary_file}")

with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = {executor.submit(save_best_hits, csv_file, best_hits) for csv_file, best_hits in best_hits_by_csv.items()}
    concurrent.futures.wait(futures)

print("🚀 All tasks completed successfully!")
