from collections import defaultdict
import os
import re
import subprocess
import csv
from matplotlib import pyplot as plt
import numpy as np
import threading
import argparse
import concurrent.futures

#############################################
#               PATHS
#############################################

# Directories containing CSV files
full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
csv_dirs = [
    os.path.join(full_genome_dir, "final_dataset/filtered_no_mito")
]

# Directory to store generated FASTA files and BLAST DBs
blast_db_dir = "/groups/itay_mayrose/alongonda/datasets/generated_blast_dbs"
os.makedirs(blast_db_dir, exist_ok=True)

parser = argparse.ArgumentParser(description="Find homolog genes for Asaph cluster")
parser.add_argument("--example_mgc", type=str, required=True, help="Path to example MGC directory")
args = parser.parse_args()

example_mgc = args.example_mgc

query_fastas = [os.path.join(example_mgc, f) for f in os.listdir(example_mgc) if f.endswith(".fasta")]

# Output directory for BLAST results
blast_results_dir = os.path.join(example_mgc,"best_hits")
os.makedirs(blast_results_dir, exist_ok=True)

#############################################
#         CSV → FASTA Conversion
#############################################

# Function to convert CSV to FASTA, skipping empty files
def csv_to_fasta(csv_path):
    fasta_filename = os.path.splitext(os.path.basename(csv_path))[0] + ".fasta"
    fasta_path = os.path.join(blast_db_dir, fasta_filename)
    valid_entries = 0
    if os.path.exists(fasta_path):
        print(f"⚠️ Skipped existing FASTA: {fasta_path}")
    else:
        with open(csv_path, 'r') as csv_in, open(fasta_path, 'w') as fasta_out:
            reader = csv.DictReader(csv_in, delimiter=',')
            for row_index, row in enumerate(reader, start=1):
                if all(k in row for k in ['id', 'sequence', 'chromosome', 'start', 'end', 'strand']) and row['sequence'].strip():
                    fasta_out.write(f">{row['id']}|{row['chromosome']}|{row['start']}|{row['end']}|{row['strand']}|{row_index}\n{row['sequence']}\n")
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

#############################################
#         Create BLAST Databases
#############################################

# Function to create BLAST database
def create_blast_db(fasta_file):
    db_name = fasta_file + "_blastdb"
    cmd = f"makeblastdb -in {fasta_file} -dbtype prot -out {db_name}"
    try:
        if os.path.exists(db_name + ".pin"):
            print(f"⚠️ Skipped existing BLAST DB: {db_name}")
        else:
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

#############################################
#       SLURM JOB SUBMISSION FOR BLAST
#############################################

MAX_RUNNING = 100

def get_running_blast_jobs():
    cmd = ["squeue", "-u", os.environ["USER"], "-h", "-n", "blast_organisms_on_single_kegg"]
    out = subprocess.check_output(cmd).decode().strip()
    if not out:
        return 0
    return len(out.splitlines())

def submit_blast_job(query_fasta, csv_file, blast_db, job_id):
    output_path = os.path.join(
        blast_results_dir, f"{os.path.basename(query_fasta)}_{csv_file}_results.txt"
    )

    job_script = f"""#!/bin/bash
#SBATCH --job-name=blast_organisms_on_single_kegg
#SBATCH --output={blast_results_dir}/job_{job_id}.OU
#SBATCH --error={blast_results_dir}/job_{job_id}.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym-pool
#SBATCH --time=6-24:00:00

module load ncbi-blast/2.13.0+

blastp -query {query_fasta} -db {blast_db} \\
  -evalue 1e-5 \\
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
  -max_target_seqs 50 -gapopen 11 -gapextend 1 \\
  -out {output_path}
"""

    script_path = f"{blast_results_dir}/job_{job_id}.sh"
    with open(script_path, "w") as f:
        f.write(job_script)

    # Submit job
    out = subprocess.check_output(["sbatch", script_path]).decode().strip()
    print(f"🧬 Submitted job {job_id}: {out}")

    # extract job number from "Submitted batch job 123456"
    jobnum = out.split()[-1]
    return jobnum


#############################################
#         SUBMISSION LOOP (max 100 jobs)
#############################################

job_id = 0
for query_fasta in query_fastas:
    for csv_file, blast_db in blast_dbs.items():

        while get_running_blast_jobs() >= MAX_RUNNING:
            print("⏳ 100 jobs running — waiting...")
            time.sleep(20)

        submit_blast_job(query_fasta, csv_file, blast_db, job_id)
        job_id += 1
        
#############################################
# WAIT FOR ALL BLAST JOBS TO FINISH
#############################################

def any_jobs_still_running(job_ids):
    if not job_ids:
        return False
    cmd = ["squeue", "-h", "-j", ",".join(job_ids)]
    out = subprocess.check_output(cmd).decode().strip()
    return bool(out)

print("⏳ Waiting for all BLAST jobs to finish...")

while any_jobs_still_running(submitted_job_ids):
    print("   …still running…")
    time.sleep(30)

print("✅ All BLAST jobs finished!")
        
#############################################
#   RESULTS PARSING & BEST HIT EXTRACTION
#############################################
            
def extract_organism_name(filename):
    match = re.search(r"fasta_(.+?)\.csv", filename)
    return match.group(1) if match else "Unknown"

used_row_indexes = defaultdict(lambda: defaultdict(set))
lock = threading.Lock()

# Function to parse BLAST results and extract the best hit for each query gene
def parse_blast_results(blast_file):
    best_hits_by_organism = defaultdict(dict)  # Organism -> Query Gene -> Best Hit

    with open(blast_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue  # Skip malformed lines
            
            query_gene = cols[0]  # The query gene ID
            subject_gene = cols[1]  # Full subject ID
            subject_info = subject_gene.split("|")
            chromosome = subject_info[1] if len(subject_info) > 1 else "Unknown"
            start = int(subject_info[2]) if len(subject_info) > 2 else 0
            end = int(subject_info[3]) if len(subject_info) > 3 else 0
            strand_value = subject_info[4] if len(subject_info) > 4 else "Unknown"
            row_index = int(subject_info[5]) if len(subject_info) > 5 else 0
            identity = float(cols[2])  # Percentage identity
            evalue = float(cols[10])  # E-value
            bit_score = float(cols[11])  # Bit score

            # Extract organism name from BLAST filename
            organism_name = os.path.basename(blast_file).replace("_results.txt", "")
            dir_name = extract_organism_name(organism_name)
            hit = [query_gene, subject_gene, chromosome, start, end, strand_value, identity, evalue, bit_score, row_index]

            with lock:
                # Skip if row_index already used in this dir_name
                if row_index in used_row_indexes[dir_name][query_gene]:
                    continue

                # Only accept hit if it's the best for this organism_name/query_gene combination
                current_hit = best_hits_by_organism[organism_name].get(query_gene)
                if current_hit is None or bit_score > current_hit[8]:  # Compare bit scores
                    # If replacing an existing hit, remove its row_index from used set
                    if current_hit is not None:
                        used_row_indexes[dir_name][query_gene].discard(current_hit[9])
                    
                    best_hits_by_organism[organism_name][query_gene] = hit
                    used_row_indexes[dir_name][query_gene].add(row_index)

    # Flatten hits for saving
    for organism in best_hits_by_organism:
        best_hits_by_organism[organism] = list(best_hits_by_organism[organism].values())

    return best_hits_by_organism

# Save best hits grouped by organism
def save_best_hits(best_hits_by_organism, output_dir):
    best_hits_dir = os.path.join(output_dir, "best_hits_by_organism")
    os.makedirs(best_hits_dir, exist_ok=True)

    for organism, hits in best_hits_by_organism.items():
        organism_name = extract_organism_name(organism)
        organism_dir = os.path.join(best_hits_dir, organism_name)
        os.makedirs(organism_dir, exist_ok=True)
        summary_file = os.path.join(organism_dir, organism)
        if os.path.exists(summary_file):
            print(f"⚠️ Skipped existing summary file: {summary_file}")
            continue
        with open(summary_file, "w", encoding="utf-8") as f:
            f.write("query_gene,subject_gene,chromosome,start,end,strand_value,identity,evalue,bit_score,row_index\n")
            for hit in hits:
                f.write(",".join(map(str, hit)) + "\n")

        print(f"✅ Best hits saved for organism: {organism}")

# Main execution loop to parse and save best hits
best_hits_by_organism = defaultdict(dict)

with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = {executor.submit(parse_blast_results, os.path.join(blast_results_dir, blast_file)): blast_file
               for blast_file in os.listdir(blast_results_dir) if blast_file.endswith(".txt")}

    for future in concurrent.futures.as_completed(futures):
        blast_file = futures[future]
        parsed_hits = future.result()

        # Merge results per organism and query gene
        for organism, hits in parsed_hits.items():
            for hit in hits:
                query_gene = hit[0]  # Query gene
                if query_gene not in best_hits_by_organism[organism] or hit[8] > best_hits_by_organism[organism][query_gene][8]:  # Compare bit scores
                    best_hits_by_organism[organism][query_gene] = hit

# Convert dictionaries to lists for final output
for organism in best_hits_by_organism:
    best_hits_by_organism[organism] = list(best_hits_by_organism[organism].values())

# Save the results
save_best_hits(best_hits_by_organism, blast_results_dir)

print("🚀 Processing completed successfully with intelligent analysis!")