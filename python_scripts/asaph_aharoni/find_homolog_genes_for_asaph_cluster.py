from collections import defaultdict
import os
import re
import subprocess
import csv
from matplotlib import pyplot as plt
import numpy as np
import threading
import concurrent.futures

# Directories containing CSV files
csv_dirs = [
    # "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test_no_chloroplast_with_sequences",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences_with_strand",
    # "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences"
]

# Directory to store generated FASTA files and BLAST DBs
blast_db_dir = "/groups/itay_mayrose/alongonda/datasets/generated_blast_dbs_without_haaap_stranded"
os.makedirs(blast_db_dir, exist_ok=True)

# Query FASTA files
example_mgc = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni"
query_fastas = [os.path.join(example_mgc, f) for f in os.listdir(example_mgc) if f.endswith(".fasta") and "HAAAP" not in f]

# Output directory for BLAST results
blast_results_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated_without_haaap_stranded"
os.makedirs(blast_results_dir, exist_ok=True)

# Function to convert CSV to FASTA, skipping empty files
def csv_to_fasta(csv_path):
    fasta_filename = os.path.splitext(os.path.basename(csv_path))[0] + ".fasta"
    fasta_path = os.path.join(blast_db_dir, fasta_filename)
    valid_entries = 0

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
    cmd = f"blastp -query {query_fasta} -db {blast_db} -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 50 -gapopen 11 -gapextend 1 -out {output_file}"
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
            
def extract_organism_name(filename):
    match = re.search(r"fasta_(.+?)\.csv", filename)
    return match.group(1) if match else "Unknown"

used_row_indexes = defaultdict(lambda: defaultdict(set))
lock = threading.Lock()

# Function to parse BLAST results and extract the best hit for each chromosome
def parse_blast_results(blast_file):
    best_hits_by_organism = defaultdict(dict)  # Organism -> Chromosome -> Best Hit

    with open(blast_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue  # Skip malformed lines
            
            query_fasta = cols[0]  # The original FASTA file name
            subject_gene = cols[1]  # Full subject ID
            subject_info = subject_gene.split("|")
            chromosome = subject_info[1] if len(subject_info) > 1 else "Unknown" # This should be the chromosome
            start = int(subject_info[2]) if len(subject_info) > 1 else "Unknown" # Correct subject start position
            end = int(subject_info[3]) if len(subject_info) > 1 else "Unknown"  # Correct subject end position
            strand_value = subject_info[4] if len(subject_info) > 1 else "Unknown"  # Strand value
            row_index = int(subject_info[5]) if len(subject_info) > 1 else "Unknown"  # Correct subject end position
            identity = float(cols[2])  # Percentage identity
            evalue = float(cols[10])  # E-value
            bit_score = float(cols[11])  # Bit score

            # Extract organism name from BLAST filename
            organism_name = os.path.basename(blast_file).replace("_results.txt", "")
            dir_name = extract_organism_name(organism_name)
            hit = [query_fasta, subject_gene, chromosome, start, end, strand_value, identity, evalue, bit_score, row_index]

            with lock:
                # skip if row_index already used in this dir_name/chromosome
                if row_index in used_row_indexes[dir_name][chromosome]:
                    continue

                # only accept hit if it's the best for this organism_name/chromosome
                current_hit = best_hits_by_organism[organism_name].get(chromosome)
                if current_hit is None or bit_score > current_hit[-2]:
                    best_hits_by_organism[organism_name][chromosome] = hit
                    used_row_indexes[dir_name][chromosome].add(row_index)

    # Flatten hits for saving
    for organism in best_hits_by_organism:
        best_hits_by_organism[organism] = list(best_hits_by_organism[organism].values())

    return best_hits_by_organism

# Save best hits grouped by chromosome
def save_best_hits(best_hits_by_organism, output_dir):
    best_hits_dir = os.path.join(output_dir, "best_hits_by_organism")
    os.makedirs(best_hits_dir, exist_ok=True)

    for organism, hits in best_hits_by_organism.items():
        organism_name = extract_organism_name(organism)
        organism_dir = os.path.join(best_hits_dir, organism_name)
        os.makedirs(organism_dir, exist_ok=True)
        summary_file = os.path.join(organism_dir, organism)
        with open(summary_file, "w", encoding="utf-8") as f:
            f.write("query_fasta,subject_gene,chromosome,start,end,strand_value,identity,evalue,bit_score,row_index\n")
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

        # Merge results per organism
        for organism, hits in parsed_hits.items():
            for hit in hits:
                chrom = hit[2]  # Chromosome
                if chrom not in best_hits_by_organism[organism] or hit[-1] > best_hits_by_organism[organism][chrom][-1]:  # Compare bit scores
                    best_hits_by_organism[organism][chrom] = hit

# Convert dictionaries to lists for final output
for organism in best_hits_by_organism:
    best_hits_by_organism[organism] = list(best_hits_by_organism[organism].values())

# Save the results
save_best_hits(best_hits_by_organism, blast_results_dir)

print("🚀 Processing completed successfully with intelligent analysis!")
