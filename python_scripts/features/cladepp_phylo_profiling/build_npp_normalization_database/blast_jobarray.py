#!/usr/bin/env python3
"""
Final end-to-end launcher:
- CSV -> FASTA (optional)
- makeblastdb for generated FASTAs (optional)
- write query_fastas.txt and blast_dbs.txt
- write & submit one SLURM array job (one task per query fasta)
- submit dependent merge job that produces combined_best_hits.csv

Usage:
    python run_blast_array_and_merge_final.py --example_mgc /groups/.../selected_genes/sly_selected_1000
"""

from collections import defaultdict
import argparse
import concurrent.futures
import csv
import glob
import os
import shlex
import subprocess
import sys
import time

# Reference doc uploaded by user (developer requested)
REFERENCE_DOC = "/mnt/data/Highlits for the Toolkit.txt"

# ---------------------------
# Arguments / defaults
# ---------------------------
parser = argparse.ArgumentParser(description="Submit BLAST array + merge best hits (TAU HPC - itay_mayrose-users_v2)")
parser.add_argument("--example_mgc", default="/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/aglu_selected_1000", help="Path to example MGC directory containing query .fasta files")
parser.add_argument("--cpus", type=int, default=8, help="CPUs per array task (passed to -num_threads)")
parser.add_argument("--mem", default="32G", help="Memory per array task")
parser.add_argument("--array_time", default="2-00:00:00", help="SLURM time for array tasks")
parser.add_argument("--array_concurrency", type=int, default=100, help="Max concurrent array tasks (SLURM %limit)")
parser.add_argument("--blast_module", default="ncbi-blast/2.13.0+", help="Module name for BLAST")
parser.add_argument("--account", default="itay_mayrose-users_v2", help="SLURM account (V2)")
parser.add_argument("--partition", default="itaym-pool", help="SLURM partition/pool")
args = parser.parse_args()

example_mgc = os.path.abspath(args.example_mgc)
if not os.path.isdir(example_mgc):
    sys.exit(f"ERROR: example_mgc not found: {example_mgc}")

# ---------------------------
# Paths
# ---------------------------
full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
csv_dirs = [os.path.join(full_genome_dir, "final_dataset/filtered_no_mito")]

blast_db_dir = "/groups/itay_mayrose/alongonda/datasets/generated_blast_dbs"
os.makedirs(blast_db_dir, exist_ok=True)

blast_results_dir = os.path.join(example_mgc, "best_hits")
os.makedirs(blast_results_dir, exist_ok=True)

# ---------------------------
# Gather query FASTAs
# ---------------------------
query_fastas = sorted([os.path.join(example_mgc, f) for f in os.listdir(example_mgc) if f.endswith(".fasta")])
if not query_fastas:
    sys.exit("ERROR: No .fasta query files found in example_mgc")

print(f"Found {len(query_fastas)} query FASTA files in {example_mgc}")

# ---------------------------
# CSV -> FASTA conversion (optional)
# ---------------------------
def csv_to_fasta(csv_path):
    fasta_filename = os.path.splitext(os.path.basename(csv_path))[0] + ".fasta"
    fasta_path = os.path.join(blast_db_dir, fasta_filename)
    valid_entries = 0
    if os.path.exists(fasta_path):
        print(f"⚠️  Skipped existing FASTA: {fasta_path}")
        return fasta_path
    try:
        with open(csv_path, "r", newline="") as csv_in, open(fasta_path, "w") as fasta_out:
            reader = csv.DictReader(csv_in, delimiter=",")
            for row_index, row in enumerate(reader, start=1):
                if all(k in row for k in ["id", "sequence", "chromosome", "start", "end", "strand"]) and row["sequence"].strip():
                    fasta_out.write(f">{row['id']}|{row['chromosome']}|{row['start']}|{row['end']}|{row['strand']}|{row_index}\n")
                    fasta_out.write(row["sequence"].strip() + "\n")
                    valid_entries += 1
        if valid_entries == 0:
            os.remove(fasta_path)
            print(f"⚠️  Skipped empty CSV: {csv_path}")
            return None
        print(f"✅  Converted {csv_path} -> {fasta_path}")
        return fasta_path
    except Exception as e:
        print(f"ERROR converting {csv_path}: {e}")
        return None

csv_to_fasta_map = {}
csv_files = []
for d in csv_dirs:
    if os.path.isdir(d):
        csv_files.extend(sorted([os.path.join(d, f) for f in os.listdir(d) if f.endswith(".csv")]))

if csv_files:
    print(f"Converting {len(csv_files)} CSV files -> FASTA (parallel)...")
    with concurrent.futures.ProcessPoolExecutor() as ex:
        futures = {ex.submit(csv_to_fasta, csv): csv for csv in csv_files}
        for fut in concurrent.futures.as_completed(futures):
            fasta = fut.result()
            if fasta:
                csv_to_fasta_map[os.path.basename(futures[fut])] = fasta
else:
    print("No CSV files found (skipping CSV->FASTA conversion).")

# ---------------------------
# Create BLAST DBs for FASTAs (optional)
# ---------------------------
def create_blast_db(fasta_file):
    db_name = fasta_file + "_blastdb"
    # avoid shell, use list args
    if os.path.exists(db_name + ".pin") or os.path.exists(db_name + ".phr"):
        print(f"⚠️  Skipped existing BLAST DB: {db_name}")
        return db_name
    cmd = ["makeblastdb", "-in", fasta_file, "-dbtype", "prot", "-out", db_name]
    try:
        subprocess.run(cmd, check=True)
        print(f"✅  Created BLAST DB: {db_name}")
        return db_name
    except subprocess.CalledProcessError as e:
        print(f"❌ makeblastdb failed for {fasta_file}: {e}")
        return None

blast_dbs = {}
if csv_to_fasta_map:
    print("Creating BLAST DBs for converted FASTAs (parallel)...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        futures = {ex.submit(create_blast_db, fasta): csvname for csvname, fasta in csv_to_fasta_map.items()}
        for fut in concurrent.futures.as_completed(futures):
            db = fut.result()
            if db:
                blast_dbs[futures[fut]] = db

print(f"Using {len(blast_dbs)} BLAST DB(s)")

# ---------------------------
# Write query_fastas.txt and blast_dbs.txt
# ---------------------------
query_list_file = os.path.join(example_mgc, "query_fastas.txt")
with open(query_list_file, "w") as fh:
    for q in query_fastas:
        fh.write(q + "\n")
print(f"Wrote queries list → {query_list_file}")

blastdb_list_file = os.path.join(example_mgc, "blast_dbs.txt")
with open(blastdb_list_file, "w") as fh:
    for key, db in sorted(blast_dbs.items()):
        fh.write(db + "\n")
print(f"Wrote BLAST DB list → {blastdb_list_file}")

# ---------------------------
# Write array sbatch
# ---------------------------
num_queries = len(query_fastas)
array_range = f"0-{num_queries-1}%{args.array_concurrency}"
array_sbatch = os.path.join(example_mgc, "blast_array.sbatch")
array_stdout = os.path.join(blast_results_dir, "blast_%A_%a.out")
array_stderr = os.path.join(blast_results_dir, "blast_%A_%a.err")

array_script = f"""#!/bin/bash
#SBATCH -J blast_array
#SBATCH -o {shlex.quote(array_stdout)}
#SBATCH -e {shlex.quote(array_stderr)}
#SBATCH -p {args.partition}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.cpus}
#SBATCH --mem={args.mem}
#SBATCH --array={array_range}
#SBATCH --time={args.array_time}

module purge
module load {args.blast_module}

QUERY=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" {shlex.quote(query_list_file)})
QBASE=$(basename "$QUERY" .fasta)

OUTDIR={shlex.quote(blast_results_dir)}/$QBASE
mkdir -p "$OUTDIR"

# iterate DBs
while read -r DB; do
  DBNAME=$(basename "$DB")
  OUTFILE="$OUTDIR/${{QBASE}}_${{DBNAME}}_results.txt"
  blastp -query "$QUERY" -db "$DB" \\
    -evalue 1e-5 \\
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
    -max_target_seqs 50 \\
    -gapopen 11 -gapextend 1 \\
    -num_threads {args.cpus} \\
    -out "$OUTFILE"
done < {shlex.quote(blastdb_list_file)}
"""

with open(array_sbatch, "w") as fh:
    fh.write(array_script)
print(f"Wrote array sbatch → {array_sbatch}")

# ---------------------------
# Submit array job
# ---------------------------
try:
    res = subprocess.check_output(["sbatch", array_sbatch]).decode().strip()
    print(res)
    array_jobid = res.split()[-1]
    print(f"Array job submitted: {array_jobid}")
except subprocess.CalledProcessError as e:
    print("ERROR submitting array job:", e.output.decode() if e.output else e)
    sys.exit(1)

# ---------------------------
# Write merge sbatch (dependent)
# ---------------------------
merge_sbatch = os.path.join(example_mgc, "merge_best_hits.sbatch")
merge_stdout = os.path.join(blast_results_dir, "merge_%j.out")
merge_stderr = os.path.join(blast_results_dir, "merge_%j.err")

# ---- Part 1: SLURM header (f-string allowed) ----
merge_script = f"""#!/bin/bash
#SBATCH -J merge_best_hits
#SBATCH -o {shlex.quote(merge_stdout)}
#SBATCH -e {shlex.quote(merge_stderr)}
#SBATCH -p {args.partition}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00

module purge
module load python/3.11

python - <<'PYCODE'
"""

# ---- Part 2: Python heredoc (NO f-string here!) ----
merge_script += f"""
import os, glob, csv
from collections import defaultdict

BLAST_ROOT = {repr(blast_results_dir)}
OUT_COMBINED = os.path.join(BLAST_ROOT, "best_hits_combined.csv")
"""

merge_script += """
# columns in blastp outfmt 6
best_per_query = {}  # key: (query_dir, qid) -> record dict

for qdir in sorted(glob.glob(os.path.join(BLAST_ROOT, "*"))):
    if not os.path.isdir(qdir):
        continue
    qbase = os.path.basename(qdir)
    for fpath in sorted(glob.glob(os.path.join(qdir, "*_results.txt"))):
        try:
            with open(fpath, "r") as fh:
                for ln in fh:
                    cols = ln.rstrip().split("\\t")
                    if len(cols) < 12:
                        continue
                    qid = cols[0]
                    sid = cols[1]
                    pident = cols[2]
                    length = cols[3]
                    mismatch = cols[4]
                    gapopen = cols[5]
                    qstart = cols[6]
                    qend = cols[7]
                    sstart = cols[8]
                    send = cols[9]
                    evalue = float(cols[10])
                    bitscore = float(cols[11])
                    key = (qbase, qid, sid)
                    rec = dict(query_dir=qbase, query_gene=qid, subject_gene=sid,
                               pident=pident, length=length, mismatch=mismatch, gapopen=gapopen,
                               qstart=qstart, qend=qend, sstart=sstart, send=send,
                               evalue=evalue, bitscore=bitscore, src=fpath)
                    if key not in best_per_query or bitscore > best_per_query[key]['bitscore'] or \
                       (bitscore == best_per_query[key]['bitscore'] and evalue < best_per_query[key]['evalue']):
                        best_per_query[key] = rec
        except Exception as e:
            print("WARN reading", fpath, e)

# write combined CSV
os.makedirs(BLAST_ROOT, exist_ok=True)
with open(OUT_COMBINED, "w", newline='', encoding="utf-8") as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow([
        "query_dir","query_gene","subject_gene","pident","length",
        "mismatch","gapopen","qstart","qend","sstart","send","evalue",
        "bitscore","source_file"
    ])
    for key in sorted(best_per_query):
        rec = best_per_query[key]
        writer.writerow([
            rec['query_dir'], rec['query_gene'], rec['subject_gene'],
            rec['pident'], rec['length'], rec['mismatch'], rec['gapopen'],
            rec['qstart'], rec['qend'], rec['sstart'], rec['send'],
            rec['evalue'], rec['bitscore'], rec['src']
        ])

print("✅ Wrote combined best-hits:", OUT_COMBINED)
PYCODE
"""

with open(merge_sbatch, "w") as fh:
    fh.write(merge_script)
print(f"Wrote merge sbatch → {merge_sbatch}")

# ---------------------------
# Submit merge job with dependency afterok
# ---------------------------
try:
    dep_opt = f"--dependency=afterok:{array_jobid}"
    res2 = subprocess.check_output(["sbatch", dep_opt, merge_sbatch]).decode().strip()
    print(res2)
    merge_jobid = res2.split()[-1]
    print(f"Merge job submitted: {merge_jobid} (will run after array {array_jobid} completes successfully)")
except subprocess.CalledProcessError as e:
    print("ERROR submitting merge job:", e.output.decode() if e.output else e)
    sys.exit(1)

# ---------------------------
# Done — print summary
# ---------------------------
print("\nSUMMARY:")
print(f"Reference doc: {REFERENCE_DOC}")
print(f"Query list: {query_list_file}")
print(f"BLAST DB list: {blastdb_list_file}")
print(f"Array script: {array_sbatch}")
print(f"Array job id: {array_jobid}")
print(f"Merge script: {merge_sbatch}")
print(f"Merge job id: {merge_jobid}")
print(f"Results root: {blast_results_dir}")
print("\nMonitor with: squeue -u $USER | grep {array_jobid}\\|{merge_jobid}")
