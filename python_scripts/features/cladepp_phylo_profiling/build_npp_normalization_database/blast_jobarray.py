#!/usr/bin/env python3
"""
Two-stage BLAST launcher for ALL MGC directories:

Stage 1 (serial):
    For each MGC directory:
        - create BLAST array
        - submit array
        - next MGC starts only after previous array completes

Stage 2 (parallel):
    For each MGC directory:
        - submit merge job
        - each merge depends on its own array job

TAU HPC (V2): itay_mayrose-users_v2, itaym-pool
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

# Reference doc
REFERENCE_DOC = "/mnt/data/Highlits for the Toolkit.txt"

# ----------------------------------------------------------------------
# Args
# ----------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--cpus", type=int, default=32)
parser.add_argument("--mem", default="32G")
parser.add_argument("--array_time", default="2-00:00:00")
parser.add_argument("--array_concurrency", type=int, default=100)
parser.add_argument("--blast_module", default="ncbi-blast/2.13.0+")
parser.add_argument("--partition", default="itaym-pool")
args = parser.parse_args()

# ----------------------------------------------------------------------
# Directories
# ----------------------------------------------------------------------
selected_genes_root = (
    "/groups/itay_mayrose/alongonda/datasets/"
    "KEGG_annotations_modules_metabolic/selected_genes"
)

mgc_dirs = sorted([
    os.path.join(selected_genes_root, d)
    for d in os.listdir(selected_genes_root)
    if os.path.isdir(os.path.join(selected_genes_root, d))
])

if not mgc_dirs:
    sys.exit("ERROR: No MGC directories found.")

print(f"Found {len(mgc_dirs)} MGC directories.")

# ----------------------------------------------------------------------
# Static dataset paths
# ----------------------------------------------------------------------
full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
csv_dirs = [os.path.join(full_genome_dir, "final_dataset/filtered_no_mito")]

blast_db_dir = "/groups/itay_mayrose/alongonda/datasets/generated_blast_dbs"
os.makedirs(blast_db_dir, exist_ok=True)

# ----------------------------------------------------------------------
# CSV → FASTA (pre-stage)
# ----------------------------------------------------------------------
def csv_to_fasta(csv_path):
    fasta_filename = os.path.splitext(os.path.basename(csv_path))[0] + ".fasta"
    fasta_path = os.path.join(blast_db_dir, fasta_filename)
    valid_entries = 0
    if os.path.exists(fasta_path):
        return fasta_path
    try:
        with open(csv_path, "r", newline="") as csv_in, open(fasta_path, "w") as out:
            reader = csv.DictReader(csv_in, delimiter=",")
            for row_index, row in enumerate(reader, 1):
                if all(k in row for k in ["id","sequence","chromosome","start","end","strand"]) \
                   and row["sequence"].strip():
                    out.write(
                        f">{row['id']}|{row['chromosome']}|{row['start']}|"
                        f"{row['end']}|{row['strand']}|{row_index}\n"
                    )
                    out.write(row["sequence"].strip() + "\n")
                    valid_entries += 1
        if valid_entries == 0:
            os.remove(fasta_path)
            return None
        return fasta_path
    except:
        return None

csv_files = []
for d in csv_dirs:
    if os.path.isdir(d):
        csv_files.extend([
            os.path.join(d, f)
            for f in os.listdir(d) if f.endswith(".csv")
        ])

csv_to_fasta_map = {}
if csv_files:
    with concurrent.futures.ProcessPoolExecutor() as ex:
        futures = {ex.submit(csv_to_fasta, csv): csv for csv in csv_files}
        for fut in concurrent.futures.as_completed(futures):
            fasta = fut.result()
            if fasta:
                csv_to_fasta_map[os.path.basename(futures[fut])] = fasta

# ----------------------------------------------------------------------
# Build BLAST DBs
# ----------------------------------------------------------------------
def create_blast_db(fasta):
    db_name = fasta + "_blastdb"
    if os.path.exists(db_name + ".pin") or os.path.exists(db_name + ".phr"):
        return db_name
    try:
        subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "prot", "-out", db_name],
                       check=True)
        return db_name
    except:
        return None

blast_dbs = {}
if csv_to_fasta_map:
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        futures = {
            ex.submit(create_blast_db, fasta): csv
            for csv, fasta in csv_to_fasta_map.items()
        }
        for fut in concurrent.futures.as_completed(futures):
            db = fut.result()
            if db:
                blast_dbs[futures[fut]] = db

# ----------------------------------------------------------------------
# Stage 1: SERIAL BLAST ARRAYS
# ----------------------------------------------------------------------
array_jobids = {}   # per MGC
previous_array_jobid = None  # serial chaining

for example_mgc in mgc_dirs:
    print(f"\n===== STAGE 1: BLAST for {example_mgc} =====")

    blast_results_dir = os.path.join(example_mgc, "best_hits_fixed")
    os.makedirs(blast_results_dir, exist_ok=True)

    # Collect queries for this MGC
    query_fastas = sorted([
        os.path.join(example_mgc, f)
        for f in os.listdir(example_mgc)
        if f.endswith(".fasta")
    ])[:100]

    if not query_fastas:
        print(f"WARNING: No FASTA files in {example_mgc}")
        continue

    # Write query list
    query_list_file = os.path.join(example_mgc, "query_fastas.txt")
    with open(query_list_file, "w") as fh:
        for q in query_fastas:
            fh.write(q + "\n")

    # Write BLAST DB list
    blastdb_list_file = os.path.join(example_mgc, "blast_dbs.txt")
    with open(blastdb_list_file, "w") as fh:
        for key, db in sorted(blast_dbs.items()):
            fh.write(db + "\n")

    # Create array script
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

    # SERIAL submission
    cmd = ["sbatch"]

    cmd.append(array_sbatch)

    res = subprocess.check_output(cmd).decode().strip()
    array_jobid = res.split()[-1]
    array_jobids[example_mgc] = array_jobid
    previous_array_jobid = array_jobid

    print(f"Submitted array for {example_mgc}: {array_jobid}")

# ----------------------------------------------------------------------
# Stage 2: PARALLEL MERGES (each depends only on its own array)
# ----------------------------------------------------------------------
for example_mgc in mgc_dirs:
    print(f"\n===== STAGE 2: MERGE for {example_mgc} =====")

    blast_results_dir = os.path.join(example_mgc, "best_hits")

    merge_sbatch = os.path.join(example_mgc, "merge_best_hits.sbatch")
    merge_stdout = os.path.join(blast_results_dir, "merge_%j.out")
    merge_stderr = os.path.join(blast_results_dir, "merge_%j.err")

    array_jobid = array_jobids.get(example_mgc)
    if not array_jobid:
        print(f"Skipping merge for {example_mgc} — no array job.")
        continue

    # ---- Part 1: SLURM header ----
    merge_script = f"""#!/bin/bash
#SBATCH -J merge_best_hits
#SBATCH -o {shlex.quote(merge_stdout)}
#SBATCH -e {shlex.quote(merge_stderr)}
#SBATCH -p {args.partition}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=5:00:00

module purge
module load python/3.11

python - <<'PYCODE'
"""

    # ---- Part 2: Python merge (no f-string) ----
    merge_script += f"""
import os, glob, csv, re
from collections import defaultdict

BLAST_ROOT = {repr(blast_results_dir)}
"""

    merge_script += r"""
def process_result_file(args):
    qbase, fpath = args
    local_best = {}
    try:
        with open(fpath, "r") as fh:
            for ln in fh:
                cols = ln.rstrip().split("\t")
                if len(cols) < 12:
                    continue

                qid    = cols[0]
                sid    = cols[1]
                pident = cols[2]
                length = cols[3]
                mismatch = cols[4]
                gapopen  = cols[5]
                qstart   = cols[6]
                qend     = cols[7]
                sstart   = cols[8]
                send     = cols[9]
                evalue   = float(cols[10])
                bitscore = float(cols[11])

                key = (qbase, qid, fpath)
                rec = dict(
                    query_dir=qbase, query_gene=qid, subject_gene=sid,
                    pident=pident, length=length, mismatch=mismatch, gapopen=gapopen,
                    qstart=qstart, qend=qend, sstart=sstart, send=send,
                    evalue=evalue, bitscore=bitscore, src=fpath
                )

                if key not in local_best or \
                   bitscore > local_best[key]['bitscore'] or \
                   (bitscore == local_best[key]['bitscore'] and evalue < local_best[key]['evalue']):
                    local_best[key] = rec
    except Exception as e:
        print("WARN reading", fpath, e)
    return local_best

# Gather all (qbase, fpath) pairs for BLAST result files
file_args = []
for qdir in sorted(glob.glob(os.path.join(BLAST_ROOT, "*"))):
    if not os.path.isdir(qdir):
        continue
    qbase = os.path.basename(qdir)
    for fpath in sorted(glob.glob(os.path.join(qdir, "*_results.txt"))):
        file_args.append((qbase, fpath))

# Use ThreadPoolExecutor for concurrency (since this is I/O bound)
best_per_query = {}
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = [executor.submit(process_result_file, arg) for arg in file_args]
    for fut in concurrent.futures.as_completed(futures):
        local_best = fut.result()
        for key, rec in local_best.items():
            if key not in best_per_query or \
               rec['bitscore'] > best_per_query[key]['bitscore'] or \
               (rec['bitscore'] == best_per_query[key]['bitscore'] and rec['evalue'] < best_per_query[key]['evalue']):
                best_per_query[key] = rec

# Create a mapping: organism -> list of records

organism_records = defaultdict(list)

for key, rec in best_per_query.items():
    source_file = os.path.basename(rec["src"])
    organism = "_".join(source_file.split("_")[2:])
    organism_records[organism].append(rec)

# Write separate output file for each organism
for organism, records in organism_records.items():
    out_path = os.path.join(BLAST_ROOT, f"best_hits_{organism}.csv")
    with open(out_path, "w", newline='', encoding="utf-8") as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow([
            "query_dir","query_gene","subject_gene","pident","length",
            "mismatch","gapopen","qstart","qend","sstart","send",
            "evalue","bitscore","source_file"
        ])
        for rec in records:
            writer.writerow([
                rec['query_dir'], rec['query_gene'], rec['subject_gene'],
                rec['pident'], rec['length'], rec['mismatch'], rec['gapopen'],
                rec['qstart'], rec['qend'], rec['sstart'], rec['send'],
                rec['evalue'], rec['bitscore'], rec['src']
            ])
    print(f"Wrote per-organism file: {out_path}")
PYCODE
"""

    with open(merge_sbatch, "w") as fh:
        fh.write(merge_script)

    # Submit merge, depends only on its MGC's array job
    res = subprocess.check_output([
        "sbatch", f"--dependency=afterok:{array_jobid}", merge_sbatch
    ]).decode().strip()

    print(f"Submitted merge for {example_mgc}: {res.split()[-1]}")

# ----------------------------------------------------------------------
# Done
# ----------------------------------------------------------------------
print("\n===== ALL JOBS SUBMITTED =====")
print("Stage 1 (arrays) ran serially.")
print("Stage 2 (merges) submitted in parallel with dependencies.")
