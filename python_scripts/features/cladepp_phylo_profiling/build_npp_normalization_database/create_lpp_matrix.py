#!/usr/bin/env python3
import os
import subprocess
import shlex
import glob

# Configuration
BASE_DIR = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/"
PARTITION = "itaym-pool"  # Update this to your partition name

print("\n===== STEP 2: CREATE LPP MATRICES =====")
print("===== FINDING DIRECTORIES WITH MERGED RESULTS =====")

# Find all directories with best_hits subdirectories that have merged CSV files
mgc_dirs = []
for item in sorted(glob.glob(os.path.join(BASE_DIR, "*"))):
    if os.path.isdir(item):
        best_hits_dir = os.path.join(item, "best_hits_fixed")
        self_hits_dir = os.path.join(item, "self_hits")
        if os.path.exists(best_hits_dir) and os.path.isdir(best_hits_dir):
            if os.path.exists(self_hits_dir) and os.path.isdir(self_hits_dir):
                # Check if merged CSV files exist
                merged_files = glob.glob(os.path.join(best_hits_dir, "best_hits_*.csv"))
                if merged_files:
                    mgc_dirs.append(item)
                    print(f"Found: {os.path.basename(item)} ({len(merged_files)} merged files)")
                else:
                    print(f"Skipping {os.path.basename(item)}: No merged CSV files found")
            else:
                print(f"Skipping {os.path.basename(item)}: No self_hits directory found")

if not mgc_dirs:
    print("No directories with both merged CSV files and self_hits found!")
    print("Make sure Step 0 (self-hits) and Step 1 (merge) have completed successfully.")
    exit(1)

print(f"\nTotal directories to process: {len(mgc_dirs)}")

# Process each directory
for example_mgc in mgc_dirs:
    mgc_name = os.path.basename(example_mgc)
    print(f"\n===== SUBMITTING LPP JOB for {mgc_name} =====")

    blast_results_dir = os.path.join(example_mgc, "best_hits_fixed")
    self_hits_dir = os.path.join(example_mgc, "self_hits")
    
    lpp_sbatch = os.path.join(example_mgc, "create_lpp_matrix.sbatch")
    lpp_stdout = os.path.join(blast_results_dir, "lpp_%j.out")
    lpp_stderr = os.path.join(blast_results_dir, "lpp_%j.err")

    # ---- Part 1: SLURM header ----
    lpp_script = f"""#!/bin/bash
#SBATCH -J lpp_{mgc_name}
#SBATCH -o {shlex.quote(lpp_stdout)}
#SBATCH -e {shlex.quote(lpp_stderr)}
#SBATCH -p {PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2:00:00

module purge
module load python/3.11

python - <<'PYCODE'
"""

    # ---- Part 2: Python LPP code (embedded, optimized for speed) ----
    lpp_script += f"""
import os, glob, csv
from collections import defaultdict
import sys

BLAST_ROOT = {repr(blast_results_dir)}
SELF_HITS_DIR = {repr(self_hits_dir)}

print("="*70)
print("LPP MATRIX GENERATION (OPTIMIZED)")
print("="*70)

# Step 1: Load self-hit bitscores from Step 0 output (FAST - single dict)
print("\\nStep 1: Loading self-hit bitscores...")
sys.stdout.flush()

self_hit_bitscores = {{}}
selfhit_files = sorted(glob.glob(os.path.join(SELF_HITS_DIR, "*_selfhit.txt")))

if not selfhit_files:
    print(f"ERROR: No self-hit files found in {{SELF_HITS_DIR}}")
    exit(1)

total_selfhit_count = 0
for selfhit_file in selfhit_files:
    filename = os.path.basename(selfhit_file)
    try:
        with open(selfhit_file, "r") as f:
            for line in f:
                # Fast split without regex
                parts = line.rstrip().split("\\t")
                if len(parts) >= 12:
                    self_hit_bitscores[parts[0]] = float(parts[11])
                    total_selfhit_count += 1
        print(f"  ✓ {{filename}}: {{len([k for k in self_hit_bitscores if k])}} unique genes")
    except Exception as e:
        print(f"  ✗ ERROR reading {{filename}}: {{e}}")

print(f"Total self-hit bitscores: {{len(self_hit_bitscores)}}")

if not self_hit_bitscores:
    print("ERROR: No self-hit bitscores loaded!")
    exit(1)

# Step 2: Load merged BLAST results (FAST - use dict readers efficiently)
print("\\nStep 2: Loading merged BLAST results...")
sys.stdout.flush()

merged_files = sorted(glob.glob(os.path.join(BLAST_ROOT, "best_hits_*.csv")))

if not merged_files:
    print(f"ERROR: No merged CSV files in {{BLAST_ROOT}}")
    exit(1)

all_records_by_organism = defaultdict(list)

for csv_file in merged_files:
    organism = os.path.basename(csv_file).replace("best_hits_", "").replace(".csv", "")
    record_count = 0
    try:
        with open(csv_file, "r", encoding="utf-8") as incsv:
            reader = csv.DictReader(incsv)
            # Collect all rows for this organism
            records = list(reader)
            all_records_by_organism[organism] = records
            record_count = len(records)
        print(f"  ✓ {{organism}}: {{record_count}} records")
    except Exception as e:
        print(f"  ✗ ERROR reading {{csv_file}}: {{e}}")

# Step 3: Create LPP matrices (FAST - vectorized writes)
print(f"\\n{{'='*70}}")
print("Step 3: Creating LPP matrices...")
print(f"{{'='*70}}\\n")
sys.stdout.flush()

for organism, records in sorted(all_records_by_organism.items()):
    lpp_out_path = os.path.join(BLAST_ROOT, f"LPP_matrix_{{organism}}")
    
    lpp_calculated = 0
    genes_without_selfhit = set()
    
    try:
        # Pre-compute all rows before writing (reduces I/O overhead)
        output_rows = []
        
        for rec in records:
            query_gene = rec['query_gene']
            bitscore = float(rec['bitscore'])
            self_bitscore = self_hit_bitscores.get(query_gene)
            
            if self_bitscore and self_bitscore > 0:
                lpp = bitscore / self_bitscore
                lpp_calculated += 1
            else:
                lpp = None
                genes_without_selfhit.add(query_gene)
            
            output_rows.append([
                rec['query_dir'], rec['query_gene'], rec['subject_gene'],
                rec['pident'], rec['length'], rec['mismatch'], rec['gapopen'],
                rec['qstart'], rec['qend'], rec['sstart'], rec['send'],
                rec['evalue'], rec['bitscore'], self_bitscore, lpp, rec['source_file']
            ])
        
        # Write all rows at once (much faster)
        with open(lpp_out_path, "w", newline='', encoding="utf-8") as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow([
                "query_dir","query_gene","subject_gene","pident","length",
                "mismatch","gapopen","qstart","qend","sstart","send",
                "evalue","bitscore","self_bitscore","LPP","source_file"
            ])
            writer.writerows(output_rows)
        
        print(f"✓ (LPP: {{lpp_calculated}}/{{len(records)}})")
        
        if genes_without_selfhit:
            print(f"    ⚠ {{len(genes_without_selfhit)}} genes without self-hit")
        
    except Exception as e:
        print(f"✗ ERROR: {{e}}")

print(f"\\n{{'='*70}}")
print("LPP MATRIX GENERATION COMPLETED!")
print(f"{{'='*70}}")
PYCODE
"""

    # Write the sbatch script
    with open(lpp_sbatch, "w") as fh:
        fh.write(lpp_script)

    print(f"Created sbatch script: {lpp_sbatch}")

    # Submit the job
    cmd = ["sbatch", lpp_sbatch]
    
    try:
        res = subprocess.check_output(cmd).decode().strip()
        job_id = res.split()[-1]
        print(f"✓ LPP job submitted successfully: {job_id}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error submitting job: {e}")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")

print("\n===== ALL LPP JOBS SUBMITTED =====")