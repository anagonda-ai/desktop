#!/usr/bin/env python3
import os
import subprocess
import shlex
import glob

# Configuration
BASE_DIR = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/"
PARTITION = "itaym-pool"  # Update this to your partition name

print("\\n===== STEP 0: GENERATE SELF-HIT BITSCORES =====")
print("===== FINDING GENE DIRECTORIES =====")

# Find all gene directories (the 100 selected genes)
gene_dirs = []
for item in sorted(glob.glob(os.path.join(BASE_DIR, "*"))):
    if os.path.isdir(item):
        gene_dirs.append(item)
        print(f"Found: {os.path.basename(item)}")

if not gene_dirs:
    print("No gene directories found!")
    exit(1)

print(f"\\nTotal gene directories to process: {len(gene_dirs)}")

# Process each gene directory
for gene_dir in gene_dirs:
    gene_name = os.path.basename(gene_dir)
    print(f"\\n===== SUBMITTING SELF-HIT JOB for {gene_name} =====")

    # Create self_hits directory if it doesn't exist
    self_hits_dir = os.path.join(gene_dir, "self_hits")
    os.makedirs(self_hits_dir, exist_ok=True)
    
    best_hits_dir = os.path.join(gene_dir, "best_hits_fixed")
    
    selected_gene_dirs = [d for d in glob.glob(os.path.join(best_hits_dir, "*")) 
                          if os.path.isdir(d)]
    
    # Find ALL FASTA files in the gene directory
    possible_fastas = glob.glob(os.path.join(gene_dir, "*.fasta")) + \
                     glob.glob(os.path.join(gene_dir, "*.fa")) + \
                     glob.glob(os.path.join(gene_dir, "*.faa"))
                     
    # Filter possible_fastas to match only those found in selected_gene_dirs
    selected_gene_names = set(os.path.basename(d) for d in selected_gene_dirs)
    filtered_fastas = []
    for fasta in possible_fastas:
        fasta_base = os.path.basename(fasta)
        fasta_name = os.path.splitext(fasta_base)[0]
        # Sometimes a fasta name might have extra suffixes; consider exact match
        if fasta_name in selected_gene_names:
            filtered_fastas.append(fasta)
    possible_fastas = filtered_fastas
    
    if not possible_fastas:
        print(f"⚠ WARNING: No FASTA files found for {gene_name}, skipping...")
        continue
    
    print(f"Found {len(possible_fastas)} FASTA files to process")
    
    selfhit_sbatch = os.path.join(gene_dir, "selfhit_blast.sbatch")
    selfhit_stdout = os.path.join(self_hits_dir, "selfhit_%j.out")
    selfhit_stderr = os.path.join(self_hits_dir, "selfhit_%j.err")

    # ---- SLURM header ----
    selfhit_script = f'''#!/bin/bash
#SBATCH -J selfhit_{gene_name}
#SBATCH -o {shlex.quote(selfhit_stdout)}
#SBATCH -e {shlex.quote(selfhit_stderr)}
#SBATCH -p {PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00

module purge
module load ncbi-blast/2.13.0+  # Adjust version as needed

echo "Starting self-hit BLAST for {gene_name}"
echo "Found {len(possible_fastas)} FASTA files to process"
echo "Start time: $(date)"
echo ""

'''

    # Add BLAST commands for each FASTA file
    for gene_fasta in possible_fastas:
        fasta_basename = os.path.basename(gene_fasta)
        fasta_name = os.path.splitext(fasta_basename)[0]
        selfhit_results = os.path.join(self_hits_dir, f"{fasta_name}_selfhit.txt")
        
        selfhit_script += f'''
echo "Processing: {fasta_basename}"
blastp \\
  -query {shlex.quote(gene_fasta)} \\
  -subject {shlex.quote(gene_fasta)} \\
  -out {shlex.quote(selfhit_results)} \\
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
  -max_target_seqs 1 \\
  -num_threads 4

if [ $? -eq 0 ]; then
    num_hits=$(grep -c "^" {shlex.quote(selfhit_results)} || echo 0)
    echo "  ✓ Completed: $num_hits self-hits found"
else
    echo "  ✗ Failed: {fasta_basename}"
fi
echo ""

'''

    selfhit_script += '''
echo "All self-hit BLASTs completed"
echo "End time: $(date)"
'''

    # Write the sbatch script
    with open(selfhit_sbatch, "w") as fh:
        fh.write(selfhit_script)

    print(f"Created sbatch script: {selfhit_sbatch}")

    # Submit the job
    cmd = ["sbatch", selfhit_sbatch]
    
    try:
        res = subprocess.check_output(cmd).decode().strip()
        job_id = res.split()[-1]
        print(f"✓ Self-hit job submitted successfully: {job_id}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error submitting job: {e}")
    except Exception as e:
        print(f"✗ Unexpected error: {e}")