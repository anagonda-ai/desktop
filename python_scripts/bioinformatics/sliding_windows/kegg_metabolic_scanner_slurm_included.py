import glob
import os
import subprocess
import time
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# ------------------------------
# Count Running Jobs
# ------------------------------

def count_running_jobs(user="alongonda", name_prefix=None):
    """Count number of running jobs for a user, optionally filtered by job name prefix."""
    result = subprocess.run(
        ["squeue", "-u", user, "-h", "-o", "%.100j"],
        capture_output=True, text=True
    )
    job_names = result.stdout.strip().splitlines()

    if name_prefix:
        job_names = [name for name in job_names if name.startswith(name_prefix)]

    return len(job_names)

# ------------------------------
# SLURM Submission
# ------------------------------

def submit_annotation_jobs(genome_files, temp_dir, annotated_dir, max_jobs=100):
    slurm_script = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/kegg_metabolic_annotation.sh"
    kegg_fasta = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/fasta/merged_metabolic_pathways.fasta"
    os.makedirs(annotated_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    for genome_file in genome_files:
        job_name = os.path.basename(genome_file).replace('.csv', '')
        output_file = os.path.join(annotated_dir, f"{job_name}_annotated.csv")
        if os.path.exists(output_file):
            print(f"✔️ Already exists: {output_file}")
            continue

        while True:
            current_jobs = count_running_jobs(user="alongonda")
            if current_jobs < max_jobs:
                break
            print(f"⏳ {current_jobs} jobs running (limit {max_jobs}). Waiting 30 seconds...")
            time.sleep(30)

        sbatch_cmd = [
            "sbatch",
            "--job-name", f"kegg_{job_name}",
            "--output", os.path.join(annotated_dir, f"{job_name}.out"),
            "--error", os.path.join(annotated_dir, f"{job_name}.err"),
            slurm_script,
            genome_file,
            annotated_dir,
            temp_dir,
            kegg_fasta
        ]
        result = subprocess.run(sbatch_cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✅ Submitted job for {job_name}: {result.stdout.strip()}")
        else:
            print(f"❌ Failed to submit job for {job_name}: {result.stderr.strip()}")

# ------------------------------
# Wait for SLURM Jobs
# ------------------------------

def wait_for_jobs(job_prefix="kegg_", user="alongonda", check_interval=60):
    print(f"⏳ Waiting for SLURM jobs starting with '{job_prefix}' by user '{user}'...")
    while True:
        result = subprocess.run(
            ["squeue", "-u", user, "-h", "-o", "%.18i %.9P %.100j %.8u %.2t %.10M %.6D %R"],
            capture_output=True, text=True
        )
        jobs = result.stdout.strip().splitlines()
        kegg_jobs = [j for j in jobs if job_prefix in j]

        if not kegg_jobs:
            print("✅ All KEGG annotation jobs are finished.")
            break
        else:
            print(f"🕒 {len(kegg_jobs)} jobs still running. Checking again in {check_interval} seconds...")
            time.sleep(check_interval)

# ------------------------------
# Subset Filtering
# ------------------------------

def remove_subset_results(output_file):
    df = pd.read_csv(output_file)
    df['metabolic_genes_set'] = df['metabolic_genes'].apply(lambda x: set(x.split(',')))
    to_remove = set()
    for i, genes_i in enumerate(df['metabolic_genes_set']):
        for j, genes_j in enumerate(df['metabolic_genes_set']):
            if i != j and genes_i < genes_j:
                to_remove.add(i)
                break
    filtered_df = df.drop(list(to_remove)).drop(columns=['metabolic_genes_set'])
    filtered_output_file = output_file.replace(".csv", "_filtered.csv")
    filtered_df.to_csv(filtered_output_file, index=False)
    print(f"🧹 Removed {len(to_remove)} subset results: {filtered_output_file}")

# ------------------------------
# Main Function
# ------------------------------

def main():
    # Directories
    full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
    genome_dirs = [
        os.path.join(full_genome_dir, "final_dataset/filtered_no_mito")
    ]
    max_jobs = 100
    min_genes = 3
    head_output_dir = f"/groups/itay_mayrose/alongonda/Plant_MGC/test_metabolic"
    output_dir = os.path.join(head_output_dir, "kegg_scanner_min_genes_based_metabolic")
    temp_dir = os.path.join(head_output_dir, "blast_temp_annotated_metabolic")
    annotated_dir = os.path.join(head_output_dir, "annotated_genomes_metabolic")
    process_annotated_file_slurm_script = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/process_annotated_file_slurm_script.sh"

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(annotated_dir, exist_ok=True)

    # Genome files
    genome_files = []
    for genome_dir in genome_dirs:
        for file in os.listdir(genome_dir):
            if file.endswith('.csv'):
                genome_files.append(os.path.join(genome_dir, file))

    print(f"Found {len(genome_files)} genome files.")

    # Submit SLURM jobs for annotation
    submit_annotation_jobs(genome_files, temp_dir, annotated_dir, max_jobs=100)
    
    # Wait until all jobs finish
    wait_for_jobs(job_prefix="kegg_", user="alongonda")

    # Sliding Window (after annotation completes manually)
    annotated_files = [os.path.join(annotated_dir, f) for f in os.listdir(annotated_dir) if f.endswith('_annotated.csv')]

    min_genes_subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    os.makedirs(min_genes_subdir, exist_ok=True)
    tmp_output_dir = os.path.join(min_genes_subdir, "tmp_outputs")
    os.makedirs(tmp_output_dir, exist_ok=True)

    for window_size in [10, 20]:
        output_file = os.path.join(min_genes_subdir, f"potential_groups_w{window_size}.csv")
        err_out_output_dir = os.path.join(min_genes_subdir, "err_out")
        os.makedirs(err_out_output_dir, exist_ok=True)
        if os.path.exists(output_file):
            os.remove(output_file)

        total_matches = 0

        for annotated_file in annotated_files:
            while True:
                current_jobs = count_running_jobs(user="alongonda")
                if current_jobs < max_jobs:
                    break
                print(f"⏳ {current_jobs} jobs running (limit {max_jobs}). Waiting 30 seconds...")
                time.sleep(30)
            
            # Submit SLURM job for processing annotated file
            job_name = os.path.basename(annotated_file).replace('.csv', '')
            tmp_output_file = os.path.join(tmp_output_dir, f"{job_name}_w{window_size}.csv")
            sbatch_cmd = [
                "sbatch",
                "--job-name", f"process_{job_name}",
                "--output", os.path.join(err_out_output_dir, f"{job_name}_{window_size}.out"),
                "--error", os.path.join(err_out_output_dir, f"{job_name}_{window_size}.err"),
                process_annotated_file_slurm_script,
                annotated_file,
                tmp_output_file,
                str(window_size),
                str(min_genes)
            ]
            result = subprocess.run(sbatch_cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✅ Submitted job for {job_name}: {result.stdout.strip()}")
            else:
                print(f"❌ Failed to submit job for {job_name}: {result.stderr.strip()}")
                
        # המתנה לסיום כל העבודות
        print("⏳ Waiting for all SLURM jobs to finish...")
        while True:
            running_jobs = count_running_jobs(user="alongonda", name_prefix="process_")
            if running_jobs == 0:
                print("✅ All jobs completed.")
                break
            print(f"⏳ {running_jobs} jobs still running. Waiting 60 seconds...")
            time.sleep(60)
        
        # איחוד כל קבצי הפלט
        merged_dfs = []
        tmp_files = sorted(glob.glob(os.path.join(tmp_output_dir, f"*w{window_size}.csv")))
        for tmp_file in tmp_files:
            try:
                df = pd.read_csv(tmp_file)
                merged_dfs.append(df)
            except Exception as e:
                print(f"⚠️ Failed to read {tmp_file}: {e}")

        if merged_dfs:
            merged = pd.concat(merged_dfs, ignore_index=True)
            merged.to_csv(output_file, index=False)
            print(f"✔️ Merged {len(tmp_files)} files into {output_file} (Total rows: {len(merged)})")
        else:
            print("⚠️ No output files to merge.")

        print(f"✔️ Total Matches for w{window_size}, min_genes={min_genes}: {total_matches}")
        remove_subset_results(output_file)


if __name__ == "__main__":
    main()
