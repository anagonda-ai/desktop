import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Config paths
root_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/test/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_fixed"
submit_find_homolog_genes_in_dir = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_homolog_genes_in_dir.sh"
submit_find_clusters_in_chromosome = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_clusters_in_chromosome.sh"
submit_multichromosome_statistics = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_multichromosome_statistics.sh"
jobs_dir = "/groups/itay_mayrose/alongonda/desktop/example_jobs"

# Get MGC dirs
all_subdirs = [
    os.path.join(root_dir, d)
    for d in os.listdir(root_dir)
    if os.path.isdir(os.path.join(root_dir, d))
]

# Helpers
def wait_for_slot(max_jobs=100):
    while True:
        result = subprocess.run(["squeue", "-u", "alongonda"], capture_output=True, text=True)
        lines = result.stdout.strip().splitlines()
        lines = [line for line in lines if "mgc_" in line]
        running_jobs = len(lines) - 1 if len(lines) > 1 else 0
        if running_jobs < max_jobs:
            break
        print(f"⏳ Waiting: {running_jobs} jobs running. Sleeping 30s...")
        time.sleep(30)

def submit_job(job_name, out_file, err_file, script_path, mgc_path):
    cmd = [
        "sbatch",
        "--job-name", job_name,
        "--output", out_file,
        "--error", err_file,
        script_path,
        mgc_path
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ Failed to submit job {job_name}: {result.stderr}")
        return None
    job_id = result.stdout.strip().split()[-1]
    print(f"✅ Submitted {job_name} with job ID {job_id}")
    return job_id

def wait_for_job_and_file(job_id, expected_output_file):
    print(f"⏳ Waiting for job {job_id} to finish...")
    while True:
        result = subprocess.run(["squeue", "-j", job_id], capture_output=True, text=True)
        if job_id not in result.stdout:
            print(f"✅ Job {job_id} completed.")
            break
        time.sleep(30)

    print(f"⏳ Waiting for output file: {expected_output_file}")
    while not os.path.exists(expected_output_file):
        time.sleep(5)

    prev_size = -1
    while True:
        curr_size = os.path.getsize(expected_output_file)
        if curr_size == prev_size:
            break
        prev_size = curr_size
        time.sleep(5)
    print(f"✅ Output file ready: {expected_output_file}")

def process_mgc_path(mgc_path):
    mgc_name = os.path.basename(mgc_path)
    try:
        # Step 1: homolog
        wait_for_slot()
        out1 = os.path.join(jobs_dir, f"{mgc_name}_homolog.out")
        err1 = os.path.join(jobs_dir, f"{mgc_name}_homolog.err")
        if not os.path.exists(out1):
            job_id1 = submit_job("mgc_" + mgc_name + "_homolog", out1, err1, submit_find_homolog_genes_in_dir, mgc_path)
            if job_id1:
                wait_for_job_and_file(job_id1, out1)

        # Step 2: cluster
        wait_for_slot()
        out2 = os.path.join(jobs_dir, f"{mgc_name}_cluster.out")
        err2 = os.path.join(jobs_dir, f"{mgc_name}_cluster.err")
        if not os.path.exists(out2):
            job_id2 = submit_job("mgc_" + mgc_name + "_cluster", out2, err2, submit_find_clusters_in_chromosome, mgc_path)
            if job_id2:
                wait_for_job_and_file(job_id2, out2)

        # Step 3: multi
        wait_for_slot()
        out3 = os.path.join(jobs_dir, f"{mgc_name}_multi.out")
        err3 = os.path.join(jobs_dir, f"{mgc_name}_multi.err")
        if not os.path.exists(out3):
            job_id3 = submit_job("mgc_" + mgc_name + "_multi", out3, err3, submit_multichromosome_statistics, mgc_path)
            if job_id3:
                wait_for_job_and_file(job_id3, out3)

        print(f"🎉 Finished all jobs for {mgc_name}")
    except Exception as e:
        print(f"❌ Error processing {mgc_name}: {e}")

# === Main execution: run multiple MGCs in parallel ===
max_workers = 30  # number of MGCs to run in parallel
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(process_mgc_path, mgc) for mgc in all_subdirs]
    for future in as_completed(futures):
        future.result()  # to raise any exceptions
