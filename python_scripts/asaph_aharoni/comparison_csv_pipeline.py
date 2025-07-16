import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Config paths
root_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir"
submit_find_homolog_genes_in_dir_array = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_homolog_genes_in_dir_array.sh"
submit_find_clusters_in_chromosome_array = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_clusters_in_chromosome_array.sh"
submit_multichromosome_statistics_array = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_multichromosome_statistics_array.sh"
jobs_dir = "/groups/itay_mayrose/alongonda/desktop/example_jobs"

# Get MGC dirs
all_subdirs = [
    os.path.join(root_dir, d)
    for d in os.listdir(root_dir)
    if os.path.isdir(os.path.join(root_dir, d))
]

# Create MGC list file for array jobs
mgc_list_file = os.path.join(jobs_dir, "mgc_list.txt")
with open(mgc_list_file, 'w') as f:
    for mgc_path in all_subdirs:
        f.write(mgc_path + '\n')

num_mgcs = len(all_subdirs)

# Helpers

def submit_array_job(job_name, script_path, mgc_list_file, num_tasks):
    max_array_index = 1000  # Maximum allowed array index on this cluster
    job_ids = []
    
    if num_tasks <= max_array_index:
        # Single array job
        cmd = [
            "sbatch",
            f"--array=1-{num_tasks}",
            "--job-name", job_name,
            "--output", os.path.join(jobs_dir, f"{job_name}_%a.out"),
            "--error", os.path.join(jobs_dir, f"{job_name}_%a.err"),
            script_path,
            mgc_list_file
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"❌ Failed to submit array job {job_name}: {result.stderr}")
            return None
        job_id = result.stdout.strip().split()[-1]
        print(f"✅ Submitted array job {job_name} with job ID {job_id}")
        return [job_id]
    else:
        # Multiple chunked array jobs with separate list files
        chunks = (num_tasks + max_array_index - 1) // max_array_index  # Ceiling division
        for chunk in range(chunks):
            start_task = chunk * max_array_index + 1
            end_task = min((chunk + 1) * max_array_index, num_tasks)
            chunk_size = end_task - start_task + 1
            chunk_name = f"{job_name}_chunk{chunk + 1}"
            
            # Create a separate list file for this chunk
            chunk_list_file = os.path.join(jobs_dir, f"mgc_list_chunk{chunk + 1}.txt")
            with open(mgc_list_file, 'r') as f:
                all_lines = f.readlines()
            
            with open(chunk_list_file, 'w') as f:
                for i in range(start_task - 1, end_task):
                    f.write(all_lines[i])
            
            cmd = [
                "sbatch",
                f"--array=1-{chunk_size}",
                "--job-name", chunk_name,
                "--output", os.path.join(jobs_dir, f"{chunk_name}_%a.out"),
                "--error", os.path.join(jobs_dir, f"{chunk_name}_%a.err"),
                script_path,
                chunk_list_file
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"❌ Failed to submit array job {chunk_name}: {result.stderr}")
                return None
            job_id = result.stdout.strip().split()[-1]
            print(f"✅ Submitted array job {chunk_name} with job ID {job_id}")
            job_ids.append(job_id)
        
        return job_ids

def wait_for_array_job(job_ids):
    if isinstance(job_ids, str):
        job_ids = [job_ids]
    
    print(f"⏳ Waiting for {len(job_ids)} array job(s) to finish...")
    remaining_jobs = set(job_ids)
    
    while remaining_jobs:
        completed_jobs = set()
        for job_id in remaining_jobs:
            result = subprocess.run(["squeue", "-j", job_id], capture_output=True, text=True)
            if job_id not in result.stdout:
                print(f"✅ Array job {job_id} completed.")
                completed_jobs.add(job_id)
        
        remaining_jobs -= completed_jobs
        if remaining_jobs:
            time.sleep(30)

# === Main execution: Submit array jobs for all MGCs ===
print(f"Processing {num_mgcs} MGC directories using array jobs...")

# Step 1: Submit homolog array job
job_ids1 = submit_array_job("mgc_homolog", submit_find_homolog_genes_in_dir_array, mgc_list_file, num_mgcs)
if job_ids1:
    wait_for_array_job(job_ids1)

# Step 2: Submit cluster array job
job_ids2 = submit_array_job("mgc_cluster", submit_find_clusters_in_chromosome_array, mgc_list_file, num_mgcs)
if job_ids2:
    wait_for_array_job(job_ids2)

# Step 3: Submit multi array job
job_ids3 = submit_array_job("mgc_multi", submit_multichromosome_statistics_array, mgc_list_file, num_mgcs)
if job_ids3:
    wait_for_array_job(job_ids3)

print("🎉 All array jobs completed!")
