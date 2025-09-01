import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Config paths
root_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/random_mgc_candidates_dir_fixed"
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

def check_active_jobs():
    """Check how many jobs are currently active for the user"""
    result = subprocess.run(["squeue", "-u", "alongonda"], capture_output=True, text=True)
    if result.returncode != 0:
        return 0
    # Count lines minus header
    lines = result.stdout.strip().split('\n')
    return max(0, len(lines) - 1)

def wait_for_job_slots(max_active_jobs=500):
    """Wait until active job count drops below threshold"""
    while True:
        active_jobs = check_active_jobs()
        if active_jobs < max_active_jobs:
            print(f"Job slots available ({active_jobs} active jobs)")
            break
        print(f"Waiting for job slots... ({active_jobs} active jobs)")
        time.sleep(60)

def submit_array_job(job_name, script_path, mgc_list_file, num_tasks, dependency_jobs=None):
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
        ]
        
        # Add dependency if specified
        if dependency_jobs:
            if isinstance(dependency_jobs, list):
                dep_string = ":".join(dependency_jobs)
            else:
                dep_string = dependency_jobs
            cmd.extend(["--dependency", f"afterok:{dep_string}"])
        
        cmd.extend([script_path, mgc_list_file])
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Failed to submit array job {job_name}: {result.stderr}")
            return None
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted array job {job_name} with job ID {job_id}")
        return [job_id]
    else:
        # Multiple chunked array jobs with dependencies between chunks
        chunks = (num_tasks + max_array_index - 1) // max_array_index  # Ceiling division
        prev_chunk_job = None
        
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
            
            # Wait for job slots before submitting (except first chunk)
            if chunk > 0:
                wait_for_job_slots(max_active_jobs=200)
            
            cmd = [
                "sbatch",
                f"--array=1-{chunk_size}",
                "--job-name", chunk_name,
                "--output", os.path.join(jobs_dir, f"{chunk_name}_%a.out"),
                "--error", os.path.join(jobs_dir, f"{chunk_name}_%a.err"),
            ]
            
            # Add dependency on previous chunk or external dependencies
            if chunk > 0 and prev_chunk_job:
                cmd.extend(["--dependency", f"afterok:{prev_chunk_job}"])
            elif dependency_jobs:
                if isinstance(dependency_jobs, list):
                    dep_string = ":".join(dependency_jobs)
                else:
                    dep_string = dependency_jobs
                cmd.extend(["--dependency", f"afterok:{dep_string}"])
            
            cmd.extend([script_path, chunk_list_file])
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Failed to submit array job {chunk_name}: {result.stderr}")
                return None
            job_id = result.stdout.strip().split()[-1]
            print(f"Submitted array job {chunk_name} with job ID {job_id}")
            job_ids.append(job_id)
            prev_chunk_job = job_id
            
            # Small delay between chunk submissions
            time.sleep(5)
        
        return job_ids

def wait_for_array_job(job_ids):
    if isinstance(job_ids, str):
        job_ids = [job_ids]
    
    print(f"Waiting for {len(job_ids)} array job(s) to finish...")
    remaining_jobs = set(job_ids)
    
    while remaining_jobs:
        completed_jobs = set()
        for job_id in remaining_jobs:
            result = subprocess.run(["squeue", "-j", job_id], capture_output=True, text=True)
            if job_id not in result.stdout:
                print(f"Array job {job_id} completed.")
                completed_jobs.add(job_id)
        
        remaining_jobs -= completed_jobs
        if remaining_jobs:
            time.sleep(30)

# === Main execution: Submit array jobs for all MGCs ===
print(f"Processing {num_mgcs} MGC directories using array jobs...")

# Wait for current jobs to clear up some space
print("Checking current job load...")
wait_for_job_slots(max_active_jobs=100)

# Step 1: Submit homolog array job
print("Submitting homolog jobs...")
job_ids1 = submit_array_job("mgc_homolog", submit_find_homolog_genes_in_dir_array, mgc_list_file, num_mgcs)
if job_ids1:
    wait_for_array_job(job_ids1)

# Step 2: Submit cluster array job (depends on homolog completion)
print("Submitting cluster jobs...")
job_ids2 = submit_array_job("mgc_cluster", submit_find_clusters_in_chromosome_array, mgc_list_file, num_mgcs, dependency_jobs=job_ids1)
if job_ids2:
    wait_for_array_job(job_ids2)

# Step 3: Submit multi array job (depends on cluster completion)
print("Submitting multi jobs...")
job_ids3 = submit_array_job("mgc_multi", submit_multichromosome_statistics_array, mgc_list_file, num_mgcs, dependency_jobs=job_ids2)
if job_ids3:
    wait_for_array_job(job_ids3)

print("All array jobs completed!")