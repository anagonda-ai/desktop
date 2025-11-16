import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob

# Config paths
root_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_candidates_dir_fixed"
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

# Initial MGC count
total_mgcs = len(all_subdirs)
print(f"Found {total_mgcs} total MGC directories")


# Helpers

def get_job_limits():
    """Get current job limits and usage"""
    try:
        # Check submitted jobs (including pending)
        result = subprocess.run(["squeue", "-u", "alongonda", "-t", "PENDING,RUNNING"], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            submitted_jobs = max(0, len(lines) - 1) if len(lines) > 1 else 0
        else:
            submitted_jobs = 0
        
        # Be much more conservative - assume lower limits based on your error
        max_submit_limit = 100  # Much more conservative based on your actual limits
        
        return submitted_jobs, max_submit_limit
    except Exception as e:
        print(f"Warning: Could not determine job limits: {e}")
        return 0, 50  # Very conservative fallback

def check_active_jobs():
    """Check how many jobs are currently active for the user"""
    result = subprocess.run(["squeue", "-u", "alongonda"], capture_output=True, text=True)
    if result.returncode != 0:
        return 0
    # Count lines minus header
    lines = result.stdout.strip().split('\n')
    return max(0, len(lines) - 1)

def wait_for_job_slots(max_active_jobs=500, max_submitted_jobs=None):
    """Wait until active job count drops below threshold and check submission limits"""
    if max_submitted_jobs is None:
        submitted_jobs, max_submit_limit = get_job_limits()
        max_submitted_jobs = max_submit_limit - 50  # Leave buffer

    while True:
        active_jobs = check_active_jobs()
        submitted_jobs, _ = get_job_limits()
        
        if active_jobs < max_active_jobs and submitted_jobs < max_submitted_jobs:
            print(f"Job slots available ({active_jobs} active, {submitted_jobs} submitted)")
            break
        print(f"Waiting for job slots... ({active_jobs} active, {submitted_jobs} submitted)")
        time.sleep(60)

def submit_array_job_rolling(job_name, script_path, mgc_list_file, num_tasks, dependency_jobs=None):
    """Submit chunks with rolling submission - new chunk starts when any previous chunk finishes"""
    
    submitted_jobs, max_submit_limit = get_job_limits()
    max_array_index = 500
    max_concurrent_chunks = 10  # Maximum chunks running simultaneously
    
    chunks = (num_tasks + max_array_index - 1) // max_array_index
    print(f"Rolling submission: {num_tasks} tasks in {chunks} chunks, max {max_concurrent_chunks} concurrent")
    
    if num_tasks <= max_array_index:
        # Single chunk - use original logic
        return submit_single_array_job(job_name, script_path, mgc_list_file, num_tasks, dependency_jobs)
    
    # Prepare all chunk data upfront
    chunk_data = []
    with open(mgc_list_file, 'r') as f:
        all_lines = f.readlines()
    
    for chunk in range(chunks):
        start_task = chunk * max_array_index + 1
        end_task = min((chunk + 1) * max_array_index, num_tasks)
        chunk_size = end_task - start_task + 1
        chunk_name = f"{job_name}_chunk{chunk + 1:03d}"
        
        # Create chunk list file
        chunk_list_file = os.path.join(jobs_dir, f"mgc_list_chunk{chunk + 1:03d}_{job_name}.txt")
        with open(chunk_list_file, 'w') as f:
            for i in range(start_task - 1, end_task):
                if i < len(all_lines):
                    f.write(all_lines[i])
        
        chunk_data.append({
            'chunk_id': chunk,
            'chunk_name': chunk_name,
            'chunk_list_file': chunk_list_file,
            'chunk_size': chunk_size,
            'submitted': False,
            'job_id': None,
            'completed': False
        })
    
    # Rolling submission logic
    job_ids = []
    chunks_submitted = 0
    chunks_completed = 0
    
    # Submit initial batch up to concurrent limit
    initial_submissions = min(max_concurrent_chunks, chunks)
    
    for i in range(initial_submissions):
        job_id = submit_single_chunk(chunk_data[i], script_path, dependency_jobs if i == 0 else None)
        if job_id:
            chunk_data[i]['submitted'] = True
            chunk_data[i]['job_id'] = job_id
            job_ids.append(job_id)
            chunks_submitted += 1
            print(f"Initial submission {i+1}/{initial_submissions}: {chunk_data[i]['chunk_name']} -> {job_id}")
        time.sleep(30)  # Pace initial submissions
    
    # Rolling submission loop
    while chunks_completed < chunks:
        # Check which jobs have completed
        newly_completed = []
        
        for chunk in chunk_data:
            if chunk['submitted'] and not chunk['completed']:
                result = subprocess.run(["squeue", "-j", chunk['job_id']], capture_output=True, text=True)
                if chunk['job_id'] not in result.stdout:
                    chunk['completed'] = True
                    chunks_completed += 1
                    newly_completed.append(chunk)
                    print(f"Chunk completed: {chunk['chunk_name']} (job {chunk['job_id']})")
        
        # Submit new chunks for each completed one
        for _ in newly_completed:
            # Find next unsubmitted chunk
            next_chunk = None
            for chunk in chunk_data:
                if not chunk['submitted']:
                    next_chunk = chunk
                    break
            
            if next_chunk:
                # Wait for submission slot
                wait_for_job_slots(max_active_jobs=30, max_submitted_jobs=max_submit_limit - 30)
                
                # Submit next chunk
                job_id = submit_single_chunk(next_chunk, script_path, dependency_jobs=None)
                if job_id:
                    next_chunk['submitted'] = True
                    next_chunk['job_id'] = job_id
                    job_ids.append(job_id)
                    chunks_submitted += 1
                    print(f"Rolling submission [{chunks_submitted}/{chunks}]: {next_chunk['chunk_name']} -> {job_id}")
                else:
                    print(f"Failed to submit {next_chunk['chunk_name']}")
        
        # Progress update
        active_chunks = sum(1 for c in chunk_data if c['submitted'] and not c['completed'])
        print(f"Progress: {chunks_completed}/{chunks} completed, {active_chunks} active, {chunks_submitted}/{chunks} submitted")
        
        # Wait before checking again
        time.sleep(120)  # Check every 2 minutes
    
    print(f"All {chunks} chunks submitted and completed!")
    return job_ids

def submit_single_chunk(chunk_data, script_path, dependency_jobs=None):
    """Submit a single chunk with error handling"""
    
    cmd = [
        "sbatch",
        f"--array=1-{chunk_data['chunk_size']}",
        "--job-name", chunk_data['chunk_name'],
        "--output", os.path.join(jobs_dir, f"{chunk_data['chunk_name']}_%a.out"),
        "--error", os.path.join(jobs_dir, f"{chunk_data['chunk_name']}_%a.err"),
    ]
    
    if dependency_jobs:
        dep_string = ":".join(dependency_jobs) if isinstance(dependency_jobs, list) else dependency_jobs
        cmd.extend(["--dependency", f"afterok:{dep_string}"])
    
    cmd.extend([script_path, chunk_data['chunk_list_file']])
    
    # Submit with retry
    for retry in range(3):
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            return job_id
        else:
            print(f"Submission failed (attempt {retry+1}): {result.stderr}")
            if "QOSMaxSubmitJobPerUserLimit" in result.stderr:
                print("Hit submission limit, waiting 10 minutes...")
                time.sleep(600)
            else:
                time.sleep(120)
    
    return None

def submit_single_array_job(job_name, script_path, mgc_list_file, num_tasks, dependency_jobs):
    """Handle single array job submission"""
    wait_for_job_slots(max_active_jobs=50, max_submitted_jobs=get_job_limits()[1] - 20)
    
    cmd = [
        "sbatch",
        f"--array=1-{num_tasks}",
        "--job-name", job_name,
        "--output", os.path.join(jobs_dir, f"{job_name}_%a.out"),
        "--error", os.path.join(jobs_dir, f"{job_name}_%a.err"),
    ]
    
    if dependency_jobs:
        dep_string = ":".join(dependency_jobs) if isinstance(dependency_jobs, list) else dependency_jobs
        cmd.extend(["--dependency", f"afterok:{dep_string}"])
    
    cmd.extend([script_path, mgc_list_file])
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Failed to submit array job {job_name}: {result.stderr}")
        return None
    
    job_id = result.stdout.strip().split()[-1]
    print(f"Submitted array job {job_name} with job ID {job_id}")
    return [job_id]

def wait_for_array_job(job_ids):
    if isinstance(job_ids, str):
        job_ids = [job_ids]
    
    if not job_ids:
        print("No jobs to wait for")
        return
    
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
            time.sleep(60)  # Check every minute

def submit_array_job_for_stage(stage, script_path, mgc_list, dependency_jobs=None):
    """Submit array job for a specific stage with filtered MGC list"""
    if not mgc_list:
        print(f"No MGCs need processing for {stage} stage - skipping")
        return []
    
    # Create stage-specific MGC list file
    stage_list_file = os.path.join(jobs_dir, f"mgc_list_{stage}.txt")
    with open(stage_list_file, 'w') as f:
        for mgc_path in mgc_list:
            f.write(mgc_path + '\n')
    
    job_name = f"mgc_{stage}"
    num_tasks = len(mgc_list)
    
    # Use the existing submit_array_job function
    return submit_array_job_rolling(job_name, script_path, stage_list_file, num_tasks, dependency_jobs)

# === Main execution: Submit array jobs for all MGCs ===
print(f"Processing {total_mgcs} MGC directories using array jobs...")

# Check current limits
submitted_jobs, max_submit_limit = get_job_limits()
print(f"Current job submission status: {submitted_jobs}/{max_submit_limit}")

# Wait for current jobs to clear up some space
print("Checking current job load...")
wait_for_job_slots(max_active_jobs=50, max_submitted_jobs=max_submit_limit - 20)

# Step 1: Submit homolog array job
print("\n=== STAGE 1: HOMOLOG PROCESSING ===")
job_ids1 = submit_array_job_for_stage('homolog', submit_find_homolog_genes_in_dir_array, 
                                       all_subdirs, dependency_jobs=None)
if job_ids1:
    wait_for_array_job(job_ids1)
elif all_subdirs:  # Only fail if there were jobs that should have been submitted
    print("Failed to submit homolog jobs. Exiting.")
    exit(1)

# Step 2: Submit cluster array job (depends on homolog completion)  
print("\n=== STAGE 2: CLUSTER PROCESSING ===")
job_ids2 = submit_array_job_for_stage('cluster', submit_find_clusters_in_chromosome_array,
                                       all_subdirs, dependency_jobs=job_ids1 if job_ids1 else None)
if job_ids2:
    wait_for_array_job(job_ids2)
elif all_subdirs:  # Only fail if there were jobs that should have been submitted
    print("Failed to submit cluster jobs. Exiting.")
    exit(1)

# Step 3: Submit multi array job (depends on cluster completion)
print("\n=== STAGE 3: MULTI PROCESSING ===")
job_ids3 = submit_array_job_for_stage('multi', submit_multichromosome_statistics_array,
                                       all_subdirs, dependency_jobs=job_ids2 if job_ids2 else None)
if job_ids3:
    wait_for_array_job(job_ids3)
elif all_subdirs:  # Only fail if there were jobs that should have been submitted
    print("Failed to submit multi jobs. Exiting.")
    exit(1)

print("All array jobs completed!")