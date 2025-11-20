import os
import subprocess
import time

# Config
DIRS_ROOT = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/"
SCRIPT = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/cladepp_phylo_profiling/build_npp_normalization_database/submit_blast_single_kegg.sh"
MAX_RUNNING = 100
USER = "alongonda"
PARTITION = "itaym-pool"

def get_immediate_subdirs(root_dir):
    return [
        os.path.join(root_dir, f)
        for f in os.listdir(root_dir)
        if os.path.isdir(os.path.join(root_dir, f))
    ]

def count_running_jobs(user, partition):
    try:
        # Use squeue to count only this user and partition, jobs with our script
        result = subprocess.run(
            [
                "squeue",
                "-u", user,
                "-p", partition,
                "-o", "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
            ],
            capture_output=True,
            text=True,
        )
        # Exclude the header
        if result.returncode != 0:
            print("Warning: squeue failed.", result.stderr)
            return 0
        lines = result.stdout.strip().split('\n')
        jobs = [line for line in lines[1:] if line.strip()]
        return len(jobs)
    except Exception as e:
        print(f"Failed to check running jobs: {e}")
        return 0

def main():
    dirs = get_immediate_subdirs(DIRS_ROOT)
    print(f"Found {len(dirs)} directories.")

    for dir_path in dirs:
        # Wait until there are fewer than MAX_RUNNING relevant jobs
        while True:
            running = count_running_jobs(USER, PARTITION)
            if running < MAX_RUNNING:
                break
            print(f"{running} jobs running on {PARTITION} (max {MAX_RUNNING}), waiting...")
            time.sleep(30)

        print(f"Submitting: {dir_path}")
        # Submit job
        result = subprocess.run(
            ["sbatch", SCRIPT, dir_path],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print(f"Submitted for {dir_path}: {result.stdout.strip()}")
        else:
            print(f"Failed to submit for {dir_path}: {result.stderr.strip()}")

if __name__ == "__main__":
    main()
