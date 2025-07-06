# Query FASTA files
import os
import subprocess
import time


root_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g3_slurm/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_fixed"

# Paths to the submission scripts
submit_find_homolog_genes_in_dir = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_homolog_genes_in_dir.sh"
submit_find_clusters_in_chromosome = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_clusters_in_chromosome.sh"
sh_script_path = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_multichromosome_statistics.sh"

all_subdirs = [
    os.path.join(root_dir, d)
    for d in os.listdir(root_dir)
    if os.path.isdir(os.path.join(root_dir, d))
]

# Submit jobs to find homolog genes in each MGC directory
for example_mgc in all_subdirs:
    while True:
        squeue_cmd = ["squeue", "-u", "alongonda"]
        squeue_result = subprocess.run(squeue_cmd, capture_output=True, text=True)
        running_jobs = len(squeue_result.stdout.strip().splitlines())
        if running_jobs < 100:
            break
        print(f"⏳ Waiting: {running_jobs} jobs running for 'alongonda'. Sleeping 30s...")
        time.sleep(30)
    out_path = os.path.join("/groups/itay_mayrose/alongonda/desktop/example_jobs",f"{os.path.basename(example_mgc)+'homolog'}.out")
    error_path = os.path.join("/groups/itay_mayrose/alongonda/desktop/example_jobs",f"{os.path.basename(example_mgc)+'homolog'}.err")
    if os.path.exists(out_path):
        print(f"⚠️ Skipped existing output file: {out_path}")
        continue
    sbatch_cmd = [
        "sbatch",
        "--job-name", os.path.basename(example_mgc)+"homolog",
        "--output", out_path,
        "--error", error_path,
        submit_find_homolog_genes_in_dir,
        example_mgc  # Pass the MGC directory as an argument to the script
    ]
    result = subprocess.run(sbatch_cmd, capture_output=True, text=True)
    print(f"Processing MGC directory: {example_mgc}")

# Wait for all jobs to finish before proceeding
while True:
    squeue_cmd = ["squeue", "-u", "alongonda"]
    squeue_result = subprocess.run(squeue_cmd, capture_output=True, text=True)
    running_jobs = len(squeue_result.stdout.strip().splitlines())
    # Check that no "homolog" jobs are running
    if all("homolog" not in line for line in squeue_result.stdout.strip().splitlines()):
        print("✅ All 'homolog' jobs completed.")
        break
    print(f"⏳ Waiting: {running_jobs} jobs running for 'alongonda'. Sleeping 30s...")
    time.sleep(30)

# Submit jobs to find clusters in each chromosome
for example_mgc in all_subdirs:
    while True:
        squeue_cmd = ["squeue", "-u", "alongonda"]
        squeue_result = subprocess.run(squeue_cmd, capture_output=True, text=True)
        running_jobs = len(squeue_result.stdout.strip().splitlines())
        if running_jobs < 100:
            break
        print(f"⏳ Waiting: {running_jobs} jobs running for 'alongonda'. Sleeping 30s...")
        time.sleep(30)
    out_path = os.path.join("/groups/itay_mayrose/alongonda/desktop/example_jobs",f"{os.path.basename(example_mgc)+'cluster'}.out")
    error_path = os.path.join("/groups/itay_mayrose/alongonda/desktop/example_jobs",f"{os.path.basename(example_mgc)+'cluster'}.err")
    if os.path.exists(out_path):
        print(f"⚠️ Skipped existing output file: {out_path}")
        continue
    sbatch_cmd = [
        "sbatch",
        "--job-name", os.path.basename(example_mgc)+"cluster",
        "--output", out_path,
        "--error", error_path,
        submit_find_clusters_in_chromosome,
        example_mgc  # Pass the MGC directory as an argument to the script
    ]
    
    result = subprocess.run(sbatch_cmd, capture_output=True, text=True)
    print(f"Processing MGC directory: {example_mgc}")
    
# Wait for all jobs to finish before proceeding
while True:
    squeue_cmd = ["squeue", "-u", "alongonda"]
    squeue_result = subprocess.run(squeue_cmd, capture_output=True, text=True)
    running_jobs = len(squeue_result.stdout.strip().splitlines())
    # Check that no "cluster" jobs are running
    if all("cluster" not in line for line in squeue_result.stdout.strip().splitlines()):
        print("✅ All 'cluster' jobs completed.")
        break
    print(f"⏳ Waiting: {running_jobs} jobs running for 'alongonda'. Sleeping 30s...")
    time.sleep(30)
    
# Submit jobs to find multi-chromosome cluster statistics
for example_mgc in all_subdirs:
    while True:
        squeue_cmd = ["squeue", "-u", "alongonda"]
        squeue_result = subprocess.run(squeue_cmd, capture_output=True, text=True)
        running_jobs = len(squeue_result.stdout.strip().splitlines())
        if running_jobs < 100:
            break
        print(f"⏳ Waiting: {running_jobs} jobs running for 'alongonda'. Sleeping 30s...")
        time.sleep(30)
    out_path = os.path.join("/groups/itay_mayrose/alongonda/desktop/example_jobs",f"{os.path.basename(example_mgc)+'multi'}.out")
    error_path = os.path.join("/groups/itay_mayrose/alongonda/desktop/example_jobs",f"{os.path.basename(example_mgc)+'multi'}.err")
    if os.path.exists(out_path):
        print(f"⚠️ Skipped existing output file: {out_path}")
        continue
    sbatch_cmd = [
        "sbatch",
        "--job-name", os.path.basename(example_mgc)+"multi",
        "--output", out_path,
        "--error", error_path,
        sh_script_path,
        example_mgc  # Pass the MGC directory as an argument to the script
    ]
    result = subprocess.run(sbatch_cmd, capture_output=True, text=True)
    print(f"Processing MGC directory: {example_mgc}")