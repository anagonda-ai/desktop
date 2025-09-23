import os
import subprocess

def analyze_tree_clades(
    tree_path,
    mgc_candidates_dir,
    mapping_file,
    compute_gain_loss_coevolution,
    max_workers=32,
    use_processes=False
):
    """
    Batch analyze tree clades for all MGC candidates.
    
    Parameters:
    - tree_path: Path to the phylogenetic tree file
    - mgc_candidates_dir: Directory containing MGC candidate subdirectories
    - mapping_file: Path to organism mapping file
    - compute_gain_loss_coevolution: Whether to compute gain/loss coevolution
    - max_workers: Number of concurrent workers for each analysis
    - use_processes: Whether to use processes instead of threads
    """
    
    # Path to the shell script
    shell_script_path = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_comparison_csv_tree_analyzer.sh"
    
    # Ensure the MGC candidates directory exists
    if not os.path.exists(mgc_candidates_dir):
        raise ValueError(f"MGC candidates directory does not exist: {mgc_candidates_dir}")
    
    # Find all subdirectories containing comparison.csv files
    subdirs_with_csv = []
    
    print(f"Scanning {mgc_candidates_dir} for subdirectories with comparison.csv...")
    
    for subdir in os.listdir(mgc_candidates_dir):
        subdir_path = os.path.join(mgc_candidates_dir, subdir)
        
        # Skip if not a directory
        if not os.path.isdir(subdir_path):
            continue
            
        # Look for comparison.csv in this subdirectory
        comparison_csv_path = os.path.join(subdir_path, "comparison_results.csv")
        
        if os.path.exists(comparison_csv_path):
            subdirs_with_csv.append((subdir, subdir_path, comparison_csv_path))
            print(f"  Found: {subdir}/comparison_results.csv")
        else:
            print(f"  Skipped: {subdir} (no comparison.csv found)")
    
    if not subdirs_with_csv:
        print("No subdirectories with comparison.csv found!")
        return
    
    print(f"\nFound {len(subdirs_with_csv)} MGC candidates to analyze")
    
    # Process each subdirectory
    successful_analyses = 0
    failed_analyses = 0
    submitted_jobs = []
    
    for i, (subdir_name, subdir_path, comparison_csv_path) in enumerate(subdirs_with_csv, 1):
        print(f"\n{'='*60}")
        print(f"Processing {i}/{len(subdirs_with_csv)}: {subdir_name}")
        print(f"{'='*60}")
        while True:
            try:
                result = subprocess.run(
                    ["squeue", "-u", "alongonda", "--noheader", "--format=%t"],
                    capture_output=True, text=True, check=True
                )
                
                if result.stdout.strip():
                    job_statuses = result.stdout.strip().split('\n')
                    pending_running_jobs = sum(1 for status in job_statuses if status.strip() in ['R', 'PD', 'CF', 'CG'])
                    
                    if pending_running_jobs >= 100:
                        print(f"  ⏳ {pending_running_jobs} jobs pending/running. Waiting before submitting more...")
                        import time
                        time.sleep(60)  # Wait 1 minute before checking again
                        continue
                    else:
                        print(f"  ✅ {pending_running_jobs} jobs pending/running. Safe to submit more.")
                        break
                else:
                    # No jobs running, safe to proceed
                    print(f"  ✅ No jobs pending/running. Safe to submit.")
                    break
                    
            except subprocess.CalledProcessError:
                print(f"  ⚠️ Could not check job queue status. Proceeding with submission.")
                break
        # Create output directory for this MGC candidate
        mgc_output_dir = os.path.join(mgc_candidates_dir, subdir_name)
        
        try:
            # Prepare arguments for the shell script
            submit_command = [
                "sbatch",
                shell_script_path,
                tree_path,
                comparison_csv_path,
                mgc_output_dir,
                mapping_file,
                str(compute_gain_loss_coevolution).lower(),
                str(max_workers),
                str(use_processes).lower()
            ]
            
            print(f"Submitting job with command: {' '.join(submit_command)}")
            
            # Submit the job
            result = subprocess.run(submit_command, capture_output=True, text=True, check=True)
            
            # Extract job ID from the output (assuming the shell script outputs it)
            job_output = result.stdout.strip()
            print(f"Job submission output: {job_output}")
            
            # Try to extract job ID (adjust this based on your shell script's output format)
            if "Submitted batch job" in job_output:
                job_id = job_output.split()[-1]
                submitted_jobs.append((subdir_name, job_id))
                print(f"✅ Successfully submitted job {job_id} for {subdir_name}")
            else:
                submitted_jobs.append((subdir_name, "unknown"))
                print(f"✅ Successfully submitted job for {subdir_name}")
            
            successful_analyses += 1
            
        except Exception as e:
            failed_analyses += 1
            print(f"❌ Failed analysis for {subdir_name}: {str(e)}")
            
            # Create error log
            error_log_path = os.path.join(mgc_output_dir, f"{subdir_name}_error.log")
            os.makedirs(os.path.dirname(error_log_path), exist_ok=True)
            
            with open(error_log_path, 'w') as f:
                f.write(f"Error analyzing {subdir_name}:\n")
                f.write(f"Error: {str(e)}\n")
                f.write(f"Comparison CSV: {comparison_csv_path}\n")
                f.write(f"Tree path: {tree_path}\n")
                f.write(f"Mapping file: {mapping_file}\n")
                
    # Wait for all jobs to complete
    if successful_analyses > 0:
        print(f"\n{'='*60}")
        print(f"WAITING FOR JOBS TO COMPLETE")
        print(f"{'='*60}")
        print(f"Submitted {successful_analyses} jobs. Waiting for completion...")
        
        # Extract valid job IDs
        valid_job_ids = [job_id for _, job_id in submitted_jobs if job_id != "unknown" and job_id.isdigit()]
        
        if valid_job_ids:
            print(f"Monitoring jobs: {', '.join(valid_job_ids)}")
            print(f"Monitor manually with: squeue -u $USER --job={','.join(valid_job_ids)}")
            
            # Wait for jobs to complete
            completed_jobs = 0
            failed_job_analyses = 0
            successful_job_analyses = 0
            
            while completed_jobs < len(valid_job_ids):
                import time
                time.sleep(30)  # Check every 30 seconds
                
                # Check job status using squeue
                try:
                    result = subprocess.run(
                        ["squeue", "-u", os.getenv("USER", ""), "--job=" + ",".join(valid_job_ids), "--noheader", "--format=%i,%t"],
                        capture_output=True, text=True, check=True
                    )
                    
                    running_jobs = []
                    if result.stdout.strip():
                        for line in result.stdout.strip().split('\n'):
                            if line.strip():
                                job_id, status = line.strip().split(',')
                                if status in ['R', 'PD', 'CF', 'CG']:  # Running, Pending, Configuring, Completing
                                    running_jobs.append(job_id)
                    
                    completed_this_round = len(valid_job_ids) - len(running_jobs)
                    if completed_this_round > completed_jobs:
                        newly_completed = completed_this_round - completed_jobs
                        completed_jobs = completed_this_round
                        print(f"  Progress: {completed_jobs}/{len(valid_job_ids)} jobs completed (+{newly_completed})")
                    
                except subprocess.CalledProcessError:
                    # If squeue fails, assume all jobs are done
                    print("  Could not check job status, assuming all jobs completed")
                    completed_jobs = len(valid_job_ids)
                    break
            
            print(f"✅ All {len(valid_job_ids)} jobs have completed!")
            
            # Check individual job results
            print(f"\nChecking individual job results...")
            for subdir_name, job_id in submitted_jobs:
                mgc_output_dir = os.path.join(mgc_candidates_dir, subdir_name)
                
                # Check if the analysis was successful by looking for output files
                summary_file = os.path.join(mgc_output_dir, "summary.csv")
                error_file = os.path.join(mgc_output_dir, f"{subdir_name}_error.log")
                
                if os.path.exists(summary_file):
                    successful_job_analyses += 1
                    print(f"  ✅ {subdir_name}: Analysis completed successfully")
                elif os.path.exists(error_file):
                    failed_job_analyses += 1
                    print(f"  ❌ {subdir_name}: Analysis failed (see error log)")
                else:
                    failed_job_analyses += 1
                    print(f"  ⚠️  {subdir_name}: Unknown status (no output files found)")
        else:
            print("No valid job IDs found for monitoring. Jobs may have been submitted but job IDs could not be extracted.")
            successful_job_analyses = 0
            failed_job_analyses = successful_analyses
    
    # Final summary
    print(f"\n{'='*60}")
    print(f"BATCH ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"Total MGC candidates processed: {len(subdirs_with_csv)}")
    print(f"Successful analyses: {successful_analyses}")
    print(f"Failed analyses: {failed_analyses}")
    
    if successful_analyses > 0:
        print(f"\nResults saved in: {mgc_candidates_dir}/")
        print("Each subdirectory contains:")
        print("  - *_summary.csv: Clade analysis results")
        print("  - *_matrix_npp_global.csv: Global normalized matrix")
        print("  - *_matrix_raw_global.csv: Global raw matrix")
        print("  - *_global_score.txt: Global CladePP score")
        print("  - *_clade_figures/: Individual clade heatmaps")
    
    if failed_analyses > 0:
        print(f"\nError logs for failed analyses saved as *_error.log files")
    
    return {
        'total': len(subdirs_with_csv),
        'successful': successful_analyses,
        'failed': failed_analyses,
        'output_dir': mgc_candidates_dir
    }
    
    
