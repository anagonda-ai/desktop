#!/usr/bin/env python3
"""
E2P2 Analysis Job Manager - Submit E2P2 analysis jobs for FASTA files to cluster
"""

import subprocess
import pandas as pd
from pathlib import Path
import os
import time

class E2P2AnalysisManager:
    def __init__(self):
        self.base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
        self.job_script_dir = self.base_dir / "random_kegg_e2p2_jobs"
        self.job_script_dir.mkdir(exist_ok=True)
        self.output_dir = self.base_dir / "random_kegg_e2p2_results"
        self.output_dir.mkdir(exist_ok=True)
        
    def get_running_job_count(self):
        """Get number of running E2P2 jobs"""
        try:
            # Check SLURM queue
            result = subprocess.run(
                ["squeue", "-u", os.environ.get("USER", "alongonda")],
                # Only count jobs with "e2p2" in the name
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=10
            )
            if result.returncode == 0:
                e2p2_jobs = [line for line in result.stdout.strip().split('\n') if "e2p2" in line]
                return len(e2p2_jobs) if result.stdout.strip() else 0
            
            
            # Fallback to PBS
            result = subprocess.run(
                ["qstat", "-u", os.environ.get("USER", "alongonda")],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=10
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                e2p2_jobs = [line for line in lines if "e2p2_tagging" in line and ("Q" in line or "R" in line)]
                return len(e2p2_jobs)
                
        except Exception:
            pass
            
        return 0
    
    def create_e2p2_job_script(self, fasta_file, fasta_id):
        """Create job script for E2P2 analysis"""
        
        job_name = f"e2p2_{fasta_id}"
        script_path = self.job_script_dir / f"{job_name}.sh"
        
        # Set up output paths
        output_subdir = self.output_dir / fasta_file.stem
        output_subdir.mkdir(parents=True, exist_ok=True)
        
        output_file = output_subdir / f"{fasta_file.stem}.pf"
        
        # Create E2P2 job script based on the template
        script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out_e2p2_{fasta_id}.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error_e2p2_{fasta_id}.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=7-00:00:00

echo "Starting E2P2 analysis for {fasta_file.name}..."
echo "Input file: {fasta_file}"
echo "Output file: {output_file}"

# Run E2P2 analysis
conda run -n deepec python /groups/itay_mayrose/alongonda/tools/E2P2/e2p2.py --input {fasta_file} --output {output_file} e2p2

# Check if E2P2 completed successfully
if [ $? -eq 0 ]; then
    echo "SUCCESS: E2P2 analysis completed for {fasta_file.name}"
    touch {output_subdir}/e2p2_complete.flag
    echo "Output saved to: {output_file}"
else
    echo "FAILED: E2P2 analysis failed for {fasta_file.name}"
    touch {output_subdir}/e2p2_failed.flag
    exit 1
fi

echo "E2P2 analysis for {fasta_file.name} finished successfully!"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_path, 0o755)
        return script_path
    
    def submit_e2p2_job(self, script_path):
        """Submit E2P2 job to scheduler"""
        try:
            # Try SLURM
            result = subprocess.run(
                ["sbatch", str(script_path)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=30
            )
            
            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                return job_id
            
            # Try PBS
            result = subprocess.run(
                ["qsub", str(script_path)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=30
            )
            
            if result.returncode == 0:
                return result.stdout.strip()
                
        except Exception as e:
            print(f"Error submitting job: {e}")
            
        return None
    
    def submit_e2p2_jobs_for_directory(self, fasta_dir):
        """Submit E2P2 jobs for all FASTA files in a directory"""
        dir_name = fasta_dir.name
        
        # Collect FASTA files
        fasta_files = list(fasta_dir.glob("*.fasta"))
        
        if not fasta_files:
            print(f"No FASTA files found in {dir_name}")
            return []
        
        print(f"Submitting E2P2 jobs for {dir_name}: {len(fasta_files)} FASTA files")
        
        submitted_jobs = []
        
        for i, fasta_file in enumerate(fasta_files):
            fasta_id = f"{dir_name}_{i}_{fasta_file.stem}"
            
            # Check if already processed
            output_subdir = self.output_dir / fasta_file.stem
            complete_flag = output_subdir / "e2p2_complete.flag"
            output_file = output_subdir / f"{fasta_file.stem}.pf"
            
            if complete_flag.exists() or output_file.exists():
                print(f"    Skipping {fasta_file.name} - already processed")
                continue
            
            # Wait if too many jobs running
            while self.get_running_job_count() >= 30:  # Limit concurrent E2P2 jobs
                print(f"  Waiting... {self.get_running_job_count()} jobs running")
                time.sleep(30)
            
            # Submit E2P2 job
            script_path = self.create_e2p2_job_script(fasta_file, fasta_id)
            job_id = self.submit_e2p2_job(script_path)
            
            if job_id:
                submitted_jobs.append({
                    'job_id': job_id,
                    'fasta_file': fasta_file,
                    'output_dir': output_subdir,
                    'fasta_id': fasta_id
                })
                print(f"    Submitted job {job_id}: {fasta_file.name}")
            else:
                print(f"    Failed to submit: {fasta_file.name}")
        
        return submitted_jobs
    
    def wait_for_all_e2p2_jobs(self, submitted_jobs):
        """Wait for all E2P2 jobs to complete"""
        if not submitted_jobs:
            return
        
        print(f"Waiting for {len(submitted_jobs)} E2P2 jobs to complete...")
        start_time = time.time()
        
        while True:
            completed = 0
            failed = 0
            
            for job_info in submitted_jobs:
                output_dir = job_info['output_dir']
                
                if (output_dir / "e2p2_complete.flag").exists():
                    completed += 1
                elif (output_dir / "e2p2_failed.flag").exists():
                    failed += 1
            
            pending = len(submitted_jobs) - completed - failed
            
            print(f"  Status: {completed} completed, {failed} failed, {pending} pending ({time.time() - start_time:.0f}s elapsed)")
            
            if pending == 0:
                break
                
            time.sleep(60)  # Check every minute
        
        print(f"All E2P2 jobs completed: {completed} successful, {failed} failed")
    
    def collect_e2p2_results(self, fasta_dir):
        """Collect E2P2 results for all FASTA files in a directory"""
        dir_name = fasta_dir.name
        
        # Collect FASTA files
        fasta_files = list(fasta_dir.glob("*.fasta"))
        
        if not fasta_files:
            print(f"No FASTA files found in {dir_name}")
            return None
        
        print(f"Collecting E2P2 results for {dir_name}")
        
        results = []
        successful = 0
        failed = 0
        
        for fasta_file in fasta_files:
            output_subdir = self.output_dir / fasta_file.stem
            complete_flag = output_subdir / "e2p2_complete.flag"
            failed_flag = output_subdir / "e2p2_failed.flag"
            output_file = output_subdir / f"{fasta_file.stem}.pf"
            
            status = "pending"
            
            if complete_flag.exists() and output_file.exists():
                status = "completed"
                successful += 1
                try:
                    # Get file size as a simple metric
                    file_size = output_file.stat().st_size
                    results.append({
                        'fasta_file': fasta_file.name,
                        'status': status,
                        'output_file': str(output_file),
                        'output_size': file_size
                    })
                    print(f"  ✓ {fasta_file.name}: {file_size} bytes")
                except Exception as e:
                    print(f"  ⚠ {fasta_file.name}: completed but error reading output - {e}")
                    
            elif failed_flag.exists():
                status = "failed"
                failed += 1
                results.append({
                    'fasta_file': fasta_file.name,
                    'status': status,
                    'output_file': None,
                    'output_size': 0
                })
                print(f"  ✗ {fasta_file.name}: failed")
            else:
                results.append({
                    'fasta_file': fasta_file.name,
                    'status': status,
                    'output_file': None,
                    'output_size': 0
                })
                print(f"  ? {fasta_file.name}: pending")
        
        # Save summary
        results_df = pd.DataFrame(results)
        summary_file = self.output_dir / f"{dir_name}_e2p2_summary.csv"
        results_df.to_csv(summary_file, index=False)
        
        print(f"  Summary: {successful} successful, {failed} failed, {len(fasta_files) - successful - failed} pending")
        print(f"  ✓ Saved summary to {summary_file}")
        return results_df
    
    
    def process_all_fasta_directories(self):
        """Process all FASTA directories with E2P2 job submission"""
        # Define the FASTA directories to process
        fasta_directories = [
            self.base_dir / "kegg_random_mgc_candidates_fasta_files"
        ]
        
        for fasta_dir in fasta_directories:
            if not fasta_dir.exists():
                print(f"Directory not found: {fasta_dir}")
                continue
            
            print(f"\n=== Processing FASTA directory: {fasta_dir.name} ===")
            
            # Phase 1: Submit E2P2 jobs
            submitted_jobs = self.submit_e2p2_jobs_for_directory(fasta_dir)
            
            # Phase 2: Wait for jobs if any were submitted
            if submitted_jobs:
                self.wait_for_all_e2p2_jobs(submitted_jobs)
            
            # Phase 3: Collect results
            self.collect_e2p2_results(fasta_dir)

def main():
    manager = E2P2AnalysisManager()
    
    print("E2P2 Analysis Job Manager")
    print("=" * 40)
    
    manager.process_all_fasta_directories()
    
    print("\nAll E2P2 analysis complete!")

if __name__ == "__main__":
    main()