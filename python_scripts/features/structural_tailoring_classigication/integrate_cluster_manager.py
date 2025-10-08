#!/usr/bin/env python3
"""
Cluster Integration Manager - Submit integration jobs to HPC cluster
Processes each cluster individually and creates separate CSV files
"""

import subprocess
import json
from pathlib import Path
import time
import os
import pandas as pd
import glob

class ClusterIntegrationManager:
    def __init__(self):
        self.base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
        self.job_script_dir = self.base_dir / "integration_jobs"
        self.job_script_dir.mkdir(exist_ok=True)
        self.output_dir = self.base_dir / "integrated_results"
        self.output_dir.mkdir(exist_ok=True)
        
        # Paths to input directories
        self.e2p2_dir = self.base_dir / "e2p2_results"
        self.domain_dir = self.base_dir / "domain_detection_results"
        self.mgc_fasta_dir = self.base_dir / "mgc_candidates_fasta_files"
        self.random_fasta_dir = self.base_dir / "random_mgc_candidates_fasta_files"
        
    def get_running_job_count(self):
        """Get number of running integration jobs"""
        try:
            result = subprocess.run(
                ["squeue", "-u", os.environ.get("USER", "alongonda"), "-h"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                return len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
        except Exception:
            pass
        return 0
    
    def find_all_clusters(self):
        """Find all clusters that have E2P2, domain, and FASTA files"""
        clusters = []
        
        # Find all E2P2 result files
        e2p2_files = list(self.e2p2_dir.glob("*/MGC_CANDIDATE_*.default.pf")) + \
                     list(self.e2p2_dir.glob("*/RANDOM_MGC_*.default.pf")) + \
                     list(self.e2p2_dir.glob("*/BGC*.default.pf"))
        
        for e2p2_file in e2p2_files:
            # Extract cluster name from file path
            cluster_name = e2p2_file.parent.name
            
            # Check if corresponding domain and FASTA files exist
            domain_file = self.domain_dir / f"{cluster_name}_domains.json"
            
            # Check both MGC and RANDOM MGC fasta directories
            fasta_file = None
            if cluster_name.startswith('MGC_CANDIDATE'):
                fasta_file = self.mgc_fasta_dir / f"{cluster_name}.fasta"
            elif cluster_name.startswith('RANDOM_MGC'):
                fasta_file = self.random_fasta_dir / f"{cluster_name}.fasta"
            elif cluster_name.startswith('BGC'):
                # BGC files might be in either directory, check both
                fasta_file = self.mgc_fasta_dir / f"{cluster_name}.fasta"
                if not fasta_file.exists():
                    fasta_file = self.random_fasta_dir / f"{cluster_name}.fasta"
            
            if domain_file.exists() and fasta_file and fasta_file.exists():
                clusters.append({
                    'cluster_name': cluster_name,
                    'e2p2_file': e2p2_file,
                    'domain_file': domain_file,
                    'fasta_file': fasta_file,
                    'output_file': self.output_dir / f"{cluster_name}_integrated_results.csv"
                })
        
        return clusters
    
    def create_integration_job_script(self, cluster_info, job_index, total_jobs):
        """Create SLURM job script for cluster integration"""
        
        cluster_name = cluster_info['cluster_name']
        job_name = f"integrate_{cluster_name}"
        script_path = self.job_script_dir / f"{job_name}.sh"
        
        # Create the job script content
        script_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            "#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/integrate_%j.out",
            "#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/integrate_%j.err",
            "#SBATCH --ntasks=1",
            "#SBATCH --time=2:00:00",
            "#SBATCH --mem=8G",
            "#SBATCH --partition=itaym",
            "#SBATCH --cpus-per-task=2",
            "",
            'echo "=========================================="',
            f'echo "Cluster Integration Job [{job_index}/{total_jobs}]"',
            f'echo "Cluster: {cluster_name}"',
            f'echo "Output: {cluster_info["output_file"]}"',
            'echo "=========================================="',
            "",
            "# Set up environment",
            "cd /groups/itay_mayrose/alongonda/desktop/python_scripts/features/structural_tailoring_classigication",
            "",
            "# Run integration for this specific cluster",
            f"python3 integrate_single_cluster.py {cluster_name} \"{cluster_info['e2p2_file']}\" \"{cluster_info['domain_file']}\" \"{cluster_info['fasta_file']}\" \"{cluster_info['output_file']}\"",
            "",
            "if [ $? -eq 0 ]; then",
            f'    echo "✓ Integration completed successfully for {cluster_name}"',
            "else",
            f'    echo "✗ Integration failed for {cluster_name}"',
            "    exit 1",
            "fi",
            "",
            'echo "=========================================="',
            f'echo "Job completed: {cluster_name}"',
            'echo "=========================================="'
        ]
        
        script_content = "\n".join(script_lines)
        
        # Write script file
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        os.chmod(script_path, 0o755)
        
        return script_path
    
    def submit_jobs(self, max_concurrent_jobs=100):
        """Submit integration jobs to SLURM cluster"""
        clusters = self.find_all_clusters()
        total_clusters = len(clusters)
        
        print(f"Found {total_clusters} clusters to process")
        
        if total_clusters == 0:
            print("No clusters found with all required files (E2P2, domain, FASTA)")
            return
        
        submitted_jobs = 0
        skipped_jobs = 0
        
        for i, cluster_info in enumerate(clusters, 1):
            cluster_name = cluster_info['cluster_name']
            output_file = cluster_info['output_file']
            
            # Check if output already exists
            if output_file.exists():
                print(f"✓ Skipping {cluster_name} - output already exists")
                skipped_jobs += 1
                continue
            
            # Check running job count
            while self.get_running_job_count() >= max_concurrent_jobs:
                print(f"Waiting for jobs to complete... ({self.get_running_job_count()}/{max_concurrent_jobs} running)")
                time.sleep(30)
            
            # Create and submit job
            script_path = self.create_integration_job_script(cluster_info, i, total_clusters)
            
            try:
                result = subprocess.run(
                    ["sbatch", str(script_path)],
                    capture_output=True, text=True, timeout=10
                )
                
                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1]
                    print(f"✓ Submitted job {job_id} for cluster {cluster_name}")
                    submitted_jobs += 1
                else:
                    print(f"✗ Failed to submit job for {cluster_name}: {result.stderr}")
                    
            except Exception as e:
                print(f"✗ Error submitting job for {cluster_name}: {e}")
            
            # Small delay between submissions
            time.sleep(1)
        
        print(f"\nJob submission complete:")
        print(f"  Total clusters: {total_clusters}")
        print(f"  Submitted jobs: {submitted_jobs}")
        print(f"  Skipped jobs: {skipped_jobs}")
        print(f"  Max concurrent: {max_concurrent_jobs}")

def main():
    """Main function"""
    manager = ClusterIntegrationManager()
    
    print("🧬 Cluster Integration Manager")
    print("=" * 50)
    
    # Check current job count
    running_jobs = manager.get_running_job_count()
    print(f"Currently running jobs: {running_jobs}")
    
    # Submit jobs
    manager.submit_jobs(max_concurrent_jobs=100)

if __name__ == "__main__":
    main()
