#!/usr/bin/env python3
"""
Domain Detection Cluster Job Manager - Submit domain detection jobs to HPC cluster
"""

import subprocess
import json
from pathlib import Path
import time
import os
import pandas as pd

class DomainDetectionClusterManager:
    def __init__(self):
        self.base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
        self.job_script_dir = self.base_dir / "domain_detection_jobs"
        self.job_script_dir.mkdir(exist_ok=True)
        
        # Tool paths
        self.hmmer_path = "hmmscan"
        self.interproscan_path = "interproscan.sh"
        self.pfam_db_path = "/groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/plantismash/antismash/generic_modules/fullhmmer/Pfam-A.hmm"
        
    def get_running_job_count(self):
        """Get number of running domain detection jobs"""
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
    
    def create_domain_detection_job_script(self, fasta_file, output_dir, file_index, total_files):
        """Create SLURM job script for domain detection"""
        
        job_name = f"domain_det_{fasta_file.stem}"
        script_path = self.job_script_dir / f"{job_name}.sh"
        
        output_json = output_dir / f"{fasta_file.stem}_domains.json"
        temp_dir = output_dir / f"temp_{fasta_file.stem}"
        
        script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/domain_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/domain_%j.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --cpus-per-task=4

echo "=========================================="
echo "Domain Detection Job [{file_index}/{total_files}]"
echo "FASTA file: {fasta_file.name}"
echo "Output: {output_json}"
echo "=========================================="

# Create temp directory
mkdir -p {temp_dir}
cd {temp_dir}

# Copy FASTA file locally
cp {fasta_file} ./input.fasta

echo "Starting Pfam analysis with HMMER..."

# Run HMMER with Pfam database
{self.hmmer_path} --domtblout pfam_output.txt --cut_ga --cpu 4 {self.pfam_db_path} ./input.fasta

if [ $? -eq 0 ]; then
    echo "✓ HMMER completed successfully"
    
    # Parse HMMER output and create JSON
    python3 << 'PYTHON_SCRIPT'
import json
from pathlib import Path

domain_hits = []

try:
    with open('pfam_output.txt', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split()
            if len(parts) >= 22:
                e_value = float(parts[6])
                score = float(parts[7])
                
                # Calculate confidence
                if e_value < 1e-10 and score > 50:
                    confidence = "high"
                elif e_value < 1e-5 and score > 30:
                    confidence = "medium"
                else:
                    confidence = "low"
                
                domain_hit = {{
                    'protein_id': parts[3],
                    'domain_id': parts[0],
                    'domain_name': parts[0],
                    'domain_description': ' '.join(parts[22:]) if len(parts) > 22 else '',
                    'start': int(parts[19]),
                    'end': int(parts[20]),
                    'e_value': e_value,
                    'score': score,
                    'database': 'Pfam',
                    'confidence': confidence
                }}
                domain_hits.append(domain_hit)
    
    # Extract unique domain names
    domain_names = list(set([hit['domain_name'] for hit in domain_hits]))
    
    # Calculate overall confidence
    if domain_hits:
        confidence_weights = {{'high': 1.0, 'medium': 0.7, 'low': 0.3}}
        weighted_scores = []
        for hit in domain_hits:
            weight = confidence_weights.get(hit['confidence'], 0.3)
            e_score = max(0, 1 - (hit['e_value'] / 1e-10) if hit['e_value'] > 0 else 1)
            weighted_scores.append(weight * e_score)
        confidence_score = sum(weighted_scores) / len(weighted_scores)
    else:
        confidence_score = 0.0
    
    # Create results dictionary
    results = {{
        'fasta_file': '{fasta_file}',
        'domain_hits': domain_hits,
        'all_domains': domain_names,
        'confidence_score': confidence_score,
        'analysis_methods': ['Pfam']
    }}
    
    # Save to JSON
    with open('{output_json}', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"✓ Found {{len(domain_hits)}} domain hits")
    print(f"✓ Confidence score: {{confidence_score:.3f}}")
    print(f"✓ Results saved to {output_json}")
    
except Exception as e:
    print(f"✗ Error parsing results: {{e}}")
    # Create empty results file
    results = {{
        'fasta_file': '{fasta_file}',
        'domain_hits': [],
        'all_domains': [],
        'confidence_score': 0.0,
        'analysis_methods': ['Pfam'],
        'error': str(e)
    }}
    with open('{output_json}', 'w') as f:
        json.dump(results, f, indent=2)

PYTHON_SCRIPT

    if [ -f "{output_json}" ]; then
        echo "✓ Job completed successfully"
        touch {output_dir}/{fasta_file.stem}_complete.flag
    else
        echo "✗ Failed to create output JSON"
        exit 1
    fi
    
else
    echo "✗ HMMER failed"
    exit 1
fi

# Cleanup temp directory
cd {output_dir}
rm -rf {temp_dir}

echo "=========================================="
echo "Domain detection complete for {fasta_file.name}"
echo "=========================================="
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_path, 0o755)
        return script_path
    
    def submit_job(self, script_path):
        """Submit job to SLURM scheduler"""
        try:
            result = subprocess.run(
                ["sbatch", str(script_path)],
                capture_output=True, text=True, timeout=30
            )
            
            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                return job_id
        except Exception as e:
            print(f"Error submitting job: {e}")
        
        return None
    
    def process_all_fasta_files(self):
        """Process all FASTA files by submitting cluster jobs"""
        
        output_dir = self.base_dir / "domain_detection_results"
        output_dir.mkdir(exist_ok=True)
        
        fasta_directories = [
            self.base_dir / "mgc_candidates_fasta_files",
            self.base_dir / "random_mgc_candidates_fasta_files"
        ]
        
        all_fasta_files = []
        
        # Collect all FASTA files
        for fasta_dir in fasta_directories:
            if not fasta_dir.exists():
                print(f"Directory not found: {fasta_dir}")
                continue
            
            print(f"\n=== Scanning FASTA directory: {fasta_dir.name} ===")
            fasta_files = list(fasta_dir.glob("*.fasta"))
            
            if not fasta_files:
                print(f"No FASTA files found in {fasta_dir.name}")
                continue
            
            print(f"Found {len(fasta_files)} FASTA files")
            all_fasta_files.extend(fasta_files)
        
        if not all_fasta_files:
            print("No FASTA files found in any directory!")
            return
        
        print(f"\n=== Submitting {len(all_fasta_files)} domain detection jobs ===")
        
        submitted_jobs = []
        
        for i, fasta_file in enumerate(all_fasta_files):
            file_index = i + 1
            
            # Check if already completed
            output_json = output_dir / f"{fasta_file.stem}_domains.json"
            complete_flag = output_dir / f"{fasta_file.stem}_complete.flag"
            
            if output_json.exists() or complete_flag.exists():
                print(f"[{file_index}/{len(all_fasta_files)}] Skipping {fasta_file.name} (already completed)")
                continue
            
            # Wait if too many jobs running
            while self.get_running_job_count() >= 100:
                running = self.get_running_job_count()
                print(f"  Waiting... {running} jobs running")
                time.sleep(30)
            
            # Create and submit job
            script_path = self.create_domain_detection_job_script(
                fasta_file, output_dir, file_index, len(all_fasta_files)
            )
            job_id = self.submit_job(script_path)
            
            if job_id:
                submitted_jobs.append({
                    'job_id': job_id,
                    'fasta_file': fasta_file,
                    'output_json': output_json
                })
                print(f"[{file_index}/{len(all_fasta_files)}] Submitted job {job_id}: {fasta_file.name}")
            else:
                print(f"[{file_index}/{len(all_fasta_files)}] Failed to submit: {fasta_file.name}")
        
        print(f"\n✓ Submitted {len(submitted_jobs)} jobs to cluster")
        return submitted_jobs, output_dir
    
    def wait_for_completion(self, submitted_jobs):
        """Wait for all jobs to complete"""
        if not submitted_jobs:
            print("No jobs to wait for")
            return
        
        print(f"\nWaiting for {len(submitted_jobs)} jobs to complete...")
        start_time = time.time()
        
        while True:
            completed = sum(1 for job in submitted_jobs if job['output_json'].exists())
            pending = len(submitted_jobs) - completed
            
            elapsed = time.time() - start_time
            print(f"  Status: {completed}/{len(submitted_jobs)} completed, {pending} pending ({elapsed:.0f}s elapsed)")
            
            if pending == 0:
                break
            
            time.sleep(60)
        
        print(f"\n✓ All jobs completed in {time.time() - start_time:.0f} seconds")
    
    def create_summary(self, output_dir):
        """Create summary CSV from all results"""
        print("\nCreating summary...")
        
        json_files = list(output_dir.glob("*_domains.json"))
        
        if not json_files:
            print("No result files found!")
            return
        
        all_results = []
        
        for json_file in json_files:
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                
                result_info = {
                    'fasta_file': data.get('fasta_file', ''),
                    'fasta_name': Path(data.get('fasta_file', '')).name,
                    'fasta_stem': Path(data.get('fasta_file', '')).stem,
                    'domain_hits': len(data.get('domain_hits', [])),
                    'domains': ', '.join(data.get('all_domains', [])),
                    'confidence_score': data.get('confidence_score', 0.0),
                    'analysis_methods': ', '.join(data.get('analysis_methods', [])),
                    'output_file': str(json_file),
                    'success': True
                }
                all_results.append(result_info)
            except Exception as e:
                print(f"  Error reading {json_file.name}: {e}")
        
        if all_results:
            df = pd.DataFrame(all_results)
            summary_file = output_dir / "domain_detection_summary.csv"
            df.to_csv(summary_file, index=False)
            print(f"✓ Summary saved: {summary_file}")
            print(f"✓ Total files processed: {len(all_results)}")
            
            # Print statistics
            total_domains = df['domain_hits'].sum()
            avg_domains = df['domain_hits'].mean()
            avg_confidence = df['confidence_score'].mean()
            
            print(f"\nStatistics:")
            print(f"  Total domain hits: {total_domains}")
            print(f"  Average domains per file: {avg_domains:.1f}")
            print(f"  Average confidence: {avg_confidence:.3f}")

def main():
    manager = DomainDetectionClusterManager()
    
    print("Domain Detection Cluster Job Manager")
    print("=" * 40)
    
    # Submit all jobs
    submitted_jobs, output_dir = manager.process_all_fasta_files()
    
    # Wait for completion
    manager.wait_for_completion(submitted_jobs)
    
    # Create summary
    manager.create_summary(output_dir)
    
    print("\nAll domain detection complete!")

if __name__ == "__main__":
    main()