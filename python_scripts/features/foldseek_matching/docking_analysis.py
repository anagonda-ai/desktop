#!/usr/bin/env python3
"""
LightDock Full Workflow Job Manager - Submit complete workflow jobs to cluster
"""

import subprocess
import pandas as pd
from pathlib import Path
import numpy as np
import os
import time
import shutil

class LightDockFullWorkflowManager:
    def __init__(self):
        self.base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
        self.job_script_dir = self.base_dir / "lightdock_full_jobs"
        self.job_script_dir.mkdir(exist_ok=True)
        
    def get_running_job_count(self):
        """Get number of running LightDock jobs"""
        try:
            # Check SLURM queue
            result = subprocess.run(
                ["squeue", "-u", os.environ.get("USER", "alongonda"), "-n", "ld_full", "-h"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                return len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
            
            # Fallback to PBS
            result = subprocess.run(
                ["qstat", "-u", os.environ.get("USER", "alongonda")],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                full_jobs = [line for line in lines if "ld_full" in line and ("Q" in line or "R" in line)]
                return len(full_jobs)
                
        except Exception:
            pass
            
        return 0
    
    def create_full_workflow_job_script(self, receptor_pdb, ligand_pdb, output_dir, pair_id):
        """Create job script for complete LightDock workflow"""
        
        job_name = f"ld_full_{pair_id}"
        script_path = self.job_script_dir / f"{job_name}.sh"
        
        # Create local copies
        local_receptor = output_dir / f"receptor_{receptor_pdb.name}"
        local_ligand = output_dir / f"ligand_{ligand_pdb.name}"
        
        output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(receptor_pdb, local_receptor)
        shutil.copy2(ligand_pdb, local_ligand)
        
        # Complete workflow job script
        script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/full_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/full_%j.err
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=64G
#SBATCH --partition=itaym
#SBATCH --cpus-per-task=2

cd {output_dir}

# Clean any existing lightdock files
rm -f lightdock_*.pdb combined_gso_100.out clustered_lightdock.out rank_by_scoring.list
rm -rf swarm_* init

# Run setup
echo "Starting LightDock setup..."
lightdock3_setup.py {local_receptor} {local_ligand} --noxt -s 5 -g 10

# Check if setup was successful
if [ -f setup.json ]; then
    echo "SUCCESS: setup.json created"
    touch setup_complete.flag
    
    # Run lightdock simulation
    echo "Starting LightDock simulation..."
    lightdock3.py setup.json 100 -c 1 -s dfire
    
    # Check if simulation completed successfully
    if [ $? -eq 0 ]; then
        echo "SUCCESS: LightDock simulation completed"
        touch lightdock_complete.flag
        
        # Generate final structures for each swarm
        echo "Generating final structures..."
        
        # Generate structures for each swarm individually
        for swarm_dir in swarm_*/; do
            if [ -d "$swarm_dir" ]; then
                swarm_id=$(basename "$swarm_dir" | cut -d'_' -f2)
                echo "Processing $swarm_dir..."
                cd "$swarm_dir"
                lgd_generate_conformations.py {local_receptor} {local_ligand} gso_100.out 10
                cd ..
            fi
        done
        
        if [ $? -eq 0 ]; then
            echo "SUCCESS: Conformations generated"
            
            # Create a combined GSO file and rank directly (skip clustering)
            echo "Creating combined GSO file for direct ranking..."
            echo "# Combined GSO output from all swarms" > combined_gso_100.out
            
            # Add all poses from all swarms with proper formatting
            for swarm_dir in swarm_*/; do
                if [ -d "$swarm_dir" ] && [ -f "$swarm_dir/gso_100.out" ]; then
                    echo "Adding poses from $swarm_dir to combined file..."
                    # Simply append the GSO file content (skip header lines starting with #)
                    grep -v '^#' "$swarm_dir/gso_100.out" >> combined_gso_100.out
                fi
            done
            
            echo "Combined GSO file created with $(grep -v '^#' combined_gso_100.out | wc -l) poses"
            
            # Skip clustering - go directly to ranking
            echo "Ranking poses directly (skipping clustering)..."
            lgd_rank.py 5 100
            
            if [ $? -eq 0 ]; then
                echo "SUCCESS: Ranking completed"
                
                # Generate top conformations from ranking results
                echo "Generating top conformations..."
                if [ -f "rank_by_scoring.list" ]; then
                    lgd_top.py {local_receptor} {local_ligand} rank_by_scoring.list 10
                else
                    echo "ERROR: No ranking file found"
                    exit 1
                fi
                
                if [ $? -eq 0 ]; then
                    echo "SUCCESS: Top conformations generated"
                    touch analysis_complete.flag
                    
                    # Extract best score for summary
                    echo "Extracting best score..."
                    if [ -f "rank_by_scoring.list" ]; then
                        BEST_SCORE=$(tail -n +2 rank_by_scoring.list | head -n1 | awk '{{print $NF}}')
                        echo "Best docking score: $BEST_SCORE" > best_score.txt
                        echo "BEST_SCORE: $BEST_SCORE"
                    else
                        echo "ERROR: Cannot extract score - no ranking file"
                        exit 1
                    fi
                    
                    # Create summary report
                    echo "Creating summary..."
                    echo "Analysis Summary:" > analysis_summary.txt
                    echo "=================" >> analysis_summary.txt
                    echo "Setup: COMPLETE" >> analysis_summary.txt
                    echo "Simulation: COMPLETE" >> analysis_summary.txt
                    echo "Structure generation: COMPLETE" >> analysis_summary.txt
                    echo "Clustering: SKIPPED (Direct ranking used)" >> analysis_summary.txt
                    echo "Ranking: COMPLETE" >> analysis_summary.txt
                    echo "Top structures: COMPLETE" >> analysis_summary.txt
                    echo "" >> analysis_summary.txt
                    echo "Output files generated in: $(pwd)" >> analysis_summary.txt
                    echo "" >> analysis_summary.txt
                    cat best_score.txt >> analysis_summary.txt
                    
                else
                    echo "FAILED: Top conformations generation failed"
                    touch analysis_failed.flag
                    exit 1
                fi
            else
                echo "FAILED: Ranking failed"
                touch analysis_failed.flag
                exit 1
            fi
        else
            echo "FAILED: Conformation generation failed"
            touch analysis_failed.flag
            exit 1
        fi
        
    else
        echo "FAILED: LightDock simulation failed"
        touch lightdock_failed.flag
        exit 1
    fi
    
else
    echo "FAILED: setup.json not created"
    touch setup_failed.flag
    exit 1
fi

echo "Complete LightDock workflow finished successfully!"
echo "Check analysis_summary.txt for details"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_path, 0o755)
        return script_path
    
    def submit_full_workflow_job(self, script_path):
        """Submit full workflow job to scheduler"""
        try:
            # Try SLURM
            result = subprocess.run(
                ["sbatch", str(script_path)],
                capture_output=True, text=True, timeout=30
            )
            
            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                return job_id
            
            # Try PBS
            result = subprocess.run(
                ["qsub", str(script_path)],
                capture_output=True, text=True, timeout=30
            )
            
            if result.returncode == 0:
                return result.stdout.strip()
                
        except Exception as e:
            print(f"Error submitting job: {e}")
            
        return None
    
    def submit_all_workflows_for_cluster(self, cluster_dir):
        """Submit full workflow jobs for all pairs in a cluster"""
        cluster_name = cluster_dir.name
        
        # Collect PDB files
        pdb_files = list(cluster_dir.glob("*.pdb"))
        for af_dir in cluster_dir.glob("*_alphafold"):
            if af_dir.is_dir():
                af_pdbs = list(af_dir.glob("*.pdb"))
                if af_pdbs:
                    pdb_files.append(af_pdbs[0])
        
        if len(pdb_files) < 2:
            return []
        
        print(f"Submitting full workflow jobs for cluster {cluster_name}: {len(pdb_files)} proteins")
        
        output_base = self.base_dir / "lightdock_results" / cluster_name
        submitted_jobs = []
        
        for i, receptor in enumerate(pdb_files):
            for j, ligand in enumerate(pdb_files):
                if i == j:  # Skip self comparisons
                    continue
                
                output_dir = output_base / f"{receptor.stem}_vs_{ligand.stem}"
                pair_id = f"{cluster_name}_{i}_{j}"
                
                # Check if already done
                analysis_complete = output_dir / "analysis_complete.flag"
                analysis_summary = output_dir / "analysis_summary.txt"
                
                if analysis_complete.exists() or analysis_summary.exists():
                    continue  # Already completed
                
                # Wait if too many jobs running
                while self.get_running_job_count() >= 50:  # Reduced concurrent jobs for full workflow
                    print(f"  Waiting... {self.get_running_job_count()} jobs running")
                    time.sleep(30)
                
                # Submit full workflow job
                script_path = self.create_full_workflow_job_script(receptor, ligand, output_dir, pair_id)
                job_id = self.submit_full_workflow_job(script_path)
                
                if job_id:
                    submitted_jobs.append({
                        'job_id': job_id,
                        'output_dir': output_dir,
                        'receptor': receptor,
                        'ligand': ligand,
                        'pair_id': pair_id
                    })
                    print(f"    Submitted job {job_id}: {receptor.stem} vs {ligand.stem}")
                else:
                    print(f"    Failed to submit: {receptor.stem} vs {ligand.stem}")
        
        return submitted_jobs
    
    def wait_for_all_workflows(self, submitted_jobs):
        """Wait for all workflow jobs to complete"""
        if not submitted_jobs:
            return
        
        print(f"Waiting for {len(submitted_jobs)} workflow jobs to complete...")
        start_time = time.time()
        
        while True:
            completed = 0
            failed = 0
            
            for job_info in submitted_jobs:
                output_dir = job_info['output_dir']
                
                if (output_dir / "analysis_complete.flag").exists():
                    completed += 1
                elif any((output_dir / f"{flag}_failed.flag").exists() 
                        for flag in ["setup", "lightdock", "analysis"]):
                    failed += 1
            
            pending = len(submitted_jobs) - completed - failed
            
            print(f"  Status: {completed} completed, {failed} failed, {pending} pending ({time.time() - start_time:.0f}s elapsed)")
            
            if pending == 0:
                break
                
            time.sleep(60)  # Check every minute
        
        print(f"All workflow jobs completed: {completed} successful, {failed} failed")
    
    def collect_results_after_workflows(self, cluster_dir):
        """Collect results after all workflow jobs complete"""
        cluster_name = cluster_dir.name
        
        # Collect PDB files
        pdb_files = list(cluster_dir.glob("*.pdb"))
        for af_dir in cluster_dir.glob("*_alphafold"):
            if af_dir.is_dir():
                af_pdbs = list(af_dir.glob("*.pdb"))
                if af_pdbs:
                    pdb_files.append(af_pdbs[0])
        
        if len(pdb_files) < 2:
            return None
        
        print(f"Collecting results for cluster {cluster_name}")
        
        output_base = self.base_dir / "lightdock_results" / cluster_name
        results = []
        
        for i, receptor in enumerate(pdb_files):
            for j, ligand in enumerate(pdb_files):
                if i == j:
                    # Self comparison
                    results.append([
                        receptor.stem, ligand.stem,
                        1.0, 0, 0, 0, 0, 0, 0, 0,
                        -999.0, 0, 1.0, 1.0, 1.0
                    ])
                    continue
                
                output_dir = output_base / f"{receptor.stem}_vs_{ligand.stem}"
                
                # Check if analysis was successful and extract score
                analysis_complete = output_dir / "analysis_complete.flag"
                best_score_file = output_dir / "best_score.txt"
                
                if analysis_complete.exists() and best_score_file.exists():
                    try:
                        with open(best_score_file, 'r') as f:
                            content = f.read()
                            # Extract score from "Best docking score: X.X" format
                            if "Best docking score:" in content:
                                score_line = [line for line in content.split('\n') if "Best docking score:" in line][0]
                                score = float(score_line.split(':')[1].strip())
                            else:
                                score = 0.0
                    except:
                        score = 0.0
                else:
                    score = 0.0
                    print(f"  No results for {receptor.stem} vs {ligand.stem}")
                
                results.append(self.create_result_row(receptor, ligand, score))
                print(f"  {receptor.stem} vs {ligand.stem}: {score:.3f}")
        
        # Save results
        df = pd.DataFrame(results, columns=[
            'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'lddt', 'tm', 'qtm'
        ])
        
        output_tsv = output_base.parent / "results_tsv" / f"{cluster_name}_all_vs_all.tsv"
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_tsv, sep='\t', index=False, header=False)
        
        print(f"  ✓ Saved results to {output_tsv}")
        return df
    
    def create_result_row(self, receptor, ligand, score):
        """Convert score to result row"""
        if score == 0.0:
            pseudo_tm = 0.1
        elif score < -100:
            pseudo_tm = 0.9
        elif score < -50:
            pseudo_tm = 0.7
        elif score < -25:
            pseudo_tm = 0.5
        elif score < 0:
            pseudo_tm = 0.3
        else:
            pseudo_tm = 0.1
        
        return [
            receptor.stem, ligand.stem,
            pseudo_tm, 0, 0, 0, 0, 0, 0, 0,
            score, 0, pseudo_tm, 0, 0
        ]
    
    def process_all_clusters(self, category_prefix):
        """Process all clusters with full workflow job submission"""
        if "RANDOM" in category_prefix:
            search_dir = self.base_dir / "random_mgc_pdb_files"
        else:
            search_dir = self.base_dir / "mgc_pdb_files"
        
        if not search_dir.exists():
            print(f"Directory not found: {search_dir}")
            return
        
        cluster_dirs = [d for d in search_dir.iterdir() 
                       if d.is_dir() and d.name.startswith(category_prefix)]
        
        print(f"Processing {len(cluster_dirs)} {category_prefix} clusters")
        
        for i, cluster_dir in enumerate(cluster_dirs):
            print(f"\n=== Cluster {i+1}/{len(cluster_dirs)}: {cluster_dir.name} ===")
            
            # Phase 1: Submit full workflow jobs
            self.submit_all_workflows_for_cluster(cluster_dir)
            
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
            
        # Collect results
        for i, cluster_dir in enumerate(cluster_dirs):
            print(f"\n=== Cluster {i+1}/{len(cluster_dirs)}: {cluster_dir.name} ===")
            try:
                result = self.collect_results_after_workflows(cluster_dir)
                if result is not None:
                    scores = [r[10] for r in result.values.tolist() if r[10] not in [0.0, -999.0]]
                    success_rate = len(scores) / len(result) if len(result) > 0 else 0
                    print(f"✓ {cluster_dir.name}: Success rate {success_rate:.1%}")
                    if scores:
                        print(f"  Score range: {min(scores):.2f} to {max(scores):.2f}")
            except Exception as e:
                print(f"✗ {cluster_dir.name}: Error - {e}")

def main():
    manager = LightDockFullWorkflowManager()
    
    print("LightDock Full Workflow Job Manager")
    print("=" * 40)
    
    for category in ['MGC_CANDIDATE', 'RANDOM', 'BGC']:
        print(f"\n=== Processing {category} clusters ===")
        manager.process_all_clusters(category)
        print(f"Completed {category}")
    
    print("\nAll analysis complete!")

if __name__ == "__main__":
    main()