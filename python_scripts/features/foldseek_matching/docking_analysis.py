#!/usr/bin/env python3
"""
Production LightDock pipeline - now that LightDock is working
"""

import subprocess
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import json
import os
import shutil

class LightDockPipeline:
    def __init__(self):
        self.base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
        
    def run_lightdock_pair(self, receptor_pdb, ligand_pdb, output_dir):
        """Run LightDock for a protein pair with better file handling"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Clear any existing results
        for old_file in output_dir.glob("*"):
            if old_file.is_file():
                old_file.unlink()
            elif old_file.is_dir():
                shutil.rmtree(old_file)
        
        # Create copies of PDB files in output directory to avoid conflicts
        local_receptor = output_dir / f"receptor_{receptor_pdb.name}"
        local_ligand = output_dir / f"ligand_{ligand_pdb.name}"
        
        try:
            # Copy PDB files locally
            shutil.copy2(receptor_pdb, local_receptor)
            shutil.copy2(ligand_pdb, local_ligand)
            
            original_cwd = Path.cwd()
            os.chdir(output_dir)
            
            # Step 1: Setup docking with local files
            setup_cmd = [
                "lightdock3_setup.py",
                str(local_receptor),
                str(local_ligand),
                "--noxt",
                "-s", "5",   # 5 swarms for speed
                "-g", "10"   # 10 glowworms per swarm
            ]
            
            setup_result = subprocess.run(
                setup_cmd, 
                capture_output=True, 
                text=True,
                timeout=60
            )
            
            if setup_result.returncode != 0:
                # Check for specific errors
                if "already exists" in setup_result.stderr:
                    # Clean up lightdock files and retry
                    for ldf in output_dir.glob("lightdock_*.pdb"):
                        ldf.unlink()
                    
                    setup_result = subprocess.run(
                        setup_cmd, 
                        capture_output=True, 
                        text=True,
                        timeout=60
                    )
                    
                    if setup_result.returncode != 0:
                        return None
                else:
                    return None
            
            setup_json = output_dir / "setup.json"
            if not setup_json.exists():
                return None
            
            # Step 2: Run docking simulation  
            dock_cmd = [
                "lightdock3.py",
                "setup.json",
                "30",  # 30 steps for reasonable accuracy
                "-c", "1",
                "-s", "dfire"
            ]
            
            dock_result = subprocess.run(
                dock_cmd, 
                capture_output=True, 
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if dock_result.returncode != 0:
                return None
            
            # Step 3: Generate conformations
            num_swarms = 5
            generated_any = False
            
            for swarm_id in range(num_swarms):
                swarm_dir = output_dir / f"swarm_{swarm_id}"
                if swarm_dir.exists():
                    try:
                        os.chdir(swarm_dir)
                        
                        gen_cmd = [
                            "lgd_generate_conformations.py",
                            str(local_receptor.resolve()),
                            str(local_ligand.resolve()),
                            "gso_30.out",
                            "15"  # Generate 15 conformations
                        ]
                        
                        gen_result = subprocess.run(
                            gen_cmd,
                            capture_output=True,
                            text=True,
                            timeout=60
                        )
                        
                        if gen_result.returncode == 0:
                            generated_any = True
                            
                    except Exception:
                        pass
                    finally:
                        os.chdir(output_dir)
            
            # Step 4: Rank solutions
            if generated_any:
                rank_cmd = [
                    "lgd_rank.py",
                    "30",
                    str(num_swarms)
                ]
                
                subprocess.run(
                    rank_cmd, 
                    capture_output=True, 
                    text=True,
                    timeout=60
                )
            
            # Step 5: Parse the best score
            score = self.parse_lightdock_score(output_dir)
            return score
            
        except subprocess.TimeoutExpired:
            return 888.0  # Timeout score
        except Exception as e:
            return 777.0  # Error score
        finally:
            try:
                os.chdir(original_cwd)
            except:
                pass
    
    def parse_lightdock_score(self, output_dir):
        """Parse LightDock output for the best score"""
        
        # Method 1: Try ranking file first (preferred)
        rank_file = output_dir / "rank_by_scoring.list"
        if rank_file.exists():
            try:
                with open(rank_file) as f:
                    first_line = f.readline().strip()
                    if first_line:
                        parts = first_line.split()
                        if len(parts) >= 2:
                            return float(parts[1])
            except Exception:
                pass
        
        # Method 2: Parse GSO files
        best_score = None
        swarm_dirs = list(output_dir.glob("swarm_*"))
        
        for swarm_dir in swarm_dirs:
            # Try different step numbers
            for step in [30, 20, 10]:
                gso_file = swarm_dir / f"gso_{step}.out"
                if gso_file.exists():
                    try:
                        with open(gso_file) as f:
                            for line in f:
                                line = line.strip()
                                if line and not line.startswith('#'):
                                    parts = line.split()
                                    if len(parts) >= 7:
                                        try:
                                            score = float(parts[-1])  # Last column is score
                                            if best_score is None or score < best_score:
                                                best_score = score
                                        except ValueError:
                                            continue
                    except Exception:
                        continue
                    break  # Found a gso file for this swarm
        
        return best_score if best_score is not None else 0.0
    
    def process_protein_pair(self, i, j, receptor, ligand, output_base, cluster_name):
        """Process a single protein pair"""
        
        if i == j:
            # Self comparison
            return [
                receptor.stem, ligand.stem,
                1.0, 0, 0, 0, 0, 0, 0, 0,
                -999.0, 0, 1.0, 1.0, 1.0
            ]
        
        output_dir = output_base / f"{receptor.stem}_vs_{ligand.stem}"
        
        # Check cache
        score_file = output_dir / "lightdock_score.txt"
        if score_file.exists():
            try:
                cached_score = float(score_file.read_text().strip())
                # Use cache if it's a reasonable score
                if cached_score not in [0.0, 777.0, 888.0]:
                    return self.create_result_row(receptor, ligand, cached_score)
            except:
                pass
        
        # Run fresh LightDock
        score = self.run_lightdock_pair(receptor, ligand, output_dir)
        
        # Cache result
        if score is not None:
            score_file.parent.mkdir(parents=True, exist_ok=True)
            score_file.write_text(str(score))
        else:
            score = 0.0
            
        return self.create_result_row(receptor, ligand, score)
    
    def create_result_row(self, receptor, ligand, score):
        """Convert score to result row"""
        
        # Convert score to pseudo-TM
        if score == 777.0:  # Error
            pseudo_tm = 0.01
        elif score == 888.0:  # Timeout
            pseudo_tm = 0.02
        elif score < -100:
            pseudo_tm = 0.9  # Very good binding
        elif score < -50:
            pseudo_tm = 0.7  # Good binding
        elif score < -25:
            pseudo_tm = 0.5  # Moderate binding
        elif score < 0:
            pseudo_tm = 0.3  # Weak binding
        else:
            pseudo_tm = 0.1  # Poor/no binding
        
        return [
            receptor.stem, ligand.stem,
            pseudo_tm, 0, 0, 0, 0, 0, 0, 0,
            score, 0, pseudo_tm, 0, 0
        ]
    
    def run_all_vs_all_cluster(self, cluster_dir):
        """Run all-vs-all docking for one cluster"""
        cluster_name = cluster_dir.name
        
        # Collect PDB files
        pdb_files = list(cluster_dir.glob("*.pdb"))
        for af_dir in cluster_dir.glob("*_alphafold"):
            if af_dir.is_dir():
                af_pdbs = list(af_dir.glob("*.pdb"))
                if af_pdbs:
                    pdb_files.append(af_pdbs[0])
        
        n = len(pdb_files)
        if n < 2:
            return None
        
        print(f"Processing cluster {cluster_name}: {n} proteins, {n*n} total pairs")
        
        output_base = self.base_dir / "lightdock_results" / cluster_name
        
        # Generate all protein pairs
        results = []
        pair_count = 0
        total_pairs = n * n
        
        for i, receptor in enumerate(pdb_files):
            for j, ligand in enumerate(pdb_files):
                pair_count += 1
                
                print(f"  [{pair_count}/{total_pairs}] {receptor.stem} vs {ligand.stem}")
                result_row = self.process_protein_pair(i, j, receptor, ligand, output_base, cluster_name)
                results.append(result_row)
                
                # Show score for non-self pairs
                if i != j:
                    score = result_row[10]  # Score is at index 10
                    if score > 500:  # Error/timeout scores
                        status = "ERROR/TIMEOUT"
                    elif score == 0.0:
                        status = "FAILED"
                    else:
                        status = f"score: {score:.3f}"
                    print(f"    → {status}")
        
        # Save results
        df = pd.DataFrame(results, columns=[
            'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'lddt', 'tm', 'qtm'
        ])
        
        output_tsv = output_base.parent / "results_tsv" / f"{cluster_name}_all_vs_all.tsv"
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_tsv, sep='\t', index=False, header=False)
        
        # Show score distribution
        scores = [r[10] for r in results if r[10] not in [0.0, 777.0, 888.0, -999.0]]
        if scores:
            print(f"  Score distribution: min={min(scores):.2f}, max={max(scores):.2f}, mean={np.mean(scores):.2f}")
        
        print(f"  ✓ Saved {len(results)} results to {output_tsv}")
        return df
    
    def process_all_clusters(self, category_prefix, max_workers=4):  # Reduced workers to avoid conflicts
        """Process clusters with controlled concurrency"""
        
        if "RANDOM" in category_prefix:
            search_dir = self.base_dir / "random_mgc_pdb_files"
        else:
            search_dir = self.base_dir / "mgc_pdb_files"
        
        if not search_dir.exists():
            print(f"Directory not found: {search_dir}")
            return
        
        cluster_dirs = [d for d in search_dir.iterdir() 
                       if d.is_dir() and d.name.startswith(category_prefix)]
        
        print(f"Processing {len(cluster_dirs)} {category_prefix} clusters with LightDock")
        
        # Start with sequential processing for stability
        for i, cluster_dir in enumerate(cluster_dirs):
            print(f"\n=== Cluster {i+1}/{len(cluster_dirs)}: {cluster_dir.name} ===")
            try:
                result = self.run_all_vs_all_cluster(cluster_dir)
                if result is not None:
                    scores = [r[10] for r in result.values.tolist() if r[10] not in [0.0, 777.0, 888.0, -999.0]]
                    success_rate = len(scores) / len(result) if len(result) > 0 else 0
                    print(f"✓ {cluster_dir.name}: Success rate {success_rate:.1%} ({len(scores)}/{len(result)} pairs)")
                    
                    if len(scores) > 0:
                        print(f"  Score range: {min(scores):.2f} to {max(scores):.2f}")
                else:
                    print(f"⚠ {cluster_dir.name}: Skipped")
            except Exception as e:
                print(f"✗ {cluster_dir.name}: Error - {e}")
                import traceback
                traceback.print_exc()

def main():
    pipeline = LightDockPipeline()
    
    print("Starting LightDock Production Pipeline")
    print("=" * 50)
    
    # Process all categories
    for category in ['MGC_CANDIDATE', 'RANDOM', 'BGC']:
        print(f"\n=== Processing {category} clusters ===")
        pipeline.process_all_clusters(category, max_workers=32)
        print(f"Completed {category}")
    
    print("\nAll docking analysis complete!")
    print("Results saved in lightdock_results/")

if __name__ == "__main__":
    main()