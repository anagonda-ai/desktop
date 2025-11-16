#!/usr/bin/env python3
"""
Recreate Lightdock CSV from existing best_score.txt files (with parallel processing)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import sys

def process_cluster(cluster_dir):
    """Process a single cluster directory and return metrics"""
    cluster_name = cluster_dir.name
    print(f"Processing cluster: {cluster_name}")
    
    # Collect all best_score.txt files
    scores = []
    for subdir in cluster_dir.iterdir():
        if subdir.is_dir():
            score_file = subdir / "best_score.txt"
            if score_file.exists():
                try:
                    with open(score_file, 'r') as f:
                        content = f.read().strip()
                        # Extract score from "Best docking score: -188.123" format
                        if "Best docking score:" in content:
                            score_str = content.split("Best docking score:")[1].strip()
                            score = float(score_str)
                            scores.append(score)
                except:
                    continue
    
    if len(scores) == 0:
        return None
        
    # Calculate metrics
    scores_array = np.array(scores)
    
    # Enhanced features calculation
    n = len(scores)
    n_squared = n * n
    n_comparisons_total = n
    n_comparisons_non_self = n  # For docking, all interactions are non-self
    
    sum_scores_all = scores_array.sum()
    sum_scores_non_self = sum_scores_all  # Same as all for docking
    
    mean_score_all = float(scores_array.mean())
    mean_score_non_self = mean_score_all  # Same as all for docking
    
    var_score_non_self = float(scores_array.var())
    std_score_non_self = float(scores_array.std())
    
    median_score_non_self = float(np.median(scores_array))
    q25_score = float(np.percentile(scores_array, 25))
    q75_score = float(np.percentile(scores_array, 75))
    min_score = float(scores_array.min())
    max_score = float(scores_array.max())
    
    # Binding fractions (assuming scores < -50 are strong, -50 to -25 medium, > -25 weak)
    fraction_strong_binders = float((scores_array < -50).mean())
    fraction_moderate_binders = float(((scores_array >= -50) & (scores_array < -25)).mean())
    fraction_weak_binders = float((scores_array >= -25).mean())
    
    # Additional metrics
    log_n = float(np.log(n)) if n > 0 else 0
    sqrt_n = float(np.sqrt(n))
    expected_random = -30.0  # Expected random docking score
    enrichment_score = mean_score_all / expected_random if expected_random != 0 else 0
    
    # Size bin
    if n <= 10:
        size_bin = "small"
    elif n <= 50:
        size_bin = "medium"
    else:
        size_bin = "large"
    
    # Z-score and effect size
    z_score = (mean_score_all - expected_random) / std_score_non_self if std_score_non_self > 0 else 0
    effect_size = (mean_score_all - expected_random) / std_score_non_self if std_score_non_self > 0 else 0
    
    metrics = {
        'name': cluster_name,
        'n': n,
        'n_squared': n_squared,
        'n_comparisons_total': n_comparisons_total,
        'n_comparisons_non_self': n_comparisons_non_self,
        'sum_scores_all': sum_scores_all,
        'sum_scores_non_self': sum_scores_non_self,
        'mean_score_all': mean_score_all,
        'mean_score_non_self': mean_score_non_self,
        'var_score_non_self': var_score_non_self,
        'std_score_non_self': std_score_non_self,
        'median_score_non_self': median_score_non_self,
        'q25_score': q25_score,
        'q75_score': q75_score,
        'min_score': min_score,
        'max_score': max_score,
        'fraction_strong_binders': fraction_strong_binders,
        'fraction_moderate_binders': fraction_moderate_binders,
        'fraction_weak_binders': fraction_weak_binders,
        'log_n': log_n,
        'sqrt_n': sqrt_n,
        'expected_random': expected_random,
        'enrichment_score': enrichment_score,
        'size_bin': size_bin,
        'z_score': z_score,
        'effect_size': effect_size,
        # Keep original simple metrics for compatibility
        'total_interactions': n,
        'mean_score': mean_score_all,
        'std_score': std_score_non_self,
        'median_score': median_score_non_self,
        'high_score_interactions': int((scores_array < -50).sum()),
        'medium_score_interactions': int(((scores_array >= -50) & (scores_array < -25)).sum()),
        'low_score_interactions': int((scores_array >= -25).sum())
    }
    
    return (cluster_name, len(scores), metrics)

def recreate_lightdock_csv():
    """Recreate Lightdock CSV from existing best_score.txt files using parallel processing"""
    print("⚡ RECREATING LIGHTDOCK CSV FROM EXISTING RESULTS (PARALLEL)...")
    
    results_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_lightdock_result")
    
    # Get all cluster directories
    cluster_dirs = [d for d in results_dir.iterdir() if d.is_dir()]
    total_clusters = len(cluster_dirs)
    
    # Determine number of workers (use all available CPUs)
    n_workers = cpu_count()
    print(f"📊 Found {total_clusters} clusters to process")
    print(f"🚀 Using {n_workers} CPU cores for parallel processing")
    
    all_metrics = []
    processed = 0
    
    # Process clusters in parallel
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all tasks
        future_to_cluster = {executor.submit(process_cluster, cluster_dir): cluster_dir 
                            for cluster_dir in cluster_dirs}
        
        # Process results as they complete
        for future in as_completed(future_to_cluster):
            cluster_dir = future_to_cluster[future]
            try:
                result = future.result()
                if result is not None:
                    cluster_name, n_scores, metrics = result
                    all_metrics.append(metrics)
                    processed += 1
                    print(f"[{processed}/{total_clusters}] ✓ {cluster_name}: {n_scores} interactions")
                else:
                    processed += 1
                    print(f"[{processed}/{total_clusters}] ✗ {cluster_dir.name}: No scores found")
            except Exception as e:
                processed += 1
                print(f"[{processed}/{total_clusters}] ✗ {cluster_dir.name}: Error - {e}")
    
    if len(all_metrics) == 0:
        print("\n❌ No valid clusters found!")
        return None
    
    # Create DataFrame and save
    df_metrics = pd.DataFrame(all_metrics)
    output_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/lightdock_cluster_metrics.csv"
    
    # Ensure output directory exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    df_metrics.to_csv(output_file, index=False)
    
    print(f"\n✅ Created Lightdock CSV with {len(df_metrics)} clusters")
    print(f"   Saved to: {output_file}")
    
    # Show MIBiG clusters
    mibig_clusters = df_metrics[df_metrics['name'].str.startswith('BGC')]
    print(f"   MIBiG clusters: {len(mibig_clusters)}")
    
    if len(mibig_clusters) > 0:
        print(f"   MIBiG cluster names: {sorted(mibig_clusters['name'].tolist())}")
    
    return df_metrics

if __name__ == "__main__":
    recreate_lightdock_csv()