#!/usr/bin/env python3
"""
Recreate Foldseek CSV from existing all_vs_all.tsv files
"""

import pandas as pd
import numpy as np
from pathlib import Path
import glob

def recreate_foldseek_csv():
    """Recreate Foldseek CSV from existing all_vs_all.tsv files"""
    print("📊 RECREATING FOLDSEEK CSV FROM EXISTING RESULTS...")
    
    # Define both directories
    foldseek_results_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/foldseek_results")
    
    all_metrics = []
    
    # Process MIBiG and KEGG-MGC files
    tsv_files = list(foldseek_results_dir.glob("*_all_vs_all.tsv"))
    print(f"Found {len(tsv_files)} MIBiG/KEGG-MGC Foldseek result files")
    
    for tsv_file in tsv_files:
        cluster_name = tsv_file.name.replace("_all_vs_all.tsv", "")
        print(f"Processing {cluster_name}...")
        
        try:
            # Read the TSV file
            df = pd.read_csv(tsv_file, sep='\t', header=None)
            
            if len(df) == 0:
                print(f"  No data in {cluster_name}")
                continue
            
            # Calculate metrics from the TSV data
            # Column 2 appears to be the score (TM-score)
            scores = df.iloc[:, 2].astype(float)
            
            # Calculate various metrics
            n = len(df)
            n_squared = n * n
            n_comparisons_total = n
            n_comparisons_non_self = len(df[df.iloc[:, 0] != df.iloc[:, 1]])  # Non-self comparisons
            
            sum_scores_all = scores.sum()
            sum_scores_non_self = df[df.iloc[:, 0] != df.iloc[:, 1]].iloc[:, 2].astype(float).sum() if n_comparisons_non_self > 0 else 0
            
            mean_score_all = scores.mean()
            mean_score_non_self = df[df.iloc[:, 0] != df.iloc[:, 1]].iloc[:, 2].astype(float).mean() if n_comparisons_non_self > 0 else 0
            
            var_score_non_self = df[df.iloc[:, 0] != df.iloc[:, 1]].iloc[:, 2].astype(float).var() if n_comparisons_non_self > 0 else 0
            std_score_non_self = np.sqrt(var_score_non_self) if var_score_non_self > 0 else 0
            
            median_score_non_self = df[df.iloc[:, 0] != df.iloc[:, 1]].iloc[:, 2].astype(float).median() if n_comparisons_non_self > 0 else 0
            q25_score = scores.quantile(0.25)
            q75_score = scores.quantile(0.75)
            min_score = scores.min()
            max_score = scores.max()
            
            # Calculate binding fractions (assuming TM-score > 0.5 is strong, > 0.3 is moderate)
            fraction_strong_binders = (scores > 0.5).mean()
            fraction_moderate_binders = ((scores > 0.3) & (scores <= 0.5)).mean()
            fraction_weak_binders = (scores <= 0.3).mean()
            
            # Additional metrics
            log_n = np.log(n) if n > 0 else 0
            sqrt_n = np.sqrt(n)
            expected_random = 0.3  # Expected random TM-score
            enrichment_score = mean_score_all / expected_random if expected_random > 0 else 0
            
            # Size bin (categorize by number of comparisons)
            if n <= 10:
                size_bin = "small"
            elif n <= 50:
                size_bin = "medium"
            else:
                size_bin = "large"
            
            # Z-score and effect size (simplified)
            z_score = (mean_score_all - expected_random) / std_score_non_self if std_score_non_self > 0 else 0
            effect_size = (mean_score_non_self - expected_random) / std_score_non_self if std_score_non_self > 0 else 0
            category = "RANDOM" if cluster_name.startswith("RANDOM") else "MIBiG" if cluster_name.startswith("BGC") else "MGC"
            # Calculate foldseek_match_coverage as the fraction of unique proteins in the cluster that have a strong structural match to at least one other member
            # "Strong" = TM-score > 0.5
            strong_pairs = df[(df.iloc[:, 0] != df.iloc[:, 1]) & (df.iloc[:, 2].astype(float) > 0.5)]
            if len(strong_pairs) > 0:
                proteins_with_strong_match = set(strong_pairs.iloc[:, 0]).union(set(strong_pairs.iloc[:, 1]))
                all_proteins = set(df.iloc[:, 0]).union(set(df.iloc[:, 1]))
                foldseek_match_coverage = len(proteins_with_strong_match) / len(all_proteins) if len(all_proteins) > 0 else 0
            else:
                foldseek_match_coverage = 0
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
                'category': category,
                'foldseek_match_coverage': foldseek_match_coverage
            }
            
            all_metrics.append(metrics)
            print(f"  Processed {n} comparisons")
            
        except Exception as e:
            print(f"  Error processing {cluster_name}: {e}")
            continue
    
    # Create DataFrame and save
    df_metrics = pd.DataFrame(all_metrics)
    output_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/foldseek_cluster_metrics_enhanced.csv"
    df_metrics.to_csv(output_file, index=False)
    
    print(f"\n✅ Created Foldseek CSV with {len(df_metrics)} clusters")
    print(f"   Saved to: {output_file}")
    
    # Show MIBiG clusters
    mibig_clusters = df_metrics[df_metrics['name'].str.startswith('BGC')]
    print(f"   MIBiG clusters: {len(mibig_clusters)}")
    
    if len(mibig_clusters) > 0:
        print(f"   MIBiG cluster names: {sorted(mibig_clusters['name'].tolist())}")
    
    return df_metrics

if __name__ == "__main__":
    recreate_foldseek_csv()
