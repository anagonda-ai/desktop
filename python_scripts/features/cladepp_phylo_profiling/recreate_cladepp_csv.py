#!/usr/bin/env python3
"""
Recreate Cladepp CSV from existing summary.csv files
"""

import pandas as pd
import numpy as np
from pathlib import Path
import glob

def recreate_cladepp_csv():
    """Recreate Cladepp CSV from existing summary.csv files"""
    print("📊 RECREATING CLADEPP CSV FROM EXISTING RESULTS...")
    
    # Define directories
    mgc_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_final")
    random_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_candidates_dir_fixed")
    
    all_metrics = []
    
    # Process MGC directory (MIBiG and KEGG-MGC)
    mgc_clusters = list(mgc_dir.glob("*"))
    print(f"Found {len(mgc_clusters)} MGC clusters")
    
    # Process RANDOM directory
    random_clusters = list(random_dir.glob("*"))
    print(f"Found {len(random_clusters)} RANDOM clusters")
    
    # Combine all clusters
    all_clusters = mgc_clusters + random_clusters
    print(f"Total clusters to process: {len(all_clusters)}")
    
    for cluster_dir in all_clusters:
        cluster_name = cluster_dir.name
        summary_file = cluster_dir / "summary.csv"
        
        if not summary_file.exists():
            print(f"  No summary.csv found for {cluster_name}")
            continue
            
        print(f"Processing {cluster_name}...")
        
        try:
            # Read the summary CSV
            df = pd.read_csv(summary_file)
            
            if len(df) == 0:
                print(f"  No data in {cluster_name}")
                continue
            
            # Calculate metrics from the summary data
            n_clades = len(df)
            n_successful = len(df[df['status'] == 'success'])
            n_failed = len(df[df['status'] == 'failed'])
            
            # Calculate various metrics
            clade_sizes = df['clade_size'].astype(float)
            mean_clade_sizes = clade_sizes.mean()
            std_clade_sizes = clade_sizes.std()
            median_clade_sizes = clade_sizes.median()
            min_clade_sizes = clade_sizes.min()
            max_clade_sizes = clade_sizes.max()
            
            # Cladepp scores
            cladepp_scores = df['cladepp_score'].astype(float)
            mean_cladepp_score = cladepp_scores.mean()
            std_cladepp_score = cladepp_scores.std()
            median_cladepp_score = cladepp_scores.median()
            min_cladepp_score = cladepp_scores.min()
            max_cladepp_score = cladepp_scores.max()
            
            # Anchor correlation metrics
            mean_anchor_corr = df['mean_anchor_corr'].astype(float).mean()
            std_anchor_corr = df['std_anchor_corr'].astype(float).mean()
            max_anchor_corr = df['max_anchor_corr'].astype(float).max()
            min_anchor_corr = df['min_anchor_corr'].astype(float).min()
            
            # Positive correlation pairs
            positive_corr_pairs = df['positive_corr_pairs'].notna().sum()
            total_corr_pairs = len(df)
            fraction_positive_corr = positive_corr_pairs / total_corr_pairs if total_corr_pairs > 0 else 0
            
            # Gain loss scores
            gain_loss_scores = df['gain_loss_score'].astype(float)
            mean_gain_loss_score = gain_loss_scores.mean()
            std_gain_loss_score = gain_loss_scores.std()
            median_gain_loss_score = gain_loss_scores.median()
            
            # Additional metrics
            success_rate = n_successful / n_clades if n_clades > 0 else 0
            total_tips = clade_sizes.sum()
            mean_tips_per_clade = total_tips / n_clades if n_clades > 0 else 0
            
            # Size bin categorization
            if n_clades <= 5:
                size_bin = "small"
            elif n_clades <= 20:
                size_bin = "medium"
            else:
                size_bin = "large"
            
            # Quality metrics
            high_quality_clades = len(df[(df['cladepp_score'] > -0.5) & (df['status'] == 'success')])
            medium_quality_clades = len(df[(df['cladepp_score'] <= -0.5) & (df['cladepp_score'] > -1.0) & (df['status'] == 'success')])
            low_quality_clades = len(df[(df['cladepp_score'] <= -1.0) & (df['status'] == 'success')])
            
            fraction_high_quality = high_quality_clades / n_successful if n_successful > 0 else 0
            fraction_medium_quality = medium_quality_clades / n_successful if n_successful > 0 else 0
            fraction_low_quality = low_quality_clades / n_successful if n_successful > 0 else 0
            
            metrics = {
                'name': cluster_name,
                'n_clades': n_clades,
                'n_successful': n_successful,
                'n_failed': n_failed,
                'success_rate': success_rate,
                'mean_clade_sizes': mean_clade_sizes,
                'std_clade_sizes': std_clade_sizes,
                'median_clade_sizes': median_clade_sizes,
                'min_clade_sizes': min_clade_sizes,
                'max_clade_sizes': max_clade_sizes,
                'total_tips': total_tips,
                'mean_tips_per_clade': mean_tips_per_clade,
                'mean_cladepp_score': mean_cladepp_score,
                'std_cladepp_score': std_cladepp_score,
                'median_cladepp_score': median_cladepp_score,
                'min_cladepp_score': min_cladepp_score,
                'max_cladepp_score': max_cladepp_score,
                'mean_anchor_corr': mean_anchor_corr,
                'std_anchor_corr': std_anchor_corr,
                'max_anchor_corr': max_anchor_corr,
                'min_anchor_corr': min_anchor_corr,
                'positive_corr_pairs': positive_corr_pairs,
                'fraction_positive_corr': fraction_positive_corr,
                'mean_gain_loss_score': mean_gain_loss_score,
                'std_gain_loss_score': std_gain_loss_score,
                'median_gain_loss_score': median_gain_loss_score,
                'high_quality_clades': high_quality_clades,
                'medium_quality_clades': medium_quality_clades,
                'low_quality_clades': low_quality_clades,
                'fraction_high_quality': fraction_high_quality,
                'fraction_medium_quality': fraction_medium_quality,
                'fraction_low_quality': fraction_low_quality,
                'size_bin': size_bin
            }
            
            all_metrics.append(metrics)
            print(f"  Processed {n_clades} clades ({n_successful} successful)")
            
        except Exception as e:
            print(f"  Error processing {cluster_name}: {e}")
            continue
    
    # Create DataFrame and save
    df_metrics = pd.DataFrame(all_metrics)
    output_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/cladepp_cluster_metrics_enhanced.csv"
    df_metrics.to_csv(output_file, index=False)
    
    print(f"\n✅ Created Cladepp CSV with {len(df_metrics)} clusters")
    print(f"   Saved to: {output_file}")
    
    # Show cluster type distribution
    mibig_clusters = df_metrics[df_metrics['name'].str.startswith('BGC')]
    mgc_clusters = df_metrics[df_metrics['name'].str.startswith('MGC_CANDIDATE')]
    random_clusters = df_metrics[df_metrics['name'].str.startswith('RANDOM_MGC')]
    
    print(f"   MIBiG clusters: {len(mibig_clusters)}")
    print(f"   KEGG-MGC clusters: {len(mgc_clusters)}")
    print(f"   RANDOM clusters: {len(random_clusters)}")
    
    if len(mibig_clusters) > 0:
        print(f"   MIBiG cluster names: {sorted(mibig_clusters['name'].tolist())}")
    
    return df_metrics

if __name__ == "__main__":
    recreate_cladepp_csv()

