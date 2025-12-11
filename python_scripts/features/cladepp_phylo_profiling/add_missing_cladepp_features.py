#!/usr/bin/env python3
"""
Add missing Cladepp features by computing them from summary.csv files
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re

def count_positive_correlations(positive_corr_pairs_str):
    """
    Count positive correlations from the positive_corr_pairs string.
    Format: "gene_1-gene_3:0.69; gene_1-gene_4:0.89; ..."
    Returns: (positive_count, total_count)
    """
    if pd.isna(positive_corr_pairs_str) or positive_corr_pairs_str == '':
        return 0, 0
    
    # Split by semicolon to get individual pairs
    pairs = str(positive_corr_pairs_str).split(';')
    positive_count = 0
    total_count = 0
    
    for pair in pairs:
        pair = pair.strip()
        if not pair:
            continue
        
        # Extract correlation value (format: "gene_1-gene_3:0.69")
        match = re.search(r':([-+]?\d*\.?\d+)', pair)
        if match:
            corr_value = float(match.group(1))
            total_count += 1
            if corr_value > 0:
                positive_count += 1
    
    return positive_count, total_count

def compute_features_from_summary(summary_df):
    """
    Compute all required Cladepp features from a summary DataFrame.
    
    Args:
        summary_df: DataFrame with columns: cladepp_score, clade_size, positive_corr_pairs, etc.
    
    Returns:
        Dictionary with all computed features
    """
    # Filter only successful clades
    if 'status' in summary_df.columns:
        summary_df = summary_df[summary_df['status'] == 'success'].copy()
    
    if len(summary_df) == 0:
        return {
            'mean_cladepp_score': np.nan,
            'weighted_cladepp_score': np.nan,
            'positive_correlation_ratio': np.nan,
            'cladepp_multi_clade_high': np.nan,
            'cladepp_multi_clade_medium': np.nan,
            'cladepp_conservation_consistency': np.nan,
            'cladepp_max_pair_score': np.nan
        }
    
    # Get cladepp scores and sizes
    cladepp_scores = summary_df['cladepp_score'].astype(float).dropna()
    clade_sizes = summary_df['clade_size'].astype(float)
    
    if len(cladepp_scores) == 0:
        return {
            'mean_cladepp_score': np.nan,
            'weighted_cladepp_score': np.nan,
            'positive_correlation_ratio': np.nan,
            'cladepp_multi_clade_high': np.nan,
            'cladepp_multi_clade_medium': np.nan,
            'cladepp_conservation_consistency': np.nan,
            'cladepp_max_pair_score': np.nan
        }
    
    # 1. mean_cladepp_score: Mean of per-clade co-evolution scores
    mean_cladepp_score = cladepp_scores.mean()
    
    # 2. weighted_cladepp_score: Clade-size weighted mean co-evolution score
    # Align sizes with scores (handle NaN scores)
    valid_mask = summary_df['cladepp_score'].notna()
    valid_scores = summary_df.loc[valid_mask, 'cladepp_score'].astype(float)
    valid_sizes = summary_df.loc[valid_mask, 'clade_size'].astype(float)
    
    if len(valid_scores) > 0 and valid_sizes.sum() > 0:
        weighted_cladepp_score = np.average(valid_scores, weights=valid_sizes)
    else:
        weighted_cladepp_score = mean_cladepp_score
    
    # 3. positive_correlation_ratio: Fraction of all anchor gene pairs showing positive correlation
    total_positive = 0
    total_pairs = 0
    
    if 'positive_corr_pairs' in summary_df.columns:
        for _, row in summary_df.iterrows():
            pos_count, pair_count = count_positive_correlations(row.get('positive_corr_pairs'))
            total_positive += pos_count
            total_pairs += pair_count
    
    # Also check mean_anchor_corr as fallback - if mean > 0, it suggests positive correlation
    if total_pairs == 0 and 'mean_anchor_corr' in summary_df.columns:
        # Use mean_anchor_corr as proxy: if mean > 0, assume positive correlation
        mean_anchor_corr = summary_df['mean_anchor_corr'].astype(float).mean()
        positive_correlation_ratio = 1.0 if mean_anchor_corr > 0 else 0.0
    else:
        positive_correlation_ratio = total_positive / total_pairs if total_pairs > 0 else np.nan
    
    # 4. cladepp_multi_clade_high: Fraction of clades with mean pairwise correlation > 0.8
    n_clades = len(cladepp_scores)
    n_high = (cladepp_scores > 0.8).sum()
    cladepp_multi_clade_high = n_high / n_clades if n_clades > 0 else 0.0
    
    # 5. cladepp_multi_clade_medium: Fraction of clades with mean pairwise correlation > 0.6
    n_medium = (cladepp_scores > 0.6).sum()
    cladepp_multi_clade_medium = n_medium / n_clades if n_clades > 0 else 0.0
    
    # 6. cladepp_conservation_consistency: 1 - (std/mean) of clade scores, clipped to [0,1]
    mean_score = cladepp_scores.mean()
    std_score = cladepp_scores.std()
    
    if mean_score > 0 and not np.isnan(mean_score) and not np.isnan(std_score):
        consistency = 1 - (std_score / mean_score)
        cladepp_conservation_consistency = max(0.0, min(1.0, consistency))  # Clip to [0, 1]
    else:
        cladepp_conservation_consistency = 0.0
    
    # 7. cladepp_max_pair_score: Maximum mean pairwise correlation observed in any single clade
    cladepp_max_pair_score = cladepp_scores.max()
    
    return {
        'mean_cladepp_score': mean_cladepp_score,
        'weighted_cladepp_score': weighted_cladepp_score,
        'positive_correlation_ratio': positive_correlation_ratio,
        'cladepp_multi_clade_high': cladepp_multi_clade_high,
        'cladepp_multi_clade_medium': cladepp_multi_clade_medium,
        'cladepp_conservation_consistency': cladepp_conservation_consistency,
        'cladepp_max_pair_score': cladepp_max_pair_score
    }

def add_missing_features():
    """Add missing Cladepp features to the enhanced CSV by computing from summary.csv files"""
    
    # Load the enhanced Cladepp CSV
    df = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/cladepp_cluster_metrics_enhanced.csv')
    
    print(f"Loaded {len(df)} clusters")
    print(f"Columns: {list(df.columns)}")
    
    # Define directories where summary.csv files are located
    mgc_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_final")
    random_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_candidates_dir_fixed")
    
    # Required features from the ML model
    required_features = [
        'mean_cladepp_score',
        'weighted_cladepp_score',
        'positive_correlation_ratio',
        'cladepp_multi_clade_high',
        'cladepp_multi_clade_medium',
        'cladepp_conservation_consistency',
        'cladepp_max_pair_score'
    ]
    
    # Initialize feature columns if they don't exist
    for feature in required_features:
        if feature not in df.columns:
            df[feature] = np.nan
    
    # Process each cluster
    print(f"\nComputing features from summary.csv files...")
    computed_count = 0
    missing_count = 0
    
    for idx, row in df.iterrows():
        cluster_name = row['name']
        
        # Try to find summary.csv in either directory
        summary_file = None
        for base_dir in [mgc_dir, random_dir]:
            potential_file = base_dir / cluster_name / "summary.csv"
            if potential_file.exists():
                summary_file = potential_file
                break
        
        if summary_file is None or not summary_file.exists():
            missing_count += 1
            continue
        
        try:
            # Load summary.csv for this cluster
            summary_df = pd.read_csv(summary_file)
            
            # Compute all features
            features = compute_features_from_summary(summary_df)
            
            # Update the dataframe
            for feature_name, feature_value in features.items():
                df.at[idx, feature_name] = feature_value
            
            computed_count += 1
            
            if (idx + 1) % 50 == 0:
                print(f"  Processed {idx + 1}/{len(df)} clusters...")
                
        except Exception as e:
            print(f"  Error processing {cluster_name}: {e}")
            missing_count += 1
            continue
    
    print(f"\n✅ Computed features for {computed_count} clusters")
    if missing_count > 0:
        print(f"⚠️  {missing_count} clusters had missing or invalid summary.csv files")
    
    # Clean up duplicate columns (remove _x and _y suffixes)
    columns_to_remove = []
    for col in df.columns:
        if col.endswith('_y'):
            columns_to_remove.append(col)
        elif col.endswith('_x'):
            # Rename _x columns to remove the suffix
            new_name = col.replace('_x', '')
            if new_name not in df.columns:
                df[new_name] = df[col]
            columns_to_remove.append(col)
    
    # Remove duplicate columns
    if columns_to_remove:
        df = df.drop(columns=columns_to_remove)
    
    # Remove duplicate cluster_name columns
    if 'cluster_name_x' in df.columns:
        df = df.drop(columns=['cluster_name_x'])
    if 'cluster_name_y' in df.columns:
        df = df.drop(columns=['cluster_name_y'])
    
    print(f"✅ Cleaned up duplicate columns")
    print(f"Final columns: {list(df.columns)}")
    
    # Save the updated CSV
    output_file = '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/cladepp_cluster_metrics_enhanced.csv'
    df.to_csv(output_file, index=False)
    
    print(f"\n✅ Saved updated Cladepp CSV with {len(df)} clusters and {len(df.columns)} columns")
    
    # Verify the required features exist and show statistics
    print(f"\n📊 Feature Verification:")
    for feature in required_features:
        if feature in df.columns:
            non_null_count = df[feature].notna().sum()
            print(f"  ✅ {feature:40s}: {non_null_count}/{len(df)} non-null values")
            if non_null_count > 0:
                print(f"     Mean: {df[feature].mean():.4f}, Std: {df[feature].std():.4f}")
        else:
            print(f"  ❌ {feature:40s}: MISSING")

if __name__ == "__main__":
    add_missing_features()
