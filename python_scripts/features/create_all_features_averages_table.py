#!/usr/bin/env python3
"""
Create comprehensive table of average values for all 26 features across MGC vs KEGG-Random groups
"""

import pandas as pd
import numpy as np
from pathlib import Path

def create_all_features_averages_table():
    """Create comprehensive table of average values for all 26 features"""
    
    BASE_DIR = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random"
    
    # Feature configurations - all 26 features across 5 categories
    feature_configs = {
        'docking': {
            'path': f'{BASE_DIR}/lightdock_cluster_metrics.csv',
            'cluster_col': 'name',
            'features': [
                'mean_score_non_self',
                'enrichment_score',
                'z_score',
                'effect_size'
            ],
            'name': 'Docking'
        },
        'foldseek': {
            'path': f'{BASE_DIR}/foldseek_cluster_metrics_enhanced.csv',
            'cluster_col': 'name',
            'features': [
                'mean_score_non_self',
                'enrichment_score',
                'z_score',
                'effect_size',
                'foldseek_match_coverage'
            ],
            'name': 'Foldseek'
        },
        'e2p2': {
            'path': f'{BASE_DIR}/e2p2_feature_extraction_detailed.csv',
            'cluster_col': 'cluster_name',
            'features': [
                'num_distinct_enzyme_classes',
                'num_distinct_enzyme_subclasses',
                'num_distinct_enzyme_families',
                'num_distinct_enzyme_subfamilies',
                'total_ec_numbers'
            ],
            'name': 'E2P2'
        },
        'promoter_similarity': {
            'path': f'{BASE_DIR}/promoter_similarity_results.csv',
            'cluster_col': 'group_name',
            'features': [
                'similarity_score',
                'normalized_similarity',
                'mean_proximal_similarity',
                'mean_distal_similarity',
                'num_tfbs_types_found'
            ],
            'name': 'Promoter Similarity'
        },
        'cladepp': {
            'path': f'{BASE_DIR}/cladepp_cluster_metrics_enhanced.csv',
            'cluster_col': 'cluster_name',
            'features': [
                'mean_cladepp_score',
                'weighted_cladepp_score',
                'positive_correlation_ratio',
                'cladepp_multi_clade_high',
                'cladepp_multi_clade_medium',
                'cladepp_conservation_consistency',
                'cladepp_max_pair_score'
            ],
            'name': 'Cladepp'
        }
    }
    
    print("="*120)
    print("COMPREHENSIVE FEATURE AVERAGES: MGC vs KEGG-RANDOM COMPARISON")
    print("All 26 Features Across 5 Categories")
    print("="*120)
    
    total_features = sum(len(c['features']) for c in feature_configs.values())
    print(f"\nTotal features: {total_features}")
    for name, config in feature_configs.items():
        print(f"  {config['name']:18s}: {len(config['features']):2d} features")
    
    # Load and merge all data
    print("\n" + "="*120)
    print("LOADING AND MERGING DATA")
    print("="*120)
    
    merged_df = None
    
    for key, config in feature_configs.items():
        try:
            df = pd.read_csv(config['path'])
            
            # Create label based on cluster name patterns
            # MIBiG: BGC*, KEGG-MGC: MGC_CANDIDATE*, RANDOM: RANDOM*
            def create_label(cluster_name):
                if str(cluster_name).startswith('BGC'):
                    return 'MIBiG'  # MIBiG
                elif str(cluster_name).startswith('MGC_CANDIDATE'):
                    return 'KEGG_MGC'  # KEGG-MGC
                elif str(cluster_name).startswith('RANDOM'):
                    return 'KEGG_Random'  # KEGG-Random
                else:
                    return None
            
            df['label'] = df[config['cluster_col']].apply(create_label)
            df = df.dropna(subset=['label'])  # Remove any unmatched clusters
            
            # Select columns
            cluster_col = config['cluster_col']
            cols_to_keep = [cluster_col, 'label'] + config['features']
            df = df[cols_to_keep]
            
            # Rename cluster column to standard name
            df = df.rename(columns={cluster_col: 'cluster_name'})
            
            # Add prefix to feature names to avoid collisions
            prefix = key + '_'
            rename_dict = {feat: prefix + feat for feat in config['features']}
            df = df.rename(columns=rename_dict)
            
            # Update config with prefixed names
            config['prefixed_features'] = [prefix + feat for feat in config['features']]
            
            print(f"\n{config['name']:18s}: {len(df):6d} clusters, {len(config['features']):2d} features")
            print(f"  MIBiG clusters: {(df['label']=='MIBiG').sum()}")
            print(f"  KEGG-MGC clusters: {(df['label']=='KEGG_MGC').sum()}")
            print(f"  KEGG-Random clusters: {(df['label']=='KEGG_Random').sum()}")
            
            # Merge
            if merged_df is None:
                merged_df = df
            else:
                merged_df = merged_df.merge(df, on=['cluster_name', 'label'], how='inner')
                
        except Exception as e:
            print(f"Error loading {config['name']}: {e}")
            continue
    
    if merged_df is None:
        print("No data loaded successfully!")
        return None
    
    print(f"\nMerged dataset: {len(merged_df)} clusters")
    print(f"  MIBiG clusters: {(merged_df['label']=='MIBiG').sum()}")
    print(f"  KEGG-MGC clusters: {(merged_df['label']=='KEGG_MGC').sum()}")
    print(f"  KEGG-Random clusters: {(merged_df['label']=='KEGG_Random').sum()}")
    
    # Collect all features
    all_features = []
    for config in feature_configs.values():
        if 'prefixed_features' in config:
            all_features.extend(config['prefixed_features'])
    
    print(f"\nTotal features in merged dataset: {len(all_features)}")
    
    # Calculate statistics for each feature
    results = []
    
    for config in feature_configs.values():
        if 'prefixed_features' not in config:
            continue
            
        for feature, original_name in zip(config['prefixed_features'], config['features']):
            if feature not in merged_df.columns:
                continue
                
            # MIBiG group
            mibig_data = merged_df[merged_df['label'] == 'MIBiG'][feature]
            mibig_mean = mibig_data.mean()
            mibig_std = mibig_data.std()
            mibig_median = mibig_data.median()
            mibig_min = mibig_data.min()
            mibig_max = mibig_data.max()
            
            # KEGG-MGC group
            kegg_mgc_data = merged_df[merged_df['label'] == 'KEGG_MGC'][feature]
            kegg_mgc_mean = kegg_mgc_data.mean()
            kegg_mgc_std = kegg_mgc_data.std()
            kegg_mgc_median = kegg_mgc_data.median()
            kegg_mgc_min = kegg_mgc_data.min()
            kegg_mgc_max = kegg_mgc_data.max()
            
            # KEGG-Random group
            random_data = merged_df[merged_df['label'] == 'KEGG_Random'][feature]
            random_mean = random_data.mean()
            random_std = random_data.std()
            random_median = random_data.median()
            random_min = random_data.min()
            random_max = random_data.max()
            
            # Calculate differences for MIBiG vs KEGG-Random
            mibig_fold_change = mibig_mean / random_mean if random_mean != 0 else np.nan
            mibig_difference = mibig_mean - random_mean
            mibig_percent_change = ((mibig_mean - random_mean) / random_mean) * 100 if random_mean != 0 else np.nan
            
            # Calculate differences for KEGG-MGC vs KEGG-Random
            kegg_mgc_fold_change = kegg_mgc_mean / random_mean if random_mean != 0 else np.nan
            kegg_mgc_difference = kegg_mgc_mean - random_mean
            kegg_mgc_percent_change = ((kegg_mgc_mean - random_mean) / random_mean) * 100 if random_mean != 0 else np.nan
            
            # Calculate differences for MIBiG vs KEGG-MGC
            mibig_vs_kegg_fold_change = mibig_mean / kegg_mgc_mean if kegg_mgc_mean != 0 else np.nan
            mibig_vs_kegg_difference = mibig_mean - kegg_mgc_mean
            mibig_vs_kegg_percent_change = ((mibig_mean - kegg_mgc_mean) / kegg_mgc_mean) * 100 if kegg_mgc_mean != 0 else np.nan
            
            results.append({
                'Category': config['name'],
                'Feature': original_name,
                'Full_Feature_Name': feature,
                'MIBiG_Mean': mibig_mean,
                'MIBiG_Std': mibig_std,
                'MIBiG_Median': mibig_median,
                'MIBiG_Min': mibig_min,
                'MIBiG_Max': mibig_max,
                'KEGG_MGC_Mean': kegg_mgc_mean,
                'KEGG_MGC_Std': kegg_mgc_std,
                'KEGG_MGC_Median': kegg_mgc_median,
                'KEGG_MGC_Min': kegg_mgc_min,
                'KEGG_MGC_Max': kegg_mgc_max,
                'KEGG_Random_Mean': random_mean,
                'KEGG_Random_Std': random_std,
                'KEGG_Random_Median': random_median,
                'KEGG_Random_Min': random_min,
                'KEGG_Random_Max': random_max,
                'MIBiG_vs_Random_Difference': mibig_difference,
                'MIBiG_vs_Random_Percent_Change': mibig_percent_change,
                'MIBiG_vs_Random_Fold_Change': mibig_fold_change,
                'KEGG_MGC_vs_Random_Difference': kegg_mgc_difference,
                'KEGG_MGC_vs_Random_Percent_Change': kegg_mgc_percent_change,
                'KEGG_MGC_vs_Random_Fold_Change': kegg_mgc_fold_change,
                'MIBiG_vs_KEGG_MGC_Difference': mibig_vs_kegg_difference,
                'MIBiG_vs_KEGG_MGC_Percent_Change': mibig_vs_kegg_percent_change,
                'MIBiG_vs_KEGG_MGC_Fold_Change': mibig_vs_kegg_fold_change
            })
    
    # Create DataFrame
    results_df = pd.DataFrame(results)
    
    # Display results by category
    print("\n" + "="*120)
    print("FEATURE AVERAGES BY CATEGORY")
    print("="*120)
    
    for category in results_df['Category'].unique():
        cat_data = results_df[results_df['Category'] == category]
        print(f"\n📊 {category.upper()} FEATURES ({len(cat_data)} features)")
        print("-" * 140)
        print(f"{'Feature':<35} {'MIBiG Mean':<12} {'KEGG-MGC Mean':<14} {'KEGG-Random Mean':<16} {'MIBiG vs Random':<16} {'KEGG-MGC vs Random':<18}")
        print("-" * 140)
        
        for _, row in cat_data.iterrows():
            print(f"{row['Feature']:<35} {row['MIBiG_Mean']:<12.4f} {row['KEGG_MGC_Mean']:<14.4f} {row['KEGG_Random_Mean']:<16.4f} {row['MIBiG_vs_Random_Difference']:<16.4f} {row['KEGG_MGC_vs_Random_Difference']:<18.4f}")
    
    # Overall summary
    print("\n\n" + "="*140)
    print("OVERALL SUMMARY - ALL 26 FEATURES")
    print("="*140)
    print(f"{'Feature':<40} {'MIBiG Mean':<12} {'KEGG-MGC Mean':<14} {'KEGG-Random Mean':<16} {'MIBiG vs Random':<16} {'KEGG-MGC vs Random':<18}")
    print("-" * 140)
    
    for _, row in results_df.iterrows():
        print(f"{row['Category'] + '_' + row['Feature']:<40} {row['MIBiG_Mean']:<12.4f} {row['KEGG_MGC_Mean']:<14.4f} {row['KEGG_Random_Mean']:<16.4f} {row['MIBiG_vs_Random_Difference']:<16.4f} {row['KEGG_MGC_vs_Random_Difference']:<18.4f}")
    
    # Statistical summary
    print("\n\n" + "="*140)
    print("STATISTICAL SUMMARY")
    print("="*140)
    
    # Count features higher/lower for each comparison
    mibig_higher = len(results_df[results_df['MIBiG_vs_Random_Difference'] > 0])
    mibig_lower = len(results_df[results_df['MIBiG_vs_Random_Difference'] < 0])
    kegg_mgc_higher = len(results_df[results_df['KEGG_MGC_vs_Random_Difference'] > 0])
    kegg_mgc_lower = len(results_df[results_df['KEGG_MGC_vs_Random_Difference'] < 0])
    
    print(f"MIBiG vs KEGG-Random:")
    print(f"  Features HIGHER in MIBiG: {mibig_higher}")
    print(f"  Features LOWER in MIBiG: {mibig_lower}")
    print(f"\nKEGG-MGC vs KEGG-Random:")
    print(f"  Features HIGHER in KEGG-MGC: {kegg_mgc_higher}")
    print(f"  Features LOWER in KEGG-MGC: {kegg_mgc_lower}")
    print(f"\nTotal features analyzed: {len(results_df)}")
    
    # Top differences for MIBiG vs Random
    print(f"\nTop 5 features with largest differences (MIBiG vs KEGG-Random):")
    top_mibig_diffs = results_df.nlargest(5, 'MIBiG_vs_Random_Difference')
    for _, row in top_mibig_diffs.iterrows():
        print(f"  {row['Category'] + '_' + row['Feature']}: {row['MIBiG_vs_Random_Difference']:+.4f} ({row['MIBiG_vs_Random_Percent_Change']:+.1f}%)")
    
    # Top differences for KEGG-MGC vs Random
    print(f"\nTop 5 features with largest differences (KEGG-MGC vs KEGG-Random):")
    top_kegg_diffs = results_df.nlargest(5, 'KEGG_MGC_vs_Random_Difference')
    for _, row in top_kegg_diffs.iterrows():
        print(f"  {row['Category'] + '_' + row['Feature']}: {row['KEGG_MGC_vs_Random_Difference']:+.4f} ({row['KEGG_MGC_vs_Random_Percent_Change']:+.1f}%)")
    
    # MIBiG vs KEGG-MGC differences
    print(f"\nTop 5 features with largest differences (MIBiG vs KEGG-MGC):")
    top_mibig_vs_kegg = results_df.nlargest(5, 'MIBiG_vs_KEGG_MGC_Difference')
    for _, row in top_mibig_vs_kegg.iterrows():
        print(f"  {row['Category'] + '_' + row['Feature']}: {row['MIBiG_vs_KEGG_MGC_Difference']:+.4f} ({row['MIBiG_vs_KEGG_MGC_Percent_Change']:+.1f}%)")
    
    # Save to CSV
    output_file = "/groups/itay_mayrose/alongonda/desktop/all_26_features_averages_table.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n💾 Detailed results saved to: {output_file}")
    
    # Create markdown table for documentation
    print("\n\n📋 MARKDOWN TABLE FOR DOCUMENTATION")
    print("-" * 120)
    print("| Category | Feature | MIBiG Mean ± Std | KEGG-MGC Mean ± Std | KEGG-Random Mean ± Std | MIBiG vs Random | KEGG-MGC vs Random |")
    print("|----------|---------|------------------|---------------------|------------------------|-----------------|-------------------|")
    
    for _, row in results_df.iterrows():
        print(f"| {row['Category']} | {row['Feature']} | {row['MIBiG_Mean']:.4f} ± {row['MIBiG_Std']:.4f} | {row['KEGG_MGC_Mean']:.4f} ± {row['KEGG_MGC_Std']:.4f} | {row['KEGG_Random_Mean']:.4f} ± {row['KEGG_Random_Std']:.4f} | {row['MIBiG_vs_Random_Difference']:+.4f} | {row['KEGG_MGC_vs_Random_Difference']:+.4f} |")
    
    return results_df

if __name__ == "__main__":
    results_df = create_all_features_averages_table()
