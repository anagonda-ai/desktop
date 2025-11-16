#!/usr/bin/env python3
"""
Create separate tables for each category (MIBiG, KEGG-MGC, KEGG-Random) with all 26 features
"""

import pandas as pd
import numpy as np
from pathlib import Path

def create_separate_category_tables():
    """Create separate tables for each category with all 26 features"""
    
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
    print("SEPARATE CATEGORY TABLES: MIBiG, KEGG-MGC, KEGG-RANDOM")
    print("All 26 Features Across 5 Categories")
    print("="*120)
    
    # Load and merge all data
    merged_df = None
    
    for key, config in feature_configs.items():
        try:
            df = pd.read_csv(config['path'])
            
            # Create label based on cluster name patterns
            def create_label(cluster_name):
                if str(cluster_name).startswith('BGC'):
                    return 'MIBiG'
                elif str(cluster_name).startswith('MGC_CANDIDATE'):
                    return 'KEGG_MGC'
                elif str(cluster_name).startswith('RANDOM'):
                    return 'KEGG_Random'
                else:
                    return None
            
            df['label'] = df[config['cluster_col']].apply(create_label)
            df = df.dropna(subset=['label'])
            
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
            
            # Merge
            if merged_df is None:
                merged_df = df
            else:
                merged_df = merged_df.merge(df, on=['cluster_name', 'label'], how='outer')
                
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
    
    # Create separate tables for each category
    categories = ['MIBiG', 'KEGG_MGC', 'KEGG_Random']
    
    for category in categories:
        print(f"\n" + "="*120)
        print(f"CREATING TABLE FOR {category.upper()}")
        print("="*120)
        
        # Filter data for this category
        category_data = merged_df[merged_df['label'] == category]
        
        if len(category_data) == 0:
            print(f"No data found for {category}")
            continue
        
        print(f"Number of clusters: {len(category_data)}")
        
        # Calculate statistics for each feature
        results = []
        
        for config in feature_configs.values():
            if 'prefixed_features' not in config:
                continue
                
            for feature, original_name in zip(config['prefixed_features'], config['features']):
                if feature not in category_data.columns:
                    continue
                
                feature_data = category_data[feature].dropna()
                
                if len(feature_data) == 0:
                    continue
                
                # Calculate statistics
                mean_val = feature_data.mean()
                std_val = feature_data.std()
                median_val = feature_data.median()
                min_val = feature_data.min()
                max_val = feature_data.max()
                q25 = feature_data.quantile(0.25)
                q75 = feature_data.quantile(0.75)
                
                results.append({
                    'Category': config['name'],
                    'Feature': original_name,
                    'Full_Feature_Name': feature,
                    'Mean': mean_val,
                    'Std': std_val,
                    'Median': median_val,
                    'Min': min_val,
                    'Max': max_val,
                    'Q25': q25,
                    'Q75': q75,
                    'Count': len(feature_data)
                })
        
        # Create DataFrame
        results_df = pd.DataFrame(results)
        
        # Display results by category
        print(f"\n📊 {category.upper()} FEATURE STATISTICS")
        print("-" * 120)
        print(f"{'Category':<18} {'Feature':<35} {'Mean':<12} {'Std':<12} {'Median':<10} {'Min':<10} {'Max':<10} {'Count':<8}")
        print("-" * 120)
        
        for _, row in results_df.iterrows():
            print(f"{row['Category']:<18} {row['Feature']:<35} {row['Mean']:<12.4f} {row['Std']:<12.4f} {row['Median']:<10.4f} {row['Min']:<10.4f} {row['Max']:<10.4f} {row['Count']:<8}")
        
        # Save to CSV
        output_file = f"/groups/itay_mayrose/alongonda/desktop/{category.lower()}_feature_statistics.csv"
        results_df.to_csv(output_file, index=False)
        print(f"\n💾 {category} statistics saved to: {output_file}")
        
        # Create markdown table for documentation
        print(f"\n📋 {category.upper()} MARKDOWN TABLE")
        print("-" * 100)
        print("| Category | Feature | Mean ± Std | Median | Min | Max | Count |")
        print("|----------|---------|------------|--------|-----|-----|-------|")
        
        for _, row in results_df.iterrows():
            print(f"| {row['Category']} | {row['Feature']} | {row['Mean']:.4f} ± {row['Std']:.4f} | {row['Median']:.4f} | {row['Min']:.4f} | {row['Max']:.4f} | {row['Count']} |")
    
    # Create comparison summary table
    print(f"\n" + "="*120)
    print("COMPARISON SUMMARY TABLE")
    print("="*120)
    
    comparison_results = []
    
    for config in feature_configs.values():
        if 'prefixed_features' not in config:
            continue
            
        for feature, original_name in zip(config['prefixed_features'], config['features']):
            if feature not in merged_df.columns:
                continue
            
            # Get means for each category (handle missing values)
            mibig_data = merged_df[merged_df['label'] == 'MIBiG'][feature].dropna()
            kegg_mgc_data = merged_df[merged_df['label'] == 'KEGG_MGC'][feature].dropna()
            random_data = merged_df[merged_df['label'] == 'KEGG_Random'][feature].dropna()
            
            mibig_mean = mibig_data.mean() if len(mibig_data) > 0 else np.nan
            kegg_mgc_mean = kegg_mgc_data.mean() if len(kegg_mgc_data) > 0 else np.nan
            random_mean = random_data.mean() if len(random_data) > 0 else np.nan
            
            comparison_results.append({
                'Category': config['name'],
                'Feature': original_name,
                'MIBiG_Mean': mibig_mean,
                'KEGG_MGC_Mean': kegg_mgc_mean,
                'KEGG_Random_Mean': random_mean,
                'MIBiG_vs_Random_Diff': mibig_mean - random_mean,
                'KEGG_MGC_vs_Random_Diff': kegg_mgc_mean - random_mean,
                'MIBiG_vs_KEGG_MGC_Diff': mibig_mean - kegg_mgc_mean
            })
    
    comparison_df = pd.DataFrame(comparison_results)
    
    print(f"{'Category':<18} {'Feature':<35} {'MIBiG':<12} {'KEGG-MGC':<12} {'KEGG-Random':<14} {'MIBiG vs R':<12} {'KEGG vs R':<12} {'MIBiG vs K':<12}")
    print("-" * 120)
    
    for _, row in comparison_df.iterrows():
        print(f"{row['Category']:<18} {row['Feature']:<35} {row['MIBiG_Mean']:<12.4f} {row['KEGG_MGC_Mean']:<12.4f} {row['KEGG_Random_Mean']:<14.4f} {row['MIBiG_vs_Random_Diff']:<12.4f} {row['KEGG_MGC_vs_Random_Diff']:<12.4f} {row['MIBiG_vs_KEGG_MGC_Diff']:<12.4f}")
    
    # Save comparison table
    comparison_output = "/groups/itay_mayrose/alongonda/desktop/comparison_summary_table.csv"
    comparison_df.to_csv(comparison_output, index=False)
    print(f"\n💾 Comparison summary saved to: {comparison_output}")
    
    return comparison_df

if __name__ == "__main__":
    comparison_df = create_separate_category_tables()
