#!/usr/bin/env python3
"""
Create comprehensive table of average E2P2 feature values for MGC vs KEGG-Random groups
"""

import pandas as pd
import numpy as np
from pathlib import Path

def create_averages_table():
    """Create comprehensive table of average values for each group"""
    
    # Load data
    feature_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/e2p2_feature_extraction_detailed.csv"
    df = pd.read_csv(feature_file)
    
    print(f"Loaded {len(df)} clusters from file")
    print(f"  MGC clusters: {df['classification_label'].sum()}")
    print(f"  KEGG-Random clusters: {(df['classification_label']==0).sum()}")
    
    # Define features
    features = [
        'num_distinct_enzyme_classes',
        'num_distinct_enzyme_subclasses', 
        'num_distinct_enzyme_families',
        'num_distinct_enzyme_subfamilies',
        'total_ec_numbers'
    ]
    
    feature_names = [
        'Enzyme Classes (X.-.-.-)',
        'Enzyme Subclasses (X.Y.-.-)',
        'Enzyme Families (X.Y.Z.-)',
        'Enzyme Subfamilies (X.Y.Z.W)',
        'Total EC Numbers'
    ]
    
    # Calculate statistics for each group
    results = []
    
    for feature, feature_name in zip(features, feature_names):
        # MGC group (classification_label = 1)
        mgc_data = df[df['classification_label'] == 1][feature]
        mgc_mean = mgc_data.mean()
        mgc_std = mgc_data.std()
        mgc_median = mgc_data.median()
        mgc_min = mgc_data.min()
        mgc_max = mgc_data.max()
        
        # KEGG-Random group (classification_label = 0)
        random_data = df[df['classification_label'] == 0][feature]
        random_mean = random_data.mean()
        random_std = random_data.std()
        random_median = random_data.median()
        random_min = random_data.min()
        random_max = random_data.max()
        
        # Calculate fold change and difference
        fold_change = mgc_mean / random_mean if random_mean != 0 else np.nan
        difference = mgc_mean - random_mean
        percent_change = ((mgc_mean - random_mean) / random_mean) * 100 if random_mean != 0 else np.nan
        
        # Determine direction
        if mgc_mean < random_mean:
            direction = "LOWER in MGCs"
            fold_reduction = random_mean / mgc_mean if mgc_mean != 0 else np.nan
        else:
            direction = "HIGHER in MGCs"
            fold_reduction = np.nan
        
        results.append({
            'Feature': feature_name,
            'MGC_Mean': mgc_mean,
            'MGC_Std': mgc_std,
            'MGC_Median': mgc_median,
            'MGC_Min': mgc_min,
            'MGC_Max': mgc_max,
            'KEGG_Random_Mean': random_mean,
            'KEGG_Random_Std': random_std,
            'KEGG_Random_Median': random_median,
            'KEGG_Random_Min': random_min,
            'KEGG_Random_Max': random_max,
            'Difference': difference,
            'Percent_Change': percent_change,
            'Fold_Change': fold_change,
            'Fold_Reduction': fold_reduction,
            'Direction': direction
        })
    
    # Create DataFrame
    results_df = pd.DataFrame(results)
    
    # Display results
    print("\n" + "="*120)
    print("E2P2 FEATURE AVERAGES: MGC vs KEGG-RANDOM COMPARISON")
    print("="*120)
    
    # Summary table
    print("\n📊 SUMMARY TABLE - AVERAGE VALUES")
    print("-" * 120)
    print(f"{'Feature':<35} {'MGC Mean':<12} {'KEGG-Random Mean':<18} {'Difference':<12} {'% Change':<10} {'Direction':<20}")
    print("-" * 120)
    
    for _, row in results_df.iterrows():
        print(f"{row['Feature']:<35} {row['MGC_Mean']:<12.4f} {row['KEGG_Random_Mean']:<18.4f} {row['Difference']:<12.4f} {row['Percent_Change']:<10.1f} {row['Direction']:<20}")
    
    # Detailed statistics table
    print("\n\n📈 DETAILED STATISTICS TABLE")
    print("-" * 140)
    print(f"{'Feature':<35} {'Group':<12} {'Mean':<8} {'Std':<8} {'Median':<8} {'Min':<6} {'Max':<6}")
    print("-" * 140)
    
    for _, row in results_df.iterrows():
        print(f"{row['Feature']:<35} {'MGC':<12} {row['MGC_Mean']:<8.4f} {row['MGC_Std']:<8.4f} {row['MGC_Median']:<8.1f} {row['MGC_Min']:<6.0f} {row['MGC_Max']:<6.0f}")
        print(f"{'':<35} {'KEGG-Random':<12} {row['KEGG_Random_Mean']:<8.4f} {row['KEGG_Random_Std']:<8.4f} {row['KEGG_Random_Median']:<8.1f} {row['KEGG_Random_Min']:<6.0f} {row['KEGG_Random_Max']:<6.0f}")
        print("-" * 140)
    
    # Fold change analysis
    print("\n\n🔄 FOLD CHANGE ANALYSIS")
    print("-" * 80)
    print(f"{'Feature':<35} {'Fold Change':<12} {'Fold Reduction':<15} {'Interpretation':<20}")
    print("-" * 80)
    
    for _, row in results_df.iterrows():
        if pd.notna(row['Fold_Reduction']):
            print(f"{row['Feature']:<35} {row['Fold_Change']:<12.3f} {row['Fold_Reduction']:<15.3f} {'MGCs are specialized':<20}")
        else:
            print(f"{row['Feature']:<35} {row['Fold_Change']:<12.3f} {'N/A':<15} {'MGCs have more enzymes':<20}")
    
    # Save to CSV
    output_file = "/groups/itay_mayrose/alongonda/desktop/e2p2_feature_averages_table.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n💾 Detailed results saved to: {output_file}")
    
    # Create markdown table for documentation
    print("\n\n📋 MARKDOWN TABLE FOR DOCUMENTATION")
    print("-" * 100)
    print("| Feature | MGC Mean ± Std | KEGG-Random Mean ± Std | Difference | % Change | Direction |")
    print("|---------|----------------|------------------------|------------|----------|-----------|")
    
    for _, row in results_df.iterrows():
        print(f"| {row['Feature']} | {row['MGC_Mean']:.4f} ± {row['MGC_Std']:.4f} | {row['KEGG_Random_Mean']:.4f} ± {row['KEGG_Random_Std']:.4f} | {row['Difference']:+.4f} | {row['Percent_Change']:+.1f}% | {row['Direction']} |")
    
    return results_df

if __name__ == "__main__":
    results_df = create_averages_table()
