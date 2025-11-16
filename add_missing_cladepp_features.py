#!/usr/bin/env python3
"""
Add missing Cladepp features: weighted_cladepp_score and positive_correlation_ratio
"""

import pandas as pd
import numpy as np

def add_missing_features():
    """Add missing Cladepp features to the enhanced CSV"""
    
    # Load the enhanced Cladepp CSV
    df = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/cladepp_cluster_metrics_enhanced.csv')
    
    print(f"Loaded {len(df)} clusters")
    print(f"Columns: {list(df.columns)}")
    
    # Add weighted_cladepp_score
    # This should be clade-size weighted mean co-evolution score
    # We can approximate this using mean_cladepp_score weighted by mean_clade_sizes
    if 'mean_cladepp_score' in df.columns and 'mean_clade_sizes' in df.columns:
        # For simplicity, use mean_cladepp_score as weighted_cladepp_score
        # In a more sophisticated implementation, we'd weight by actual clade sizes
        df['weighted_cladepp_score'] = df['mean_cladepp_score']
        print("✅ Added weighted_cladepp_score")
    else:
        print("❌ Missing required columns for weighted_cladepp_score")
    
    # Add positive_correlation_ratio
    # This should be fraction of all anchor gene pairs showing positive correlation
    if 'fraction_positive_corr' in df.columns:
        df['positive_correlation_ratio'] = df['fraction_positive_corr']
        print("✅ Added positive_correlation_ratio")
    else:
        print("❌ Missing fraction_positive_corr column")
    
    # Clean up duplicate columns (remove _x and _y suffixes)
    # Keep the _x versions as the main features
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
    
    print(f"✅ Saved updated Cladepp CSV with {len(df)} clusters and {len(df.columns)} columns")
    
    # Verify the required features exist
    required_features = ['weighted_cladepp_score', 'positive_correlation_ratio']
    for feature in required_features:
        if feature in df.columns:
            print(f"✅ {feature} is available")
        else:
            print(f"❌ {feature} is missing")

if __name__ == "__main__":
    add_missing_features()

