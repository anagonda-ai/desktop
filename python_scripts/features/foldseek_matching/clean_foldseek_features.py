#!/usr/bin/env python3
"""
Clean up Foldseek enhanced features by removing duplicate columns and renaming _x columns
"""

import pandas as pd
import numpy as np

def clean_foldseek_features():
    """Clean up Foldseek enhanced features"""
    
    # Load the enhanced Foldseek CSV
    df = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/foldseek_cluster_metrics_enhanced.csv')
    
    print(f"Loaded {len(df)} clusters")
    print(f"Original columns: {len(df.columns)}")
    
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
    print(f"Final columns: {len(df.columns)}")
    
    # Save the updated CSV
    output_file = '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/foldseek_cluster_metrics_enhanced.csv'
    df.to_csv(output_file, index=False)
    
    print(f"✅ Saved updated Foldseek CSV with {len(df)} clusters and {len(df.columns)} columns")
    
    # Show the final column names
    print(f"Final column names: {list(df.columns)}")

if __name__ == "__main__":
    clean_foldseek_features()

