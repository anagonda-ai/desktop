#!/usr/bin/env python3
"""
Fix category names in Foldseek CSV to use 'RANDOM' instead of 'KEGG_Random'
"""

import pandas as pd

def fix_categories():
    """Fix category names in Foldseek CSV"""
    
    # Load the Foldseek CSV
    df = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/foldseek_cluster_metrics_enhanced.csv')
    
    print(f"Loaded {len(df)} clusters")
    print("Current category distribution:")
    print(df['category'].value_counts())
    
    # Fix category names
    df['category'] = df['category'].replace({
        'KEGG_Random': 'RANDOM',
        'KEGG_MGC': 'MGC',
        'MIBiG': 'MIBiG'
    })
    
    print("\nUpdated category distribution:")
    print(df['category'].value_counts())
    
    # Save the updated CSV
    output_file = '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/foldseek_cluster_metrics_enhanced.csv'
    df.to_csv(output_file, index=False)
    
    print(f"✅ Updated Foldseek CSV with corrected category names")

if __name__ == "__main__":
    fix_categories()

