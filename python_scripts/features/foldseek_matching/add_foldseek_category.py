#!/usr/bin/env python3
"""
Add category column to Foldseek CSV
"""

import pandas as pd

def add_category_to_foldseek():
    """Add category column to Foldseek CSV"""
    print("📊 Adding category column to Foldseek CSV...")
    
    # Load the CSV file
    df = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/foldseek_cluster_metrics_enhanced.csv')
    
    print(f"Loaded {len(df)} rows")
    
    # Create category column
    def create_category(name):
        if name.startswith('BGC'):
            return 'MIBiG'
        elif name.startswith('MGC_CANDIDATE'):
            return 'KEGG_MGC'
        elif name.startswith('RANDOM_MGC'):
            return 'KEGG_Random'
        else:
            return 'Unknown'
    
    df['category'] = df['name'].apply(create_category)
    
    # Show distribution
    print("\nCategory distribution:")
    print(df['category'].value_counts())
    
    # Save the updated CSV
    output_file = '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/foldseek_cluster_metrics_enhanced.csv'
    df.to_csv(output_file, index=False)
    
    print(f"\n✅ Updated Foldseek CSV saved to: {output_file}")

if __name__ == "__main__":
    add_category_to_foldseek()

