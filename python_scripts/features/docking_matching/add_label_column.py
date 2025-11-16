#!/usr/bin/env python3
"""
Add label column to Lightdock CSV
"""

import pandas as pd

def add_label_column():
    """Add label column to Lightdock CSV based on cluster name"""
    
    # Load the Lightdock CSV
    df = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/lightdock_cluster_metrics.csv')
    
    print(f"Loaded {len(df)} clusters")
    
    # Create label column based on cluster name
    def create_label(name):
        if name.startswith('BGC'):
            return 'MIBiG'
        elif name.startswith('MGC_CANDIDATE'):
            return 'KEGG_MGC'
        elif name.startswith('RANDOM'):
            return 'KEGG_Random'
        else:
            return 'Unknown'
    
    df['category'] = df['name'].apply(create_label)
    
    # Show distribution
    print("\nCategory distribution:")
    print(df['category'].value_counts())
    
    # Save updated CSV
    output_file = '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/lightdock_cluster_metrics.csv'
    df.to_csv(output_file, index=False)
    
    print(f"\n✅ Updated Lightdock CSV with label column")
    print(f"   Saved to: {output_file}")

if __name__ == "__main__":
    add_label_column()
