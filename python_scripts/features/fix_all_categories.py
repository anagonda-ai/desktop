#!/usr/bin/env python3
"""
Fix category names in all CSV files to use 'RANDOM' instead of 'KEGG_Random'
"""

import pandas as pd

def fix_all_categories():
    """Fix category names in all CSV files"""
    
    files_to_fix = [
        '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/lightdock_cluster_metrics.csv',
        '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/cladepp_cluster_metrics_enhanced.csv',
        '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/promoter_similarity_results.csv',
        '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/e2p2_feature_extraction_detailed.csv'
    ]
    
    for file_path in files_to_fix:
        try:
            print(f"\nProcessing: {file_path}")
            df = pd.read_csv(file_path)
            
            if 'category' in df.columns:
                print("Current category distribution:")
                print(df['category'].value_counts())
                
                # Fix category names
                df['category'] = df['category'].replace({
                    'KEGG_Random': 'RANDOM',
                    'KEGG_MGC': 'MGC',
                    'MIBiG': 'MIBiG'
                })
                
                print("Updated category distribution:")
                print(df['category'].value_counts())
                
                # Save the updated CSV
                df.to_csv(file_path, index=False)
                print(f"✅ Updated {file_path}")
            else:
                print(f"❌ No 'category' column found in {file_path}")
                
        except Exception as e:
            print(f"❌ Error processing {file_path}: {e}")
    
    print("\n✅ All category names fixed!")

if __name__ == "__main__":
    fix_all_categories()

