#!/usr/bin/env python3
"""
Recreate CSV Files from Existing Results

This script recreates the feature metrics CSV files from the existing analysis results
that are already present in the directories.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import glob

def recreate_lightdock_csv():
    """Recreate Lightdock CSV from existing results"""
    print("⚡ RECREATING LIGHTDOCK CSV FROM EXISTING RESULTS...")
    
    results_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/lightdock_results")
    
    all_metrics = []
    
    # Process all cluster directories
    for cluster_dir in results_dir.iterdir():
        if not cluster_dir.is_dir():
            continue
            
        cluster_name = cluster_dir.name
        
        # Look for results_tsv directory
        results_tsv_dir = cluster_dir / "results_tsv"
        if not results_tsv_dir.exists():
            continue
            
        # Look for the main results file
        results_file = results_tsv_dir / f"{cluster_name}_all_vs_all.tsv"
        if not results_file.exists():
            continue
            
        try:
            # Read the results file
            df = pd.read_csv(results_file, sep='\t', header=None)
            
            if len(df) == 0:
                continue
                
            # Calculate metrics
            scores = df.iloc[:, 10] if len(df.columns) > 10 else df.iloc[:, -1]  # Score column
            
            metrics = {
                'name': cluster_name,
                'total_interactions': len(df),
                'mean_score': float(scores.mean()) if len(scores) > 0 else 0.0,
                'std_score': float(scores.std()) if len(scores) > 0 else 0.0,
                'min_score': float(scores.min()) if len(scores) > 0 else 0.0,
                'max_score': float(scores.max()) if len(scores) > 0 else 0.0,
                'median_score': float(scores.median()) if len(scores) > 0 else 0.0,
                'high_score_interactions': int((scores < -50).sum()) if len(scores) > 0 else 0,
                'medium_score_interactions': int(((scores >= -50) & (scores < -25)).sum()) if len(scores) > 0 else 0,
                'low_score_interactions': int((scores >= -25).sum()) if len(scores) > 0 else 0
            }
            
            all_metrics.append(metrics)
            
        except Exception as e:
            print(f"Error processing {cluster_name}: {e}")
            continue
    
    # Create DataFrame and save
    df_metrics = pd.DataFrame(all_metrics)
    output_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/lightdock_cluster_metrics.csv"
    df_metrics.to_csv(output_file, index=False)
    
    print(f"✅ Created Lightdock CSV with {len(df_metrics)} clusters")
    print(f"   Saved to: {output_file}")
    
    # Show MIBiG clusters
    mibig_clusters = df_metrics[df_metrics['name'].str.startswith('BGC')]
    print(f"   MIBiG clusters: {len(mibig_clusters)}")
    
    return df_metrics

def recreate_foldseek_csv():
    """Recreate Foldseek CSV from existing results"""
    print("📊 RECREATING FOLDSEEK CSV FROM EXISTING RESULTS...")
    
    # Check if enhanced features exist
    enhanced_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/foldseek_cluster_metrics_enhanced.csv"
    
    if os.path.exists(enhanced_file):
        df = pd.read_csv(enhanced_file)
        print(f"✅ Using existing Foldseek enhanced CSV with {len(df)} clusters")
        
        # Show MIBiG clusters
        mibig_clusters = df[df['name'].str.startswith('BGC')]
        print(f"   MIBiG clusters: {len(mibig_clusters)}")
        
        return df
    else:
        print("❌ No existing Foldseek results found")
        return None

def recreate_cladepp_csv():
    """Recreate Cladepp CSV from existing results"""
    print("🌳 RECREATING CLADEPP CSV FROM EXISTING RESULTS...")
    
    # Check if enhanced features exist
    enhanced_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/cladepp_cluster_metrics_enhanced.csv"
    
    if os.path.exists(enhanced_file):
        df = pd.read_csv(enhanced_file)
        print(f"✅ Using existing Cladepp enhanced CSV with {len(df)} clusters")
        
        # Show MIBiG clusters
        mibig_clusters = df[df['cluster_name'].str.startswith('BGC')]
        print(f"   MIBiG clusters: {len(mibig_clusters)}")
        
        return df
    else:
        print("❌ No existing Cladepp results found")
        return None

def recreate_promoter_csv():
    """Recreate Promoter CSV from existing results"""
    print("🧬 RECREATING PROMOTER CSV FROM EXISTING RESULTS...")
    
    # Check if results exist
    results_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/promoter_similarity_results.csv"
    
    if os.path.exists(results_file):
        df = pd.read_csv(results_file)
        print(f"✅ Using existing Promoter CSV with {len(df)} clusters")
        
        # Show MIBiG clusters
        mibig_clusters = df[df['group_name'].str.startswith('BGC')]
        print(f"   MIBiG clusters: {len(mibig_clusters)}")
        
        return df
    else:
        print("❌ No existing Promoter results found")
        return None

def recreate_e2p2_csv():
    """Recreate E2P2 CSV from existing results"""
    print("🧪 RECREATING E2P2 CSV FROM EXISTING RESULTS...")
    
    # Check if results exist
    results_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/e2p2_feature_extraction_detailed.csv"
    
    if os.path.exists(results_file):
        df = pd.read_csv(results_file)
        print(f"✅ Using existing E2P2 CSV with {len(df)} clusters")
        
        # Show MIBiG clusters
        mibig_clusters = df[df['cluster_name'].str.startswith('BGC')]
        print(f"   MIBiG clusters: {len(mibig_clusters)}")
        
        return df
    else:
        print("❌ No existing E2P2 results found")
        return None

def main():
    print("="*80)
    print("🔄 RECREATING CSV FILES FROM EXISTING RESULTS")
    print("="*80)
    print()
    print("This script will recreate the feature metrics CSV files from")
    print("the existing analysis results that are already present.")
    print()
    
    # Recreate each CSV file
    lightdock_df = recreate_lightdock_csv()
    print()
    
    foldseek_df = recreate_foldseek_csv()
    print()
    
    cladepp_df = recreate_cladepp_csv()
    print()
    
    promoter_df = recreate_promoter_csv()
    print()
    
    e2p2_df = recreate_e2p2_csv()
    print()
    
    print("="*80)
    print("📋 SUMMARY OF RECREATED CSV FILES:")
    print("="*80)
    print()
    
    if lightdock_df is not None:
        print(f"✅ Lightdock: {len(lightdock_df)} clusters")
    else:
        print("❌ Lightdock: Failed")
        
    if foldseek_df is not None:
        print(f"✅ Foldseek: {len(foldseek_df)} clusters")
    else:
        print("❌ Foldseek: Failed")
        
    if cladepp_df is not None:
        print(f"✅ Cladepp: {len(cladepp_df)} clusters")
    else:
        print("❌ Cladepp: Failed")
        
    if promoter_df is not None:
        print(f"✅ Promoter: {len(promoter_df)} clusters")
    else:
        print("❌ Promoter: Failed")
        
    if e2p2_df is not None:
        print(f"✅ E2P2: {len(e2p2_df)} clusters")
    else:
        print("❌ E2P2: Failed")
    
    print()
    print("🎉 CSV RECREATION COMPLETE!")
    print("Now you can run the comprehensive analysis with the recreated data.")

if __name__ == "__main__":
    main()

