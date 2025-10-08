#!/usr/bin/env python3
"""
Runner script for E2P2 feature extraction and statistical analysis
"""

import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.append(str(Path(__file__).parent))

from e2p2_feature_extractor import E2P2FeatureExtractor
from e2p2_statistical_analysis import main as run_statistical_analysis

def main():
    """
    Complete E2P2 analysis pipeline:
    1. Extract features from E2P2 results
    2. Perform statistical analysis and threshold optimization
    """
    print("E2P2 Complete Analysis Pipeline")
    print("=" * 50)
    
    # Step 1: Feature extraction
    print("\nStep 1: Extracting E2P2 features...")
    print("-" * 40)
    
    # Initialize extractor
    base_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    extractor = E2P2FeatureExtractor(base_dir=base_dir)
    
    # Process all clusters
    results_df = extractor.process_all_clusters()
    
    if results_df is None:
        print("No E2P2 results found. Please ensure E2P2 analysis has been completed.")
        return
    
    print(f"\nFeature extraction completed:")
    print(f"  Total clusters processed: {len(results_df)}")
    print(f"  Results saved to: {extractor.output_dir}")
    
    # Step 2: Statistical analysis
    print(f"\nStep 2: Performing statistical analysis...")
    print("-" * 40)
    
    try:
        # Run statistical analysis
        df, feature_stats, group_results, threshold_results = run_statistical_analysis()
        
        print(f"\nStatistical analysis completed successfully!")
        print(f"  Key findings:")
        
        if threshold_results:
            # Find best performing feature
            best_feature = None
            best_auc = 0
            
            for feature, results in threshold_results.items():
                if results['robust'] and 'test_roc_auc' in results['robust']:
                    if results['robust']['test_roc_auc'] > best_auc:
                        best_auc = results['robust']['test_roc_auc']
                        best_feature = feature
            
            if best_feature:
                print(f"    Best performing feature: {best_feature}")
                print(f"    Test ROC AUC: {best_auc:.4f}")
                
                # Get recommended threshold
                if 'f1' in threshold_results[best_feature]['robust']['test_results']:
                    best_threshold = threshold_results[best_feature]['robust']['test_results']['f1']['threshold']
                    best_f1 = threshold_results[best_feature]['robust']['test_results']['f1']['f1_score']
                    print(f"    Recommended threshold: {best_threshold:.4f}")
                    print(f"    Expected F1 score: {best_f1:.4f}")
        
        print(f"\nAnalysis complete! Check the output files in:")
        print(f"  {extractor.output_dir}")
        
    except Exception as e:
        print(f"Error in statistical analysis: {e}")
        print("Feature extraction completed, but statistical analysis failed.")
        print("You can run the statistical analysis separately using:")
        print("  python e2p2_statistical_analysis.py")

if __name__ == "__main__":
    main()
