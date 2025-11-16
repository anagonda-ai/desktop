#!/usr/bin/env python3
"""
Extract Enhanced Features for CladePP Only

This script computes additional features based on clade-level information.
(Foldseek processing has been commented out)
"""

import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

class EnhancedFeatureExtractor:
    def __init__(self):
        # Base paths
        self.mgc_base = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
        self.mgc_dir = self.mgc_base / "mgc_candidates_dir_final"
        self.random_dir = self.mgc_base / "kegg_random_mgc_candidates_dir_fixed"
        
        # Existing metrics
        # self.foldseek_metrics = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/foldseek_cluster_metrics_enhanced.csv')
        self.cladepp_metrics = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/cladepp_cluster_metrics_enhanced.csv')
        
        print("="*80)
        print("🔬 ENHANCED FEATURE EXTRACTION - CladePP Only")
        print("="*80)
        
    # def compute_foldseek_enhanced_features(self):
    #     """
    #     Compute enhanced Foldseek features:
    #     1. Match coverage ratio
    #     2. Core module fraction
    #     3. Network density
    #     4. Multi-match fraction
    #     """
    #     print("\n📊 Computing Enhanced Foldseek Features...")
    #     print("   Note: Using fraction_strong_binders as proxy for match coverage")
    #     print("   (Actual pairwise data would require re-processing raw Foldseek output)")
    #     
    #     enhanced_features = []
    #     
    #     for idx, row in self.foldseek_metrics.iterrows():
    #         cluster_name = row['name']
    #         n_proteins = row['n']
    #         
    #         # Feature 1: Match Coverage Ratio
    #         # Use fraction_strong_binders as proxy (proteins with strong matches)
    #         match_coverage = row.get('fraction_strong_binders', 0.0)
    #         
    #         # Feature 2: Core Module Fraction
    #         # Approximate using: if high enrichment + high mean score → large core
    #         # Real implementation would need pairwise matrix
    #         mean_score = row['mean_score_non_self']
    #         enrichment = row['enrichment_score']
    #         
    #         if mean_score > 0.8 and enrichment > 100:
    #             core_module_fraction = 0.9  # Very tight module
    #         elif mean_score > 0.6 and enrichment > 10:
    #             core_module_fraction = 0.7  # Moderate module
    #         elif mean_score > 0.4:
    #             core_module_fraction = 0.5  # Loose module
    #         else:
    #             core_module_fraction = 0.3  # Weak module
    #         
    #         # Feature 3: Network Density
    #         # Approximate: strong_binders^2 (connectivity scales with match fraction)
    #         network_density = match_coverage ** 2
    #         
    #         # Feature 4: Multi-Match Fraction
    #         # Use moderate + strong binders as proxy
    #         multi_match_fraction = row.get('fraction_strong_binders', 0) + \
    #                              row.get('fraction_moderate_binders', 0)
    #         
    #         enhanced_features.append({
    #             'cluster_name': cluster_name,
    #             'foldseek_match_coverage': match_coverage,
    #             'foldseek_core_module_fraction': core_module_fraction,
    #             'foldseek_network_density': network_density,
    #             'foldseek_multi_match_fraction': min(multi_match_fraction, 1.0)
    #         })
    #     
    #     df_foldseek = pd.DataFrame(enhanced_features)
    #     print(f"   ✅ Computed enhanced features for {len(df_foldseek)} clusters")
    #     return df_foldseek
    
    def load_cladepp_summary(self, cluster_dir):
        """Load CladePP summary.csv for a cluster"""
        summary_file = cluster_dir / "summary.csv"
        if not summary_file.exists():
            return None
        
        try:
            df = pd.read_csv(summary_file)
            return df
        except Exception as e:
            return None
    
    def compute_cladepp_enhanced_features_for_cluster(self, cluster_name, cluster_dir):
        """
        Compute enhanced CladePP features for a single cluster:
        1. Largest high-scoring clade fraction
        2. Multi-clade conservation (high/medium/low thresholds)
        3. Conservation consistency
        4. Max pair co-evolution score
        """
        summary_df = self.load_cladepp_summary(cluster_dir)
        
        default_result = {
            'cluster_name': cluster_name,
            'cladepp_largest_clade_fraction': 0.0,
            'cladepp_multi_clade_high': 0.0,
            'cladepp_multi_clade_medium': 0.0,
            'cladepp_multi_clade_low': 0.0,
            'cladepp_conservation_consistency': 0.0,
            'cladepp_max_pair_score': 0.0
        }
        
        if summary_df is None or len(summary_df) == 0:
            return default_result
        
        # Filter only successful rows
        if 'status' in summary_df.columns:
            summary_df = summary_df[summary_df['status'] == 'success'].copy()
        
        if len(summary_df) == 0:
            return default_result
        
        # Use cladepp_score column (actual co-evolution score)
        if 'cladepp_score' not in summary_df.columns:
            return default_result
        
        scores = summary_df['cladepp_score'].dropna()
        
        if len(scores) == 0:
            return default_result
        
        # Feature 1: Largest high-scoring clade fraction
        if 'clade_size' in summary_df.columns:
            # Filter high-scoring clades (score > 0.7)
            high_scoring_df = summary_df[summary_df['cladepp_score'] > 0.7]
            
            if len(high_scoring_df) > 0:
                largest_size = high_scoring_df['clade_size'].max()
                total_size = summary_df['clade_size'].sum()
                largest_clade_fraction = largest_size / total_size if total_size > 0 else 0.0
            else:
                largest_clade_fraction = 0.0
        else:
            # No size info: use fraction of clades with high scores
            largest_clade_fraction = (scores > 0.7).mean()
        
        # Feature 2-4: Multi-clade conservation at different thresholds
        multi_clade_high = (scores > 0.8).mean()
        multi_clade_medium = (scores > 0.6).mean()
        multi_clade_low = (scores > 0.4).mean()
        
        # Feature 5: Conservation consistency (inverse of coefficient of variation)
        mean_score = scores.mean()
        std_score = scores.std()
        if mean_score > 0:
            consistency = 1 - (std_score / mean_score)
            consistency = max(0, min(1, consistency))  # Clip to [0, 1]
        else:
            consistency = 0.0
        
        # Feature 6: Max pair co-evolution score
        max_pair_score = scores.max()
        
        return {
            'cluster_name': cluster_name,
            'cladepp_largest_clade_fraction': largest_clade_fraction,
            'cladepp_multi_clade_high': multi_clade_high,
            'cladepp_multi_clade_medium': multi_clade_medium,
            'cladepp_multi_clade_low': multi_clade_low,
            'cladepp_conservation_consistency': consistency,
            'cladepp_max_pair_score': max_pair_score
        }
    
    def compute_cladepp_enhanced_features(self):
        """Compute enhanced CladePP features for all clusters"""
        print("\n🌳 Computing Enhanced CladePP Features...")
        print(f"   MGC directory: {self.mgc_dir}")
        print(f"   Random directory: {self.random_dir}")
        
        enhanced_features = []
        
        # Process MGC clusters
        if self.mgc_dir.exists():
            cluster_dirs = [d for d in self.mgc_dir.iterdir() if d.is_dir()]
            print(f"   Processing {len(cluster_dirs)} MGC clusters...")
            
            for cluster_dir in cluster_dirs:
                cluster_name = cluster_dir.name
                features = self.compute_cladepp_enhanced_features_for_cluster(cluster_name, cluster_dir)
                enhanced_features.append(features)
        
        # Process Random clusters
        if self.random_dir.exists():
            cluster_dirs = [d for d in self.random_dir.iterdir() if d.is_dir()]
            print(f"   Processing {len(cluster_dirs)} Random clusters...")
            
            for cluster_dir in cluster_dirs:
                cluster_name = cluster_dir.name
                features = self.compute_cladepp_enhanced_features_for_cluster(cluster_name, cluster_dir)
                enhanced_features.append(features)
        
        df_cladepp = pd.DataFrame(enhanced_features)
        print(f"   ✅ Computed enhanced features for {len(df_cladepp)} clusters")
        return df_cladepp
    
    def merge_with_existing(self, df_enhanced_cladepp):
        """Merge enhanced features with existing metrics"""
        print("\n🔗 Merging enhanced features with existing metrics...")
        
        # # Merge with existing Foldseek
        # foldseek_merged = pd.merge(
        #     self.foldseek_metrics,
        #     df_enhanced_foldseek,
        #     left_on='name',
        #     right_on='cluster_name',
        #     how='left'
        # )
        
        # Merge with existing CladePP
        cladepp_merged = pd.merge(
            self.cladepp_metrics,
            df_enhanced_cladepp,
            left_on='name',
            right_on='cluster_name',
            how='left'
        )
        
        # # Fill NaN values with 0
        # for col in df_enhanced_foldseek.columns:
        #     if col != 'cluster_name' and col in foldseek_merged.columns:
        #         foldseek_merged[col] = foldseek_merged[col].fillna(0)
        
        for col in df_enhanced_cladepp.columns:
            if col != 'cluster_name' and col in cladepp_merged.columns:
                cladepp_merged[col] = cladepp_merged[col].fillna(0)
        
        # print(f"   ✅ Foldseek: {len(foldseek_merged)} rows, {len(foldseek_merged.columns)} columns")
        print(f"   ✅ CladePP: {len(cladepp_merged)} rows, {len(cladepp_merged.columns)} columns")
        
        return cladepp_merged
    
    def save_results(self, df_cladepp):
        """Save enhanced feature sets"""
        output_dir = Path('/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/')
        
        # foldseek_out = output_dir / 'foldseek_cluster_metrics_enhanced.csv'
        cladepp_out = output_dir / 'cladepp_cluster_metrics_enhanced.csv'
        
        # df_foldseek.to_csv(foldseek_out, index=False)
        df_cladepp.to_csv(cladepp_out, index=False)
        
        print("\n💾 Saved enhanced features:")
        # print(f"   📊 Foldseek: {foldseek_out}")
        print(f"   🌳 CladePP: {cladepp_out}")
        
        # Print feature counts
        # foldseek_new_features = [col for col in df_foldseek.columns if col.startswith('foldseek_') and col not in self.foldseek_metrics.columns]
        cladepp_new_features = [col for col in df_cladepp.columns if col.startswith('cladepp_') and col not in self.cladepp_metrics.columns]
        
        print(f"\n📈 New Features Added:")
        # print(f"   Foldseek: {len(foldseek_new_features)} new features")
        # for feat in foldseek_new_features:
        #     print(f"      • {feat}")
        print(f"   CladePP: {len(cladepp_new_features)} new features")
        for feat in cladepp_new_features:
            print(f"      • {feat}")
        
        # total_old = len(self.foldseek_metrics.columns) - 1 + len(self.cladepp_metrics.columns) - 1  # Subtract cluster name cols
        total_old = len(self.cladepp_metrics.columns) - 1  # Subtract cluster name cols
        # total_new = len(foldseek_new_features) + len(cladepp_new_features)
        total_new = len(cladepp_new_features)
        
        print(f"\n🎯 Total Features: {total_old} → {total_old + total_new} (+{total_new})")
        
        return cladepp_out

def main():
    print("="*80)
    print("🚀 ENHANCED FEATURE EXTRACTION - CladePP Only")
    print("="*80)
    print()
    print("This script extracts additional features:")
    # print("  📊 Foldseek: Coverage, core module, network topology")
    print("  🌳 CladePP: Clade size, multi-clade conservation, consistency")
    print()
    
    extractor = EnhancedFeatureExtractor()
    
    # # Extract enhanced Foldseek features
    # df_foldseek_enhanced = extractor.compute_foldseek_enhanced_features()
    
    # Extract enhanced CladePP features
    df_cladepp_enhanced = extractor.compute_cladepp_enhanced_features()
    
    # Merge with existing
    df_cladepp_merged = extractor.merge_with_existing(
        df_cladepp_enhanced
    )
    
    # Save results
    cladepp_file = extractor.save_results(
        df_cladepp_merged
    )
    
    print("\n" + "="*80)
    print("✅ ENHANCED FEATURE EXTRACTION COMPLETE!")
    print("="*80)
    print()
    print("Next steps:")
    print("  1. Update comprehensive_incremental_analysis.py to use enhanced files")
    print("  2. Re-run the incremental analysis with new features")
    print("  3. Compare performance with and without enhanced features")

if __name__ == "__main__":
    main()

