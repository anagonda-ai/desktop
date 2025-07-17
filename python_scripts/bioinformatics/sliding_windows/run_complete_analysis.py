#!/usr/bin/env python3
"""
Complete MGC Analysis Pipeline

This script runs the complete analysis pipeline for KEGG-MGC research,
including genome statistics extraction, advanced visualization, and
comprehensive statistical analysis.
"""

import os
import sys
import time
from pathlib import Path

# Add the current directory to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from genome_statistics_extractor import GenomeStatisticsExtractor
from advanced_mgc_visualization import AdvancedMGCVisualizer
from statistical_analysis import MGCStatisticalAnalyzer


def main():
    """Main function to run complete analysis pipeline."""
    
    print("🧬 KEGG-MGC COMPREHENSIVE ANALYSIS PIPELINE")
    print("=" * 60)
    
    # Define paths
    base_dir = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3'
    mgc_results_file = os.path.join(base_dir, 'potential_groups_w10_filtered.csv')
    genome_stats_file = os.path.join(base_dir, 'genome_statistics.csv')
    
    # Check if MGC results file exists
    if not os.path.exists(mgc_results_file):
        print(f"❌ MGC results file not found: {mgc_results_file}")
        return
    
    print(f"📁 Working directory: {base_dir}")
    print(f"📊 MGC results file: {os.path.basename(mgc_results_file)}")
    
    # Step 1: Extract genome statistics
    print("\n" + "=" * 60)
    print("STEP 1: GENOME STATISTICS EXTRACTION")
    print("=" * 60)
    
    start_time = time.time()
    
    extractor = GenomeStatisticsExtractor(mgc_results_file)
    stats_file = extractor.save_statistics()
    
    step1_time = time.time() - start_time
    print(f"⏱️  Step 1 completed in {step1_time:.2f} seconds")
    
    # Step 2: Advanced visualization
    print("\n" + "=" * 60)
    print("STEP 2: ADVANCED VISUALIZATION")
    print("=" * 60)
    
    start_time = time.time()
    
    visualizer = AdvancedMGCVisualizer(genome_stats_file, mgc_results_file)
    visualization_files = visualizer.run_full_analysis()
    
    step2_time = time.time() - start_time
    print(f"⏱️  Step 2 completed in {step2_time:.2f} seconds")
    
    # Step 3: Statistical analysis
    print("\n" + "=" * 60)
    print("STEP 3: STATISTICAL ANALYSIS")
    print("=" * 60)
    
    start_time = time.time()
    
    analyzer = MGCStatisticalAnalyzer(genome_stats_file, mgc_results_file)
    statistical_results = analyzer.run_complete_analysis()
    
    step3_time = time.time() - start_time
    print(f"⏱️  Step 3 completed in {step3_time:.2f} seconds")
    
    # Summary
    print("\n" + "=" * 60)
    print("ANALYSIS PIPELINE COMPLETE")
    print("=" * 60)
    
    total_time = step1_time + step2_time + step3_time
    print(f"⏱️  Total analysis time: {total_time:.2f} seconds")
    
    print("\n📋 GENERATED FILES:")
    print(f"   • Genome statistics: {os.path.basename(stats_file)}")
    
    for viz_file in visualization_files:
        print(f"   • Visualization: {os.path.basename(viz_file)}")
    
    print(f"   • Statistical results: {os.path.basename(statistical_results['results_file'])}")
    
    print("\n✅ All analysis files generated successfully!")
    print(f"📁 Check output directory: {base_dir}")
    
    # Generate research text suggestions
    print("\n" + "=" * 60)
    print("RESEARCH TEXT SUGGESTIONS")
    print("=" * 60)
    
    # Load final genome statistics
    import pandas as pd
    final_stats = pd.read_csv(genome_stats_file)
    
    # Find extreme examples
    max_mgc_genome = final_stats.loc[final_stats['mgc_count'].idxmax()]
    min_mgc_genome = final_stats.loc[final_stats['mgc_count'].idxmin()]
    
    print(f"""
📝 SUGGESTED TEXT FOR RESEARCH PAPER:

"The number of identified KEGG-MGCs varied dramatically across analyzed genomes 
(Figure [MGC_DISTRIBUTION]). For example, {max_mgc_genome['genome_name']} contained 
{int(max_mgc_genome['mgc_count'])} clusters while {min_mgc_genome['genome_name']} 
had only {int(min_mgc_genome['mgc_count'])}. This difference may be explained by 
the total number of annotated genes (r = [CORRELATION_GENES], p = [P_VALUE]), 
gene density (genome size / number of genes) (r = [CORRELATION_DENSITY], p = [P_VALUE]), 
and the number of species-specific KEGG pathways (r = [CORRELATION_PATHWAYS], p = [P_VALUE])."

📊 KEY STATISTICS:
   • Total genomes analyzed: {len(final_stats)}
   • MGC range: {int(final_stats['mgc_count'].min())} - {int(final_stats['mgc_count'].max())}
   • Mean MGCs per genome: {final_stats['mgc_count'].mean():.1f} ± {final_stats['mgc_count'].std():.1f}
   • Median MGCs per genome: {final_stats['mgc_count'].median():.1f}

🔍 FIGURES TO REFERENCE:
   • mgc_distribution_analysis.png - Shows dramatic variation across genomes
   • correlation_analysis.png - Shows relationships with genome characteristics
   • correlation_heatmap.png - Comprehensive correlation matrix
   • pathway_analysis.png - KEGG pathway distribution analysis

📊 STATISTICAL TESTS:
   • Normality tests performed for all variables
   • Correlation analysis with multiple testing correction
   • Group comparisons using non-parametric tests
   • Regression analysis for predictive modeling
   • Pathway enrichment analysis
""")
    
    print("\n🎯 ANALYSIS COMPLETE - READY FOR PUBLICATION!")


if __name__ == "__main__":
    main()