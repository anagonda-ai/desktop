#!/usr/bin/env python3
"""
Combined MGC Analysis Figure for Thesis
Creates a comprehensive figure combining all MGC analysis results.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Set publication-ready style with larger, readable fonts
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 20
})

def clean_genome_name(source_file):
    """Clean genome name from source file path."""
    name = os.path.basename(source_file).lower()
    suffixes = [
        '_updated_annotated.csv', '_annotated.csv', '.filtered.csv',
        '.filtered', '.csv', '_metabolic.csv'
    ]
    
    for suffix in suffixes:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    
    name = name.strip('._')
    name = name.replace('..', '.')
    name = name.replace('_', ' ')
    name = name.title()
    
    return name

def load_and_process_data(file_path):
    """Load and process MGC data."""
    df = pd.read_csv(file_path)
    
    # Clean genome names
    df['genome_name'] = df['source_file'].apply(clean_genome_name)
    
    # Count genes in each MGC
    df['gene_count'] = df['genes'].apply(lambda x: len(x.split(',')) if pd.notna(x) else 0)
    
    # Create organism-specific pathway identifiers
    df['organism_pathway'] = df['source_file'].apply(
        lambda x: os.path.basename(x).split('_')[0].lower()
    ) + '_' + df['pathway']
    
    return df

def get_genome_statistics(df):
    """Get genome-level statistics."""
    genome_stats = df.groupby('genome_name').agg({
        'pathway': 'count',
        'organism_pathway': 'nunique',
        'length': ['mean', 'median', 'std', 'min', 'max'],
        'gene_count': ['mean', 'median', 'std', 'min', 'max']
    }).round(2)
    
    genome_stats.columns = ['mgc_count', 'unique_pathways', 'avg_length', 'median_length', 'std_length', 
                           'min_length', 'max_length', 'avg_genes', 'median_genes', 
                           'std_genes', 'min_genes', 'max_genes']
    
    return genome_stats.reset_index()

def get_source_file_stats(df):
    """Get source file statistics for annotation analysis."""
    source_stats = []
    
    for source_file in df['source_file'].unique():
        if os.path.exists(source_file):
            try:
                source_df = pd.read_csv(source_file)
                total_genes = len(source_df)
                
                # Count annotated genes
                if 'annotation' in source_df.columns:
                    annotated_genes = source_df['annotation'].notna().sum()
                    annotation_rate = annotated_genes / total_genes
                else:
                    annotated_genes = total_genes
                    annotation_rate = 1.0
                
                # Count KEGG annotations
                kegg_cols = [col for col in source_df.columns if 'kegg' in col.lower() or 'ko' in col.lower()]
                if kegg_cols:
                    kegg_annotated = source_df[kegg_cols[0]].notna().sum()
                    kegg_rate = kegg_annotated / total_genes
                else:
                    kegg_annotated = 0
                    kegg_rate = 0.0
                
                source_stats.append({
                    'source_file': source_file,
                    'genome_name': clean_genome_name(source_file),
                    'total_genes': total_genes,
                    'annotated_genes': annotated_genes,
                    'annotation_rate': annotation_rate,
                    'kegg_annotated': kegg_annotated,
                    'kegg_annotation_rate': kegg_rate
                })
            except Exception as e:
                print(f"Error processing {source_file}: {e}")
    
    return pd.DataFrame(source_stats)

def create_combined_figure(df, genome_stats, source_stats, output_dir):
    """Create combined figure for thesis."""
    # Create figure with subplots - larger for readability
    fig = plt.figure(figsize=(24, 16))
    
    # Create main title
    fig.suptitle('KEGG-MGC Analysis Across Plant Genomes', fontsize=24, fontweight='bold', y=0.97)
    
    # Define grid layout with more spacing
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1], width_ratios=[1, 1, 1, 1], 
                          hspace=0.4, wspace=0.35, top=0.92, bottom=0.08, left=0.06, right=0.98)
    
    # ===== PANEL A: MGC Distribution =====
    ax1 = fig.add_subplot(gs[0, :2])
    
    # Sort by MGC count
    sorted_stats = genome_stats.sort_values('mgc_count', ascending=True)
    
    # Show only top 10 and bottom 10 genomes for readability
    top_10 = sorted_stats.tail(10)
    bottom_10 = sorted_stats.head(10)
    
    # Combine with a gap indicator
    display_genomes = pd.concat([bottom_10, top_10])
    display_positions = list(range(10)) + list(range(12, 22))  # Leave gap at position 10-11
    
    # Create horizontal bar chart
    bars = ax1.barh(display_positions, display_genomes['mgc_count'], 
                   color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Set up y-axis
    ax1.set_yticks(display_positions)
    genome_labels = [name[:15] + '...' if len(name) > 15 else name 
                    for name in display_genomes['genome_name']]
    ax1.set_yticklabels(genome_labels, fontsize=11)
    
    # Add gap indicator
    ax1.text(0, 10.5, '... (middle genomes omitted) ...', ha='left', va='center', 
             fontsize=10, style='italic', color='gray')
    
    ax1.set_xlabel('Number of KEGG-MGCs', fontweight='bold', fontsize=15)
    ax1.set_title('A. Top 10 and Bottom 10 Genomes by MGC Count', fontweight='bold', pad=15, fontsize=18)
    ax1.grid(axis='x', alpha=0.3)
    
    # Add value labels for all values now that there are fewer
    for i, (pos, value) in enumerate(zip(display_positions, display_genomes['mgc_count'])):
        ax1.text(value + 0.5, pos, str(int(value)), va='center', fontsize=12, fontweight='bold')
    
    # ===== PANEL B: MGC Count Histogram =====
    ax2 = fig.add_subplot(gs[0, 2:])
    
    ax2.hist(genome_stats['mgc_count'], bins=15, color='lightcoral', alpha=0.7, 
             edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Number of KEGG-MGCs', fontweight='bold', fontsize=15)
    ax2.set_ylabel('Number of Genomes', fontweight='bold', fontsize=15)
    ax2.set_title('B. Distribution of MGC Counts', fontweight='bold', pad=15, fontsize=18)
    ax2.grid(alpha=0.3)
    
    # Add statistics
    mean_mgc = genome_stats['mgc_count'].mean()
    median_mgc = genome_stats['mgc_count'].median()
    ax2.axvline(mean_mgc, color='red', linestyle='--', alpha=0.8, linewidth=3, 
                label=f'Mean: {mean_mgc:.1f}')
    ax2.axvline(median_mgc, color='blue', linestyle='--', alpha=0.8, linewidth=3, 
                label=f'Median: {median_mgc:.1f}')
    ax2.legend(fontsize=13)
    
    # ===== PANEL C: Key Correlations =====
    # Merge genome stats with source stats
    merged_df = pd.merge(genome_stats, source_stats, on='genome_name', how='inner')
    
    # C1: Total genes vs MGC count
    ax3 = fig.add_subplot(gs[1, 0])
    plot_data = merged_df[['total_genes', 'mgc_count']].dropna()
    ax3.scatter(plot_data['total_genes'], plot_data['mgc_count'], 
               alpha=0.7, s=60, color='darkgreen', edgecolors='black', linewidth=0.8)
    
    pearson_r, pearson_p = pearsonr(plot_data['total_genes'], plot_data['mgc_count'])
    if abs(pearson_r) > 0.1:
        z = np.polyfit(plot_data['total_genes'], plot_data['mgc_count'], 1)
        p = np.poly1d(z)
        ax3.plot(plot_data['total_genes'], p(plot_data['total_genes']), 
                "r--", alpha=0.8, linewidth=3)
    
    ax3.set_xlabel('Total Genes', fontweight='bold', fontsize=14)
    ax3.set_ylabel('MGC Count', fontweight='bold', fontsize=14)
    ax3.set_title(f'C1. Total Genes vs MGC Count\nr = {pearson_r:.3f}, p = {pearson_p:.2e}', 
                  fontweight='bold', pad=12, fontsize=16)
    ax3.grid(alpha=0.3)
    
    # C2: Unique pathways vs MGC count
    ax4 = fig.add_subplot(gs[1, 1])
    plot_data = merged_df[['unique_pathways', 'mgc_count']].dropna()
    ax4.scatter(plot_data['unique_pathways'], plot_data['mgc_count'], 
               alpha=0.7, s=60, color='purple', edgecolors='black', linewidth=0.8)
    
    pearson_r, pearson_p = pearsonr(plot_data['unique_pathways'], plot_data['mgc_count'])
    if abs(pearson_r) > 0.1:
        z = np.polyfit(plot_data['unique_pathways'], plot_data['mgc_count'], 1)
        p = np.poly1d(z)
        ax4.plot(plot_data['unique_pathways'], p(plot_data['unique_pathways']), 
                "r--", alpha=0.8, linewidth=3)
    
    ax4.set_xlabel('Unique Organism-Pathways', fontweight='bold', fontsize=14)
    ax4.set_ylabel('MGC Count', fontweight='bold', fontsize=14)
    ax4.set_title(f'C2. Unique Pathways vs MGC Count\nr = {pearson_r:.3f}, p = {pearson_p:.2e}', 
                  fontweight='bold', pad=12, fontsize=16)
    ax4.grid(alpha=0.3)
    
    # C3: KEGG annotation rate vs MGC count
    ax5 = fig.add_subplot(gs[1, 2])
    plot_data = merged_df[['kegg_annotation_rate', 'mgc_count']].dropna()
    ax5.scatter(plot_data['kegg_annotation_rate'], plot_data['mgc_count'], 
               alpha=0.7, s=60, color='orange', edgecolors='black', linewidth=0.8)
    
    pearson_r, pearson_p = pearsonr(plot_data['kegg_annotation_rate'], plot_data['mgc_count'])
    if abs(pearson_r) > 0.1:
        z = np.polyfit(plot_data['kegg_annotation_rate'], plot_data['mgc_count'], 1)
        p = np.poly1d(z)
        ax5.plot(plot_data['kegg_annotation_rate'], p(plot_data['kegg_annotation_rate']), 
                "r--", alpha=0.8, linewidth=3)
    
    ax5.set_xlabel('KEGG Annotation Rate', fontweight='bold', fontsize=14)
    ax5.set_ylabel('MGC Count', fontweight='bold', fontsize=14)
    ax5.set_title(f'C3. KEGG Annotation vs MGC Count\nr = {pearson_r:.3f}, p = {pearson_p:.2e}', 
                  fontweight='bold', pad=12, fontsize=16)
    ax5.grid(alpha=0.3)
    
    # C4: Average MGC length vs count
    ax6 = fig.add_subplot(gs[1, 3])
    plot_data = merged_df[['avg_length', 'mgc_count']].dropna()
    ax6.scatter(plot_data['avg_length'], plot_data['mgc_count'], 
               alpha=0.7, s=60, color='brown', edgecolors='black', linewidth=0.8)
    
    pearson_r, pearson_p = pearsonr(plot_data['avg_length'], plot_data['mgc_count'])
    if abs(pearson_r) > 0.1:
        z = np.polyfit(plot_data['avg_length'], plot_data['mgc_count'], 1)
        p = np.poly1d(z)
        ax6.plot(plot_data['avg_length'], p(plot_data['avg_length']), 
                "r--", alpha=0.8, linewidth=3)
    
    ax6.set_xlabel('Average MGC Length (bp)', fontweight='bold', fontsize=14)
    ax6.set_ylabel('MGC Count', fontweight='bold', fontsize=14)
    ax6.set_title(f'C4. MGC Length vs Count\nr = {pearson_r:.3f}, p = {pearson_p:.2e}', 
                  fontweight='bold', pad=12, fontsize=16)
    ax6.grid(alpha=0.3)
    
    # ===== PANEL D: Pathway Analysis =====
    # D1: Top organism-pathways
    ax7 = fig.add_subplot(gs[2, :2])
    
    pathway_counts = df['organism_pathway'].value_counts()
    top_pathways = pathway_counts.head(15)
    
    bars = ax7.barh(range(len(top_pathways)), top_pathways.values, 
                   color='lightgreen', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax7.set_yticks(range(len(top_pathways)))
    ax7.set_yticklabels([path[:25] + '...' if len(path) > 25 else path 
                        for path in top_pathways.index], fontsize=12)
    ax7.set_xlabel('Number of MGCs', fontweight='bold', fontsize=15)
    ax7.set_title('D1. Top 15 Organism-Specific KEGG Pathways', fontweight='bold', pad=15, fontsize=18)
    ax7.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, top_pathways.values)):
        ax7.text(value + 0.3, i, str(int(value)), va='center', fontsize=12, fontweight='bold')
    
    # D2: Pathway distribution
    ax8 = fig.add_subplot(gs[2, 2:])
    
    ax8.hist(pathway_counts.values, bins=20, color='orange', 
             alpha=0.7, edgecolor='black', linewidth=0.5)
    ax8.set_xlabel('Number of MGCs per Organism-Pathway', fontweight='bold', fontsize=15)
    ax8.set_ylabel('Number of Organism-Pathways', fontweight='bold', fontsize=15)
    ax8.set_title('D2. Distribution of MGCs per Organism-Pathway', fontweight='bold', pad=15, fontsize=18)
    ax8.grid(alpha=0.3)
    
    # Add statistics text box
    stats_text = f"""Key Statistics:
• Total genomes: {len(genome_stats)}
• Total MGCs: {len(df)}
• MGCs per genome: {genome_stats['mgc_count'].mean():.1f} ± {genome_stats['mgc_count'].std():.1f}
• Range: {genome_stats['mgc_count'].min()}-{genome_stats['mgc_count'].max()}
• Unique organism-pathways: {df['organism_pathway'].nunique()}"""
    
    ax8.text(0.98, 0.98, stats_text, transform=ax8.transAxes, fontsize=12,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=0.8))
    
    # Save figure
    output_file = os.path.join(output_dir, 'combined_mgc_analysis_thesis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return output_file

def main():
    """Main function to create combined figure."""
    # File paths
    mgc_file = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/potential_groups_w10_filtered_with_length.csv'
    output_dir = os.path.dirname(mgc_file)
    
    print("🧬 CREATING COMBINED MGC ANALYSIS FIGURE FOR THESIS")
    print("=" * 55)
    
    # Load and process data
    print("📊 Loading and processing data...")
    df = load_and_process_data(mgc_file)
    
    # Get genome statistics
    print("📈 Calculating genome statistics...")
    genome_stats = get_genome_statistics(df)
    
    # Get source file statistics
    print("🔍 Analyzing source files...")
    source_stats = get_source_file_stats(df)
    
    # Create combined figure
    print("🎨 Creating combined figure...")
    output_file = create_combined_figure(df, genome_stats, source_stats, output_dir)
    
    print("\n✅ Combined figure created successfully!")
    print(f"📁 Output file: {output_file}")
    print("\n📋 Figure panels:")
    print("   A. Distribution of KEGG-MGCs across plant genomes")
    print("   B. Distribution of MGC counts")
    print("   C1-C4. Correlation analyses")
    print("   D1. Top organism-specific KEGG pathways")
    print("   D2. Distribution of MGCs per organism-pathway")

if __name__ == "__main__":
    main()