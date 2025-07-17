#!/usr/bin/env python3
"""
Direct MGC Analysis for Research Paper

Direct analysis of MGC results to generate figures and statistics for the research paper.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from collections import Counter
import warnings
import csv
warnings.filterwarnings('ignore')

# Set publication-ready style
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 14
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
    
    return df

def get_genome_statistics(df):
    """Get genome-level statistics."""
    # Create organism-specific pathway identifiers
    df['organism_pathway'] = df['source_file'].apply(
        lambda x: os.path.basename(x).split('_')[0].lower()
    ) + '_' + df['pathway']
    
    genome_stats = df.groupby('genome_name').agg({
        'pathway': 'count',
        'organism_pathway': 'nunique',  # Count unique organism-pathways
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
                
                # Count annotated genes (assuming 'annotation' column exists)
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

def check_missing_files_from_mapping(mgc_df, mapping_file_path):
    """Check for files with KEGG annotations that are missing from MGC results."""
    print("🔍 Checking for missing files from organism mapping...")
    
    # Get files referenced in MGC results
    mgc_files = set()
    for source_file in mgc_df['source_file'].unique():
        if 'annotated_genomes_metabolic' in source_file:
            filename = os.path.basename(source_file).replace('_annotated.csv', '.csv')
            mgc_files.add(filename)
    
    # Get files with KEGG annotations from mapping
    mapping_files_with_kegg = set()
    try:
        with open(mapping_file_path, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                if len(row) >= 3 and row[2].strip():  # non-empty kegg_fasta
                    filename = os.path.basename(row[1])
                    mapping_files_with_kegg.add(filename)
    except Exception as e:
        print(f"Warning: Could not read mapping file {mapping_file_path}: {e}")
        return []
    
    # Find missing files
    missing_files = mapping_files_with_kegg - mgc_files
    
    print(f"📊 Files with KEGG annotations: {len(mapping_files_with_kegg)}")
    print(f"📊 Files found in MGC results: {len(mgc_files)}")
    print(f"📊 Missing files: {len(missing_files)}")
    
    if missing_files:
        print("\n⚠️  Files with KEGG annotations missing from MGC results:")
        for filename in sorted(missing_files):
            organism_name = clean_genome_name(filename)
            print(f"   - {organism_name} ({filename})")
    
    return sorted(missing_files)

def create_mgc_distribution_figure(genome_stats, output_dir, missing_files=None):
    """Create MGC distribution figure."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # Create combined dataset including missing organisms
    all_organisms = []
    mgc_counts = []
    colors = []
    
    # Add existing organisms
    for _, row in genome_stats.iterrows():
        all_organisms.append(row['genome_name'])
        mgc_counts.append(row['mgc_count'])
        colors.append('steelblue')
    
    # Add missing organisms if they exist
    if missing_files:
        for filename in missing_files:
            organism_name = clean_genome_name(filename)
            all_organisms.append(organism_name)
            mgc_counts.append(0)
            colors.append('red')
    
    # Sort by MGC count
    combined_data = list(zip(all_organisms, mgc_counts, colors))
    combined_data.sort(key=lambda x: x[1])
    
    organisms_sorted, counts_sorted, colors_sorted = zip(*combined_data)
    
    # 1. Horizontal bar chart - show all organisms including missing ones
    bars = ax1.barh(range(len(organisms_sorted)), counts_sorted, 
                   color=colors_sorted, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(organisms_sorted)))
    ax1.set_yticklabels([name[:20] + '...' if len(name) > 20 else name 
                        for name in organisms_sorted], fontsize=8)
    ax1.set_xlabel('Number of KEGG-MGCs')
    ax1.set_title('KEGG-MGC Distribution Across Genomes (Red = Missing)')
    ax1.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for i, (value, color) in enumerate(zip(counts_sorted, colors_sorted)):
        if value > 0:
            ax1.text(value + 0.5, i, str(int(value)), va='center', fontsize=8)
        elif color == 'red':
            ax1.text(0.5, i, '0', va='center', fontsize=8, 
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.7))
    
    # 2. Histogram - include missing organisms as 0 count
    all_counts = list(genome_stats['mgc_count'])
    if missing_files:
        all_counts.extend([0] * len(missing_files))
    
    ax2.hist(all_counts, bins=min(15, len(all_counts)//2), 
            color='lightcoral', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Number of KEGG-MGCs')
    ax2.set_ylabel('Number of Genomes')
    ax2.set_title('Distribution of KEGG-MGC Counts (Including Missing)')
    ax2.grid(alpha=0.3)
    
    # Add statistics
    mean_mgc = genome_stats['mgc_count'].mean()
    median_mgc = genome_stats['mgc_count'].median()
    ax2.axvline(mean_mgc, color='red', linestyle='--', alpha=0.7, label=f'Mean: {mean_mgc:.1f}')
    ax2.axvline(median_mgc, color='blue', linestyle='--', alpha=0.7, label=f'Median: {median_mgc:.1f}')
    if missing_files:
        ax2.axvline(0, color='red', linestyle='-', alpha=0.7, linewidth=2, label=f'Missing: {len(missing_files)}')
    ax2.legend()
    
    # 3. MGC length distribution
    ax3.hist(genome_stats['avg_length'], bins=min(15, len(genome_stats)//2), 
            color='lightgreen', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax3.set_xlabel('Average MGC Length (bp)')
    ax3.set_ylabel('Number of Genomes')
    ax3.set_title('Distribution of Average MGC Lengths')
    ax3.grid(alpha=0.3)
    
    # 4. Gene count distribution
    ax4.hist(genome_stats['avg_genes'], bins=min(15, len(genome_stats)//2), 
            color='orange', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax4.set_xlabel('Average Genes per MGC')
    ax4.set_ylabel('Number of Genomes')
    ax4.set_title('Distribution of Average Genes per MGC')
    ax4.grid(alpha=0.3)
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'combined_mgc_analysis_thesis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return output_file

def create_comprehensive_thesis_figure(genome_stats, source_stats, df, output_dir, missing_files=None):
    """Create comprehensive 6-panel thesis figure with missing organisms."""
    fig = plt.figure(figsize=(20, 12))
    
    # Create subplot layout: 2 rows, 3 columns
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # Top row: Distribution analysis
    ax1 = fig.add_subplot(gs[0, 0])  # MGC distribution bar chart
    ax2 = fig.add_subplot(gs[0, 1])  # MGC count histogram
    ax3 = fig.add_subplot(gs[0, 2])  # Correlation scatter
    
    # Bottom row: Pathway analysis
    ax4 = fig.add_subplot(gs[1, 0])  # Pathway bar chart
    ax5 = fig.add_subplot(gs[1, 1])  # Pathway histogram
    ax6 = fig.add_subplot(gs[1, 2])  # Additional analysis
    
    # Panel A: MGC Distribution with missing organisms
    all_organisms = []
    mgc_counts = []
    colors = []
    
    # Add existing organisms
    for _, row in genome_stats.iterrows():
        all_organisms.append(row['genome_name'])
        mgc_counts.append(row['mgc_count'])
        colors.append('steelblue')
    
    # Add missing organisms
    if missing_files:
        for filename in missing_files:
            organism_name = clean_genome_name(filename)
            all_organisms.append(organism_name)
            mgc_counts.append(0)
            colors.append('red')
    
    # Sort by MGC count and show top/bottom 10
    combined_data = list(zip(all_organisms, mgc_counts, colors))
    combined_data.sort(key=lambda x: x[1])
    
    # Show bottom 10 and top 10
    bottom_10 = combined_data[:10]
    top_10 = combined_data[-10:]
    display_data = bottom_10 + top_10
    
    orgs_display, counts_display, colors_display = zip(*display_data)
    
    bars = ax1.barh(range(len(orgs_display)), counts_display, 
                   color=colors_display, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(orgs_display)))
    ax1.set_yticklabels([name[:20] + '...' if len(name) > 20 else name 
                        for name in orgs_display], fontsize=8)
    ax1.set_xlabel('Number of KEGG-MGCs')
    ax1.set_title('A. Top 10 and Bottom 10 Genomes by MGC Count')
    ax1.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for i, (value, color) in enumerate(zip(counts_display, colors_display)):
        if value > 0:
            ax1.text(value + 0.5, i, str(int(value)), va='center', fontsize=8)
        elif color == 'red':
            ax1.text(0.5, i, '0', va='center', fontsize=8, 
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.7))
    
    # Panel B: Distribution histogram with missing organisms
    all_counts = list(genome_stats['mgc_count'])
    if missing_files:
        all_counts.extend([0] * len(missing_files))
    
    ax2.hist(all_counts, bins=min(15, len(all_counts)//2), 
            color='lightcoral', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Number of KEGG-MGCs')
    ax2.set_ylabel('Number of Genomes')
    ax2.set_title('B. Distribution of MGC Counts')
    ax2.grid(alpha=0.3)
    
    # Add statistics
    mean_mgc = genome_stats['mgc_count'].mean()
    median_mgc = genome_stats['mgc_count'].median()
    ax2.axvline(mean_mgc, color='red', linestyle='--', alpha=0.7, label=f'Mean: {mean_mgc:.1f}')
    ax2.axvline(median_mgc, color='blue', linestyle='--', alpha=0.7, label=f'Median: {median_mgc:.1f}')
    if missing_files:
        ax2.axvline(0, color='red', linestyle='-', alpha=0.7, linewidth=2, label=f'Missing: {len(missing_files)}')
    ax2.legend()
    
    # Panel C: Correlation analysis (using existing correlation logic)
    merged_df = pd.merge(genome_stats, source_stats, on='genome_name', how='inner')
    if len(merged_df) > 2:
        ax3.scatter(merged_df['unique_pathways'], merged_df['mgc_count'], 
                   alpha=0.6, s=50, color='darkblue', edgecolors='black', linewidth=0.5)
        
        # Calculate correlation
        from scipy.stats import pearsonr
        pearson_r, pearson_p = pearsonr(merged_df['unique_pathways'], merged_df['mgc_count'])
        
        # Add trend line
        if abs(pearson_r) > 0.1:
            z = np.polyfit(merged_df['unique_pathways'], merged_df['mgc_count'], 1)
            p = np.poly1d(z)
            ax3.plot(merged_df['unique_pathways'], p(merged_df['unique_pathways']), 
                    "r--", alpha=0.8, linewidth=2)
        
        ax3.set_xlabel('Unique Organism-Pathways')
        ax3.set_ylabel('MGC Count')
        ax3.set_title(f'C2. Unique Pathways vs MGC Count\nr = {pearson_r:.3f}, p = {pearson_p:.2e}')
        ax3.grid(alpha=0.3)
    
    # Panel D: Pathway analysis
    df['organism_pathway'] = df['source_file'].apply(
        lambda x: os.path.basename(x).split('_')[0].lower()
    ) + '_' + df['pathway']
    
    pathway_counts = df['organism_pathway'].value_counts()
    top_pathways = pathway_counts.head(15)
    
    ax4.barh(range(len(top_pathways)), top_pathways.values, 
            color='lightgreen', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax4.set_yticks(range(len(top_pathways)))
    ax4.set_yticklabels(top_pathways.index, fontsize=8)
    ax4.set_xlabel('Number of MGCs')
    ax4.set_title('D1. Top 15 Organism-Specific KEGG Pathways')
    ax4.grid(axis='x', alpha=0.3)
    
    # Panel E: Pathway distribution
    ax5.hist(pathway_counts.values, bins=20, color='orange', 
            alpha=0.7, edgecolor='black', linewidth=0.5)
    ax5.set_xlabel('Number of MGCs per Organism-Pathway')
    ax5.set_ylabel('Number of Organism-Pathways')
    ax5.set_title('D2. Distribution of MGCs per Organism-Pathway')
    ax5.grid(alpha=0.3)
    
    # Panel F: Key statistics summary
    ax6.axis('off')
    stats_text = f"""Key Statistics:
• Total genomes: {len(genome_stats) + len(missing_files if missing_files else [])}
• Total MGCs: {df.shape[0]}
• MGCs per genome: {mean_mgc:.1f} ± {genome_stats['mgc_count'].std():.1f}
• Range: {genome_stats['mgc_count'].min()} - {genome_stats['mgc_count'].max()}
• Unique organism-pathways: {df['organism_pathway'].nunique()}"""
    
    if missing_files:
        stats_text += f"\n• Missing organisms: {len(missing_files)}"
    
    ax6.text(0.1, 0.9, stats_text, transform=ax6.transAxes, fontsize=12, 
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.5", facecolor='wheat', alpha=0.8))
    
    plt.suptitle('KEGG-MGC Analysis Across Plant Genomes', fontsize=16, fontweight='bold')
    
    output_file = os.path.join(output_dir, 'combined_mgc_analysis_thesis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return output_file

def create_correlation_analysis(genome_stats, source_stats, output_dir):
    """Create correlation analysis between genome characteristics and MGC count."""
    # Merge genome stats with source stats
    merged_df = pd.merge(genome_stats, source_stats, on='genome_name', how='inner')
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    # Define relationships to analyze
    relationships = [
        ('total_genes', 'mgc_count', 'Total Genes vs MGC Count'),
        ('kegg_annotation_rate', 'mgc_count', 'KEGG Annotation Rate vs MGC Count'),
        ('unique_pathways', 'mgc_count', 'Unique Organism-Pathways vs MGC Count'),
        ('avg_length', 'mgc_count', 'Average MGC Length vs MGC Count'),
        ('avg_genes', 'mgc_count', 'Average Genes per MGC vs MGC Count'),
        ('annotation_rate', 'mgc_count', 'Annotation Rate vs MGC Count')
    ]
    
    correlation_results = {}
    
    for i, (x_var, y_var, title) in enumerate(relationships):
        if x_var in merged_df.columns and y_var in merged_df.columns:
            plot_data = merged_df[[x_var, y_var]].dropna()
            
            if len(plot_data) > 2:
                # Scatter plot
                axes[i].scatter(plot_data[x_var], plot_data[y_var], 
                              alpha=0.6, s=50, color='darkblue', edgecolors='black', linewidth=0.5)
                
                # Calculate correlation
                pearson_r, pearson_p = pearsonr(plot_data[x_var], plot_data[y_var])
                spearman_r, spearman_p = spearmanr(plot_data[x_var], plot_data[y_var])
                
                # Add trend line if correlation is meaningful
                if abs(pearson_r) > 0.1:
                    z = np.polyfit(plot_data[x_var], plot_data[y_var], 1)
                    p = np.poly1d(z)
                    axes[i].plot(plot_data[x_var], p(plot_data[x_var]), 
                               "r--", alpha=0.8, linewidth=2)
                
                axes[i].set_xlabel(x_var.replace('_', ' ').title())
                axes[i].set_ylabel(y_var.replace('_', ' ').title())
                axes[i].set_title(f'{title}\\nPearson r={pearson_r:.3f}, p={pearson_p:.2e}')
                axes[i].grid(alpha=0.3)
                
                correlation_results[f'{x_var}_vs_{y_var}'] = {
                    'pearson_r': pearson_r,
                    'pearson_p': pearson_p,
                    'spearman_r': spearman_r,
                    'spearman_p': spearman_p
                }
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'correlation_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return output_file, correlation_results

def create_pathway_analysis(df, output_dir):
    """Create pathway analysis figure with organism-specific pathway identifiers."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Create organism-specific pathway identifiers
    df['organism_pathway'] = df['source_file'].apply(
        lambda x: os.path.basename(x).split('_')[0].lower()
    ) + '_' + df['pathway']
    
    # Pathway frequency using organism-specific identifiers
    pathway_counts = df['organism_pathway'].value_counts()
    
    # Top pathways
    top_pathways = pathway_counts.head(15)
    ax1.barh(range(len(top_pathways)), top_pathways.values, 
            color='lightgreen', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(top_pathways)))
    ax1.set_yticklabels(top_pathways.index, fontsize=8)
    ax1.set_xlabel('Number of MGCs')
    ax1.set_title('Top 15 Organism-Specific KEGG Pathways in MGCs')
    ax1.grid(axis='x', alpha=0.3)
    
    # Pathway distribution
    ax2.hist(pathway_counts.values, bins=20, color='orange', 
            alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Number of MGCs per Organism-Pathway')
    ax2.set_ylabel('Number of Organism-Pathways')
    ax2.set_title('Distribution of MGCs per Organism-Pathway')
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'pathway_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return output_file


def generate_research_summary(genome_stats, correlation_results, df, output_dir, missing_files=None):
    """Generate research summary and text suggestions."""
    # Find extreme examples
    max_mgc_genome = genome_stats.loc[genome_stats['mgc_count'].idxmax()]
    min_mgc_genome = genome_stats.loc[genome_stats['mgc_count'].idxmin()]
    
    # Get correlation values
    corr_genes = correlation_results.get('total_genes_vs_mgc_count', {})
    corr_pathways = correlation_results.get('unique_pathways_vs_mgc_count', {})
    
    missing_section = ""
    if missing_files:
        missing_section = f"""
MISSING FILES ANALYSIS:
- {len(missing_files)} organisms with KEGG annotations are missing from MGC results
- This may indicate these genomes had no MGCs that met the filtering criteria
- Missing organisms: {', '.join([clean_genome_name(f) for f in missing_files[:5]])}{'...' if len(missing_files) > 5 else ''}
"""
    
    summary = f"""
KEGG-MGC ANALYSIS SUMMARY
========================

BASIC STATISTICS:
- Total genomes analyzed: {len(genome_stats)}
- Total MGCs identified: {df.shape[0]}
- MGCs per genome: {genome_stats['mgc_count'].mean():.1f} ± {genome_stats['mgc_count'].std():.1f}
- Range: {genome_stats['mgc_count'].min()} - {genome_stats['mgc_count'].max()}

EXTREME EXAMPLES:
- Highest MGC count: {max_mgc_genome['genome_name']} ({int(max_mgc_genome['mgc_count'])} MGCs)
- Lowest MGC count: {min_mgc_genome['genome_name']} ({int(min_mgc_genome['mgc_count'])} MGCs)

PATHWAY ANALYSIS:
- Total unique organism-pathways: {df['organism_pathway'].nunique()}
- Most common organism-pathway: {df['organism_pathway'].value_counts().index[0]} ({df['organism_pathway'].value_counts().iloc[0]} MGCs)
{missing_section}
RESEARCH TEXT SUGGESTION:
"The number of identified KEGG-MGCs varied dramatically across analyzed genomes 
(Figure XXX). For example, {max_mgc_genome['genome_name']} contained 
{int(max_mgc_genome['mgc_count'])} clusters while {min_mgc_genome['genome_name']} 
had only {int(min_mgc_genome['mgc_count'])}. This difference may be explained by 
the total number of annotated genes (r = {corr_genes.get('pearson_r', 0):.3f}, 
p = {corr_genes.get('pearson_p', 1):.2e}) and the number of unique organism-specific 
KEGG pathways (r = {corr_pathways.get('pearson_r', 0):.3f}, 
p = {corr_pathways.get('pearson_p', 1):.2e})."

FIGURES GENERATED:
- mgc_distribution_analysis.png: Shows dramatic variation across genomes
- correlation_analysis.png: Shows relationships with genome characteristics  
- pathway_analysis.png: KEGG pathway distribution analysis
"""
    
    # Save summary
    summary_file = os.path.join(output_dir, 'research_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(summary)
    
    print(summary)
    return summary_file

def main():
    """Main analysis function."""
    # File paths
    mgc_file = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/potential_groups_w10_filtered_with_length.csv'
    mapping_file = '/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping_with_fasta.csv'
    output_dir = os.path.dirname(mgc_file)
    
    print("🧬 DIRECT MGC ANALYSIS FOR RESEARCH PAPER")
    print("=" * 50)
    
    # Load and process data
    print("📊 Loading and processing data...")
    df = load_and_process_data(mgc_file)
    
    # Check for missing files from organism mapping
    missing_files = check_missing_files_from_mapping(df, mapping_file)
    
    # Get genome statistics
    print("📈 Calculating genome statistics...")
    genome_stats = get_genome_statistics(df)
    
    # Get source file statistics
    print("🔍 Analyzing source files...")
    source_stats = get_source_file_stats(df)
    
    # Create comprehensive thesis figure
    print("🎨 Creating comprehensive thesis figure...")
    thesis_file = create_comprehensive_thesis_figure(genome_stats, source_stats, df, output_dir, missing_files)
    
    print("📊 Creating correlation analysis...")
    corr_file, corr_results = create_correlation_analysis(genome_stats, source_stats, output_dir)
    
    print("🛤️  Creating pathway analysis...")
    path_file = create_pathway_analysis(df, output_dir)
    
    # Generate research summary
    print("📝 Generating research summary...")
    summary_file = generate_research_summary(genome_stats, corr_results, df, output_dir, missing_files)
    
    print("\n✅ Analysis complete!")
    print(f"📁 Output files saved to: {output_dir}")
    print(f"   - {os.path.basename(thesis_file)}")
    print(f"   - {os.path.basename(corr_file)}")
    print(f"   - {os.path.basename(path_file)}")
    print(f"   - {os.path.basename(summary_file)}")
    
    if missing_files:
        print(f"\n⚠️  Note: {len(missing_files)} organisms with KEGG annotations had no MGCs in results")
        print(f"   Missing organisms are shown in red in the comprehensive thesis figure")

if __name__ == "__main__":
    main()