#!/usr/bin/env python3
"""
MGC Length Comparison Analysis
Compares length statistics between MiBIG clusters and KEGG verified clusters
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_data():
    """Load both datasets"""
    mibig_df = pd.read_csv('/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/mgc_stats.csv')
    kegg_df = pd.read_csv('/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/potential_groups_w10_filtered_with_length.csv')
    
    # Calculate gene counts for KEGG data using metabolic_genes column
    kegg_df['gene_count'] = kegg_df['metabolic_genes'].apply(lambda x: len(x.split(',')))
    kegg_df['length_per_gene'] = kegg_df['length'] / kegg_df['gene_count']
    
    # Calculate gene counts for MiBIG data by reading CSV files
    csv_files_dir = '/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/csv_files'
    gene_counts = []
    
    for _, row in mibig_df.iterrows():
        mgc_name = row['mgc_name']
        csv_file_path = f"{csv_files_dir}/{mgc_name}.csv"
        try:
            # Read the CSV file and count data rows (excluding header)
            gene_df = pd.read_csv(csv_file_path)
            gene_count = len(gene_df)
            gene_counts.append(gene_count)
        except FileNotFoundError:
            print(f"Warning: CSV file not found for {mgc_name}")
            gene_counts.append(None)
    
    mibig_df['gene_count'] = gene_counts
    mibig_df['length_per_gene'] = mibig_df['length'] / mibig_df['gene_count']
    
    return mibig_df, kegg_df

def calculate_statistics(mibig_df, kegg_df):
    """Calculate descriptive statistics for both datasets"""
    
    mibig_lengths = mibig_df['length']
    kegg_lengths = kegg_df['length']
    mibig_length_per_gene = mibig_df['length_per_gene']
    kegg_length_per_gene = kegg_df['length_per_gene']
    
    print('=== DATASET COMPARISON ===')
    print(f'MiBIG clusters: {len(mibig_lengths)} samples')
    print(f'KEGG verified clusters: {len(kegg_lengths)} samples')
    print()

    print('=== DESCRIPTIVE STATISTICS ===')
    print('MiBIG Clusters (absolute length):')
    print(f'  Mean: {mibig_lengths.mean():.2f}')
    print(f'  Median: {mibig_lengths.median():.2f}')
    print(f'  Std Dev: {mibig_lengths.std():.2f}')
    print(f'  Min: {mibig_lengths.min()}')
    print(f'  Max: {mibig_lengths.max()}')
    print(f'  Q1: {mibig_lengths.quantile(0.25):.2f}')
    print(f'  Q3: {mibig_lengths.quantile(0.75):.2f}')
    print()

    print('KEGG Verified Clusters (absolute length):')
    print(f'  Mean: {kegg_lengths.mean():.2f}')
    print(f'  Median: {kegg_lengths.median():.2f}')
    print(f'  Std Dev: {kegg_lengths.std():.2f}')
    print(f'  Min: {kegg_lengths.min()}')
    print(f'  Max: {kegg_lengths.max()}')
    print(f'  Q1: {kegg_lengths.quantile(0.25):.2f}')
    print(f'  Q3: {kegg_lengths.quantile(0.75):.2f}')
    print()

    print('KEGG Verified Clusters (length per gene):')
    print(f'  Mean: {kegg_length_per_gene.mean():.2f}')
    print(f'  Median: {kegg_length_per_gene.median():.2f}')
    print(f'  Std Dev: {kegg_length_per_gene.std():.2f}')
    print(f'  Min: {kegg_length_per_gene.min():.2f}')
    print(f'  Max: {kegg_length_per_gene.max():.2f}')
    print(f'  Q1: {kegg_length_per_gene.quantile(0.25):.2f}')
    print(f'  Q3: {kegg_length_per_gene.quantile(0.75):.2f}')
    print()

    print('MiBIG Clusters (length per gene):')
    print(f'  Mean: {mibig_length_per_gene.mean():.2f}')
    print(f'  Median: {mibig_length_per_gene.median():.2f}')
    print(f'  Std Dev: {mibig_length_per_gene.std():.2f}')
    print(f'  Min: {mibig_length_per_gene.min():.2f}')
    print(f'  Max: {mibig_length_per_gene.max():.2f}')
    print(f'  Q1: {mibig_length_per_gene.quantile(0.25):.2f}')
    print(f'  Q3: {mibig_length_per_gene.quantile(0.75):.2f}')
    print()

    print('=== COMPARISON ===')
    print('Absolute Length Differences:')
    print(f'  Mean difference: {mibig_lengths.mean() - kegg_lengths.mean():.2f} bp')
    print(f'  Median difference: {mibig_lengths.median() - kegg_lengths.median():.2f} bp')
    print(f'  MiBIG/KEGG mean ratio: {mibig_lengths.mean() / kegg_lengths.mean():.2f}x')
    print(f'  MiBIG/KEGG median ratio: {mibig_lengths.median() / kegg_lengths.median():.2f}x')
    print()

    print('Length Per Gene Differences:')
    print(f'  Mean difference: {mibig_length_per_gene.mean() - kegg_length_per_gene.mean():.2f} bp/gene')
    print(f'  Median difference: {mibig_length_per_gene.median() - kegg_length_per_gene.median():.2f} bp/gene')
    print(f'  MiBIG/KEGG mean ratio: {mibig_length_per_gene.mean() / kegg_length_per_gene.mean():.2f}x')
    print(f'  MiBIG/KEGG median ratio: {mibig_length_per_gene.median() / kegg_length_per_gene.median():.2f}x')
    print()

    print('=== DISTRIBUTION ANALYSIS ===')
    # Percentile analysis
    print('Length percentiles:')
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    for p in percentiles:
        mibig_p = np.percentile(mibig_lengths, p)
        kegg_p = np.percentile(kegg_lengths, p)
        print(f'  {p}th percentile - MiBIG: {mibig_p:.0f}, KEGG: {kegg_p:.0f}, Ratio: {mibig_p/kegg_p:.2f}x')
    print()

    # Range analysis
    print('Range analysis:')
    print(f'  MiBIG range: {mibig_lengths.max() - mibig_lengths.min():.0f} bp')
    print(f'  KEGG range: {kegg_lengths.max() - kegg_lengths.min():.0f} bp')
    print(f'  MiBIG IQR: {mibig_lengths.quantile(0.75) - mibig_lengths.quantile(0.25):.0f} bp')
    print(f'  KEGG IQR: {kegg_lengths.quantile(0.75) - kegg_lengths.quantile(0.25):.0f} bp')
    print()

    print('=== SUMMARY ===')
    print('Key findings:')
    if mibig_lengths.mean() > kegg_lengths.mean():
        print(f'• MiBIG clusters are significantly longer on average ({mibig_lengths.mean():.0f} vs {kegg_lengths.mean():.0f} bp)')
    else:
        print(f'• KEGG clusters are significantly longer on average ({kegg_lengths.mean():.0f} vs {mibig_lengths.mean():.0f} bp)')
        
    print(f'• MiBIG has higher variability (std: {mibig_lengths.std():.0f} vs {kegg_lengths.std():.0f})')
    print(f'• Largest MiBIG: {mibig_lengths.max():.0f} bp, Largest KEGG: {kegg_lengths.max():.0f} bp')
    print(f'• Smallest MiBIG: {mibig_lengths.min():.0f} bp, Smallest KEGG: {kegg_lengths.min():.0f} bp')
    
def create_figures_length(mibig_df, kegg_df, output_dir='./'):
    """Create comparison figures"""
    
    mibig_length = mibig_df['length']
    kegg_length = kegg_df['length']
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Histogram comparison
    ax1 = plt.subplot(2, 3, 1)
    plt.hist(mibig_length, bins=20, alpha=0.7, label='MiBIG', color='skyblue', edgecolor='black')
    plt.hist(kegg_length, bins=50, alpha=0.7, label='KEGG', color='lightcoral', edgecolor='black')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Length Distribution Comparison')
    plt.legend()
    plt.yscale('log')
    
    # 2. Box plot comparison
    ax2 = plt.subplot(2, 3, 2)
    data_to_plot = [mibig_length, kegg_length]
    labels = ['MiBIG', 'KEGG']
    plt.boxplot(data_to_plot, labels=labels, patch_artist=True,
                boxprops=dict(facecolor='lightblue', alpha=0.7),
                medianprops=dict(color='red', linewidth=2))
    plt.ylabel('Length (bp)')
    plt.title('Length Distribution Box Plot')
    plt.yscale('log')
    
    # 3. Violin plot
    ax3 = plt.subplot(2, 3, 3)
    combined_data = pd.DataFrame({
        'Length': pd.concat([mibig_length, kegg_length]),
        'Dataset': ['MiBIG']*len(mibig_length) + ['KEGG']*len(kegg_length)
    })
    sns.violinplot(data=combined_data, x='Dataset', y='Length', ax=ax3)
    plt.title('Length Distribution Violin Plot')
    plt.yscale('log')
    
    # 4. Cumulative distribution
    ax4 = plt.subplot(2, 3, 4)
    mibig_sorted = np.sort(mibig_length)
    kegg_sorted = np.sort(kegg_length)
    mibig_cumulative = np.arange(1, len(mibig_sorted) + 1) / len(mibig_sorted)
    kegg_cumulative = np.arange(1, len(kegg_sorted) + 1) / len(kegg_sorted)
    
    plt.plot(mibig_sorted, mibig_cumulative, label='MiBIG', linewidth=2)
    plt.plot(kegg_sorted, kegg_cumulative, label='KEGG', linewidth=2)
    plt.xlabel('Length (bp)')
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Distribution Function')
    plt.legend()
    plt.xscale('log')
    
    # 5. Percentile comparison
    ax5 = plt.subplot(2, 3, 5)
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    mibig_percentiles = [np.percentile(mibig_length, p) for p in percentiles]
    kegg_percentiles = [np.percentile(kegg_length, p) for p in percentiles]
    
    x = np.arange(len(percentiles))
    width = 0.35
    
    plt.bar(x - width/2, mibig_percentiles, width, label='MiBIG', alpha=0.7)
    plt.bar(x + width/2, kegg_percentiles, width, label='KEGG', alpha=0.7)
    plt.xlabel('Percentile')
    plt.ylabel('Length (bp)')
    plt.title('Percentile Comparison')
    plt.xticks(x, [f'{p}th' for p in percentiles])
    plt.legend()
    plt.yscale('log')
    
    # 6. Statistics summary table
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('tight')
    ax6.axis('off')
    
    stats_data = [
        ['Statistic', 'MiBIG', 'KEGG', 'Ratio'],
        ['Count', f'{len(mibig_length)}', f'{len(kegg_length)}', f'{len(mibig_length)/len(kegg_length):.2f}'],
        ['Mean', f'{mibig_length.mean():.0f}', f'{kegg_length.mean():.0f}', f'{mibig_length.mean()/kegg_length.mean():.2f}'],
        ['Median', f'{mibig_length.median():.0f}', f'{kegg_length.median():.0f}', f'{mibig_length.median()/kegg_length.median():.2f}'],
        ['Std Dev', f'{mibig_length.std():.0f}', f'{kegg_length.std():.0f}', f'{mibig_length.std()/kegg_length.std():.2f}'],
        ['Min', f'{mibig_length.min():.0f}', f'{kegg_length.min():.0f}', f'{mibig_length.min()/kegg_length.min():.2f}'],
        ['Max', f'{mibig_length.max():.0f}', f'{kegg_length.max():.0f}', f'{mibig_length.max()/kegg_length.max():.2f}']
    ]
    
    table = ax6.table(cellText=stats_data, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    
    # Style the header row
    for i in range(4):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    plt.title('Statistical Summary (Length)', pad=20)
    
    plt.tight_layout()
    
    # Save the figure
    output_path = Path(output_dir) / 'mgc_length_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f'Figure saved to: {output_path}')
    
    # Show the plot
    plt.show()

def create_figures_length_per_gene(mibig_df, kegg_df, output_dir='./'):
    """Create comparison figures"""
    
    mibig_length_per_gene = mibig_df['length_per_gene']
    kegg_length_per_gene = kegg_df['length_per_gene']
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Histogram comparison
    ax1 = plt.subplot(2, 3, 1)
    plt.hist(mibig_length_per_gene, bins=20, alpha=0.7, label='MiBIG', color='skyblue', edgecolor='black')
    plt.hist(kegg_length_per_gene, bins=50, alpha=0.7, label='KEGG', color='lightcoral', edgecolor='black')
    plt.xlabel('Length per Gene (bp/gene)')
    plt.ylabel('Frequency')
    plt.title('Length per Gene Distribution Comparison')
    plt.legend()
    plt.yscale('log')
    
    # 2. Box plot comparison
    ax2 = plt.subplot(2, 3, 2)
    data_to_plot = [mibig_length_per_gene, kegg_length_per_gene]
    labels = ['MiBIG', 'KEGG']
    plt.boxplot(data_to_plot, labels=labels, patch_artist=True,
                boxprops=dict(facecolor='lightblue', alpha=0.7),
                medianprops=dict(color='red', linewidth=2))
    plt.ylabel('Length per Gene (bp/gene)')
    plt.title('Length per Gene Distribution Box Plot')
    plt.yscale('log')
    
    # 3. Violin plot
    ax3 = plt.subplot(2, 3, 3)
    combined_data = pd.DataFrame({
        'Length_per_Gene': pd.concat([mibig_length_per_gene, kegg_length_per_gene]),
        'Dataset': ['MiBIG']*len(mibig_length_per_gene) + ['KEGG']*len(kegg_length_per_gene)
    })
    sns.violinplot(data=combined_data, x='Dataset', y='Length_per_Gene', ax=ax3)
    plt.title('Length per Gene Distribution Violin Plot')
    plt.yscale('log')
    
    # 4. Cumulative distribution
    ax4 = plt.subplot(2, 3, 4)
    mibig_sorted = np.sort(mibig_length_per_gene)
    kegg_sorted = np.sort(kegg_length_per_gene)
    mibig_cumulative = np.arange(1, len(mibig_sorted) + 1) / len(mibig_sorted)
    kegg_cumulative = np.arange(1, len(kegg_sorted) + 1) / len(kegg_sorted)
    
    plt.plot(mibig_sorted, mibig_cumulative, label='MiBIG', linewidth=2)
    plt.plot(kegg_sorted, kegg_cumulative, label='KEGG', linewidth=2)
    plt.xlabel('Length per Gene (bp/gene)')
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Distribution Function')
    plt.legend()
    plt.xscale('log')
    
    # 5. Percentile comparison
    ax5 = plt.subplot(2, 3, 5)
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    mibig_percentiles = [np.percentile(mibig_length_per_gene, p) for p in percentiles]
    kegg_percentiles = [np.percentile(kegg_length_per_gene, p) for p in percentiles]
    
    x = np.arange(len(percentiles))
    width = 0.35
    
    plt.bar(x - width/2, mibig_percentiles, width, label='MiBIG', alpha=0.7)
    plt.bar(x + width/2, kegg_percentiles, width, label='KEGG', alpha=0.7)
    plt.xlabel('Percentile')
    plt.ylabel('Length per Gene (bp/gene)')
    plt.title('Percentile Comparison')
    plt.xticks(x, [f'{p}th' for p in percentiles])
    plt.legend()
    plt.yscale('log')
    
    # 6. Statistics summary table
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('tight')
    ax6.axis('off')
    
    stats_data = [
        ['Statistic', 'MiBIG', 'KEGG', 'Ratio'],
        ['Count', f'{len(mibig_length_per_gene)}', f'{len(kegg_length_per_gene)}', f'{len(mibig_length_per_gene)/len(kegg_length_per_gene):.2f}'],
        ['Mean', f'{mibig_length_per_gene.mean():.0f}', f'{kegg_length_per_gene.mean():.0f}', f'{mibig_length_per_gene.mean()/kegg_length_per_gene.mean():.2f}'],
        ['Median', f'{mibig_length_per_gene.median():.0f}', f'{kegg_length_per_gene.median():.0f}', f'{mibig_length_per_gene.median()/kegg_length_per_gene.median():.2f}'],
        ['Std Dev', f'{mibig_length_per_gene.std():.0f}', f'{kegg_length_per_gene.std():.0f}', f'{mibig_length_per_gene.std()/kegg_length_per_gene.std():.2f}'],
        ['Min', f'{mibig_length_per_gene.min():.0f}', f'{kegg_length_per_gene.min():.0f}', f'{mibig_length_per_gene.min()/kegg_length_per_gene.min():.2f}'],
        ['Max', f'{mibig_length_per_gene.max():.0f}', f'{kegg_length_per_gene.max():.0f}', f'{mibig_length_per_gene.max()/kegg_length_per_gene.max():.2f}']
    ]
    
    table = ax6.table(cellText=stats_data, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    
    # Style the header row
    for i in range(4):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    plt.title('Statistical Summary (Length per Gene)', pad=20)
    
    plt.tight_layout()
    
    # Save the figure
    output_path = Path(output_dir) / 'mgc_length_per_gene_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f'Figure saved to: {output_path}')
    
    # Show the plot
    plt.show()

def main():
    """Main function to run the analysis"""
    print("MGC Length Comparison Analysis")
    print("=" * 50)
    
    # Load data
    mibig_df, kegg_df = load_data()
    
    # Calculate and print statistics
    calculate_statistics(mibig_df, kegg_df)
    
    create_figures_length(mibig_df, kegg_df, output_dir='/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/')
    
    # Create figures
    create_figures_length_per_gene(mibig_df, kegg_df, output_dir='/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/')

if __name__ == "__main__":
    main()