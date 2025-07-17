"""
Advanced MGC Visualization Suite

This module provides comprehensive visualization and statistical analysis
for KEGG-MGC identification results and genome characteristics.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# Set publication-ready style
plt.style.use('default')
sns.set_palette("husl")


class AdvancedMGCVisualizer:
    """Advanced visualization and analysis for MGC results."""
    
    def __init__(self, genome_stats_file: str, mgc_results_file: str):
        """
        Initialize with genome statistics and MGC results.
        
        Args:
            genome_stats_file: Path to genome statistics CSV
            mgc_results_file: Path to MGC results CSV
        """
        self.genome_stats_file = genome_stats_file
        self.mgc_results_file = mgc_results_file
        self.output_dir = os.path.dirname(genome_stats_file)
        
        # Load data
        self.genome_stats = pd.read_csv(genome_stats_file)
        self.mgc_results = pd.read_csv(mgc_results_file)
        
        # Set up matplotlib for publication quality
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
    
    def create_mgc_distribution_plot(self, figsize: Tuple[int, int] = (12, 8)) -> str:
        """
        Create comprehensive MGC distribution visualization.
        
        Args:
            figsize: Figure size tuple
            
        Returns:
            Path to saved figure
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        
        # Sort genomes by MGC count
        sorted_stats = self.genome_stats.sort_values('mgc_count', ascending=True)
        
        # 1. Horizontal bar chart of MGC counts
        ax1.barh(range(len(sorted_stats)), sorted_stats['mgc_count'], 
                color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
        ax1.set_yticks(range(len(sorted_stats)))
        ax1.set_yticklabels([name[:15] + '...' if len(name) > 15 else name 
                           for name in sorted_stats['genome_name']], fontsize=8)
        ax1.set_xlabel('Number of KEGG-MGCs')
        ax1.set_title('KEGG-MGC Distribution Across Genomes')
        ax1.grid(axis='x', alpha=0.3)
        
        # Add value labels on bars
        for i, v in enumerate(sorted_stats['mgc_count']):
            if v > 0:
                ax1.text(v + 0.5, i, str(int(v)), va='center', fontsize=8)
        
        # 2. Histogram of MGC counts
        ax2.hist(self.genome_stats['mgc_count'], bins=15, color='lightcoral', 
                alpha=0.7, edgecolor='black', linewidth=0.5)
        ax2.set_xlabel('Number of KEGG-MGCs')
        ax2.set_ylabel('Number of Genomes')
        ax2.set_title('Distribution of KEGG-MGC Counts')
        ax2.grid(alpha=0.3)
        
        # Add statistics text
        mean_mgc = self.genome_stats['mgc_count'].mean()
        median_mgc = self.genome_stats['mgc_count'].median()
        std_mgc = self.genome_stats['mgc_count'].std()
        ax2.axvline(mean_mgc, color='red', linestyle='--', alpha=0.7, label=f'Mean: {mean_mgc:.1f}')
        ax2.axvline(median_mgc, color='blue', linestyle='--', alpha=0.7, label=f'Median: {median_mgc:.1f}')
        ax2.legend()
        
        # 3. Box plot by genome categories (if applicable)
        if 'total_genes' in self.genome_stats.columns:
            # Create genome size categories
            genome_sizes = pd.cut(self.genome_stats['total_genes'], 
                                bins=3, labels=['Small', 'Medium', 'Large'])
            plot_data = pd.DataFrame({
                'MGC_Count': self.genome_stats['mgc_count'],
                'Genome_Size': genome_sizes
            })
            
            sns.boxplot(data=plot_data, x='Genome_Size', y='MGC_Count', ax=ax3)
            ax3.set_title('MGC Count by Genome Size Category')
            ax3.set_xlabel('Genome Size Category')
            ax3.set_ylabel('Number of KEGG-MGCs')
        
        # 4. Cumulative distribution
        sorted_mgc = np.sort(self.genome_stats['mgc_count'])
        cumulative = np.arange(1, len(sorted_mgc) + 1) / len(sorted_mgc)
        ax4.plot(sorted_mgc, cumulative, marker='o', markersize=4, 
                color='darkgreen', linewidth=2)
        ax4.set_xlabel('Number of KEGG-MGCs')
        ax4.set_ylabel('Cumulative Probability')
        ax4.set_title('Cumulative Distribution of KEGG-MGCs')
        ax4.grid(alpha=0.3)
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, 'mgc_distribution_analysis.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✅ MGC distribution analysis saved to: {output_file}")
        return output_file
    
    def create_correlation_analysis(self, figsize: Tuple[int, int] = (14, 10)) -> str:
        """
        Create comprehensive correlation analysis visualization.
        
        Args:
            figsize: Figure size tuple
            
        Returns:
            Path to saved figure
        """
        # Select numeric columns for correlation
        numeric_cols = self.genome_stats.select_dtypes(include=[np.number]).columns
        correlation_vars = [col for col in numeric_cols if col not in ['file_size_mb']]
        
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
        
        # Define key relationships to analyze
        relationships = [
            ('total_genes', 'mgc_count', 'Total Genes vs MGC Count'),
            ('gene_density', 'mgc_count', 'Gene Density vs MGC Count'),
            ('kegg_annotation_rate', 'mgc_count', 'KEGG Annotation Rate vs MGC Count'),
            ('unique_kegg_pathways', 'mgc_count', 'Unique KEGG Pathways vs MGC Count'),
            ('total_kegg_annotations', 'mgc_count', 'Total KEGG Annotations vs MGC Count'),
            ('annotated_genes', 'mgc_count', 'Annotated Genes vs MGC Count')
        ]
        
        correlation_results = {}
        
        for i, (x_var, y_var, title) in enumerate(relationships):
            if x_var in self.genome_stats.columns and y_var in self.genome_stats.columns:
                # Remove rows with missing values
                plot_data = self.genome_stats[[x_var, y_var]].dropna()
                
                if len(plot_data) > 3:  # Need at least 4 points for meaningful correlation
                    # Scatter plot
                    axes[i].scatter(plot_data[x_var], plot_data[y_var], 
                                  alpha=0.6, s=50, color='darkblue', edgecolors='black', linewidth=0.5)
                    
                    # Calculate correlations
                    pearson_r, pearson_p = pearsonr(plot_data[x_var], plot_data[y_var])
                    spearman_r, spearman_p = spearmanr(plot_data[x_var], plot_data[y_var])
                    
                    # Add trend line
                    if abs(pearson_r) > 0.1:  # Only add trend line if correlation is meaningful
                        z = np.polyfit(plot_data[x_var], plot_data[y_var], 1)
                        p = np.poly1d(z)
                        axes[i].plot(plot_data[x_var], p(plot_data[x_var]), 
                                   "r--", alpha=0.8, linewidth=2)
                    
                    # Format axes
                    axes[i].set_xlabel(x_var.replace('_', ' ').title())
                    axes[i].set_ylabel(y_var.replace('_', ' ').title())
                    axes[i].set_title(f'{title}\\nr={pearson_r:.3f}, p={pearson_p:.2e}')
                    axes[i].grid(alpha=0.3)
                    
                    # Store correlation results
                    correlation_results[f'{x_var}_vs_{y_var}'] = {
                        'pearson_r': pearson_r,
                        'pearson_p': pearson_p,
                        'spearman_r': spearman_r,
                        'spearman_p': spearman_p,
                        'n_samples': len(plot_data)
                    }
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, 'correlation_analysis.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        # Save correlation results
        correlation_df = pd.DataFrame(correlation_results).T
        correlation_file = os.path.join(self.output_dir, 'correlation_results.csv')
        correlation_df.to_csv(correlation_file)
        
        print(f"✅ Correlation analysis saved to: {output_file}")
        print(f"✅ Correlation results saved to: {correlation_file}")
        
        return output_file
    
    def create_heatmap_analysis(self, figsize: Tuple[int, int] = (10, 8)) -> str:
        """
        Create correlation heatmap of genome characteristics.
        
        Args:
            figsize: Figure size tuple
            
        Returns:
            Path to saved figure
        """
        # Select numeric columns
        numeric_cols = self.genome_stats.select_dtypes(include=[np.number]).columns
        heatmap_vars = [col for col in numeric_cols if col not in ['file_size_mb']]
        
        # Create correlation matrix
        correlation_matrix = self.genome_stats[heatmap_vars].corr()
        
        # Create heatmap
        plt.figure(figsize=figsize)
        mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
        
        sns.heatmap(correlation_matrix, mask=mask, annot=True, cmap='RdBu_r', 
                   center=0, square=True, fmt='.3f', cbar_kws={'shrink': 0.8})
        
        plt.title('Correlation Matrix: Genome Characteristics vs MGC Count')
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, 'correlation_heatmap.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✅ Correlation heatmap saved to: {output_file}")
        return output_file
    
    def create_pathway_analysis(self, figsize: Tuple[int, int] = (12, 8)) -> str:
        """
        Create KEGG pathway analysis visualization.
        
        Args:
            figsize: Figure size tuple
            
        Returns:
            Path to saved figure
        """
        # Analyze pathway distribution in MGC results
        pathway_counts = self.mgc_results['pathway'].value_counts()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # 1. Top pathways bar chart
        top_pathways = pathway_counts.head(15)
        ax1.barh(range(len(top_pathways)), top_pathways.values, 
                color='lightgreen', alpha=0.7, edgecolor='black', linewidth=0.5)
        ax1.set_yticks(range(len(top_pathways)))
        ax1.set_yticklabels(top_pathways.index, fontsize=8)
        ax1.set_xlabel('Number of MGCs')
        ax1.set_title('Top 15 KEGG Pathways in MGCs')
        ax1.grid(axis='x', alpha=0.3)
        
        # 2. Pathway distribution histogram
        ax2.hist(pathway_counts.values, bins=20, color='orange', 
                alpha=0.7, edgecolor='black', linewidth=0.5)
        ax2.set_xlabel('Number of MGCs per Pathway')
        ax2.set_ylabel('Number of Pathways')
        ax2.set_title('Distribution of MGCs per Pathway')
        ax2.grid(alpha=0.3)
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, 'pathway_analysis.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✅ Pathway analysis saved to: {output_file}")
        return output_file
    
    def create_pca_analysis(self, figsize: Tuple[int, int] = (10, 8)) -> str:
        """
        Create PCA analysis of genome characteristics.
        
        Args:
            figsize: Figure size tuple
            
        Returns:
            Path to saved figure
        """
        # Select numeric columns for PCA
        numeric_cols = self.genome_stats.select_dtypes(include=[np.number]).columns
        pca_vars = [col for col in numeric_cols if col not in ['file_size_mb']]
        
        # Prepare data
        pca_data = self.genome_stats[pca_vars].dropna()
        
        if len(pca_data) < 3:
            print("Insufficient data for PCA analysis")
            return None
        
        # Standardize data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(pca_data)
        
        # Perform PCA
        pca = PCA()
        pca_result = pca.fit_transform(scaled_data)
        
        # Create plots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # 1. PCA scatter plot
        scatter = ax1.scatter(pca_result[:, 0], pca_result[:, 1], 
                            c=self.genome_stats.loc[pca_data.index, 'mgc_count'],
                            cmap='viridis', alpha=0.7, s=60, edgecolors='black', linewidth=0.5)
        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax1.set_title('PCA of Genome Characteristics')
        ax1.grid(alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax1)
        cbar.set_label('MGC Count')
        
        # 2. Explained variance plot
        ax2.bar(range(1, len(pca.explained_variance_ratio_) + 1), 
               pca.explained_variance_ratio_, alpha=0.7, color='coral', edgecolor='black', linewidth=0.5)
        ax2.set_xlabel('Principal Component')
        ax2.set_ylabel('Explained Variance Ratio')
        ax2.set_title('PCA Explained Variance')
        ax2.grid(alpha=0.3)
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, 'pca_analysis.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"✅ PCA analysis saved to: {output_file}")
        return output_file
    
    def generate_summary_report(self) -> str:
        """
        Generate comprehensive summary report.
        
        Returns:
            Path to saved report
        """
        report_lines = []
        report_lines.append("# KEGG-MGC Analysis Summary Report")
        report_lines.append("=" * 50)
        report_lines.append("")
        
        # Basic statistics
        report_lines.append("## Basic Statistics")
        report_lines.append(f"- Total genomes analyzed: {len(self.genome_stats)}")
        report_lines.append(f"- Total MGCs identified: {self.genome_stats['mgc_count'].sum()}")
        report_lines.append(f"- Average MGCs per genome: {self.genome_stats['mgc_count'].mean():.2f} ± {self.genome_stats['mgc_count'].std():.2f}")
        report_lines.append(f"- Median MGCs per genome: {self.genome_stats['mgc_count'].median():.2f}")
        report_lines.append(f"- Range: {self.genome_stats['mgc_count'].min()} - {self.genome_stats['mgc_count'].max()}")
        report_lines.append("")
        
        # Top and bottom genomes
        sorted_genomes = self.genome_stats.sort_values('mgc_count', ascending=False)
        report_lines.append("## Top 5 Genomes by MGC Count")
        for i, (_, row) in enumerate(sorted_genomes.head().iterrows()):
            report_lines.append(f"{i+1}. {row['genome_name']}: {row['mgc_count']} MGCs")
        report_lines.append("")
        
        report_lines.append("## Bottom 5 Genomes by MGC Count")
        for i, (_, row) in enumerate(sorted_genomes.tail().iterrows()):
            report_lines.append(f"{i+1}. {row['genome_name']}: {row['mgc_count']} MGCs")
        report_lines.append("")
        
        # Correlation analysis
        if 'total_genes' in self.genome_stats.columns:
            correlation_vars = ['total_genes', 'gene_density', 'kegg_annotation_rate', 
                              'unique_kegg_pathways', 'total_kegg_annotations']
            
            report_lines.append("## Correlation with MGC Count")
            for var in correlation_vars:
                if var in self.genome_stats.columns:
                    valid_data = self.genome_stats[[var, 'mgc_count']].dropna()
                    if len(valid_data) > 3:
                        r, p = pearsonr(valid_data[var], valid_data['mgc_count'])
                        report_lines.append(f"- {var.replace('_', ' ').title()}: r = {r:.3f}, p = {p:.2e}")
            report_lines.append("")
        
        # Pathway analysis
        report_lines.append("## KEGG Pathway Analysis")
        pathway_counts = self.mgc_results['pathway'].value_counts()
        report_lines.append(f"- Total unique pathways: {len(pathway_counts)}")
        report_lines.append(f"- Most common pathway: {pathway_counts.index[0]} ({pathway_counts.iloc[0]} MGCs)")
        report_lines.append(f"- Average MGCs per pathway: {pathway_counts.mean():.2f}")
        report_lines.append("")
        
        # Save report
        report_file = os.path.join(self.output_dir, 'mgc_analysis_report.txt')
        with open(report_file, 'w') as f:
            f.write('\n'.join(report_lines))
        
        print(f"✅ Summary report saved to: {report_file}")
        return report_file
    
    def run_full_analysis(self) -> List[str]:
        """
        Run complete analysis pipeline.
        
        Returns:
            List of output file paths
        """
        print("Running comprehensive MGC analysis...")
        
        output_files = []
        
        # Create all visualizations
        output_files.append(self.create_mgc_distribution_plot())
        output_files.append(self.create_correlation_analysis())
        output_files.append(self.create_heatmap_analysis())
        output_files.append(self.create_pathway_analysis())
        
        # PCA analysis (if enough data)
        pca_file = self.create_pca_analysis()
        if pca_file:
            output_files.append(pca_file)
        
        # Generate summary report
        output_files.append(self.generate_summary_report())
        
        print(f"\n✅ Complete analysis finished. Generated {len(output_files)} output files.")
        return output_files


def main():
    """Main function to run advanced MGC visualization."""
    # Define paths
    base_dir = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3'
    genome_stats_file = os.path.join(base_dir, 'genome_statistics.csv')
    mgc_results_file = os.path.join(base_dir, 'potential_groups_w10_filtered.csv')
    
    # Check if genome statistics file exists
    if not os.path.exists(genome_stats_file):
        print(f"Genome statistics file not found: {genome_stats_file}")
        print("Please run genome_statistics_extractor.py first.")
        return
    
    # Initialize visualizer
    visualizer = AdvancedMGCVisualizer(genome_stats_file, mgc_results_file)
    
    # Run full analysis
    output_files = visualizer.run_full_analysis()
    
    return output_files


if __name__ == "__main__":
    main()