"""
Genome Statistics Extraction Module

This module extracts comprehensive genome statistics from annotated genome files
to analyze relationships between genome characteristics and KEGG-MGC identification.
"""

import os
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import re
from collections import defaultdict


class GenomeStatisticsExtractor:
    """Extract and analyze genome statistics from annotated genome files."""
    
    def __init__(self, mgc_results_file: str):
        """
        Initialize with MGC results file path.
        
        Args:
            mgc_results_file: Path to the MGC results CSV file
        """
        self.mgc_results_file = mgc_results_file
        self.mgc_df = pd.read_csv(mgc_results_file)
        self.genome_stats = {}
        
    def extract_genome_name(self, file_path: str) -> str:
        """
        Extract clean genome name from file path.
        
        Args:
            file_path: Path to genome annotation file
            
        Returns:
            Clean genome name
        """
        name = os.path.basename(file_path).lower()
        
        # Remove common suffixes
        suffixes = [
            '_updated_annotated.csv', '_annotated.csv', '.filtered.csv',
            '.filtered', '.csv', '_metabolic.csv'
        ]
        
        for suffix in suffixes:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        
        # Clean up the name
        name = name.strip('._')
        name = name.replace('..', '.')
        name = name.replace('_', ' ')
        name = name.title()
        
        return name
    
    def get_genome_file_info(self, file_path: str) -> Dict:
        """
        Extract basic file information and gene statistics.
        
        Args:
            file_path: Path to genome annotation file
            
        Returns:
            Dictionary with file statistics
        """
        try:
            # Read the annotation file
            df = pd.read_csv(file_path)
            
            # Basic statistics
            stats = {
                'file_path': file_path,
                'genome_name': self.extract_genome_name(file_path),
                'total_genes': len(df),
                'file_size_mb': os.path.getsize(file_path) / (1024 * 1024),
                'columns': list(df.columns)
            }
            
            # Gene annotation statistics
            if 'annotation' in df.columns:
                annotated_genes = df['annotation'].notna().sum()
                stats['annotated_genes'] = annotated_genes
                stats['annotation_rate'] = annotated_genes / len(df) if len(df) > 0 else 0
            
            # Coordinate-based statistics (if available)
            if 'start' in df.columns and 'end' in df.columns:
                df_coords = df.dropna(subset=['start', 'end'])
                if not df_coords.empty:
                    stats['genome_span'] = df_coords['end'].max() - df_coords['start'].min()
                    stats['average_gene_length'] = (df_coords['end'] - df_coords['start']).mean()
                    stats['median_gene_length'] = (df_coords['end'] - df_coords['start']).median()
            
            # Chromosome/scaffold information
            if 'chromosome' in df.columns or 'scaffold' in df.columns:
                chrom_col = 'chromosome' if 'chromosome' in df.columns else 'scaffold'
                stats['num_chromosomes'] = df[chrom_col].nunique()
                stats['genes_per_chromosome'] = df.groupby(chrom_col).size().to_dict()
            
            return stats
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            return {
                'file_path': file_path,
                'genome_name': self.extract_genome_name(file_path),
                'error': str(e)
            }
    
    def calculate_gene_density(self, genome_stats: Dict) -> Dict:
        """
        Calculate gene density metrics.
        
        Args:
            genome_stats: Dictionary with genome statistics
            
        Returns:
            Updated dictionary with density metrics
        """
        if 'genome_span' in genome_stats and 'total_genes' in genome_stats:
            genome_span = genome_stats['genome_span']
            total_genes = genome_stats['total_genes']
            
            if genome_span > 0:
                genome_stats['gene_density'] = total_genes / genome_span
                genome_stats['genes_per_mb'] = total_genes / (genome_span / 1_000_000)
            
        return genome_stats
    
    def extract_kegg_pathways(self, file_path: str) -> Dict:
        """
        Extract KEGG pathway information from annotation file.
        
        Args:
            file_path: Path to genome annotation file
            
        Returns:
            Dictionary with KEGG pathway statistics
        """
        try:
            df = pd.read_csv(file_path)
            
            kegg_stats = {
                'total_kegg_annotations': 0,
                'unique_kegg_pathways': 0,
                'kegg_pathway_distribution': {},
                'kegg_annotation_rate': 0
            }
            
            # Look for KEGG annotations in various columns
            kegg_columns = [col for col in df.columns if 'kegg' in col.lower() or 'ko' in col.lower()]
            
            if kegg_columns:
                kegg_col = kegg_columns[0]  # Use first KEGG column found
                
                # Count KEGG annotations
                kegg_annotated = df[kegg_col].notna()
                kegg_stats['total_kegg_annotations'] = kegg_annotated.sum()
                kegg_stats['kegg_annotation_rate'] = kegg_annotated.sum() / len(df) if len(df) > 0 else 0
                
                # Extract unique pathways
                kegg_entries = df[kegg_col].dropna().astype(str)
                
                # Parse KEGG identifiers (K numbers)
                k_numbers = set()
                for entry in kegg_entries:
                    k_matches = re.findall(r'K\d{5}', entry)
                    k_numbers.update(k_matches)
                
                kegg_stats['unique_kegg_pathways'] = len(k_numbers)
                
                # Create pathway distribution
                pathway_counts = defaultdict(int)
                for entry in kegg_entries:
                    k_matches = re.findall(r'K\d{5}', entry)
                    for k_num in k_matches:
                        pathway_counts[k_num] += 1
                
                kegg_stats['kegg_pathway_distribution'] = dict(pathway_counts)
            
            return kegg_stats
            
        except Exception as e:
            print(f"Error extracting KEGG pathways from {file_path}: {e}")
            return {
                'total_kegg_annotations': 0,
                'unique_kegg_pathways': 0,
                'kegg_pathway_distribution': {},
                'kegg_annotation_rate': 0,
                'error': str(e)
            }
    
    def count_mgc_per_genome(self) -> Dict:
        """
        Count MGCs identified per genome.
        
        Returns:
            Dictionary mapping genome names to MGC counts
        """
        mgc_counts = {}
        
        if 'source_file' in self.mgc_df.columns:
            # Count MGCs per source file
            source_counts = self.mgc_df['source_file'].value_counts()
            
            for source_file, count in source_counts.items():
                genome_name = self.extract_genome_name(source_file)
                mgc_counts[genome_name] = count
        
        return mgc_counts
    
    def process_all_genomes(self) -> pd.DataFrame:
        """
        Process all genomes and extract comprehensive statistics.
        
        Returns:
            DataFrame with all genome statistics
        """
        print("Processing genome statistics...")
        
        # Get unique source files from MGC results
        source_files = self.mgc_df['source_file'].unique()
        
        all_stats = []
        
        for source_file in source_files:
            print(f"Processing: {source_file}")
            
            # Handle relative paths
            if not os.path.isabs(source_file):
                source_file = os.path.join(os.path.dirname(self.mgc_results_file), source_file)
            
            if not os.path.exists(source_file):
                print(f"Warning: File not found: {source_file}")
                continue
            
            # Extract basic genome statistics
            genome_stats = self.get_genome_file_info(source_file)
            
            # Calculate gene density
            genome_stats = self.calculate_gene_density(genome_stats)
            
            # Extract KEGG pathway information
            kegg_stats = self.extract_kegg_pathways(source_file)
            genome_stats.update(kegg_stats)
            
            all_stats.append(genome_stats)
        
        # Convert to DataFrame
        stats_df = pd.DataFrame(all_stats)
        
        # Add MGC counts
        mgc_counts = self.count_mgc_per_genome()
        stats_df['mgc_count'] = stats_df['genome_name'].map(mgc_counts).fillna(0)
        
        return stats_df
    
    def save_statistics(self, output_dir: str = None) -> str:
        """
        Save genome statistics to CSV file.
        
        Args:
            output_dir: Directory to save output (default: same as MGC results)
            
        Returns:
            Path to saved statistics file
        """
        if output_dir is None:
            output_dir = os.path.dirname(self.mgc_results_file)
        
        # Process all genomes
        stats_df = self.process_all_genomes()
        
        # Save to CSV
        output_file = os.path.join(output_dir, 'genome_statistics.csv')
        stats_df.to_csv(output_file, index=False)
        
        print(f"✅ Genome statistics saved to: {output_file}")
        
        # Print summary
        print(f"\nSummary:")
        print(f"- Total genomes processed: {len(stats_df)}")
        print(f"- Total MGCs identified: {stats_df['mgc_count'].sum()}")
        print(f"- Average MGCs per genome: {stats_df['mgc_count'].mean():.2f}")
        print(f"- Genome with most MGCs: {stats_df.loc[stats_df['mgc_count'].idxmax(), 'genome_name']} ({stats_df['mgc_count'].max()} MGCs)")
        
        return output_file


def main():
    """Main function to run genome statistics extraction."""
    # Define paths
    mgc_results_file = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/potential_groups_w10_filtered.csv'
    
    # Initialize extractor
    extractor = GenomeStatisticsExtractor(mgc_results_file)
    
    # Extract and save statistics
    stats_file = extractor.save_statistics()
    
    return stats_file


if __name__ == "__main__":
    main()