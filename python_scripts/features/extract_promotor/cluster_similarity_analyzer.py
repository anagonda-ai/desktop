#!/usr/bin/env python3
"""
Single MGC Cluster Analysis Using Existing Computed Files
Leverages all the pre-computed CSV files for comprehensive analysis
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from itertools import combinations
import time
import argparse
from collections import Counter

class ExistingDataClusterAnalyzer:
    """Analyzer that uses pre-computed CSV files plus FASTA sequences"""
    
    def __init__(self, cluster_id, promoter_dir):
        self.cluster_id = cluster_id
        self.promoter_dir = promoter_dir
        self.sequences = {}
        self.csv_data = {}
        
        print(f"=== Analyzing cluster: {cluster_id} ===")
        print(f"Using existing files from: {promoter_dir}")
    
    def determine_cluster_group(self):
        """Determine cluster group from ID"""
        if self.cluster_id.startswith('BGC'):
            return 'MIBIG'
        elif self.cluster_id.startswith('RAND'):
            return 'Random'
        else:
            return 'KEGG'
    
    def load_all_data(self):
        """Load FASTA and all existing CSV files"""
        print("\nLoading all available data...")
        
        # 1. Load FASTA sequences
        fasta_path = f"{self.promoter_dir}/{self.cluster_id}_promoters.fasta"
        if not os.path.exists(fasta_path):
            print(f"ERROR: FASTA file not found: {fasta_path}")
            return False
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            gene_id = record.id.split('|')[0]
            sequence = str(record.seq).upper()
            self.sequences[gene_id] = sequence
        
        print(f"✓ Loaded {len(self.sequences)} sequences")
        
        # 2. Load all CSV files
        csv_files = {
            'kmer_motifs': f"{self.cluster_id}_kmer_motifs.csv",
            'regulatory_motifs': f"{self.cluster_id}_regulatory_motifs.csv",
            'promoters_analysis': f"{self.cluster_id}_promoters_analysis.csv", 
            'palindromic_motifs': f"{self.cluster_id}_palindromic_motifs.csv",
            'consensus_motifs': f"{self.cluster_id}_consensus_motifs.csv",
            'motif_clusters': f"{self.cluster_id}_motif_clusters.csv"
        }
        
        loaded_files = []
        for key, filename in csv_files.items():
            filepath = f"{self.promoter_dir}/{filename}"
            if os.path.exists(filepath):
                try:
                    df = pd.read_csv(filepath)
                    df.columns = df.columns.str.strip()  # Clean column names
                    self.csv_data[key] = df
                    loaded_files.append(key)
                    print(f"✓ Loaded {key}: {len(df)} rows")
                except Exception as e:
                    print(f"✗ Failed to load {key}: {e}")
                    self.csv_data[key] = None
            else:
                print(f"- Missing: {key}")
                self.csv_data[key] = None
        
        print(f"Successfully loaded {len(loaded_files)}/6 CSV files")
        return len(loaded_files) > 0
    
    def compute_sequence_similarity(self, k=6):
        """Compute pairwise sequence similarities"""
        sequences = list(self.sequences.values())
        gene_ids = list(self.sequences.keys())
        n_seqs = len(sequences)
        
        if n_seqs < 2:
            return None
        
        def kmer_jaccard(seq1, seq2, k):
            if len(seq1) < k or len(seq2) < k:
                return 0.0
            
            kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
            kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
            
            intersection = len(kmers1 & kmers2)
            union = len(kmers1 | kmers2)
            
            return intersection / union if union > 0 else 0.0
        
        # Calculate all pairwise similarities
        similarities = []
        pairwise_results = []
        
        for i, j in combinations(range(n_seqs), 2):
            sim = kmer_jaccard(sequences[i], sequences[j], k)
            similarities.append(sim)
            
            pairwise_results.append({
                'gene1': gene_ids[i],
                'gene2': gene_ids[j],
                'similarity': sim
            })
        
        return {
            'mean_similarity': np.mean(similarities),
            'std_similarity': np.std(similarities),
            'min_similarity': np.min(similarities),
            'max_similarity': np.max(similarities),
            'pairwise_similarities': similarities,
            'pairwise_details': pairwise_results
        }
    
    def analyze_existing_kmer_data(self):
        """Analyze pre-computed k-mer motif data"""
        if self.csv_data.get('kmer_motifs') is None:
            return None
        
        kmer_df = self.csv_data['kmer_motifs']
        
        # Conservation analysis
        highly_conserved = kmer_df[kmer_df['sequence_coverage_pct'] > 80]
        perfectly_conserved = kmer_df[kmer_df['sequence_coverage_pct'] == 100]
        
        # GC content analysis of conserved vs non-conserved
        if 'gc_content' in kmer_df.columns:
            conserved_gc = highly_conserved['gc_content'].mean()
            all_gc = kmer_df['gc_content'].mean()
        else:
            conserved_gc = all_gc = None
        
        return {
            'total_kmers': len(kmer_df),
            'highly_conserved_count': len(highly_conserved),
            'perfectly_conserved_count': len(perfectly_conserved),
            'conservation_ratio': len(highly_conserved) / len(kmer_df),
            'avg_coverage': kmer_df['sequence_coverage_pct'].mean(),
            'max_coverage': kmer_df['sequence_coverage_pct'].max(),
            'conserved_gc_content': conserved_gc,
            'overall_gc_content': all_gc,
            'top_10_conserved': [
                {'rank': row['rank'], 'kmer': row['kmer'], 'coverage': row['sequence_coverage_pct']}
                for _, row in kmer_df.nlargest(10, 'sequence_coverage_pct').iterrows()
            ]
        }
    
    def analyze_existing_regulatory_data(self):
        """Analyze pre-computed regulatory motif data"""
        if self.csv_data.get('regulatory_motifs') is None:
            return None
        
        reg_df = self.csv_data['regulatory_motifs']
        
        # Family distribution
        family_counts = reg_df['motif_family'].value_counts()
        
        # Shannon diversity
        total = len(reg_df)
        shannon_div = -sum((count/total) * np.log2(count/total) for count in family_counts)
        
        # Positional analysis
        positions = reg_df['position'].values
        
        return {
            'total_instances': len(reg_df),
            'unique_families': len(family_counts),
            'family_distribution': family_counts.to_dict(),
            'shannon_diversity': shannon_div,
            'dominant_family': family_counts.index[0],
            'dominant_count': family_counts.iloc[0],
            'avg_position': np.mean(positions),
            'position_std': np.std(positions),
            'sequences_with_motifs': len(reg_df['sequence_id'].unique())
        }
    
    def analyze_existing_palindromic_data(self):
        """Analyze pre-computed palindromic motif data"""
        if self.csv_data.get('palindromic_motifs') is None:
            return None
        
        pal_df = self.csv_data['palindromic_motifs']
        
        if len(pal_df) == 0:
            return {'total_palindromes': 0}
        
        return {
            'total_palindromes': len(pal_df),
            'perfect_conservation': len(pal_df[pal_df['sequence_coverage_pct'] == 100]),
            'avg_coverage': pal_df['sequence_coverage_pct'].mean(),
            'avg_length': pal_df['length'].mean(),
            'length_distribution': pal_df['length'].value_counts().to_dict(),
            'total_occurrences': pal_df['total_occurrences'].sum(),
            'top_palindromes': [
                {'motif': row['motif'], 'coverage': row['sequence_coverage_pct'], 'length': row['length']}
                for _, row in pal_df.nlargest(5, 'sequence_coverage_pct').iterrows()
            ]
        }
    
    def analyze_existing_consensus_data(self):
        """Analyze pre-computed consensus motif data"""
        if self.csv_data.get('consensus_motifs') is None:
            return None
        
        cons_df = self.csv_data['consensus_motifs']
        
        if len(cons_df) == 0:
            return {'total_consensus': 0}
        
        return {
            'total_consensus': len(cons_df),
            'avg_conservation_score': cons_df['conservation_score'].mean(),
            'max_conservation_score': cons_df['conservation_score'].max(),
            'avg_length': cons_df['length'].mean(),
            'consensus_motifs': [
                {'rank': row['rank'], 'motif': row['consensus_motif'], 
                 'score': row['conservation_score'], 'length': row['length'], 'type': row['motif_type']}
                for _, row in cons_df.iterrows()
            ]
        }
    
    def analyze_existing_promoter_data(self):
        """Analyze pre-computed promoter characteristics"""
        if self.csv_data.get('promoters_analysis') is None:
            return None
        
        prom_df = self.csv_data['promoters_analysis']
        
        # Basic statistics for all columns
        stats = {}
        for col in ['length', 'GC_content', 'AT_content', 'TATA_count', 'CAAT_count', 'CpG_count', 'poly_A_tracts', 'poly_T_tracts']:
            if col in prom_df.columns:
                values = prom_df[col]
                stats[col] = {
                    'mean': values.mean(),
                    'std': values.std(),
                    'min': values.min(),
                    'max': values.max(),
                    'values': values.tolist()
                }
        
        # Regulatory element classification
        tata_count = len(prom_df[prom_df['TATA_count'] > 0])
        caat_count = len(prom_df[prom_df['CAAT_count'] > 0])
        both_count = len(prom_df[(prom_df['TATA_count'] > 0) & (prom_df['CAAT_count'] > 0)])
        
        # Density calculations (per kb)
        total_length_kb = prom_df['length'].sum() / 1000
        densities = {
            'tata_per_kb': prom_df['TATA_count'].sum() / total_length_kb,
            'caat_per_kb': prom_df['CAAT_count'].sum() / total_length_kb,
            'cpg_per_kb': prom_df['CpG_count'].sum() / total_length_kb
        }
        
        return {
            'basic_statistics': stats,
            'regulatory_classification': {
                'total_sequences': len(prom_df),
                'tata_containing': tata_count,
                'caat_containing': caat_count,
                'both_elements': both_count,
                'tata_percentage': 100 * tata_count / len(prom_df),
                'caat_percentage': 100 * caat_count / len(prom_df)
            },
            'element_densities': densities,
            'per_sequence_data': prom_df.to_dict('records')
        }
    
    def run_comprehensive_analysis(self):
        """Run comprehensive analysis using all existing data"""
        start_time = time.time()
        
        # Load all data
        if not self.load_all_data():
            print("Failed to load cluster data")
            return None
        
        print(f"\n=== COMPREHENSIVE ANALYSIS ===")
        
        # 1. Sequence similarity (computed from FASTA)
        print("\n1. SEQUENCE SIMILARITY")
        print("-" * 25)
        similarity_results = {}
        
        # Test multiple k-mer sizes
        for k in [4, 5, 6, 7, 8]:
            sim_result = self.compute_sequence_similarity(k)
            if sim_result:
                similarity_results[f'k{k}'] = sim_result
                print(f"  k={k}: Mean similarity = {sim_result['mean_similarity']:.4f}")
        
        # 2. K-mer motif analysis (from existing CSV)
        print("\n2. K-MER MOTIF ANALYSIS (from existing data)")
        print("-" * 50)
        kmer_analysis = self.analyze_existing_kmer_data()
        if kmer_analysis:
            print(f"  Total k-mers: {kmer_analysis['total_kmers']}")
            print(f"  Highly conserved (>80%): {kmer_analysis['highly_conserved_count']}")
            print(f"  Perfect conservation (100%): {kmer_analysis['perfectly_conserved_count']}")
            print(f"  Conservation ratio: {kmer_analysis['conservation_ratio']:.3f}")
        
        # 3. Regulatory motif analysis (from existing CSV)
        print("\n3. REGULATORY MOTIF ANALYSIS (from existing data)")
        print("-" * 55)
        regulatory_analysis = self.analyze_existing_regulatory_data()
        if regulatory_analysis:
            print(f"  Total regulatory instances: {regulatory_analysis['total_instances']}")
            print(f"  Unique motif families: {regulatory_analysis['unique_families']}")
            print(f"  Shannon diversity: {regulatory_analysis['shannon_diversity']:.3f}")
            print(f"  Dominant family: {regulatory_analysis['dominant_family']} ({regulatory_analysis['dominant_count']} instances)")
        
        # 4. Palindromic motif analysis (from existing CSV)
        print("\n4. PALINDROMIC MOTIF ANALYSIS (from existing data)")
        print("-" * 55)
        palindromic_analysis = self.analyze_existing_palindromic_data()
        if palindromic_analysis and palindromic_analysis['total_palindromes'] > 0:
            print(f"  Total palindromes: {palindromic_analysis['total_palindromes']}")
            print(f"  Perfect conservation: {palindromic_analysis['perfect_conservation']}")
            print(f"  Average coverage: {palindromic_analysis['avg_coverage']:.1f}%")
        
        # 5. Consensus motif analysis (from existing CSV)
        print("\n5. CONSENSUS MOTIF ANALYSIS (from existing data)")
        print("-" * 52)
        consensus_analysis = self.analyze_existing_consensus_data()
        if consensus_analysis and consensus_analysis['total_consensus'] > 0:
            print(f"  Total consensus motifs: {consensus_analysis['total_consensus']}")
            print(f"  Average conservation score: {consensus_analysis['avg_conservation_score']:.3f}")
            print(f"  Top motif: {consensus_analysis['consensus_motifs'][0]['motif']}")
        
        # 6. Promoter characteristics (from existing CSV)
        print("\n6. PROMOTER CHARACTERISTICS (from existing data)")
        print("-" * 53)
        promoter_analysis = self.analyze_existing_promoter_data()
        if promoter_analysis:
            stats = promoter_analysis['basic_statistics']
            classification = promoter_analysis['regulatory_classification']
            
            print(f"  GC content: {stats['GC_content']['mean']:.1f}%±{stats['GC_content']['std']:.1f}%")
            print(f"  TATA boxes per sequence: {stats['TATA_count']['mean']:.1f}±{stats['TATA_count']['std']:.1f}")
            print(f"  CAAT boxes per sequence: {stats['CAAT_count']['mean']:.1f}±{stats['CAAT_count']['std']:.1f}")
            print(f"  Sequences with TATA: {classification['tata_containing']}/{classification['total_sequences']} ({classification['tata_percentage']:.1f}%)")
            print(f"  Sequences with CAAT: {classification['caat_containing']}/{classification['total_sequences']} ({classification['caat_percentage']:.1f}%)")
        
        # Combine all results
        comprehensive_results = {
            'cluster_id': self.cluster_id,
            'group': self.determine_cluster_group(),
            'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'processing_time_seconds': time.time() - start_time,
            'data_files_available': {
                'sequences': len(self.sequences),
                'csv_files_loaded': len([k for k, v in self.csv_data.items() if v is not None]),
                'available_files': [k for k, v in self.csv_data.items() if v is not None]
            },
            'sequence_similarity': similarity_results,
            'kmer_motif_analysis': kmer_analysis,
            'regulatory_motif_analysis': regulatory_analysis,
            'palindromic_motif_analysis': palindromic_analysis,
            'consensus_motif_analysis': consensus_analysis,
            'promoter_characteristics': promoter_analysis
        }
        
        # Generate biological summary
        bio_summary = self.generate_biological_summary(comprehensive_results)
        comprehensive_results['biological_summary'] = bio_summary
        
        return comprehensive_results
    
    def generate_biological_summary(self, results):
        """Generate biological interpretation summary"""
        summary = {
            'cluster_characteristics': [],
            'regulatory_features': [],
            'conservation_patterns': [],
            'functional_implications': []
        }
        
        # Sequence similarity insights
        if 'sequence_similarity' in results and 'k6' in results['sequence_similarity']:
            k6_sim = results['sequence_similarity']['k6']['mean_similarity']
            
            if k6_sim > 0.3:
                summary['cluster_characteristics'].append(f"High sequence conservation (similarity: {k6_sim:.3f})")
            elif k6_sim > 0.15:
                summary['cluster_characteristics'].append(f"Moderate sequence conservation (similarity: {k6_sim:.3f})")
            else:
                summary['cluster_characteristics'].append(f"Low sequence conservation (similarity: {k6_sim:.3f})")
        
        # Motif conservation insights
        if 'kmer_motif_analysis' in results and results['kmer_motif_analysis']:
            kmer_data = results['kmer_motif_analysis']
            conservation_ratio = kmer_data['conservation_ratio']
            
            if conservation_ratio > 0.4:
                summary['conservation_patterns'].append("Strong motif conservation across sequences")
            elif conservation_ratio > 0.2:
                summary['conservation_patterns'].append("Moderate motif conservation")
            else:
                summary['conservation_patterns'].append("Weak motif conservation")
            
            if kmer_data['perfectly_conserved_count'] > 50:
                summary['conservation_patterns'].append(f"Many perfectly conserved motifs ({kmer_data['perfectly_conserved_count']})")
        
        # Regulatory element insights
        if 'promoter_characteristics' in results and results['promoter_characteristics']:
            prom_data = results['promoter_characteristics']
            stats = prom_data['basic_statistics']
            classification = prom_data['regulatory_classification']
            
            # GC content classification
            gc_mean = stats['GC_content']['mean']
            if gc_mean < 30:
                summary['regulatory_features'].append(f"Extremely AT-rich promoters ({gc_mean:.1f}% GC)")
            elif gc_mean < 45:
                summary['regulatory_features'].append(f"AT-rich promoters ({gc_mean:.1f}% GC)")
            else:
                summary['regulatory_features'].append(f"Balanced nucleotide composition ({gc_mean:.1f}% GC)")
            
            # TATA box analysis
            tata_mean = stats['TATA_count']['mean']
            if tata_mean > 8:
                summary['regulatory_features'].append(f"Very high TATA box density ({tata_mean:.1f} per sequence)")
            elif tata_mean > 3:
                summary['regulatory_features'].append(f"High TATA box density ({tata_mean:.1f} per sequence)")
            
            # Coverage analysis
            tata_pct = classification['tata_percentage']
            if tata_pct == 100:
                summary['regulatory_features'].append("All sequences contain TATA boxes")
            elif tata_pct > 80:
                summary['regulatory_features'].append("Most sequences contain TATA boxes")
        
        # Functional implications
        if 'regulatory_motif_analysis' in results and results['regulatory_motif_analysis']:
            reg_data = results['regulatory_motif_analysis']
            diversity = reg_data['shannon_diversity']
            
            if diversity < 1.5:
                summary['functional_implications'].append("Low motif diversity suggests specialized function")
            elif diversity > 2.5:
                summary['functional_implications'].append("High motif diversity suggests complex regulation")
            
            # Family-specific insights
            families = reg_data['family_distribution']
            if 'TATA_box' in families and families['TATA_box'] > 50:
                summary['functional_implications'].append("TATA-dominated regulation suggests strong transcriptional control")
        
        return summary
    
    def save_results(self, results, output_dir):
        """Save comprehensive results"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Main comprehensive results
        json_file = f"{output_dir}/{self.cluster_id}_comprehensive_analysis.json"
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        # Summary CSV for easy analysis
        summary_row = {
            'cluster_id': self.cluster_id,
            'group': results['group'],
            'n_sequences': results['data_files_available']['sequences'],
            'csv_files_loaded': results['data_files_available']['csv_files_loaded']
        }
        
        # Add similarity metrics
        if 'sequence_similarity' in results and 'k6' in results['sequence_similarity']:
            k6_data = results['sequence_similarity']['k6']
            summary_row.update({
                'mean_similarity_k6': k6_data['mean_similarity'],
                'std_similarity_k6': k6_data['std_similarity'],
                'max_similarity_k6': k6_data['max_similarity']
            })
        
        # Add motif metrics
        if 'kmer_motif_analysis' in results and results['kmer_motif_analysis']:
            kmer_data = results['kmer_motif_analysis']
            summary_row.update({
                'total_kmers': kmer_data['total_kmers'],
                'highly_conserved_kmers': kmer_data['highly_conserved_count'],
                'conservation_ratio': kmer_data['conservation_ratio']
            })
        
        # Add regulatory metrics
        if 'promoter_characteristics' in results and results['promoter_characteristics']:
            prom_data = results['promoter_characteristics']['basic_statistics']
            summary_row.update({
                'mean_gc_content': prom_data['GC_content']['mean'],
                'mean_tata_count': prom_data['TATA_count']['mean'],
                'mean_caat_count': prom_data['CAAT_count']['mean']
            })
        
        # Save summary
        summary_df = pd.DataFrame([summary_row])
        summary_file = f"{output_dir}/{self.cluster_id}_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        
        print(f"\nResults saved:")
        print(f"  Comprehensive: {json_file}")
        print(f"  Summary: {summary_file}")
        
        return json_file

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Comprehensive MGC cluster analysis using existing files')
    parser.add_argument('cluster_id', help='Cluster ID to analyze (e.g., BGC0000669)')
    parser.add_argument('--promoter-dir', 
                       default='/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_promotor_analysis',
                       help='Directory containing promoter analysis files')
    parser.add_argument('--output-dir',
                       default='/groups/itay_mayrose/alongonda/desktop/mgc_analysis_results',
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    print(f"=== MGC Comprehensive Cluster Analysis ===")
    print(f"Cluster: {args.cluster_id}")
    print(f"Data directory: {args.promoter_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 45)
    
    try:
        # Initialize analyzer
        analyzer = ExistingDataClusterAnalyzer(
            cluster_id=args.cluster_id,
            promoter_dir=args.promoter_dir
        )
        
        # Run comprehensive analysis
        results = analyzer.run_comprehensive_analysis()
        
        if results is None:
            print("Analysis failed")
            sys.exit(1)
        
        # Save results
        output_file = analyzer.save_results(results, args.output_dir)
        
        # Display final summary
        print(f"\n=== FINAL SUMMARY ===")
        print(f"Cluster: {results['cluster_id']} ({results['group']} group)")
        print(f"Processing time: {results['processing_time_seconds']:.2f} seconds")
        print(f"Files used: {results['data_files_available']['csv_files_loaded']}/6 CSV files + FASTA")
        
        # Key metrics
        if 'sequence_similarity' in results and 'k6' in results['sequence_similarity']:
            sim = results['sequence_similarity']['k6']['mean_similarity']
            print(f"Mean sequence similarity (k=6): {sim:.4f}")
        
        if 'kmer_motif_analysis' in results:
            kmer_data = results['kmer_motif_analysis']
            print(f"Motif conservation ratio: {kmer_data['conservation_ratio']:.3f}")
        
        if 'biological_summary' in results:
            print(f"\nBiological insights:")
            bio = results['biological_summary']
            for category, insights in bio.items():
                if insights:
                    print(f"  {category.replace('_', ' ').title()}:")
                    for insight in insights[:2]:  # Show top 2 insights per category
                        print(f"    • {insight}")
        
        print(f"\nDetailed results: {output_file}")
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()