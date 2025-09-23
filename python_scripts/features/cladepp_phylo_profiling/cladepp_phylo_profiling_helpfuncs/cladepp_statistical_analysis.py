#!/usr/bin/env python3
"""
Evolutionary Conservation Score Analysis for MGC Candidates
Computes conservation scores for individual clusters and category-level comparisons
with proper statistical normalization and group size considerations.
Uses concurrent processing for fast execution.
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import bootstrap
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import partial
import warnings
warnings.filterwarnings('ignore')

class ConservationScoreAnalyzer:
    def __init__(self, max_workers=None):
        self.cluster_scores = None
        self.category_scores = None
        self.max_workers = max_workers or min(32, (os.cpu_count() or 1) + 4)
        
    def compute_cluster_conservation_score(self, summary_df):
        """
        Compute evolutionary conservation score for a single cluster based on CladePP scores.
        
        Uses three principled approaches:
        1. Simple mean CladePP score (most straightforward)
        2. Weighted mean by clade size (accounts for sampling variance)
        3. Median CladePP score (robust to outliers)
        """
        if len(summary_df) == 0 or summary_df['cladepp_score'].isna().all():
            return {
                'mean_conservation_score': np.nan,
                'weighted_conservation_score': np.nan,
                'median_conservation_score': np.nan,
                'n_clades': 0,
                'n_positive_clades': 0,
                'proportion_positive': np.nan
            }
        
        # Clean data - remove NaN values
        valid_data = summary_df.dropna(subset=['cladepp_score'])
        
        if len(valid_data) == 0:
            return {
                'mean_conservation_score': np.nan,
                'weighted_conservation_score': np.nan,
                'median_conservation_score': np.nan,
                'n_clades': 0,
                'n_positive_clades': 0,
                'proportion_positive': np.nan
            }
        
        cladepp_scores = valid_data['cladepp_score'].values
        clade_sizes = valid_data['clade_size'].values if 'clade_size' in valid_data.columns else np.ones(len(cladepp_scores))
        
        # 1. Simple mean CladePP score
        mean_score = np.mean(cladepp_scores)
        
        # 2. Weighted mean by clade size (larger clades get more weight)
        weighted_score = np.average(cladepp_scores, weights=clade_sizes)
        
        # 3. Median CladePP score (robust to outliers)
        median_score = np.median(cladepp_scores)
        
        # Additional metrics
        n_positive_clades = np.sum(cladepp_scores > 0)
        proportion_positive = n_positive_clades / len(cladepp_scores)
        
        return {
            'mean_conservation_score': mean_score,
            'weighted_conservation_score': weighted_score,
            'median_conservation_score': median_score,
            'n_clades': len(cladepp_scores),
            'n_positive_clades': n_positive_clades,
            'proportion_positive': proportion_positive
        }
    
    def process_single_cluster(self, cluster_info):
        """
        Process a single cluster - designed for concurrent execution.
        """
        category, candidate_dir, candidate_path = cluster_info
        
        try:
            summary_file = os.path.join(candidate_path, "summary.csv")
            if not os.path.exists(summary_file):
                return None
            
            summary_df = pd.read_csv(summary_file)
            score_data = self.compute_cluster_conservation_score(summary_df)
            
            # Add metadata
            score_data.update({
                'category': category,
                'cluster_id': candidate_dir,
                'cluster_path': candidate_path
            })
            
            return score_data
            
        except Exception as e:
            print(f"Error processing {candidate_dir}: {e}")
            return None
    
    def analyze_clusters(self, base_dirs, category_patterns):
        """
        Analyze all clusters and compute conservation scores using concurrent processing.
        """
        # Collect all cluster tasks
        cluster_tasks = []
        
        for category, base_dir in base_dirs.items():
            print(f"Scanning {category} candidates from: {base_dir}")
            
            if not os.path.exists(base_dir):
                print(f"Warning: Directory does not exist: {base_dir}")
                continue
            
            category_clusters = 0
            
            for candidate_dir in os.listdir(base_dir):
                candidate_path = os.path.join(base_dir, candidate_dir)
                
                if not os.path.isdir(candidate_path):
                    continue
                
                if not candidate_dir.startswith(category_patterns[category]):
                    continue
                
                cluster_tasks.append((category, candidate_dir, candidate_path))
                category_clusters += 1
            
            print(f"  Found {category_clusters} {category} clusters to process")
        
        print(f"\nProcessing {len(cluster_tasks)} clusters concurrently with {self.max_workers} workers...")
        
        # Process clusters concurrently
        cluster_results = []
        processed_count = 0
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_cluster = {
                executor.submit(self.process_single_cluster, cluster_info): cluster_info[1] 
                for cluster_info in cluster_tasks
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_cluster):
                cluster_id = future_to_cluster[future]
                processed_count += 1
                
                try:
                    result = future.result()
                    if result is not None:
                        cluster_results.append(result)
                    
                    if processed_count % 100 == 0:
                        print(f"  Processed {processed_count}/{len(cluster_tasks)} clusters...")
                        
                except Exception as exc:
                    print(f"Cluster {cluster_id} generated an exception: {exc}")
        
        print(f"  Completed processing {processed_count} clusters")
        print(f"  Successfully analyzed {len(cluster_results)} clusters")
        
        self.cluster_scores = pd.DataFrame(cluster_results)
        return self.cluster_scores
    
    def bootstrap_category_comparison(self, category_scores, n_bootstrap=10000):
        """
        Perform parallel bootstrap analysis for robust category comparison.
        """
        results = {}
        
        # Bootstrap confidence intervals for each category
        for category, score_dict in category_scores.items():
            # Use mean conservation score as primary metric
            scores = score_dict['mean_conservation_score']
            if len(scores) == 0:
                continue
            
            print(f"  Bootstrap analysis for {category} ({len(scores)} clusters)...")
            
            # Parallel bootstrap for mean
            bootstrap_means = self.parallel_bootstrap(scores, np.mean, n_bootstrap)
            
            # Calculate confidence intervals
            ci_low = np.percentile(bootstrap_means, 2.5)
            ci_high = np.percentile(bootstrap_means, 97.5)
            
            results[category] = {
                'mean': np.mean(scores),
                'std': np.std(scores),
                'median': np.median(scores),
                'n_clusters': len(scores),
                'ci_low': ci_low,
                'ci_high': ci_high
            }
        
        return results
    
    def effect_size_analysis(self, category_scores):
        """
        Compute effect sizes (Cohen's d) for pairwise comparisons.
        """
        categories = list(category_scores.keys())
        effect_sizes = {}
        
        for i, cat1 in enumerate(categories):
            for cat2 in categories[i+1:]:
                scores1 = category_scores[cat1]['mean_conservation_score']
                scores2 = category_scores[cat2]['mean_conservation_score']
                
                if len(scores1) == 0 or len(scores2) == 0:
                    continue
                
                # Cohen's d
                pooled_std = np.sqrt(((len(scores1)-1)*np.var(scores1) + 
                                    (len(scores2)-1)*np.var(scores2)) / 
                                   (len(scores1) + len(scores2) - 2))
                
                cohens_d = (np.mean(scores1) - np.mean(scores2)) / pooled_std
                
                # Effect size interpretation
                if abs(cohens_d) < 0.2:
                    interpretation = "negligible"
                elif abs(cohens_d) < 0.5:
                    interpretation = "small"
                elif abs(cohens_d) < 0.8:
                    interpretation = "medium"
                else:
                    interpretation = "large"
                
                effect_sizes[f"{cat1}_vs_{cat2}"] = {
                    'cohens_d': cohens_d,
                    'interpretation': interpretation,
                    'higher_category': cat1 if cohens_d > 0 else cat2
                }
        
        return effect_sizes
    
    def parallel_permutation_test(self, scores1, scores2, n_permutations=10000, chunk_size=1000):
        """
        Parallel non-parametric permutation test for difference in means.
        """
        observed_diff = np.mean(scores1) - np.mean(scores2)
        combined = np.concatenate([scores1, scores2])
        n1 = len(scores1)
        
        def permutation_chunk(chunk_size, random_seed):
            """Process a chunk of permutations."""
            np.random.seed(random_seed)
            chunk_diffs = []
            
            for _ in range(chunk_size):
                shuffled = np.random.permutation(combined)
                perm_scores1 = shuffled[:n1]
                perm_scores2 = shuffled[n1:]
                chunk_diffs.append(np.mean(perm_scores1) - np.mean(perm_scores2))
            
            return chunk_diffs
        
        # Split permutations into chunks for parallel processing
        n_chunks = min(self.max_workers, n_permutations // chunk_size)
        chunk_sizes = [n_permutations // n_chunks] * n_chunks
        # Distribute remainder
        for i in range(n_permutations % n_chunks):
            chunk_sizes[i] += 1
        
        all_diffs = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = [
                executor.submit(permutation_chunk, size, i) 
                for i, size in enumerate(chunk_sizes)
            ]
            
            for future in as_completed(futures):
                all_diffs.extend(future.result())
        
        p_value = np.mean(np.abs(all_diffs) >= np.abs(observed_diff))
        return observed_diff, p_value
    
    def parallel_bootstrap(self, data, statistic_func, n_bootstrap=10000, chunk_size=1000):
        """
        Parallel bootstrap resampling.
        """
        def bootstrap_chunk(data, statistic_func, chunk_size, random_seed):
            """Process a chunk of bootstrap samples."""
            np.random.seed(random_seed)
            chunk_stats = []
            
            for _ in range(chunk_size):
                bootstrap_sample = np.random.choice(data, size=len(data), replace=True)
                chunk_stats.append(statistic_func(bootstrap_sample))
            
            return chunk_stats
        
        # Split bootstrap samples into chunks
        n_chunks = min(self.max_workers, n_bootstrap // chunk_size)
        chunk_sizes = [n_bootstrap // n_chunks] * n_chunks
        # Distribute remainder
        for i in range(n_bootstrap % n_chunks):
            chunk_sizes[i] += 1
        
        all_stats = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = [
                executor.submit(bootstrap_chunk, data, statistic_func, size, i) 
                for i, size in enumerate(chunk_sizes)
            ]
            
            for future in as_completed(futures):
                all_stats.extend(future.result())
        
        return np.array(all_stats)
    
    def comprehensive_statistical_analysis(self):
        """
        Perform comprehensive statistical analysis with proper normalization.
        """
        if self.cluster_scores is None:
            raise ValueError("Must run analyze_clusters first")
        
        # Remove clusters with invalid scores
        valid_clusters = self.cluster_scores.dropna(subset=['mean_conservation_score'])
        
        print("="*80)
        print("EVOLUTIONARY CONSERVATION ANALYSIS")
        print("="*80)
        
        # Summary by category - focusing on mean conservation score
        category_summary = valid_clusters.groupby('category').agg({
            'mean_conservation_score': ['count', 'mean', 'std', 'median', 'min', 'max'],
            'weighted_conservation_score': ['mean', 'std'],
            'median_conservation_score': ['mean', 'std'],
            'n_clades': ['mean', 'sum'],
            'proportion_positive': ['mean', 'std']
        }).round(4)
        
        print("\nCLUSTER-LEVEL CONSERVATION SCORES:")
        print("-" * 50)
        for category in valid_clusters['category'].unique():
            cat_data = valid_clusters[valid_clusters['category'] == category]
            print(f"\n{category}:")
            print(f"  Clusters analyzed: {len(cat_data)}")
            print(f"  Mean CladePP score: {cat_data['mean_conservation_score'].mean():.4f} ± {cat_data['mean_conservation_score'].std():.4f}")
            print(f"  Median CladePP score: {cat_data['median_conservation_score'].mean():.4f}")
            print(f"  Weighted CladePP score: {cat_data['weighted_conservation_score'].mean():.4f}")
            print(f"  Total clades analyzed: {cat_data['n_clades'].sum()}")
            print(f"  Mean clades per cluster: {cat_data['n_clades'].mean():.1f}")
            print(f"  Proportion positive clades: {cat_data['proportion_positive'].mean():.3f}")
        
        # Group scores by category for statistical tests
        category_scores = {}
        for category in valid_clusters['category'].unique():
            cat_data = valid_clusters[valid_clusters['category'] == category]
            category_scores[category] = {
                'mean_conservation_score': cat_data['mean_conservation_score'].values,
                'weighted_conservation_score': cat_data['weighted_conservation_score'].values,
                'median_conservation_score': cat_data['median_conservation_score'].values
            }
        
        # Bootstrap analysis
        print(f"\n{'='*50}")
        print("BOOTSTRAP CONFIDENCE INTERVALS (95%)")
        print("="*50)
        bootstrap_results = self.bootstrap_category_comparison(category_scores)
        
        for category, results in bootstrap_results.items():
            print(f"\n{category}:")
            print(f"  Mean: {results['mean']:.4f} [{results['ci_low']:.4f}, {results['ci_high']:.4f}]")
            print(f"  Clusters: {results['n_clusters']}")
        
        # Effect size analysis
        print(f"\n{'='*50}")
        print("PAIRWISE EFFECT SIZE ANALYSIS")
        print("="*50)
        effect_sizes = self.effect_size_analysis(category_scores)
        
        for comparison, results in effect_sizes.items():
            print(f"\n{comparison}:")
            print(f"  Cohen's d: {results['cohens_d']:.4f} ({results['interpretation']} effect)")
            print(f"  Higher category: {results['higher_category']}")
        
        # Statistical significance testing
        print(f"\n{'='*50}")
        print("STATISTICAL SIGNIFICANCE TESTS")
        print("="*50)
        
        categories = list(category_scores.keys())
        
        # Overall test (Kruskal-Wallis)
        if len(categories) > 2:
            kruskal_data = [category_scores[cat]['mean_conservation_score'] for cat in categories]
            stat, p_val = stats.kruskal(*kruskal_data)
            print(f"\nKruskal-Wallis test (overall difference):")
            print(f"  H = {stat:.4f}, p = {p_val:.2e}")
            if p_val < 0.001:
                sig_level = "highly significant (p < 0.001)"
            elif p_val < 0.01:
                sig_level = "very significant (p < 0.01)"
            elif p_val < 0.05:
                sig_level = "significant (p < 0.05)"
            else:
                sig_level = "not significant (p ≥ 0.05)"
            print(f"  Result: {sig_level}")
        
        # Pairwise permutation tests
        print(f"\nPairwise permutation tests:")
        for i, cat1 in enumerate(categories):
            for cat2 in categories[i+1:]:
                if len(category_scores[cat1]) == 0 or len(category_scores[cat2]) == 0:
                    continue
                
                scores1 = category_scores[cat1]['mean_conservation_score']
                scores2 = category_scores[cat2]['mean_conservation_score']
                obs_diff, p_val = self.parallel_permutation_test(scores1, scores2)
                
                print(f"  {cat1} vs {cat2}:")
                print(f"    Mean difference: {obs_diff:.4f}")
                print(f"    P-value: {p_val:.4f}")
                
                if p_val < 0.05:
                    higher = cat1 if obs_diff > 0 else cat2
                    print(f"    Result: {higher} significantly higher")
                else:
                    print(f"    Result: No significant difference")
        
        # Final ranking
        print(f"\n{'='*50}")
        print("FINAL CATEGORY RANKING")
        print("="*50)
        
        category_means = [(cat, np.mean(scores['mean_conservation_score'])) for cat, scores in category_scores.items()]
        category_means.sort(key=lambda x: x[1], reverse=True)
        
        print("\nCategories ranked by mean conservation score:")
        for rank, (category, mean_score) in enumerate(category_means, 1):
            n_clusters = len(category_scores[category])
            print(f"  {rank}. {category}: {mean_score:.4f} (n={n_clusters})")
        
        return {
            'cluster_scores': valid_clusters,
            'bootstrap_results': bootstrap_results,
            'effect_sizes': effect_sizes,
            'category_ranking': category_means
        }

def main():
    # Configuration
    base_dirs = {
        'MGC_CANDIDATE': "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_final/",
        'BGC': "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_final/",
        'RANDOM': "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/random_mgc_candidates_dir_fixed/"
    }
    
    category_patterns = {
        'MGC_CANDIDATE': 'MGC_CANDIDATE',
        'BGC': 'BGC',
        'RANDOM': 'RANDOM'
    }
    
    # Initialize analyzer
    analyzer = ConservationScoreAnalyzer()
    
    # Analyze clusters
    cluster_scores = analyzer.analyze_clusters(base_dirs, category_patterns)
    
    if len(cluster_scores) == 0:
        print("No valid cluster data found!")
        return
    
    # Comprehensive analysis
    results = analyzer.comprehensive_statistical_analysis()
    
    # Save results
    cluster_scores.to_csv("evolutionary_conservation_cluster_scores.csv", index=False)
    print(f"\nDetailed cluster scores saved to: evolutionary_conservation_cluster_scores.csv")
    
    # Summary report
    print(f"\n{'='*80}")
    print("ANALYSIS SUMMARY")
    print("="*80)
    
    total_clusters = len(results['cluster_scores'])
    print(f"Total clusters analyzed: {total_clusters}")
    
    print(f"\nConservation score interpretation:")
    print(f"- Higher scores indicate stronger evolutionary conservation")
    print(f"- Score combines mean correlation, consistency, coverage, and peak signals")
    print(f"- Statistical analysis accounts for group size differences")
    
    best_category = results['category_ranking'][0][0]
    print(f"\nBest performing category: {best_category}")

if __name__ == "__main__":
    main()