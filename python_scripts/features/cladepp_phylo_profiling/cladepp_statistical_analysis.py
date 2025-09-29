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
from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score
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
    
    def find_optimal_threshold_binary(self, scores, labels, target_category, metric='youden'):
        """
        Find optimal threshold for binary classification using various metrics.
        
        Parameters:
        -----------
        scores : array-like
            Conservation scores
        labels : array-like 
            Category labels
        target_category : str
            Category to classify as positive class
        metric : str
            Optimization metric ('youden', 'f1', 'accuracy', 'sensitivity', 'specificity')
            
        Returns:
        --------
        dict with optimal threshold and performance metrics
        """
        # Convert to binary classification problem
        binary_labels = (labels == target_category).astype(int)
        
        # Calculate ROC curve
        fpr, tpr, thresholds = roc_curve(binary_labels, scores)
        roc_auc = auc(fpr, tpr)
        
        # Calculate precision-recall curve
        precision, recall, pr_thresholds = precision_recall_curve(binary_labels, scores)
        pr_auc = auc(recall, precision)
        
        # Find optimal threshold based on different criteria
        optimal_results = {}
        
        if metric == 'youden' or metric == 'all':
            # Youden's J statistic (sensitivity + specificity - 1)
            youden_scores = tpr - fpr
            optimal_idx = np.argmax(youden_scores)
            optimal_threshold = thresholds[optimal_idx]
            
            optimal_results['youden'] = {
                'threshold': optimal_threshold,
                'sensitivity': tpr[optimal_idx],
                'specificity': 1 - fpr[optimal_idx],
                'youden_j': youden_scores[optimal_idx],
                'fpr': fpr[optimal_idx],
                'tpr': tpr[optimal_idx]
            }
        
        if metric == 'f1' or metric == 'all':
            # F1 score optimization
            f1_scores = []
            for thresh in thresholds:
                pred_binary = (scores >= thresh).astype(int)
                if len(np.unique(pred_binary)) == 1:
                    f1_scores.append(0)
                else:
                    f1_scores.append(f1_score(binary_labels, pred_binary))
            
            f1_scores = np.array(f1_scores)
            optimal_idx = np.argmax(f1_scores)
            optimal_threshold = thresholds[optimal_idx]
            
            pred_binary = (scores >= optimal_threshold).astype(int)
            tn = np.sum((pred_binary == 0) & (binary_labels == 0))
            fp = np.sum((pred_binary == 1) & (binary_labels == 0))
            fn = np.sum((pred_binary == 0) & (binary_labels == 1))
            tp = np.sum((pred_binary == 1) & (binary_labels == 1))
            
            optimal_results['f1'] = {
                'threshold': optimal_threshold,
                'f1_score': f1_scores[optimal_idx],
                'precision': tp / (tp + fp) if (tp + fp) > 0 else 0,
                'recall': tp / (tp + fn) if (tp + fn) > 0 else 0,
                'accuracy': (tp + tn) / len(binary_labels)
            }
        
        if metric == 'accuracy' or metric == 'all':
            # Accuracy optimization
            accuracies = []
            for thresh in thresholds:
                pred_binary = (scores >= thresh).astype(int)
                accuracy = np.mean(pred_binary == binary_labels)
                accuracies.append(accuracy)
            
            accuracies = np.array(accuracies)
            optimal_idx = np.argmax(accuracies)
            optimal_threshold = thresholds[optimal_idx]
            
            optimal_results['accuracy'] = {
                'threshold': optimal_threshold,
                'accuracy': accuracies[optimal_idx],
                'sensitivity': tpr[optimal_idx],
                'specificity': 1 - fpr[optimal_idx]
            }
        
        # Add ROC and PR AUC to results
        base_results = {
            'roc_auc': roc_auc,
            'pr_auc': pr_auc,
            'n_positive': np.sum(binary_labels),
            'n_negative': np.sum(1 - binary_labels),
            'target_category': target_category,
            'roc_curve': {'fpr': fpr, 'tpr': tpr, 'thresholds': thresholds},
            'pr_curve': {'precision': precision, 'recall': recall, 'thresholds': pr_thresholds}
        }
        
        if metric == 'all':
            return {**base_results, **optimal_results}
        else:
            return {**base_results, metric: optimal_results[metric]}
    
    def multiclass_threshold_optimization(self, scores, labels, metrics=['youden', 'f1']):
        """
        Perform threshold optimization for multiclass classification using one-vs-rest approach.
        """
        unique_categories = np.unique(labels)
        results = {}
        
        print(f"\n{'='*60}")
        print("THRESHOLD OPTIMIZATION FOR CLASSIFICATION")
        print("="*60)
        
        for category in unique_categories:
            print(f"\nOptimizing thresholds for {category} vs Others:")
            print("-" * 40)
            
            category_results = {}
            
            for metric in metrics:
                result = self.find_optimal_threshold_binary(scores, labels, category, metric)
                category_results[metric] = result[metric]
                
                print(f"\n{metric.upper()} optimization:")
                if metric == 'youden':
                    thresh_info = result[metric]
                    print(f"  Optimal threshold: {thresh_info['threshold']:.4f}")
                    print(f"  Sensitivity: {thresh_info['sensitivity']:.4f}")
                    print(f"  Specificity: {thresh_info['specificity']:.4f}")
                    print(f"  Youden's J: {thresh_info['youden_j']:.4f}")
                elif metric == 'f1':
                    thresh_info = result[metric]
                    print(f"  Optimal threshold: {thresh_info['threshold']:.4f}")
                    print(f"  F1 score: {thresh_info['f1_score']:.4f}")
                    print(f"  Precision: {thresh_info['precision']:.4f}")
                    print(f"  Recall: {thresh_info['recall']:.4f}")
                    print(f"  Accuracy: {thresh_info['accuracy']:.4f}")
                elif metric == 'accuracy':
                    thresh_info = result[metric]
                    print(f"  Optimal threshold: {thresh_info['threshold']:.4f}")
                    print(f"  Accuracy: {thresh_info['accuracy']:.4f}")
                    print(f"  Sensitivity: {thresh_info['sensitivity']:.4f}")
                    print(f"  Specificity: {thresh_info['specificity']:.4f}")
            
            # Add base metrics
            category_results['roc_auc'] = result['roc_auc']
            category_results['pr_auc'] = result['pr_auc']
            category_results['n_positive'] = result['n_positive']
            category_results['n_negative'] = result['n_negative']
            
            print(f"\nDiscriminative Performance:")
            print(f"  ROC AUC: {result['roc_auc']:.4f}")
            print(f"  PR AUC: {result['pr_auc']:.4f}")
            print(f"  Positive samples: {result['n_positive']}")
            print(f"  Negative samples: {result['n_negative']}")
            
            results[category] = category_results
        
        # Summary of best thresholds
        print(f"\n{'='*60}")
        print("OPTIMAL THRESHOLD SUMMARY")
        print("="*60)
        
        for category in unique_categories:
            print(f"\n{category}:")
            for metric in metrics:
                thresh = results[category][metric]['threshold']
                print(f"  Best {metric} threshold: {thresh:.4f}")
        
        return results
    
    def generate_threshold_recommendations(self, threshold_results):
        """
        Generate practical recommendations for threshold selection.
        """
        print(f"\n{'='*60}")
        print("THRESHOLD SELECTION RECOMMENDATIONS")
        print("="*60)
        
        for category, results in threshold_results.items():
            print(f"\n{category}:")
            print("-" * 30)
            
            # Get ROC AUC for discriminative power assessment
            roc_auc = results['roc_auc']
            
            if roc_auc < 0.6:
                discriminative_power = "Poor"
                recommendation = "Not suitable for classification"
            elif roc_auc < 0.7:
                discriminative_power = "Fair"
                recommendation = "Use with caution, consider additional features"
            elif roc_auc < 0.8:
                discriminative_power = "Good"
                recommendation = "Suitable for classification"
            elif roc_auc < 0.9:
                discriminative_power = "Excellent"
                recommendation = "Strong classifier"
            else:
                discriminative_power = "Outstanding"
                recommendation = "Exceptional classifier"
            
            print(f"  Discriminative Power: {discriminative_power} (AUC = {roc_auc:.3f})")
            print(f"  Recommendation: {recommendation}")
            
            # Recommend best threshold based on use case
            if 'youden' in results and 'f1' in results:
                youden_thresh = results['youden']['threshold']
                f1_thresh = results['f1']['threshold']
                
                print(f"\n  Threshold Recommendations:")
                print(f"    Balanced classification: {youden_thresh:.4f} (Youden's J)")
                print(f"    Precision-recall balance: {f1_thresh:.4f} (F1 optimization)")
                
                if roc_auc >= 0.7:
                    if abs(youden_thresh - f1_thresh) < 0.1:
                        print(f"    Consensus threshold: ~{np.mean([youden_thresh, f1_thresh]):.4f}")
                    else:
                        print(f"    Choose based on cost of false positives vs false negatives")
        
        return threshold_results
    
    
    def create_threshold_summary_table(self, threshold_results):
        """
        Create a summary table of optimal thresholds for easy reference.
        """
        summary_data = []
        
        for category, results in threshold_results.items():
            row = {
                'Category': category,
                'ROC_AUC': results['roc_auc'],
                'PR_AUC': results['pr_auc'],
                'N_Positive': results['n_positive'],
                'N_Negative': results['n_negative']
            }
            
            # Add threshold information
            for metric in ['youden', 'f1', 'accuracy']:
                if metric in results:
                    row[f'{metric}_threshold'] = results[metric]['threshold']
                    
                    if metric == 'youden':
                        row[f'{metric}_sensitivity'] = results[metric]['sensitivity']
                        row[f'{metric}_specificity'] = results[metric]['specificity']
                        row[f'{metric}_youden_j'] = results[metric]['youden_j']
                    elif metric == 'f1':
                        row[f'{metric}_f1_score'] = results[metric]['f1_score']
                        row[f'{metric}_precision'] = results[metric]['precision']
                        row[f'{metric}_recall'] = results[metric]['recall']
                    elif metric == 'accuracy':
                        row[f'{metric}_accuracy'] = results[metric]['accuracy']
            
            summary_data.append(row)
        
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.round(4)
        
        return summary_df
    
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
        
        # Threshold optimization for classification
        print(f"\n{'='*50}")
        print("THRESHOLD OPTIMIZATION FOR ML CLASSIFICATION")
        print("="*50)
        
        # Prepare data for binary classification: BGC/MGC_CANDIDATE = positive, RANDOM = negative
        scores = valid_clusters['mean_conservation_score'].values
        original_labels = valid_clusters['category'].values
        
        # Convert to binary classification
        binary_labels = np.where(
            (original_labels == 'BGC') | (original_labels == 'MGC_CANDIDATE'), 
            'POSITIVE', 
            'NEGATIVE'
        )
        
        print("Binary Classification Setup:")
        print(f"  Positive class: BGC + MGC_CANDIDATE")
        print(f"  Negative class: RANDOM")
        print(f"  Total positive samples: {np.sum(binary_labels == 'POSITIVE')}")
        print(f"  Total negative samples: {np.sum(binary_labels == 'NEGATIVE')}")
        
        # Perform binary threshold optimization
        threshold_results = self.find_optimal_threshold_binary(
            scores, binary_labels, 'POSITIVE', metric='all'
        )
        
        # Print detailed binary classification results
        print(f"\n{'='*60}")
        print("BINARY CLASSIFICATION THRESHOLD OPTIMIZATION")
        print("="*60)
        
        print(f"\nROC Analysis:")
        print(f"  ROC AUC: {threshold_results['roc_auc']:.4f}")
        print(f"  PR AUC: {threshold_results['pr_auc']:.4f}")
        
        # Youden's J optimization
        if 'youden' in threshold_results:
            youden = threshold_results['youden']
            print(f"\nYouden's J Optimization (Balanced Classification):")
            print(f"  Optimal threshold: {youden['threshold']:.4f}")
            print(f"  Sensitivity (True Positive Rate): {youden['sensitivity']:.4f}")
            print(f"  Specificity (True Negative Rate): {youden['specificity']:.4f}")
            print(f"  Youden's J statistic: {youden['youden_j']:.4f}")
            print(f"  False Positive Rate: {youden['fpr']:.4f}")
        
        # F1 optimization
        if 'f1' in threshold_results:
            f1 = threshold_results['f1']
            print(f"\nF1 Score Optimization (Precision-Recall Balance):")
            print(f"  Optimal threshold: {f1['threshold']:.4f}")
            print(f"  F1 score: {f1['f1_score']:.4f}")
            print(f"  Precision: {f1['precision']:.4f}")
            print(f"  Recall (Sensitivity): {f1['recall']:.4f}")
            print(f"  Accuracy: {f1['accuracy']:.4f}")
        
        # Accuracy optimization
        if 'accuracy' in threshold_results:
            acc = threshold_results['accuracy']
            print(f"\nAccuracy Optimization:")
            print(f"  Optimal threshold: {acc['threshold']:.4f}")
            print(f"  Accuracy: {acc['accuracy']:.4f}")
            print(f"  Sensitivity: {acc['sensitivity']:.4f}")
            print(f"  Specificity: {acc['specificity']:.4f}")
        
        # Practical recommendations
        print(f"\n{'='*60}")
        print("CLASSIFICATION RECOMMENDATIONS")
        print("="*60)
        
        roc_auc = threshold_results['roc_auc']
        if roc_auc < 0.6:
            discriminative_power = "Poor"
            recommendation = "Not suitable for classification"
        elif roc_auc < 0.7:
            discriminative_power = "Fair" 
            recommendation = "Use with caution, consider additional features"
        elif roc_auc < 0.8:
            discriminative_power = "Good"
            recommendation = "Suitable for classification"
        elif roc_auc < 0.9:
            discriminative_power = "Excellent"
            recommendation = "Strong classifier"
        else:
            discriminative_power = "Outstanding"
            recommendation = "Exceptional classifier"
        
        print(f"\nDiscriminative Power: {discriminative_power} (AUC = {roc_auc:.3f})")
        print(f"Recommendation: {recommendation}")
        
        if roc_auc >= 0.7:
            print(f"\nRecommended Thresholds:")
            if 'youden' in threshold_results and 'f1' in threshold_results:
                youden_thresh = threshold_results['youden']['threshold']
                f1_thresh = threshold_results['f1']['threshold']
                print(f"  Balanced classification: {youden_thresh:.4f} (Youden's J)")
                print(f"  Precision-recall balance: {f1_thresh:.4f} (F1 optimization)")
                
                if abs(youden_thresh - f1_thresh) < 0.1:
                    consensus_thresh = np.mean([youden_thresh, f1_thresh])
                    print(f"  Consensus threshold: ~{consensus_thresh:.4f}")
                else:
                    print(f"  Choose based on cost of false positives vs false negatives")
        
        # Create simplified summary table
        threshold_summary_table = pd.DataFrame([{
            'Metric': 'Youden_J',
            'Threshold': threshold_results.get('youden', {}).get('threshold', np.nan),
            'Sensitivity': threshold_results.get('youden', {}).get('sensitivity', np.nan),
            'Specificity': threshold_results.get('youden', {}).get('specificity', np.nan),
            'ROC_AUC': threshold_results['roc_auc']
        }, {
            'Metric': 'F1_Score',
            'Threshold': threshold_results.get('f1', {}).get('threshold', np.nan),
            'F1_Score': threshold_results.get('f1', {}).get('f1_score', np.nan),
            'Precision': threshold_results.get('f1', {}).get('precision', np.nan),
            'Recall': threshold_results.get('f1', {}).get('recall', np.nan),
            'ROC_AUC': threshold_results['roc_auc']
        }, {
            'Metric': 'Accuracy',
            'Threshold': threshold_results.get('accuracy', {}).get('threshold', np.nan),
            'Accuracy': threshold_results.get('accuracy', {}).get('accuracy', np.nan),
            'Sensitivity': threshold_results.get('accuracy', {}).get('sensitivity', np.nan),
            'Specificity': threshold_results.get('accuracy', {}).get('specificity', np.nan),
            'ROC_AUC': threshold_results['roc_auc']
        }]).round(4)
        
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
            'category_ranking': category_means,
            'threshold_optimization': threshold_results,
            'threshold_summary_table': threshold_summary_table
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
    
    # ML Classification Summary
    print(f"\n{'='*80}")
    print("BINARY CLASSIFICATION SUMMARY")
    print("="*80)
    
    threshold_results = results['threshold_optimization']
    print(f"\nBinary Classification Performance:")
    print(f"  ROC AUC: {threshold_results['roc_auc']:.4f}")
    print(f"  PR AUC: {threshold_results['pr_auc']:.4f}")
    print(f"  Positive samples (BGC + MGC_CANDIDATE): {threshold_results['n_positive']}")
    print(f"  Negative samples (RANDOM): {threshold_results['n_negative']}")
    
    print(f"\nDetailed Threshold Performance:")
    print("-" * 60)
    
    if 'youden' in threshold_results:
        youden = threshold_results['youden']
        print(f"\n1. BALANCED CLASSIFICATION (Youden's J):")
        print(f"   Threshold: {youden['threshold']:.4f}")
        print(f"   Sensitivity (True Positive Rate): {youden['sensitivity']:.4f}")
        print(f"   Specificity (True Negative Rate): {youden['specificity']:.4f}")
        print(f"   False Positive Rate: {youden['fpr']:.4f}")
        print(f"   Youden's J statistic: {youden['youden_j']:.4f}")
    
    if 'f1' in threshold_results:
        f1 = threshold_results['f1']
        print(f"\n2. PRECISION-RECALL OPTIMIZATION (F1 Score):")
        print(f"   Threshold: {f1['threshold']:.4f}")
        print(f"   F1 Score: {f1['f1_score']:.4f}")
        print(f"   Precision: {f1['precision']:.4f}")
        print(f"   Recall (Sensitivity): {f1['recall']:.4f}")
        print(f"   Accuracy: {f1['accuracy']:.4f}")
    
    if 'accuracy' in threshold_results:
        acc = threshold_results['accuracy']
        print(f"\n3. ACCURACY OPTIMIZATION:")
        print(f"   Threshold: {acc['threshold']:.4f}")
        print(f"   Accuracy: {acc['accuracy']:.4f}")
        print(f"   Sensitivity: {acc['sensitivity']:.4f}")
        print(f"   Specificity: {acc['specificity']:.4f}")
    
    # Best precision-recall summary
    print(f"\n{'='*60}")
    print("PRECISION-RECALL SUMMARY")
    print("="*60)
    
    if 'f1' in threshold_results:
        f1 = threshold_results['f1']
        print(f"\nBest Precision-Recall Performance:")
        print(f"  Optimal Threshold: {f1['threshold']:.4f}")
        print(f"  At this threshold:")
        print(f"    - Precision: {f1['precision']:.4f} ({f1['precision']*100:.1f}% of positive predictions are correct)")
        print(f"    - Recall: {f1['recall']:.4f} ({f1['recall']*100:.1f}% of actual positives are found)")
        print(f"    - F1 Score: {f1['f1_score']:.4f} (harmonic mean of precision and recall)")
        
        # Calculate additional metrics for interpretation
        total_samples = threshold_results['n_positive'] + threshold_results['n_negative']
        print(f"\n  Performance Interpretation:")
        print(f"    - Out of {threshold_results['n_positive']} true positives, you'll find {int(f1['recall'] * threshold_results['n_positive'])}")
        print(f"    - When you predict positive, you'll be right {f1['precision']*100:.1f}% of the time")
        print(f"    - Overall accuracy: {f1['accuracy']*100:.1f}% of all {total_samples} samples")

if __name__ == "__main__":
    main()