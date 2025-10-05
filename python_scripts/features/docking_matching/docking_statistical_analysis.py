#!/usr/bin/env python3
"""
Statistically rigorous analysis of LightDock all-vs-all docking results.
Accounts for n×n matrix structure and cluster size effects.
Text-only output with concurrent processing.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from multiprocessing import Pool, cpu_count
from sklearn.linear_model import HuberRegressor
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, VotingClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, classification_report, precision_score, recall_score, f1_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
try:
    from imblearn.over_sampling import SMOTE
    from imblearn.under_sampling import RandomUnderSampler
    from imblearn.pipeline import Pipeline as ImbPipeline
    IMBLEARN_AVAILABLE = True
except ImportError:
    print("Warning: imblearn not available. Using basic class balancing.")
    IMBLEARN_AVAILABLE = False
    # Create sklearn-compatible dummy classes
    class SMOTE:
        def __init__(self, random_state=None, **kwargs):
            self.random_state = random_state
        def fit_resample(self, X, y):
            return X, y
        def get_params(self, deep=True):
            return {'random_state': self.random_state}
        def set_params(self, **params):
            for key, value in params.items():
                setattr(self, key, value)
            return self
    
    class RandomUnderSampler:
        def __init__(self, random_state=None, **kwargs):
            self.random_state = random_state
        def fit_resample(self, X, y):
            return X, y
        def get_params(self, deep=True):
            return {'random_state': self.random_state}
        def set_params(self, **params):
            for key, value in params.items():
                setattr(self, key, value)
            return self
    
    class ImbPipeline:
        def __init__(self, steps):
            self.steps = steps
            self.named_steps = dict(steps)
        def fit(self, X, y):
            return self
        def predict(self, X):
            return np.zeros(len(X))
        def predict_proba(self, X):
            return np.column_stack([np.ones(len(X)), np.zeros(len(X))])
        def get_params(self, deep=True):
            params = {}
            for name, step in self.steps:
                if hasattr(step, 'get_params'):
                    params[f'{name}__'] = step.get_params(deep=deep)
                else:
                    params[f'{name}__'] = {}
            return params
        def set_params(self, **params):
            for param_name, param_value in params.items():
                step_name = param_name.split('__')[0]
                if step_name in self.named_steps:
                    step = self.named_steps[step_name]
                    if hasattr(step, 'set_params'):
                        step.set_params(**{param_name.split('__', 1)[1]: param_value})
            return self
import warnings
warnings.filterwarnings('ignore')

def collect_lightdock_results_from_cluster(args):
    """Collect all LightDock results from a cluster directory"""
    cluster_dir, category = args
    
    try:
        cluster_name = cluster_dir.name
        
        # Find all completed docking results
        results = []
        result_dirs = list(cluster_dir.glob("*_vs_*"))
        
        if not result_dirs:
            return None
            
        for result_dir in result_dirs:
            # Check if analysis completed successfully
            analysis_complete = result_dir / "analysis_complete.flag"
            best_score_file = result_dir / "best_score.txt"
            
            if analysis_complete.exists() and best_score_file.exists():
                try:
                    with open(best_score_file, 'r') as f:
                        content = f.read()
                        # Extract score from "Best docking score: X.X" format
                        if "Best docking score:" in content:
                            score_line = [line for line in content.split('\n') if "Best docking score:" in line][0]
                            score = float(score_line.split(':')[1].strip())
                        else:
                            continue  # Skip if no valid score
                    
                    # Parse protein names from directory name
                    parts = result_dir.name.split('_vs_')
                    if len(parts) == 2:
                        receptor = parts[0]
                        ligand = parts[1]
                        
                        results.append({
                            'query': receptor,
                            'target': ligand,
                            'docking_score': score,
                            'is_self': receptor == ligand
                        })
                        
                except Exception:
                    continue  # Skip failed extractions
        
        if not results:
            return None
            
        return {
            'cluster_name': cluster_name,
            'category': category,
            'results': results
        }
        
    except Exception as e:
        print(f"Error processing {cluster_dir}: {e}")
        return None

def process_cluster_results(cluster_data):
    """Process collected results into statistical metrics"""
    if not cluster_data:
        return None
        
    try:
        results_df = pd.DataFrame(cluster_data['results'])
        
        # Get unique proteins
        unique_proteins = set(results_df['query'].unique()) | set(results_df['target'].unique())
        n = len(unique_proteins)
        
        # All scores (including self-comparisons if they exist)
        all_scores = results_df['docking_score'].values
        
        # Non-self comparisons
        non_self_mask = ~results_df['is_self']
        non_self_scores = results_df.loc[non_self_mask, 'docking_score'].values
        
        if len(non_self_scores) == 0:
            return None
            
        # Expected number of comparisons
        expected_total = n * n
        expected_non_self = n * (n - 1)
        
        return {
            'category': cluster_data['category'],
            'name': cluster_data['cluster_name'],
            'n': n,
            'n_squared': n * n,
            'n_comparisons_total': len(results_df),
            'n_comparisons_non_self': non_self_mask.sum(),
            
            # Raw sums (for accurate group averaging)
            'sum_scores_all': np.sum(all_scores),
            'sum_scores_non_self': np.sum(non_self_scores),
            
            # Means (note: more negative = better binding)
            'mean_score_all': np.mean(all_scores),
            'mean_score_non_self': np.mean(non_self_scores),
            
            # Variance components
            'var_score_non_self': np.var(non_self_scores),
            'std_score_non_self': np.std(non_self_scores),
            
            # Distribution characteristics
            'median_score_non_self': np.median(non_self_scores),
            'q25_score': np.percentile(non_self_scores, 25),
            'q75_score': np.percentile(non_self_scores, 75),
            'min_score': np.min(non_self_scores),  # Best (most negative) score
            'max_score': np.max(non_self_scores),  # Worst (least negative) score
            
            # Binding strength metrics (more negative = stronger)
            'fraction_strong_binders': np.mean(non_self_scores < -50),
            'fraction_moderate_binders': np.mean((non_self_scores < -25) & (non_self_scores >= -50)),
            'fraction_weak_binders': np.mean(non_self_scores >= -25),
            
            # For size effect modeling
            'log_n': np.log(n),
            'sqrt_n': np.sqrt(n)
        }
        
    except Exception as e:
        print(f"Error processing cluster {cluster_data['cluster_name']}: {e}")
        return None

def calculate_size_normalized_scores(df):
    """Apply statistical normalization for cluster size effects"""
    
    # 1. Model size effect on random clusters to establish null distribution
    random_df = df[df['category'] == 'RANDOM'].copy()
    
    if len(random_df) > 10:  # Need enough data for modeling
        # Fit regression model: score ~ f(n)
        X = random_df[['n', 'n_squared']].values
        y = random_df['mean_score_non_self'].values
        
        # Robust regression to handle outliers
        model = HuberRegressor()
        model.fit(X, y)
        
        # Predict expected random score for each cluster size
        df['expected_random'] = model.predict(df[['n', 'n_squared']].values)
    else:
        # Fallback: use mean of random clusters
        df['expected_random'] = random_df['mean_score_non_self'].mean() if len(random_df) > 0 else 0.0
    
    # 2. Calculate enrichment scores (for docking: more negative is better)
    # We want to see how much MORE negative (better) than random
    df['enrichment_score'] = df['mean_score_non_self'] - df['expected_random']
    
    # 3. Z-score normalization within size bins
    df['size_bin'] = pd.qcut(df['n'], q=min(5, len(df)//10), labels=False, duplicates='drop')
    
    def z_normalize(group):
        if len(group) > 1:
            group['z_score'] = (group['mean_score_non_self'] - group['mean_score_non_self'].mean()) / group['mean_score_non_self'].std()
        else:
            group['z_score'] = 0
        return group
    
    df = df.groupby('size_bin').apply(z_normalize)
    
    # 4. Calculate effect size (Cohen's d) relative to random clusters of similar size
    def calculate_effect_size(row):
        similar_random = random_df[np.abs(random_df['n'] - row['n']) <= 1]
        if len(similar_random) > 0:
            random_mean = similar_random['mean_score_non_self'].mean()
            random_std = similar_random['mean_score_non_self'].std()
            if random_std > 0:
                return (row['mean_score_non_self'] - random_mean) / random_std
        return 0
    
    df['effect_size'] = df.apply(calculate_effect_size, axis=1)
    
    return df

def weighted_group_statistics(df):
    """Calculate properly weighted group statistics"""
    results = {}
    
    for category in df['category'].unique():
        cat_df = df[df['category'] == category]
        
        # Total comparisons in this category
        total_comparisons = cat_df['n_comparisons_non_self'].sum()
        total_clusters = len(cat_df)
        
        # Weighted mean (by number of comparisons)
        if total_comparisons > 0:
            weights = cat_df['n_comparisons_non_self'] / total_comparisons
            weighted_mean_score = np.sum(weights * cat_df['mean_score_non_self'])
            
            # Alternative: Direct calculation from sums
            direct_mean_score = cat_df['sum_scores_non_self'].sum() / total_comparisons
            
            # Weighted variance
            weighted_var = np.sum(weights * (cat_df['mean_score_non_self'] - weighted_mean_score)**2)
            
            # Bootstrap confidence intervals
            bootstrap_means = []
            for _ in range(1000):
                sample = cat_df.sample(n=len(cat_df), replace=True)
                if sample['n_comparisons_non_self'].sum() > 0:
                    boot_mean = sample['sum_scores_non_self'].sum() / sample['n_comparisons_non_self'].sum()
                    bootstrap_means.append(boot_mean)
            
            if bootstrap_means:
                ci_lower = np.percentile(bootstrap_means, 2.5)
                ci_upper = np.percentile(bootstrap_means, 97.5)
            else:
                ci_lower = ci_upper = weighted_mean_score
        else:
            weighted_mean_score = direct_mean_score = 0
            weighted_var = 0
            ci_lower = ci_upper = 0
        
        results[category] = {
            'n_clusters': total_clusters,
            'total_proteins': cat_df['n'].sum(),
            'total_comparisons': total_comparisons,
            'weighted_mean_score': weighted_mean_score,
            'direct_mean_score': direct_mean_score,
            'weighted_std': np.sqrt(weighted_var) if weighted_var > 0 else 0,
            'ci_95_lower': ci_lower,
            'ci_95_upper': ci_upper,
            'mean_enrichment': cat_df['enrichment_score'].mean(),
            'median_enrichment': cat_df['enrichment_score'].median(),
            'mean_cluster_size': cat_df['n'].mean(),
            'std_cluster_size': cat_df['n'].std(),
            'fraction_strong_binders': cat_df['fraction_strong_binders'].mean(),
            'fraction_moderate_binders': cat_df['fraction_moderate_binders'].mean(),
            'fraction_weak_binders': cat_df['fraction_weak_binders'].mean(),
            'best_score_overall': cat_df['min_score'].min()  # Most negative = best
        }
    
    return results

def statistical_tests(df, group_stats):
    """Perform rigorous statistical tests with size correction"""
    
    test_results = {}
    
    # 1. Permutation test (most robust for this data structure)
    def permutation_test(group1_df, group2_df, n_permutations=10000):
        # Observed difference in weighted means
        obs_diff = (group1_df['sum_scores_non_self'].sum() / group1_df['n_comparisons_non_self'].sum() - 
                   group2_df['sum_scores_non_self'].sum() / group2_df['n_comparisons_non_self'].sum())
        
        # Combine data
        combined = pd.concat([group1_df, group2_df])
        n1 = len(group1_df)
        
        # Permutation
        perm_diffs = []
        for _ in range(n_permutations):
            shuffled = combined.sample(frac=1, replace=False)
            perm_g1 = shuffled.iloc[:n1]
            perm_g2 = shuffled.iloc[n1:]
            
            sum1 = perm_g1['n_comparisons_non_self'].sum()
            sum2 = perm_g2['n_comparisons_non_self'].sum()
            
            if sum1 > 0 and sum2 > 0:
                perm_diff = (perm_g1['sum_scores_non_self'].sum() / sum1 - 
                            perm_g2['sum_scores_non_self'].sum() / sum2)
                perm_diffs.append(perm_diff)
        
        if perm_diffs:
            p_value = np.mean(np.array(perm_diffs) >= obs_diff)
        else:
            p_value = 1.0
            
        return obs_diff, p_value
    
    # 2. ANOVA-style test (controlling for cluster size)
    from scipy.stats import f_oneway
    
    categories = df['category'].unique()
    
    # Test enrichment scores (size-normalized)
    if len(categories) >= 2:
        groups_enrichment = [df[df['category'] == cat]['enrichment_score'].values for cat in categories]
        # Remove empty groups
        groups_enrichment = [g for g in groups_enrichment if len(g) > 0]
        if len(groups_enrichment) >= 2:
            f_stat, p_anova = f_oneway(*groups_enrichment)
            test_results['ANOVA_enrichment'] = {'F': f_stat, 'p': p_anova}
    
    # Pairwise comparisons
    comparisons = []
    for cat1 in categories:
        for cat2 in categories:
            if cat1 != cat2:
                g1 = df[df['category'] == cat1]
                g2 = df[df['category'] == cat2]
                
                if len(g1) > 0 and len(g2) > 0:
                    # Permutation test
                    diff, p_perm = permutation_test(g1, g2)
                    
                    # Mann-Whitney on enrichment scores
                    try:
                        u_stat, p_mw = stats.mannwhitneyu(
                            g1['enrichment_score'], 
                            g2['enrichment_score'],
                            alternative='two-sided'
                        )
                    except:
                        p_mw = 1.0
                    
                    comparisons.append({
                        'comparison': f"{cat1} vs {cat2}",
                        'mean_diff': group_stats[cat1]['weighted_mean_score'] - group_stats[cat2]['weighted_mean_score'],
                        'enrichment_diff': g1['enrichment_score'].mean() - g2['enrichment_score'].mean(),
                        'p_permutation': p_perm,
                        'p_mannwhitney': p_mw,
                        'significant': p_perm < 0.05
                    })
    
    test_results['pairwise'] = comparisons
    return test_results

def print_detailed_results(df, group_stats, test_results):
    """Print comprehensive textual results"""
    
    print("=" * 80)
    print("LIGHTDOCK DOCKING SCORES - COMPREHENSIVE STATISTICAL ANALYSIS")
    print("=" * 80)
    
    # Dataset overview
    print("\nDATASET OVERVIEW")
    print("-" * 40)
    print(f"Total clusters analyzed: {len(df)}")
    print(f"Total protein comparisons: {df['n_comparisons_non_self'].sum():,}")
    print(f"Total unique proteins: {df['n'].sum()}")
    
    print("\nCluster distribution by category:")
    for cat in sorted(df['category'].unique()):
        count = len(df[df['category'] == cat])
        total_proteins = df[df['category'] == cat]['n'].sum()
        total_comparisons = df[df['category'] == cat]['n_comparisons_non_self'].sum()
        print(f"  {cat}: {count} clusters, {total_proteins} proteins, {total_comparisons:,} comparisons")
    
    # Cluster size analysis
    print("\nCLUSTER SIZE DISTRIBUTION")
    print("-" * 40)
    for cat in sorted(df['category'].unique()):
        cat_df = df[df['category'] == cat]
        sizes = cat_df['n'].values
        print(f"{cat}:")
        print(f"  Mean cluster size: {np.mean(sizes):.1f} ± {np.std(sizes):.1f}")
        print(f"  Median cluster size: {np.median(sizes):.1f}")
        print(f"  Size range: {np.min(sizes)} - {np.max(sizes)} proteins")
        print(f"  Quartiles (Q1, Q3): {np.percentile(sizes, 25):.1f}, {np.percentile(sizes, 75):.1f}")
    
    # Main statistical results
    print("\nGROUP STATISTICS (SIZE-CORRECTED)")
    print("-" * 40)
    
    for cat in ['RANDOM', 'MGC (KEGG)', 'BGC (MIBiG)']:
        if cat in group_stats:
            stats = group_stats[cat]
            print(f"\n{cat.upper()}:")
            print(f"  Number of clusters: {stats['n_clusters']}")
            print(f"  Total proteins: {stats['total_proteins']}")
            print(f"  Total pairwise comparisons: {stats['total_comparisons']:,}")
            print(f"  Mean cluster size: {stats['mean_cluster_size']:.1f} ± {stats['std_cluster_size']:.1f}")
            print(f"  Weighted mean docking score: {stats['weighted_mean_score']:.3f}")
            print(f"  95% Bootstrap CI: [{stats['ci_95_lower']:.3f}, {stats['ci_95_upper']:.3f}]")
            print(f"  Standard deviation: {stats['weighted_std']:.3f}")
            print(f"  Best score overall: {stats['best_score_overall']:.2f}")
            print(f"  Mean enrichment vs random: {stats['mean_enrichment']:.3f}")
            print(f"  Median enrichment vs random: {stats['median_enrichment']:.3f}")
            
            print(f"  Binding strength distribution:")
            print(f"    Strong binders (< -50): {stats['fraction_strong_binders']:.1%}")
            print(f"    Moderate binders (-50 to -25): {stats['fraction_moderate_binders']:.1%}")
            print(f"    Weak binders (≥ -25): {stats['fraction_weak_binders']:.1%}")
    
    # Statistical significance tests
    print("\nSTATISTICAL SIGNIFICANCE TESTS")
    print("-" * 40)
    
    if 'ANOVA_enrichment' in test_results:
        anova = test_results['ANOVA_enrichment']
        print(f"\nOne-way ANOVA on enrichment scores:")
        print(f"  F-statistic: {anova['F']:.3f}")
        print(f"  p-value: {anova['p']:.2e}")
        print(f"  Significance: {'***' if anova['p'] < 0.001 else '**' if anova['p'] < 0.01 else '*' if anova['p'] < 0.05 else 'ns'}")
    
    print(f"\nPairwise comparisons (Permutation tests with {10000:,} iterations):")
    print("  Comparison                    | Mean Diff | Enrichment Diff | p-value | Significance")
    print("  " + "-" * 75)
    
    for comp in test_results['pairwise']:
        significance = "***" if comp['p_permutation'] < 0.001 else "**" if comp['p_permutation'] < 0.01 else "*" if comp['p_permutation'] < 0.05 else "ns"
        print(f"  {comp['comparison']:<28} | {comp['mean_diff']:8.3f} | {comp['enrichment_diff']:14.3f} | {comp['p_permutation']:7.3f} | {significance:>12}")
    
    # Effect sizes
    print("\nEFFECT SIZES (Cohen's d relative to RANDOM)")
    print("-" * 40)
    
    for cat in ['MGC (KEGG)', 'BGC (MIBiG)']:
        if cat in df['category'].unique():
            cat_df = df[df['category'] == cat]
            effect_sizes = cat_df['effect_size'].values
            print(f"{cat}:")
            print(f"  Mean effect size: {np.mean(effect_sizes):.3f}")
            print(f"  Median effect size: {np.median(effect_sizes):.3f}")
            print(f"  Effect size range: {np.min(effect_sizes):.3f} to {np.max(effect_sizes):.3f}")
            
            # Interpret effect sizes
            large_effects = np.sum(np.abs(effect_sizes) > 0.8)
            medium_effects = np.sum((np.abs(effect_sizes) > 0.5) & (np.abs(effect_sizes) <= 0.8))
            small_effects = np.sum((np.abs(effect_sizes) > 0.2) & (np.abs(effect_sizes) <= 0.5))
            
            print(f"  Large effects (|d| > 0.8): {large_effects}/{len(effect_sizes)} ({large_effects/len(effect_sizes):.1%})")
            print(f"  Medium effects (0.5 < |d| ≤ 0.8): {medium_effects}/{len(effect_sizes)} ({medium_effects/len(effect_sizes):.1%})")
            print(f"  Small effects (0.2 < |d| ≤ 0.5): {small_effects}/{len(effect_sizes)} ({small_effects/len(effect_sizes):.1%})")
    
    # Detailed score distributions
    print("\nSCORE DISTRIBUTION ANALYSIS")
    print("-" * 40)
    
    for cat in sorted(df['category'].unique()):
        cat_df = df[df['category'] == cat]
        scores = []
        for _, row in cat_df.iterrows():
            # We don't have individual scores, so use cluster means as proxy
            scores.extend([row['mean_score_non_self']] * int(row['n_comparisons_non_self']))
        
        if scores:
            scores = np.array(scores)
            print(f"\n{cat} - Score distribution:")
            print(f"  Mean: {np.mean(scores):.3f}")
            print(f"  Median: {np.median(scores):.3f}")
            print(f"  Standard deviation: {np.std(scores):.3f}")
            print(f"  Range: {np.min(scores):.2f} to {np.max(scores):.2f}")
            print(f"  Quartiles (Q1, Q3): {np.percentile(scores, 25):.3f}, {np.percentile(scores, 75):.3f}")
            print(f"  Fraction < -50 (strong): {np.mean(scores < -50):.1%}")
            print(f"  Fraction < -25 (moderate+): {np.mean(scores < -25):.1%}")
    
    # Summary and interpretation
    print("\nSUMMARY AND INTERPRETATION")
    print("-" * 40)
    
    # Find best performing category
    if len(group_stats) > 1:
        best_cat = min(group_stats.keys(), key=lambda k: group_stats[k]['weighted_mean_score'])
        worst_cat = max(group_stats.keys(), key=lambda k: group_stats[k]['weighted_mean_score'])
        
        best_score = group_stats[best_cat]['weighted_mean_score']
        worst_score = group_stats[worst_cat]['weighted_mean_score']
        
        print(f"Best performing category: {best_cat} (mean score: {best_score:.3f})")
        print(f"Worst performing category: {worst_cat} (mean score: {worst_score:.3f})")
        print(f"Difference: {worst_score - best_score:.3f} units")
        
        # Statistical significance of best vs worst
        best_vs_worst = [c for c in test_results['pairwise'] 
                        if (best_cat in c['comparison'] and worst_cat in c['comparison'])]
        if best_vs_worst:
            comp = best_vs_worst[0]
            print(f"Statistical significance: p = {comp['p_permutation']:.3f}")
    
    print(f"\nNote: More negative docking scores indicate stronger predicted binding affinity.")
    print(f"Enrichment scores show deviation from random expectation (negative = worse than random).")
    
    print("\n" + "=" * 80)

def advanced_cross_validation_analysis(scores, binary_labels, n_splits=5, random_state=42):
    """
    Perform advanced k-fold cross-validation with multiple algorithms and class balancing.
    """
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    
    # Create feature matrix with additional engineered features
    X = np.column_stack([
        scores,
        scores**2,  # Quadratic term
        np.log(np.abs(scores) + 1),  # Log transformation
        np.sqrt(np.abs(scores))  # Square root transformation
    ])
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Define multiple algorithms with hyperparameter tuning
    algorithms = {
        'LogisticRegression': Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', LogisticRegression(random_state=random_state, max_iter=1000))
        ]),
        'RandomForest': Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', RandomForestClassifier(
                n_estimators=100, 
                max_depth=10, 
                min_samples_split=5,
                min_samples_leaf=2,
                random_state=random_state
            ))
        ]),
        'GradientBoosting': Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', GradientBoostingClassifier(
                n_estimators=100,
                max_depth=6,
                learning_rate=0.1,
                random_state=random_state
            ))
        ]),
        'SVM': Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', SVC(
                kernel='rbf',
                C=1.0,
                gamma='scale',
                probability=True,
                random_state=random_state
            ))
        ])
    }
    
    # Class balancing pipeline (only if imblearn is available)
    all_algorithms = algorithms.copy()
    
    if IMBLEARN_AVAILABLE:
        balanced_algorithms = {}
        for name, pipeline in algorithms.items():
            balanced_algorithms[f'{name}_Balanced'] = ImbPipeline([
                ('sampler', SMOTE(random_state=random_state)),
                ('classifier', pipeline)
            ])
        all_algorithms.update(balanced_algorithms)
    else:
        print("Skipping balanced algorithms - imblearn not available")
    
    results = {}
    
    for name, model in all_algorithms.items():
        try:
            # Cross-validation scores
            cv_accuracy = cross_val_score(model, X_scaled, binary_labels, 
                                        cv=skf, scoring='accuracy')
            cv_precision = cross_val_score(model, X_scaled, binary_labels, 
                                         cv=skf, scoring='precision')
            cv_recall = cross_val_score(model, X_scaled, binary_labels, 
                                      cv=skf, scoring='recall')
            cv_f1 = cross_val_score(model, X_scaled, binary_labels, 
                                   cv=skf, scoring='f1')
            cv_roc_auc = cross_val_score(model, X_scaled, binary_labels, 
                                       cv=skf, scoring='roc_auc')
            
            results[name] = {
                'accuracy': {'mean': cv_accuracy.mean(), 'std': cv_accuracy.std()},
                'precision': {'mean': cv_precision.mean(), 'std': cv_precision.std()},
                'recall': {'mean': cv_recall.mean(), 'std': cv_recall.std()},
                'f1': {'mean': cv_f1.mean(), 'std': cv_f1.std()},
                'roc_auc': {'mean': cv_roc_auc.mean(), 'std': cv_roc_auc.std()}
            }
            
        except Exception as e:
            print(f"Error with {name}: {e}")
            continue
    
    # Print results
    print(f"\n{n_splits}-Fold Cross-Validation Results (Multiple Algorithms):")
    print(f"{'Algorithm':<20} {'Accuracy':<12} {'Precision':<12} {'Recall':<12} {'F1':<12} {'ROC AUC':<12}")
    print("-" * 80)
    
    for name, metrics in results.items():
        print(f"{name:<20} {metrics['accuracy']['mean']:.4f}±{metrics['accuracy']['std']:.3f} "
              f"{metrics['precision']['mean']:.4f}±{metrics['precision']['std']:.3f} "
              f"{metrics['recall']['mean']:.4f}±{metrics['recall']['std']:.3f} "
              f"{metrics['f1']['mean']:.4f}±{metrics['f1']['std']:.3f} "
              f"{metrics['roc_auc']['mean']:.4f}±{metrics['roc_auc']['std']:.3f}")
    
    # Find best algorithm
    best_algorithm = max(results.keys(), 
                        key=lambda x: results[x]['roc_auc']['mean'])
    
    print(f"\nBest Algorithm: {best_algorithm}")
    print(f"  ROC AUC: {results[best_algorithm]['roc_auc']['mean']:.4f} ± {results[best_algorithm]['roc_auc']['std']:.4f}")
    print(f"  F1 Score: {results[best_algorithm]['f1']['mean']:.4f} ± {results[best_algorithm]['f1']['std']:.4f}")
    
    return results, best_algorithm

def advanced_threshold_optimization_with_validation(df, score_column='mean_score_non_self', 
                                                test_size=0.3, n_splits=5, random_state=42):
    """
    Perform advanced threshold optimization with ensemble methods, hyperparameter tuning,
    and class balancing to prevent overfitting/underfitting.
    """
    print(f"\n{'='*70}")
    print("ADVANCED THRESHOLD OPTIMIZATION WITH VALIDATION")
    print("="*70)
    
    # Create binary labels: 1 for positives (MIBiG + MGC), 0 for negatives (RANDOM)
    binary_labels = df['category'].isin(['BGC (MIBiG)', 'MGC (KEGG)']).astype(int)
    scores = df[score_column].values
    
    # Remove any NaN values
    valid_mask = ~np.isnan(scores)
    scores = scores[valid_mask]
    binary_labels = binary_labels[valid_mask]
    
    if len(np.unique(binary_labels)) < 2:
        print(f"Warning: Not enough classes for {score_column}")
        return None
    
    print(f"\nBinary Classification Setup:")
    print(f"  Positive class: BGC (MIBiG) + MGC (KEGG)")
    print(f"  Negative class: RANDOM")
    print(f"  Total positive samples: {np.sum(binary_labels)}")
    print(f"  Total negative samples: {len(binary_labels) - np.sum(binary_labels)}")
    print(f"  Class balance ratio: {np.sum(binary_labels) / len(binary_labels):.3f}")
    
    # Create enhanced feature matrix
    X = np.column_stack([
        scores,
        scores**2,  # Quadratic term
        np.log(np.abs(scores) + 1),  # Log transformation
        np.sqrt(np.abs(scores))  # Square root transformation
    ])
    
    # Train/Test Split with stratification
    X_train, X_test, y_train, y_test = train_test_split(
        X, binary_labels, 
        test_size=test_size, stratify=binary_labels, 
        random_state=random_state
    )
    
    print(f"\nDataset Split:")
    print(f"  Training set: {len(X_train)} samples ({np.sum(y_train)} positive, {len(y_train)-np.sum(y_train)} negative)")
    print(f"  Test set: {len(X_test)} samples ({np.sum(y_test)} positive, {len(y_test)-np.sum(y_test)} negative)")
    
    # 1. Advanced Cross-Validation with Multiple Algorithms
    print(f"\n1. ADVANCED CROSS-VALIDATION ANALYSIS")
    print("-" * 50)
    
    cv_results, best_algorithm = advanced_cross_validation_analysis(scores, binary_labels, n_splits, random_state)
    
    # 2. Hyperparameter Tuning for Best Algorithm
    print(f"\n2. HYPERPARAMETER TUNING")
    print("-" * 50)
    
    # Define parameter grids for different algorithms
    param_grids = {
        'LogisticRegression': {
            'classifier__C': [0.1, 1.0, 10.0, 100.0],
            'classifier__penalty': ['l1', 'l2'],
            'classifier__solver': ['liblinear', 'saga']
        },
        'RandomForest': {
            'classifier__n_estimators': [50, 100, 200],
            'classifier__max_depth': [5, 10, 15, None],
            'classifier__min_samples_split': [2, 5, 10],
            'classifier__min_samples_leaf': [1, 2, 4]
        },
        'GradientBoosting': {
            'classifier__n_estimators': [50, 100, 200],
            'classifier__max_depth': [3, 6, 10],
            'classifier__learning_rate': [0.01, 0.1, 0.2]
        },
        'SVM': {
            'classifier__C': [0.1, 1.0, 10.0, 100.0],
            'classifier__gamma': ['scale', 'auto', 0.001, 0.01, 0.1, 1.0],
            'classifier__kernel': ['rbf', 'poly', 'sigmoid']
        }
    }
    
    # Get base algorithm name (remove _Balanced suffix)
    base_algorithm = best_algorithm.replace('_Balanced', '')
    
    if base_algorithm in param_grids:
        print(f"Tuning hyperparameters for {base_algorithm}...")
        
        # Create the base pipeline
        if 'Balanced' in best_algorithm:
            base_pipeline = ImbPipeline([
                ('sampler', SMOTE(random_state=random_state)),
                ('scaler', StandardScaler()),
                ('classifier', None)  # Will be set by GridSearchCV
            ])
        else:
            base_pipeline = Pipeline([
                ('scaler', StandardScaler()),
                ('classifier', None)  # Will be set by GridSearchCV
            ])
        
        # Set the classifier
        if base_algorithm == 'LogisticRegression':
            base_pipeline.set_params(classifier=LogisticRegression(random_state=random_state, max_iter=1000))
        elif base_algorithm == 'RandomForest':
            base_pipeline.set_params(classifier=RandomForestClassifier(random_state=random_state))
        elif base_algorithm == 'GradientBoosting':
            base_pipeline.set_params(classifier=GradientBoostingClassifier(random_state=random_state))
        elif base_algorithm == 'SVM':
            base_pipeline.set_params(classifier=SVC(probability=True, random_state=random_state))
        
        # Grid search with cross-validation
        grid_search = GridSearchCV(
            base_pipeline, 
            param_grids[base_algorithm], 
            cv=StratifiedKFold(n_splits=3, shuffle=True, random_state=random_state),
            scoring='roc_auc',
            n_jobs=-1,
            verbose=0
        )
        
        grid_search.fit(X_train, y_train)
        
        print(f"Best parameters: {grid_search.best_params_}")
        print(f"Best CV score: {grid_search.best_score_:.4f}")
        
        best_model = grid_search.best_estimator_
    else:
        print(f"Using default parameters for {best_algorithm}")
        # Use the best algorithm from CV without tuning
        best_model = None  # Will be handled in the next section
    
    # 3. Ensemble Method
    print(f"\n3. ENSEMBLE METHOD")
    print("-" * 50)
    
    # Create ensemble of best performing algorithms
    ensemble_models = []
    for name, metrics in cv_results.items():
        if metrics['roc_auc']['mean'] > 0.6:  # Only include decent performers
            if 'Balanced' in name and IMBLEARN_AVAILABLE:
                base_name = name.replace('_Balanced', '')
                if base_name == 'LogisticRegression':
                    model = ImbPipeline([
                        ('sampler', SMOTE(random_state=random_state)),
                        ('scaler', StandardScaler()),
                        ('classifier', LogisticRegression(random_state=random_state, max_iter=1000))
                    ])
                elif base_name == 'RandomForest':
                    model = ImbPipeline([
                        ('sampler', SMOTE(random_state=random_state)),
                        ('scaler', StandardScaler()),
                        ('classifier', RandomForestClassifier(random_state=random_state))
                    ])
                elif base_name == 'GradientBoosting':
                    model = ImbPipeline([
                        ('sampler', SMOTE(random_state=random_state)),
                        ('scaler', StandardScaler()),
                        ('classifier', GradientBoostingClassifier(random_state=random_state))
                    ])
                elif base_name == 'SVM':
                    model = ImbPipeline([
                        ('sampler', SMOTE(random_state=random_state)),
                        ('scaler', StandardScaler()),
                        ('classifier', SVC(probability=True, random_state=random_state))
                    ])
            else:
                if name == 'LogisticRegression':
                    model = Pipeline([
                        ('scaler', StandardScaler()),
                        ('classifier', LogisticRegression(random_state=random_state, max_iter=1000))
                    ])
                elif name == 'RandomForest':
                    model = Pipeline([
                        ('scaler', StandardScaler()),
                        ('classifier', RandomForestClassifier(random_state=random_state))
                    ])
                elif name == 'GradientBoosting':
                    model = Pipeline([
                        ('scaler', StandardScaler()),
                        ('classifier', GradientBoostingClassifier(random_state=random_state))
                    ])
                elif name == 'SVM':
                    model = Pipeline([
                        ('scaler', StandardScaler()),
                        ('classifier', SVC(probability=True, random_state=random_state))
                    ])
            
            ensemble_models.append((name, model))
    
    if len(ensemble_models) > 1:
        # Create voting classifier
        voting_classifier = VotingClassifier(
            estimators=ensemble_models,
            voting='soft'  # Use predicted probabilities
        )
        
        print(f"Created ensemble with {len(ensemble_models)} models")
        
        # Train ensemble
        voting_classifier.fit(X_train, y_train)
        
        # Evaluate ensemble
        ensemble_pred = voting_classifier.predict(X_test)
        ensemble_pred_proba = voting_classifier.predict_proba(X_test)[:, 1]
        
        # Calculate ensemble metrics
        tn = np.sum((ensemble_pred == 0) & (y_test == 0))
        fp = np.sum((ensemble_pred == 1) & (y_test == 0))
        fn = np.sum((ensemble_pred == 0) & (y_test == 1))
        tp = np.sum((ensemble_pred == 1) & (y_test == 1))
        
        ensemble_accuracy = (tp + tn) / len(y_test)
        ensemble_precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        ensemble_recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        ensemble_f1 = 2 * (ensemble_precision * ensemble_recall) / (ensemble_precision + ensemble_recall) if (ensemble_precision + ensemble_recall) > 0 else 0
        
        # ROC AUC for ensemble
        fpr_ensemble, tpr_ensemble, _ = roc_curve(y_test, ensemble_pred_proba)
        ensemble_roc_auc = auc(fpr_ensemble, tpr_ensemble)
        
        print(f"\nEnsemble Performance:")
        print(f"  Accuracy: {ensemble_accuracy:.4f}")
        print(f"  Precision: {ensemble_precision:.4f}")
        print(f"  Recall: {ensemble_recall:.4f}")
        print(f"  F1 Score: {ensemble_f1:.4f}")
        print(f"  ROC AUC: {ensemble_roc_auc:.4f}")
        
        # Use ensemble as final model
        final_model = voting_classifier
        final_pred_proba = ensemble_pred_proba
    else:
        # Use single best model
        if best_model is not None:
            final_model = best_model
        else:
            # Fallback to simple logistic regression
            final_model = Pipeline([
                ('scaler', StandardScaler()),
                ('classifier', LogisticRegression(random_state=random_state, max_iter=1000))
            ])
            final_model.fit(X_train, y_train)
        
        final_pred_proba = final_model.predict_proba(X_test)[:, 1]
    
    # 4. Threshold Optimization on Final Model
    print(f"\n4. THRESHOLD OPTIMIZATION")
    print("-" * 50)
    
    # Get training probabilities for threshold optimization
    train_pred_proba = final_model.predict_proba(X_train)[:, 1]
    
    # Calculate ROC curve on training data
    fpr_train, tpr_train, thresholds_train = roc_curve(y_train, train_pred_proba)
    roc_auc_train = auc(fpr_train, tpr_train)
    
    # Calculate precision-recall curve on training data
    precision_train, recall_train, pr_thresholds_train = precision_recall_curve(y_train, train_pred_proba)
    pr_auc_train = auc(recall_train, precision_train)
    
    # Find optimal thresholds
    # Youden's J
    youden_j = tpr_train - fpr_train
    optimal_idx_youden = np.argmax(youden_j)
    optimal_threshold_youden = thresholds_train[optimal_idx_youden]
    
    # F1 score
    f1_scores = 2 * (precision_train * recall_train) / (precision_train + recall_train + 1e-8)
    optimal_idx_f1 = np.argmax(f1_scores)
    optimal_threshold_f1 = pr_thresholds_train[optimal_idx_f1]
    
    # Balanced accuracy
    tnr_train = 1 - fpr_train
    balanced_acc = (tpr_train + tnr_train) / 2
    optimal_idx_balanced = np.argmax(balanced_acc)
    optimal_threshold_balanced = thresholds_train[optimal_idx_balanced]
    
    print(f"\nTraining Performance:")
    print(f"  ROC AUC: {roc_auc_train:.4f}")
    print(f"  PR AUC: {pr_auc_train:.4f}")
    print(f"\nOptimal thresholds:")
    print(f"  Youden's J: {optimal_threshold_youden:.4f}")
    print(f"  F1-optimized: {optimal_threshold_f1:.4f}")
    print(f"  Accuracy-optimized: {optimal_threshold_balanced:.4f}")
    
    # 5. Test Set Evaluation
    print(f"\n5. TEST SET EVALUATION")
    print("-" * 50)
    
    # Test set ROC analysis
    fpr_test, tpr_test, _ = roc_curve(y_test, final_pred_proba)
    roc_auc_test = auc(fpr_test, tpr_test)
    
    precision_test, recall_test, _ = precision_recall_curve(y_test, final_pred_proba)
    pr_auc_test = auc(recall_test, precision_test)
    
    print(f"Test Set Performance:")
    print(f"  ROC AUC: {roc_auc_test:.4f}")
    print(f"  PR AUC: {pr_auc_test:.4f}")
    
    # Evaluate different thresholds on test set
    test_results = {}
    thresholds_to_test = {
        'youden': optimal_threshold_youden,
        'f1': optimal_threshold_f1,
        'accuracy': optimal_threshold_balanced
    }
    
    for metric_name, threshold in thresholds_to_test.items():
        y_pred_test = (final_pred_proba >= threshold).astype(int)
        
        # Calculate metrics
        tn = np.sum((y_pred_test == 0) & (y_test == 0))
        fp = np.sum((y_pred_test == 1) & (y_test == 0))
        fn = np.sum((y_pred_test == 0) & (y_test == 1))
        tp = np.sum((y_pred_test == 1) & (y_test == 1))
        
        test_accuracy = (tp + tn) / len(y_test)
        test_precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        test_recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        test_specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        test_f1 = 2 * (test_precision * test_recall) / (test_precision + test_recall) if (test_precision + test_recall) > 0 else 0
        
        test_results[metric_name] = {
            'threshold': threshold,
            'accuracy': test_accuracy,
            'precision': test_precision,
            'recall': test_recall,
            'specificity': test_specificity,
            'f1_score': test_f1
        }
        
        print(f"\n{metric_name.upper()} threshold:")
        print(f"  Threshold: {threshold:.4f}")
        print(f"  Accuracy: {test_accuracy:.4f}")
        print(f"  Precision: {test_precision:.4f}")
        print(f"  Recall: {test_recall:.4f}")
        print(f"  Specificity: {test_specificity:.4f}")
        print(f"  F1 Score: {test_f1:.4f}")
    
    # 6. Overfitting/Underfitting Assessment
    print(f"\n6. OVERFITTING/UNDERFITTING ASSESSMENT")
    print("-" * 50)
    
    # Calculate training performance
    train_pred = final_model.predict(X_train)
    train_pred_proba = final_model.predict_proba(X_train)[:, 1]
    
    # Training metrics
    train_accuracy = (train_pred == y_train).mean()
    train_precision = precision_score(y_train, train_pred)
    train_recall = recall_score(y_train, train_pred)
    train_f1 = f1_score(y_train, train_pred)
    
    # Test metrics (using F1 threshold)
    test_pred = (final_pred_proba >= optimal_threshold_f1).astype(int)
    test_accuracy = (test_pred == y_test).mean()
    test_precision = precision_score(y_test, test_pred)
    test_recall = recall_score(y_test, test_pred)
    test_f1 = f1_score(y_test, test_pred)
    
    print(f"Training vs Test Performance:")
    print(f"  Accuracy:  {train_accuracy:.4f} vs {test_accuracy:.4f} (diff: {abs(train_accuracy - test_accuracy):.4f})")
    print(f"  Precision: {train_precision:.4f} vs {test_precision:.4f} (diff: {abs(train_precision - test_precision):.4f})")
    print(f"  Recall:    {train_recall:.4f} vs {test_recall:.4f} (diff: {abs(train_recall - test_recall):.4f})")
    print(f"  F1 Score:  {train_f1:.4f} vs {test_f1:.4f} (diff: {abs(train_f1 - test_f1):.4f})")
    
    # Assess overfitting/underfitting
    accuracy_diff = abs(train_accuracy - test_accuracy)
    f1_diff = abs(train_f1 - test_f1)
    
    if accuracy_diff < 0.05 and f1_diff < 0.05:
        fit_assessment = "Well-fitted model - no significant overfitting or underfitting"
    elif train_accuracy > test_accuracy + 0.1 or train_f1 > test_f1 + 0.1:
        fit_assessment = "Potential overfitting - training performance much better than test"
    elif test_accuracy < 0.6 and test_f1 < 0.6:
        fit_assessment = "Potential underfitting - both training and test performance are low"
    else:
        fit_assessment = "Moderate fit - some performance difference between train and test"
    
    print(f"\nFit Assessment: {fit_assessment}")
    
    # 7. Final Recommendations
    print(f"\n7. FINAL RECOMMENDATIONS")
    print("-" * 50)
    
    # Determine recommendation based on performance and fit
    if roc_auc_test >= 0.8 and accuracy_diff < 0.05:
        recommendation = "Excellent model - highly recommended for production"
    elif roc_auc_test >= 0.7 and accuracy_diff < 0.1:
        recommendation = "Good model - suitable for production with monitoring"
    elif roc_auc_test >= 0.6:
        recommendation = "Moderate model - consider additional features or more data"
    else:
        recommendation = "Poor model - requires significant improvement"
    
    print(f"Overall Recommendation: {recommendation}")
    
    if 'f1' in test_results:
        best_threshold = test_results['f1']['threshold']
        best_f1 = test_results['f1']['f1_score']
        print(f"\nRecommended Production Threshold: {best_threshold:.4f}")
        print(f"  Expected F1 Score: {best_f1:.4f}")
        print(f"  Expected Precision: {test_results['f1']['precision']:.4f}")
        print(f"  Expected Recall: {test_results['f1']['recall']:.4f}")
    
    return {
        'cv_results': cv_results,
        'best_algorithm': best_algorithm,
        'test_results': test_results,
        'test_roc_auc': roc_auc_test,
        'test_pr_auc': pr_auc_test,
        'fit_assessment': fit_assessment,
        'recommendation': recommendation,
        'final_model': final_model
    }

def evaluate_all_scores_with_validation(df):
    """
    Evaluate advanced threshold optimization for all available score columns.
    """
    score_columns = ['mean_score_non_self', 'enrichment_score', 'z_score', 'effect_size']
    all_results = {}
    
    print("\n" + "="*80)
    print("ADVANCED THRESHOLD OPTIMIZATION WITH VALIDATION")
    print("="*80)
    
    for score_col in score_columns:
        if score_col not in df.columns:
            print(f"Skipping {score_col} - not available in data")
            continue
            
        print(f"\n{'='*60}")
        print(f"ANALYZING SCORE: {score_col.upper()}")
        print("="*60)
        
        try:
            result = advanced_threshold_optimization_with_validation(df, score_col)
            if result:
                all_results[score_col] = result
            else:
                print(f"Failed to analyze {score_col}")
        except Exception as e:
            print(f"Error analyzing {score_col}: {e}")
    
    # Summary comparison
    print(f"\n{'='*80}")
    print("SCORE COMPARISON SUMMARY")
    print("="*80)
    
    if all_results:
        print(f"{'Score':<20} {'Test ROC AUC':<12} {'Test F1':<10} {'Fit Assessment':<30} {'Recommendation':<25}")
        print("-" * 100)
        
        for score_col, result in all_results.items():
            if 'test_results' in result and 'f1' in result['test_results']:
                roc_auc = result['test_roc_auc']
                f1 = result['test_results']['f1']['f1_score']
                fit = result['fit_assessment'][:30] + "..." if len(result['fit_assessment']) > 30 else result['fit_assessment']
                rec = result['recommendation'][:25] + "..." if len(result['recommendation']) > 25 else result['recommendation']
                print(f"{score_col:<20} {roc_auc:<12.4f} {f1:<10.4f} {fit:<30} {rec:<25}")
    
    return all_results

def main():
    base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
    
    # Collect all cluster directories with LightDock results
    lightdock_results_dir = base_dir / "lightdock_results"
    
    if not lightdock_results_dir.exists():
        print(f"LightDock results directory not found: {lightdock_results_dir}")
        return
    
    cluster_dirs = []
    for cluster_dir in lightdock_results_dir.iterdir():
        if cluster_dir.is_dir():
            # Determine category based on name
            if cluster_dir.name.startswith('RANDOM_'):
                category = 'RANDOM'
            elif cluster_dir.name.startswith('BGC'):
                category = 'BGC (MIBiG)'
            elif cluster_dir.name.startswith('MGC_CANDIDATE_'):
                category = 'MGC (KEGG)'
            else:
                continue
            
            cluster_dirs.append((cluster_dir, category))
    
    print(f"Found {len(cluster_dirs)} cluster directories to process...")
    
    # Use concurrent processing to collect results
    print("Collecting LightDock results using concurrent processing...")
    with Pool(cpu_count()) as pool:
        all_cluster_data = pool.map(collect_lightdock_results_from_cluster, cluster_dirs)
    
    # Filter out None results
    all_cluster_data = [data for data in all_cluster_data if data is not None]
    print(f"Successfully collected results from {len(all_cluster_data)} clusters")
    
    # Use concurrent processing to analyze results
    print("Processing cluster statistics using concurrent processing...")
    with Pool(cpu_count()) as pool:
        processed_results = pool.map(process_cluster_results, all_cluster_data)
    
    # Filter out None results
    processed_results = [result for result in processed_results if result is not None]
    
    if not processed_results:
        print("No valid processed results found!")
        return
    
    # Create DataFrame
    df = pd.DataFrame(processed_results)
    
    print(f"Final dataset: {len(df)} clusters")
    for cat in sorted(df['category'].unique()):
        count = len(df[df['category'] == cat])
        print(f"  {cat}: {count} clusters")
    
    # Apply size normalization
    print("Applying statistical normalizations...")
    df = calculate_size_normalized_scores(df)
    
    # Calculate group statistics with proper weighting
    print("Calculating group statistics...")
    group_stats = weighted_group_statistics(df)
    
    # Statistical tests
    print("Performing statistical tests...")
    test_results = statistical_tests(df, group_stats)
    
    # Print comprehensive results
    print_detailed_results(df, group_stats, test_results)
    
    # ML Threshold Analysis
    print("\n" + "="*80)
    print("MACHINE LEARNING THRESHOLD ANALYSIS")
    print("="*80)
    
    # Perform ML threshold optimization with validation
    ml_results = evaluate_all_scores_with_validation(df)
    
    # Save results to files
    print("\nSaving results to files...")
    
    # Export data
    df.to_csv(base_dir / 'lightdock_cluster_metrics.csv', index=False)
    print(f"Saved: lightdock_cluster_metrics.csv ({len(df)} clusters)")
    
    group_stats_df = pd.DataFrame(group_stats).T
    group_stats_df.to_csv(base_dir / 'lightdock_group_statistics.csv')
    print(f"Saved: lightdock_group_statistics.csv ({len(group_stats)} categories)")
    
    # Export detailed comparison results
    if test_results.get('pairwise'):
        comparisons_df = pd.DataFrame(test_results['pairwise'])
        comparisons_df.to_csv(base_dir / 'lightdock_statistical_comparisons.csv', index=False)
        print(f"Saved: lightdock_statistical_comparisons.csv ({len(comparisons_df)} comparisons)")
    
    # Export ML results
    if ml_results:
        ml_summary = {}
        for score_col, result in ml_results.items():
            if 'test_results' in result and 'f1' in result['test_results']:
                ml_summary[score_col] = {
                    'threshold': result['test_results']['f1']['threshold'],
                    'test_roc_auc': result['test_roc_auc'],
                    'test_f1': result['test_results']['f1']['f1_score'],
                    'test_precision': result['test_results']['f1']['precision'],
                    'test_recall': result['test_results']['f1']['recall'],
                    'recommendation': result['recommendation']
                }
        
        if ml_summary:
            ml_df = pd.DataFrame(ml_summary).T
            ml_df.to_csv(base_dir / 'lightdock_ml_thresholds.csv')
            print(f"Saved: lightdock_ml_thresholds.csv ({len(ml_df)} score types)")
    
    print(f"\nAll files saved to: {base_dir}")
    
    return df, group_stats, test_results, ml_results

if __name__ == "__main__":
    df, group_stats, test_results, ml_results = main()