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
    
    print(f"\nAll files saved to: {base_dir}")
    
    return df, group_stats, test_results

if __name__ == "__main__":
    df, group_stats, test_results = main()