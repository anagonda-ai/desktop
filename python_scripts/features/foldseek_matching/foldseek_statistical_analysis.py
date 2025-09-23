#!/usr/bin/env python3
"""
Statistically rigorous normalization of Foldseek all-vs-all results.
Accounts for n×n matrix structure and cluster size effects.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.optimize import curve_fit
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

def process_complete_tsv(file_path):
    """Process entire TSV file - use ALL data, no sampling"""
    try:
        name = file_path.stem.replace('_all_vs_all', '')
        
        # Determine category
        if name.startswith('RANDOM_'):
            category = 'RANDOM'
        elif name.startswith('BGC'):
            category = 'BGC (MIBiG)'
        elif name.startswith('MGC_CANDIDATE_'):
            category = 'MGC (KEGG)'
        else:
            return None
        
        # Read complete TSV
        df = pd.read_csv(file_path, sep='\t', header=None)
        
        if df.empty:
            return None
            
        # Assign columns based on count
        if len(df.columns) >= 13:
            col_names = ['query', 'target', 'fident', 'alnlen', 'mismatch', 
                        'gapopen', 'qstart', 'qend', 'tstart', 'tend', 
                        'evalue', 'bits', 'tmscore']
            df.columns = col_names + [f'col_{i}' for i in range(len(df.columns) - len(col_names))]
        else:
            return None
        
        # Get cluster size
        n = len(df['query'].unique())
        
        # ALL pairwise comparisons (including self)
        all_tm_scores = df['tmscore'].values
        all_fident = df['fident'].values
        
        # Non-self comparisons
        non_self_mask = df['query'] != df['target']
        non_self_tm = df.loc[non_self_mask, 'tmscore'].values
        non_self_fident = df.loc[non_self_mask, 'fident'].values
        
        # Expected number of comparisons
        expected_total = n * n  # Including self
        expected_non_self = n * (n - 1)  # Excluding self
        
        # Calculate raw sums and means
        sum_tm_all = np.sum(all_tm_scores)
        sum_tm_non_self = np.sum(non_self_tm)
        
        return {
            'category': category,
            'name': name,
            'n': n,
            'n_squared': n * n,
            'n_comparisons_total': len(df),
            'n_comparisons_non_self': non_self_mask.sum(),
            
            # Raw sums (for accurate group averaging)
            'sum_tm_all': sum_tm_all,
            'sum_tm_non_self': sum_tm_non_self,
            'sum_fident_all': np.sum(all_fident),
            'sum_fident_non_self': np.sum(non_self_fident),
            
            # Means
            'mean_tm_all': sum_tm_all / expected_total if expected_total > 0 else 0,
            'mean_tm_non_self': sum_tm_non_self / expected_non_self if expected_non_self > 0 else 0,
            'mean_fident_non_self': np.mean(non_self_fident) if len(non_self_fident) > 0 else 0,
            
            # Variance components
            'var_tm_non_self': np.var(non_self_tm) if len(non_self_tm) > 0 else 0,
            'std_tm_non_self': np.std(non_self_tm) if len(non_self_tm) > 0 else 0,
            
            # Distribution characteristics
            'median_tm_non_self': np.median(non_self_tm) if len(non_self_tm) > 0 else 0,
            'q25_tm': np.percentile(non_self_tm, 25) if len(non_self_tm) > 0 else 0,
            'q75_tm': np.percentile(non_self_tm, 75) if len(non_self_tm) > 0 else 0,
            'max_tm': np.max(non_self_tm) if len(non_self_tm) > 0 else 0,
            
            # For size effect modeling
            'log_n': np.log(n),
            'sqrt_n': np.sqrt(n)
        }
        
    except Exception as e:
        print(f"Error in {file_path.name}: {e}")
        return None

def calculate_size_normalized_scores(df):
    """Apply statistical normalization for cluster size effects"""
    
    # 1. Model size effect on random clusters to establish null distribution
    random_df = df[df['category'] == 'RANDOM'].copy()
    
    if len(random_df) > 10:  # Need enough data for modeling
        # Fit regression model: score ~ f(n)
        X = random_df[['n', 'n_squared']].values
        y = random_df['mean_tm_non_self'].values
        
        # Robust regression to handle outliers
        from sklearn.linear_model import HuberRegressor
        model = HuberRegressor()
        model.fit(X, y)
        
        # Predict expected random score for each cluster size
        df['expected_random'] = model.predict(df[['n', 'n_squared']].values)
        df['expected_random'] = np.maximum(df['expected_random'], 1e-6)  # Avoid division by zero
    else:
        # Fallback: use mean of random clusters
        df['expected_random'] = random_df['mean_tm_non_self'].mean() if len(random_df) > 0 else 0.001
    
    # 2. Calculate normalized scores (observed / expected under null)
    df['enrichment_score'] = df['mean_tm_non_self'] / df['expected_random']
    
    # 3. Z-score normalization within size bins
    df['size_bin'] = pd.qcut(df['n'], q=min(5, len(df)//10), labels=False, duplicates='drop')
    
    def z_normalize(group):
        if len(group) > 1:
            group['z_score'] = (group['mean_tm_non_self'] - group['mean_tm_non_self'].mean()) / group['mean_tm_non_self'].std()
        else:
            group['z_score'] = 0
        return group
    
    df = df.groupby('size_bin').apply(z_normalize)
    
    # 4. Calculate effect size (Cohen's d) relative to random clusters of similar size
    def calculate_effect_size(row):
        similar_random = random_df[np.abs(random_df['n'] - row['n']) <= 1]
        if len(similar_random) > 0:
            random_mean = similar_random['mean_tm_non_self'].mean()
            random_std = similar_random['mean_tm_non_self'].std()
            if random_std > 0:
                return (row['mean_tm_non_self'] - random_mean) / random_std
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
        weights = cat_df['n_comparisons_non_self'] / total_comparisons
        weighted_mean_tm = np.sum(weights * cat_df['mean_tm_non_self'])
        
        # Alternative: Direct calculation from sums
        direct_mean_tm = cat_df['sum_tm_non_self'].sum() / total_comparisons
        
        # Weighted variance
        weighted_var = np.sum(weights * (cat_df['mean_tm_non_self'] - weighted_mean_tm)**2)
        
        # Bootstrap confidence intervals
        bootstrap_means = []
        for _ in range(1000):
            sample = cat_df.sample(n=len(cat_df), replace=True)
            boot_mean = sample['sum_tm_non_self'].sum() / sample['n_comparisons_non_self'].sum()
            bootstrap_means.append(boot_mean)
        
        ci_lower = np.percentile(bootstrap_means, 2.5)
        ci_upper = np.percentile(bootstrap_means, 97.5)
        
        results[category] = {
            'n_clusters': total_clusters,
            'total_proteins': cat_df['n'].sum(),
            'total_comparisons': total_comparisons,
            'weighted_mean_tm': weighted_mean_tm,
            'direct_mean_tm': direct_mean_tm,
            'weighted_std': np.sqrt(weighted_var),
            'ci_95_lower': ci_lower,
            'ci_95_upper': ci_upper,
            'mean_enrichment': cat_df['enrichment_score'].mean(),
            'median_enrichment': cat_df['enrichment_score'].median(),
            'mean_cluster_size': cat_df['n'].mean(),
            'std_cluster_size': cat_df['n'].std()
        }
    
    return results

def statistical_tests(df, group_stats):
    """Perform rigorous statistical tests with size correction"""
    
    test_results = {}
    
    # 1. Permutation test (most robust for this data structure)
    def permutation_test(group1_df, group2_df, n_permutations=10000):
        # Observed difference in weighted means
        obs_diff = (group1_df['sum_tm_non_self'].sum() / group1_df['n_comparisons_non_self'].sum() - 
                   group2_df['sum_tm_non_self'].sum() / group2_df['n_comparisons_non_self'].sum())
        
        # Combine data
        combined = pd.concat([group1_df, group2_df])
        n1 = len(group1_df)
        
        # Permutation
        perm_diffs = []
        for _ in range(n_permutations):
            shuffled = combined.sample(frac=1, replace=False)
            perm_g1 = shuffled.iloc[:n1]
            perm_g2 = shuffled.iloc[n1:]
            perm_diff = (perm_g1['sum_tm_non_self'].sum() / perm_g1['n_comparisons_non_self'].sum() - 
                        perm_g2['sum_tm_non_self'].sum() / perm_g2['n_comparisons_non_self'].sum())
            perm_diffs.append(perm_diff)
        
        p_value = np.mean(np.array(perm_diffs) >= obs_diff)
        return obs_diff, p_value
    
    # 2. ANCOVA-style test (controlling for cluster size)
    from scipy.stats import f_oneway
    
    categories = df['category'].unique()
    
    # Test enrichment scores (size-normalized)
    if len(categories) >= 2:
        groups_enrichment = [df[df['category'] == cat]['enrichment_score'].values for cat in categories]
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
                    u_stat, p_mw = stats.mannwhitneyu(
                        g1['enrichment_score'], 
                        g2['enrichment_score'],
                        alternative='two-sided'
                    )
                    
                    comparisons.append({
                        'comparison': f"{cat1} vs {cat2}",
                        'mean_diff': group_stats[cat1]['weighted_mean_tm'] - group_stats[cat2]['weighted_mean_tm'],
                        'enrichment_diff': g1['enrichment_score'].mean() - g2['enrichment_score'].mean(),
                        'p_permutation': p_perm,
                        'p_mannwhitney': p_mw,
                        'significant': p_perm < 0.05
                    })
    
    test_results['pairwise'] = comparisons
    return test_results

def main():
    base_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test")
    
    # Collect ALL files
    files = []
    for subdir in ["mibig_kegg_foldseek_predictions", "random_foldseek_predictions"]:
        results_dir = base_dir / subdir / "foldseek_results"
        if results_dir.exists():
            files.extend(results_dir.glob("*.tsv"))
    
    print(f"Processing {len(files)} complete TSV files...")
    
    # Process all files completely
    with Pool(cpu_count()) as pool:
        results = pool.map(process_complete_tsv, files)
    
    # Create DataFrame
    df = pd.DataFrame([r for r in results if r is not None])
    
    print(f"Loaded {len(df)} clusters")
    
    # Apply size normalization
    df = calculate_size_normalized_scores(df)
    
    # Calculate group statistics with proper weighting
    group_stats = weighted_group_statistics(df)
    
    # Statistical tests
    test_results = statistical_tests(df, group_stats)
    
    # Print results
    print("\n" + "="*80)
    print("STATISTICALLY NORMALIZED RESULTS")
    print("="*80)
    
    print("\n1. GROUP STATISTICS (Size-Corrected)")
    print("-"*50)
    
    for cat in ['RANDOM', 'MGC (KEGG)', 'BGC (MIBiG)']:
        if cat in group_stats:
            stats = group_stats[cat]
            print(f"\n{cat}:")
            print(f"  Clusters: {stats['n_clusters']}")
            print(f"  Total proteins: {stats['total_proteins']}")
            print(f"  Total comparisons: {stats['total_comparisons']:,}")
            print(f"  Mean cluster size: {stats['mean_cluster_size']:.1f} ± {stats['std_cluster_size']:.1f}")
            print(f"  Weighted mean TM-score: {stats['weighted_mean_tm']:.4f}")
            print(f"  95% CI: [{stats['ci_95_lower']:.4f}, {stats['ci_95_upper']:.4f}]")
            print(f"  Mean enrichment over random: {stats['mean_enrichment']:.2f}x")
    
    print("\n2. STATISTICAL TESTS")
    print("-"*50)
    
    if 'ANOVA_enrichment' in test_results:
        print(f"\nANOVA on enrichment scores: F={test_results['ANOVA_enrichment']['F']:.2f}, p={test_results['ANOVA_enrichment']['p']:.2e}")
    
    print("\nPairwise comparisons (permutation test):")
    for comp in test_results['pairwise']:
        print(f"  {comp['comparison']}: Δ={comp['mean_diff']:.4f}, p={comp['p_permutation']:.3f} {('***' if comp['significant'] else '')}")
    
    # Visualization
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. Raw means by cluster size
    ax = axes[0, 0]
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        ax.scatter(cat_df['n'], cat_df['mean_tm_non_self'], label=cat, alpha=0.6)
    ax.set_xlabel('Cluster Size (n proteins)')
    ax.set_ylabel('Mean TM-score')
    ax.set_title('Raw Scores vs Cluster Size')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Enrichment scores (size-normalized)
    ax = axes[0, 1]
    data_by_cat = [df[df['category'] == cat]['enrichment_score'].values 
                   for cat in ['RANDOM', 'MGC (KEGG)', 'BGC (MIBiG)'] 
                   if cat in df['category'].unique()]
    labels = [cat for cat in ['RANDOM', 'MGC (KEGG)', 'BGC (MIBiG)'] 
              if cat in df['category'].unique()]
    bp = ax.boxplot(data_by_cat, labels=labels, patch_artist=True)
    ax.set_ylabel('Enrichment Score')
    ax.set_title('Size-Normalized Enrichment')
    ax.axhline(y=1, color='r', linestyle='--', alpha=0.5)
    ax.grid(True, alpha=0.3)
    
    # 3. Z-scores within size bins
    ax = axes[0, 2]
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        ax.scatter(cat_df['n'], cat_df['z_score'], label=cat, alpha=0.6)
    ax.set_xlabel('Cluster Size')
    ax.set_ylabel('Z-score')
    ax.set_title('Z-scores Within Size Bins')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. Weighted means with CI
    ax = axes[1, 0]
    cats = []
    means = []
    ci_lower = []
    ci_upper = []
    for cat in ['RANDOM', 'MGC (KEGG)', 'BGC (MIBiG)']:
        if cat in group_stats:
            cats.append(cat)
            means.append(group_stats[cat]['weighted_mean_tm'])
            ci_lower.append(group_stats[cat]['ci_95_lower'])
            ci_upper.append(group_stats[cat]['ci_95_upper'])
    
    if cats:
        x_pos = np.arange(len(cats))
        # Calculate error bars properly
        errors = [[means[i] - ci_lower[i] for i in range(len(means))],
                  [ci_upper[i] - means[i] for i in range(len(means))]]
        # Ensure no negative values
        errors = [[max(0, e) for e in errors[0]], 
                  [max(0, e) for e in errors[1]]]
        
        ax.bar(x_pos, means, yerr=errors, capsize=5, color=['#FF6B6B', '#4ECDC4', '#45B7D1'][:len(cats)])
        ax.set_xticks(x_pos)
        ax.set_xticklabels(cats)
        ax.set_ylabel('Weighted Mean TM-score')
        ax.set_title('Group Means with 95% Bootstrap CI')
        ax.grid(True, alpha=0.3, axis='y')
    
    # 5. Distribution of cluster sizes
    ax = axes[1, 1]
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        ax.hist(cat_df['n'], alpha=0.5, label=cat, bins=20)
    ax.set_xlabel('Cluster Size')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Cluster Sizes')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 6. Effect sizes
    ax = axes[1, 2]
    for cat in df['category'].unique():
        if cat != 'RANDOM':
            cat_df = df[df['category'] == cat]
            ax.scatter(cat_df['n'], cat_df['effect_size'], label=cat, alpha=0.6)
    ax.set_xlabel('Cluster Size')
    ax.set_ylabel("Cohen's d vs Random")
    ax.set_title('Effect Size Relative to Random')
    ax.axhline(y=0, color='r', linestyle='--', alpha=0.5)
    ax.axhline(y=0.8, color='g', linestyle='--', alpha=0.3)  # Large effect threshold
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Statistically Normalized Foldseek Analysis', fontsize=14)
    plt.tight_layout()
    
    output_path = base_dir / 'statistically_normalized_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Export results
    df.to_csv(base_dir / 'normalized_cluster_metrics.csv', index=False)
    pd.DataFrame(group_stats).T.to_csv(base_dir / 'group_statistics.csv')
    
    print(f"\n\nResults saved to {base_dir}")
    
    plt.show()
    
    return df, group_stats, test_results

if __name__ == "__main__":
    df, group_stats, test_results = main()