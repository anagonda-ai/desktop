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
from sklearn.linear_model import  LogisticRegression
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
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

def find_optimal_threshold(df, score_column='mean_tm_non_self'):
    """
    Find optimal threshold for separating positives (MIBiG + MGC) from negatives (RANDOM)
    using ROC curve analysis and multiple optimization criteria.
    """
    # Create binary labels: 1 for positives (MIBiG + MGC), 0 for negatives (RANDOM)
    df['is_positive'] = df['category'].isin(['BGC (MIBiG)', 'MGC (KEGG)']).astype(int)
    
    # Get scores and labels
    scores = df[score_column].values
    labels = df['is_positive'].values
    
    # Remove any NaN values
    valid_mask = ~np.isnan(scores)
    scores = scores[valid_mask]
    labels = labels[valid_mask]
    
    if len(np.unique(labels)) < 2:
        print(f"Warning: Not enough classes for {score_column}")
        return None
    
    # Calculate ROC curve
    fpr, tpr, thresholds = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)
    
    # Calculate precision-recall curve
    precision, recall, pr_thresholds = precision_recall_curve(labels, scores)
    pr_auc = auc(recall, precision)
    
    # Find optimal thresholds using different criteria
    results = {}
    
    # 1. Youden's J statistic (maximizes TPR - FPR)
    youden_j = tpr - fpr
    optimal_idx_youden = np.argmax(youden_j)
    optimal_threshold_youden = thresholds[optimal_idx_youden]
    
    # 2. Maximum F1 score
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
    optimal_idx_f1 = np.argmax(f1_scores)
    optimal_threshold_f1 = pr_thresholds[optimal_idx_f1]
    
    # 3. Maximum precision (for high-confidence predictions)
    optimal_idx_precision = np.argmax(precision)
    optimal_threshold_precision = pr_thresholds[optimal_idx_precision]
    
    # 4. Balanced accuracy (maximizes (TPR + TNR) / 2)
    tnr = 1 - fpr  # True Negative Rate
    balanced_acc = (tpr + tnr) / 2
    optimal_idx_balanced = np.argmax(balanced_acc)
    optimal_threshold_balanced = thresholds[optimal_idx_balanced]
    
    # 5. Maximum Matthews Correlation Coefficient (MCC)
    mcc_scores = []
    for thresh in thresholds:
        pred = (scores >= thresh).astype(int)
        tn, fp, fn, tp = confusion_matrix(labels, pred).ravel()
        mcc = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-8)
        mcc_scores.append(mcc)
    
    optimal_idx_mcc = np.argmax(mcc_scores)
    optimal_threshold_mcc = thresholds[optimal_idx_mcc]
    
    # Calculate metrics for each optimal threshold
    def calculate_metrics_at_threshold(threshold):
        pred = (scores >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(labels, pred).ravel()
        
        precision_val = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall_val = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity_val = tn / (tn + fp) if (tn + fp) > 0 else 0
        f1_val = 2 * (precision_val * recall_val) / (precision_val + recall_val + 1e-8)
        mcc_val = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-8)
        accuracy_val = (tp + tn) / (tp + tn + fp + fn)
        balanced_acc_val = (recall_val + specificity_val) / 2
        
        return {
            'threshold': threshold,
            'precision': precision_val,
            'recall': recall_val,
            'specificity': specificity_val,
            'f1_score': f1_val,
            'mcc': mcc_val,
            'accuracy': accuracy_val,
            'balanced_accuracy': balanced_acc_val,
            'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn
        }
    
    # Store results for each method
    methods = {
        'youden': (optimal_threshold_youden, optimal_idx_youden),
        'f1': (optimal_threshold_f1, optimal_idx_f1),
        'precision': (optimal_threshold_precision, optimal_idx_precision),
        'balanced': (optimal_threshold_balanced, optimal_idx_balanced),
        'mcc': (optimal_threshold_mcc, optimal_idx_mcc)
    }
    
    for method_name, (threshold, idx) in methods.items():
        results[method_name] = calculate_metrics_at_threshold(threshold)
    
    # Overall results
    results['roc_auc'] = roc_auc
    results['pr_auc'] = pr_auc
    results['score_column'] = score_column
    results['n_positives'] = np.sum(labels)
    results['n_negatives'] = len(labels) - np.sum(labels)
    results['roc_curve'] = (fpr, tpr)
    results['pr_curve'] = (precision, recall)
    
    return results

def evaluate_all_scores(df):
    """
    Evaluate threshold optimization for all available score columns.
    """
    score_columns = ['mean_tm_non_self', 'enrichment_score', 'z_score', 'effect_size']
    all_results = {}
    
    print("\n" + "="*80)
    print("THRESHOLD OPTIMIZATION RESULTS")
    print("="*80)
    
    for score_col in score_columns:
        if score_col not in df.columns:
            print(f"Skipping {score_col} - not available in data")
            continue
            
        print(f"\nAnalyzing {score_col}:")
        print("-" * 50)
        
        results = find_optimal_threshold(df, score_col)
        if results is None:
            continue
            
        all_results[score_col] = results
        
        # Print results for each optimization method
        print(f"ROC AUC: {results['roc_auc']:.3f}")
        print(f"PR AUC: {results['pr_auc']:.3f}")
        print(f"Data: {results['n_positives']} positives, {results['n_negatives']} negatives")
        
        print("\nOptimal thresholds by method:")
        for method in ['youden', 'f1', 'precision', 'balanced', 'mcc']:
            if method in results:
                m = results[method]
                print(f"  {method.upper()}: threshold={m['threshold']:.4f}, "
                      f"F1={m['f1_score']:.3f}, MCC={m['mcc']:.3f}, "
                      f"Precision={m['precision']:.3f}, Recall={m['recall']:.3f}")
    
    return all_results

def robust_threshold_optimization_with_validation(df, score_column='mean_tm_non_self', 
                                                test_size=0.3, n_splits=5, random_state=42):
    """
    Perform robust threshold optimization with proper train/test splits and cross-validation
    to prevent overfitting and ensure generalizability.
    """
    print(f"\n{'='*70}")
    print("ROBUST THRESHOLD OPTIMIZATION WITH VALIDATION")
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
    
    # Train/Test Split
    X_train, X_test, y_train, y_test = train_test_split(
        scores.reshape(-1, 1), binary_labels, 
        test_size=test_size, stratify=binary_labels, 
        random_state=random_state
    )
    
    print(f"\nDataset Split:")
    print(f"  Training set: {len(X_train)} samples ({np.sum(y_train)} positive, {len(y_train)-np.sum(y_train)} negative)")
    print(f"  Test set: {len(X_test)} samples ({np.sum(y_test)} positive, {len(y_test)-np.sum(y_test)} negative)")
    
    # 1. Find optimal thresholds on training set
    print(f"\n1. THRESHOLD OPTIMIZATION ON TRAINING SET")
    print("-" * 50)
    
    # Calculate ROC curve on training data
    fpr_train, tpr_train, thresholds_train = roc_curve(y_train, X_train.flatten())
    roc_auc_train = auc(fpr_train, tpr_train)
    
    # Calculate precision-recall curve on training data
    precision_train, recall_train, pr_thresholds_train = precision_recall_curve(y_train, X_train.flatten())
    pr_auc_train = auc(recall_train, precision_train)
    
    # Find optimal thresholds on training set
    train_results = {}
    
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
    
    train_results = {
        'youden': {'threshold': optimal_threshold_youden},
        'f1': {'threshold': optimal_threshold_f1},
        'accuracy': {'threshold': optimal_threshold_balanced},
        'roc_auc': roc_auc_train,
        'pr_auc': pr_auc_train
    }
    
    # Display training results
    print(f"\nROC AUC (Training): {roc_auc_train:.4f}")
    print(f"PR AUC (Training): {pr_auc_train:.4f}")
    print(f"\nOptimal thresholds on training set:")
    print(f"  Youden's J: {optimal_threshold_youden:.4f}")
    print(f"  F1-optimized: {optimal_threshold_f1:.4f}")
    print(f"  Accuracy-optimized: {optimal_threshold_balanced:.4f}")
    
    # 2. Evaluate on test set using training thresholds
    print(f"\n2. VALIDATION ON INDEPENDENT TEST SET")
    print("-" * 50)
    
    test_results = {}
    for metric_name, threshold_data in train_results.items():
        if metric_name in ['youden', 'f1', 'accuracy']:
            threshold = threshold_data['threshold']
            
            # Apply threshold to test set
            y_pred_test = (X_test.flatten() >= threshold).astype(int)
            
            # Calculate test metrics
            tn = np.sum((y_pred_test == 0) & (y_test == 0))
            fp = np.sum((y_pred_test == 1) & (y_test == 0))
            fn = np.sum((y_pred_test == 0) & (y_test == 1))
            tp = np.sum((y_pred_test == 1) & (y_test == 1))
            
            test_sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
            test_specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            test_precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            test_recall = test_sensitivity
            test_accuracy = (tp + tn) / len(y_test)
            test_f1 = 2 * (test_precision * test_recall) / (test_precision + test_recall) if (test_precision + test_recall) > 0 else 0
            
            test_results[metric_name] = {
                'threshold': threshold,
                'sensitivity': test_sensitivity,
                'specificity': test_specificity,
                'precision': test_precision,
                'recall': test_recall,
                'accuracy': test_accuracy,
                'f1_score': test_f1
            }
            
            print(f"\n{metric_name.upper()} threshold validation:")
            print(f"  Threshold: {threshold:.4f}")
            print(f"  Test Accuracy: {test_accuracy:.4f}")
            print(f"  Test Sensitivity: {test_sensitivity:.4f}")
            print(f"  Test Specificity: {test_specificity:.4f}")
            print(f"  Test Precision: {test_precision:.4f}")
            print(f"  Test F1 Score: {test_f1:.4f}")
    
    # 3. Cross-validation analysis
    print(f"\n3. CROSS-VALIDATION ANALYSIS")
    print("-" * 50)
    
    cv_results = cross_validation_analysis(scores, binary_labels, n_splits, random_state)
    
    # 4. Test set ROC analysis
    fpr_test, tpr_test, _ = roc_curve(y_test, X_test.flatten())
    roc_auc_test = auc(fpr_test, tpr_test)
    
    precision_test, recall_test, _ = precision_recall_curve(y_test, X_test.flatten())
    pr_auc_test = auc(recall_test, precision_test)
    
    print(f"\n4. TEST SET PERFORMANCE")
    print("-" * 50)
    print(f"  Test ROC AUC: {roc_auc_test:.4f}")
    print(f"  Test PR AUC: {pr_auc_test:.4f}")
    
    # Generalizability assessment
    print(f"\n{'='*60}")
    print("GENERALIZABILITY ASSESSMENT")
    print("="*60)
    
    cv_roc_auc_mean = cv_results['roc_auc']['mean']
    cv_roc_auc_std = cv_results['roc_auc']['std']
    
    if cv_roc_auc_std < 0.05:
        stability = "Highly stable"
    elif cv_roc_auc_std < 0.1:
        stability = "Stable"
    else:
        stability = "Variable - consider more data or feature engineering"
    
    print(f"\nModel Stability: {stability}")
    print(f"  CV ROC AUC variability: {cv_roc_auc_std:.4f}")
    
    if abs(roc_auc_test - cv_roc_auc_mean) < 0.05:
        generalization = "Excellent generalization"
    elif abs(roc_auc_test - cv_roc_auc_mean) < 0.1:
        generalization = "Good generalization"
    else:
        generalization = "Potential overfitting - validate on additional data"
    
    print(f"\nGeneralization Quality: {generalization}")
    print(f"  Test vs CV AUC difference: {abs(roc_auc_test - cv_roc_auc_mean):.4f}")
    
    # Final recommendations
    print(f"\n{'='*60}")
    print("ROBUST THRESHOLD RECOMMENDATIONS")
    print("="*60)
    
    if roc_auc_test >= 0.8 and cv_roc_auc_std < 0.1:
        recommendation = "Highly recommended for production use"
    elif roc_auc_test >= 0.7 and cv_roc_auc_std < 0.15:
        recommendation = "Suitable for use with monitoring"
    else:
        recommendation = "Requires improvement before production use"
    
    print(f"\nOverall Recommendation: {recommendation}")
    
    if 'f1' in test_results:
        best_threshold = test_results['f1']['threshold']
        best_f1 = test_results['f1']['f1_score']
        print(f"\nRecommended Production Threshold: {best_threshold:.4f}")
        print(f"  Expected F1 Score: {best_f1:.4f}")
        print(f"  Expected Precision: {test_results['f1']['precision']:.4f}")
        print(f"  Expected Recall: {test_results['f1']['recall']:.4f}")
    
    return {
        'train_results': train_results,
        'test_results': test_results,
        'cv_results': cv_results,
        'test_roc_auc': roc_auc_test,
        'test_pr_auc': pr_auc_test,
        'train_size': len(X_train),
        'test_size': len(X_test),
        'score_column': score_column
    }

def cross_validation_analysis(scores, binary_labels, n_splits=5, random_state=42):
    """
    Perform k-fold cross-validation to assess model generalizability.
    """
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    
    # Simple logistic regression model
    lr_model = LogisticRegression(random_state=random_state, max_iter=1000)
    
    # Cross-validation scores
    cv_accuracy = cross_val_score(lr_model, scores.reshape(-1, 1), binary_labels, 
                                cv=skf, scoring='accuracy')
    cv_precision = cross_val_score(lr_model, scores.reshape(-1, 1), binary_labels, 
                                 cv=skf, scoring='precision')
    cv_recall = cross_val_score(lr_model, scores.reshape(-1, 1), binary_labels, 
                              cv=skf, scoring='recall')
    cv_f1 = cross_val_score(lr_model, scores.reshape(-1, 1), binary_labels, 
                           cv=skf, scoring='f1')
    cv_roc_auc = cross_val_score(lr_model, scores.reshape(-1, 1), binary_labels, 
                               cv=skf, scoring='roc_auc')
    
    print(f"\n{n_splits}-Fold Cross-Validation Results:")
    print(f"  Accuracy:  {cv_accuracy.mean():.4f} ± {cv_accuracy.std():.4f}")
    print(f"  Precision: {cv_precision.mean():.4f} ± {cv_precision.std():.4f}")
    print(f"  Recall:    {cv_recall.mean():.4f} ± {cv_recall.std():.4f}")
    print(f"  F1 Score:  {cv_f1.mean():.4f} ± {cv_f1.std():.4f}")
    print(f"  ROC AUC:   {cv_roc_auc.mean():.4f} ± {cv_roc_auc.std():.4f}")
    
    return {
        'accuracy': {'mean': cv_accuracy.mean(), 'std': cv_accuracy.std(), 'scores': cv_accuracy},
        'precision': {'mean': cv_precision.mean(), 'std': cv_precision.std(), 'scores': cv_precision},
        'recall': {'mean': cv_recall.mean(), 'std': cv_recall.std(), 'scores': cv_recall},
        'f1': {'mean': cv_f1.mean(), 'std': cv_f1.std(), 'scores': cv_f1},
        'roc_auc': {'mean': cv_roc_auc.mean(), 'std': cv_roc_auc.std(), 'scores': cv_roc_auc}
    }

def evaluate_all_scores_with_validation(df):
    """
    Evaluate robust threshold optimization for all available score columns.
    """
    score_columns = ['mean_tm_non_self', 'enrichment_score', 'z_score', 'effect_size']
    all_results = {}
    
    print("\n" + "="*80)
    print("ROBUST THRESHOLD OPTIMIZATION WITH VALIDATION")
    print("="*80)
    
    for score_col in score_columns:
        if score_col not in df.columns:
            print(f"Skipping {score_col} - not available in data")
            continue
            
        print(f"\n{'='*70}")
        print(f"ANALYZING {score_col.upper()}")
        print("="*70)
        
        results = robust_threshold_optimization_with_validation(df, score_col)
        if results is None:
            continue
            
        all_results[score_col] = results
    
    return all_results

def plot_threshold_analysis(df, all_results, output_dir):
    """
    Create comprehensive visualizations for threshold analysis.
    """
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    
    # 1. Score distributions by category
    ax = axes[0, 0]
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        ax.hist(cat_df['mean_tm_non_self'], alpha=0.6, label=cat, bins=30)
    ax.set_xlabel('Mean TM-score')
    ax.set_ylabel('Count')
    ax.set_title('TM-score Distribution by Category')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. ROC curves for all scores
    ax = axes[0, 1]
    colors = ['blue', 'red', 'green', 'orange']
    for i, (score_col, results) in enumerate(all_results.items()):
        if 'roc_curve' in results:
            fpr, tpr = results['roc_curve']
            ax.plot(fpr, tpr, color=colors[i % len(colors)], 
                   label=f'{score_col} (AUC={results["roc_auc"]:.3f})')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curves')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Precision-Recall curves
    ax = axes[0, 2]
    for i, (score_col, results) in enumerate(all_results.items()):
        if 'pr_curve' in results:
            precision, recall = results['pr_curve']
            ax.plot(recall, precision, color=colors[i % len(colors)], 
                   label=f'{score_col} (AUC={results["pr_auc"]:.3f})')
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Precision-Recall Curves')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. Threshold comparison heatmap
    ax = axes[1, 0]
    methods = ['youden', 'f1', 'precision', 'balanced', 'mcc']
    metrics = ['f1_score', 'mcc', 'precision', 'recall', 'balanced_accuracy']
    
    heatmap_data = []
    score_cols = []
    for score_col, results in all_results.items():
        row = []
        for metric in metrics:
            if 'youden' in results:
                row.append(results['youden'][metric])
        heatmap_data.append(row)
        score_cols.append(score_col)
    
    if heatmap_data:
        im = ax.imshow(heatmap_data, cmap='viridis', aspect='auto')
        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels(metrics, rotation=45)
        ax.set_yticks(range(len(score_cols)))
        ax.set_yticklabels(score_cols)
        ax.set_title('Performance Metrics Heatmap')
        
        # Add colorbar
        plt.colorbar(im, ax=ax)
    
    # 5. Score vs cluster size with thresholds
    ax = axes[1, 1]
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        ax.scatter(cat_df['n'], cat_df['mean_tm_non_self'], 
                  label=cat, alpha=0.6, s=20)
    
    # Add threshold lines
    if 'mean_tm_non_self' in all_results:
        best_threshold = all_results['mean_tm_non_self']['youden']['threshold']
        ax.axhline(y=best_threshold, color='red', linestyle='--', 
                  label=f'Optimal threshold: {best_threshold:.3f}')
    
    ax.set_xlabel('Cluster Size')
    ax.set_ylabel('Mean TM-score')
    ax.set_title('TM-score vs Cluster Size with Threshold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 6. Enrichment score distribution
    ax = axes[1, 2]
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        ax.hist(cat_df['enrichment_score'], alpha=0.6, label=cat, bins=30)
    ax.set_xlabel('Enrichment Score')
    ax.set_ylabel('Count')
    ax.set_title('Enrichment Score Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 7. Confusion matrix for best threshold
    ax = axes[2, 0]
    if 'mean_tm_non_self' in all_results:
        best_results = all_results['mean_tm_non_self']['youden']
        cm_data = [[best_results['tn'], best_results['fp']],
                   [best_results['fn'], best_results['tp']]]
        im = ax.imshow(cm_data, cmap='Blues')
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(['Predicted Negative', 'Predicted Positive'])
        ax.set_yticklabels(['Actual Negative', 'Actual Positive'])
        ax.set_title('Confusion Matrix (Youden\'s J)')
        
        # Add text annotations
        for i in range(2):
            for j in range(2):
                ax.text(j, i, str(cm_data[i][j]), ha='center', va='center', 
                       color='white' if cm_data[i][j] > max(map(max, cm_data))/2 else 'black')
    
    # 8. Threshold sensitivity analysis
    ax = axes[2, 1]
    if 'mean_tm_non_self' in all_results:
        # Create threshold range
        scores = df['mean_tm_non_self'].values
        thresholds = np.linspace(scores.min(), scores.max(), 100)
        
        f1_scores = []
        mcc_scores = []
        
        for thresh in thresholds:
            pred = (scores >= thresh).astype(int)
            labels = df['category'].isin(['BGC (MIBiG)', 'MGC (KEGG)']).astype(int)
            
            if len(np.unique(pred)) > 1 and len(np.unique(labels)) > 1:
                tn, fp, fn, tp = confusion_matrix(labels, pred).ravel()
                precision = tp / (tp + fp) if (tp + fp) > 0 else 0
                recall = tp / (tp + fn) if (tp + fn) > 0 else 0
                f1 = 2 * (precision * recall) / (precision + recall + 1e-8)
                mcc = (tp * tn - fp * fn) / np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) + 1e-8)
                
                f1_scores.append(f1)
                mcc_scores.append(mcc)
            else:
                f1_scores.append(0)
                mcc_scores.append(0)
        
        ax.plot(thresholds, f1_scores, label='F1 Score', color='blue')
        ax.plot(thresholds, mcc_scores, label='MCC', color='red')
        ax.set_xlabel('Threshold')
        ax.set_ylabel('Score')
        ax.set_title('Threshold Sensitivity')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # 9. Performance comparison
    ax = axes[2, 2]
    if all_results:
        score_names = list(all_results.keys())
        f1_values = [all_results[score]['youden']['f1_score'] for score in score_names]
        mcc_values = [all_results[score]['youden']['mcc'] for score in score_names]
        
        x = np.arange(len(score_names))
        width = 0.35
        
        ax.bar(x - width/2, f1_values, width, label='F1 Score', alpha=0.8)
        ax.bar(x + width/2, mcc_values, width, label='MCC', alpha=0.8)
        
        ax.set_xlabel('Score Type')
        ax.set_ylabel('Performance')
        ax.set_title('Performance Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels(score_names, rotation=45)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('Comprehensive Threshold Analysis', fontsize=16)
    plt.tight_layout()
    
    # Save the plot
    output_path = output_dir / 'threshold_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Threshold analysis plot saved to: {output_path}")
    
    plt.show()

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
    
    # Robust threshold optimization with validation
    print("\n" + "="*80)
    print("ROBUST THRESHOLD OPTIMIZATION WITH VALIDATION")
    print("="*80)
    
    # First run the old method for comparison
    print("\nBasic threshold optimization (for reference):")
    basic_threshold_results = evaluate_all_scores(df)
    
    # Then run the robust validation
    robust_threshold_results = evaluate_all_scores_with_validation(df)
    
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
    
    # Print robust threshold validation results
    print("\n" + "="*80)
    print("ROBUST THRESHOLD VALIDATION SUMMARY")
    print("="*80)
    
    if robust_threshold_results:
        for score_col, results in robust_threshold_results.items():
            print(f"\n{'='*60}")
            print(f"VALIDATED RESULTS FOR {score_col.upper()}")
            print("="*60)
            
            print(f"\nModel Generalizability Assessment:")
            print(f"  Training samples: {results['train_size']}")
            print(f"  Test samples: {results['test_size']}")
            print(f"  Test ROC AUC: {results['test_roc_auc']:.4f}")
            print(f"  Test PR AUC: {results['test_pr_auc']:.4f}")
            
            # Cross-validation stability
            cv_results = results['cv_results']
            print(f"\nCross-Validation Stability:")
            print(f"  CV ROC AUC: {cv_results['roc_auc']['mean']:.4f} ± {cv_results['roc_auc']['std']:.4f}")
            print(f"  CV F1 Score: {cv_results['f1']['mean']:.4f} ± {cv_results['f1']['std']:.4f}")
            print(f"  CV Precision: {cv_results['precision']['mean']:.4f} ± {cv_results['precision']['std']:.4f}")
            print(f"  CV Recall: {cv_results['recall']['mean']:.4f} ± {cv_results['recall']['std']:.4f}")
            
            # Test set performance for each threshold
            test_results_data = results['test_results']
            print(f"\nTest Set Performance by Threshold:")
            for metric_name in ['youden', 'f1', 'accuracy']:
                if metric_name in test_results_data:
                    result = test_results_data[metric_name]
                    print(f"\n{metric_name.upper()} threshold ({result['threshold']:.4f}):")
                    print(f"  Test Accuracy: {result['accuracy']:.4f}")
                    print(f"  Test Precision: {result['precision']:.4f}")
                    print(f"  Test Recall: {result['recall']:.4f}")
                    print(f"  Test F1 Score: {result['f1_score']:.4f}")
            
            # Model generalizability assessment
            test_roc_auc = results['test_roc_auc']
            cv_roc_auc_mean = cv_results['roc_auc']['mean']
            cv_roc_auc_std = cv_results['roc_auc']['std']
            
            print(f"\nGeneralizability Assessment:")
            
            if cv_roc_auc_std < 0.05:
                stability = "Highly stable"
            elif cv_roc_auc_std < 0.1:
                stability = "Stable"
            else:
                stability = "Variable - consider more data or feature engineering"
            
            print(f"  Model Stability: {stability}")
            print(f"  CV ROC AUC variability: {cv_roc_auc_std:.4f}")
            
            if abs(test_roc_auc - cv_roc_auc_mean) < 0.05:
                generalization = "Excellent generalization"
            elif abs(test_roc_auc - cv_roc_auc_mean) < 0.1:
                generalization = "Good generalization"
            else:
                generalization = "Potential overfitting - validate on additional data"
            
            print(f"  Generalization Quality: {generalization}")
            print(f"  Test vs CV AUC difference: {abs(test_roc_auc - cv_roc_auc_mean):.4f}")
            
            # Final recommendations
            print(f"\nProduction Recommendations:")
            
            if test_roc_auc >= 0.8 and cv_roc_auc_std < 0.1:
                recommendation = "Highly recommended for production use"
            elif test_roc_auc >= 0.7 and cv_roc_auc_std < 0.15:
                recommendation = "Suitable for use with monitoring"
            else:
                recommendation = "Requires improvement before production use"
            
            print(f"  Overall Recommendation: {recommendation}")
            
            if 'f1' in test_results_data:
                best_threshold = test_results_data['f1']['threshold']
                best_f1 = test_results_data['f1']['f1_score']
                print(f"\nRecommended Production Threshold: {best_threshold:.4f}")
                print(f"  Expected F1 Score: {best_f1:.4f}")
                print(f"  Expected Precision: {test_results_data['f1']['precision']:.4f}")
                print(f"  Expected Recall: {test_results_data['f1']['recall']:.4f}")
                print(f"  Expected Accuracy: {test_results_data['f1']['accuracy']:.4f}")
    
    else:
        print("\nNo robust threshold validation results available.")
    
    # Final summary comparison
    print(f"\n{'='*80}")
    print("FINAL SUMMARY - BEST PERFORMING FEATURES")
    print("="*80)
    
    if robust_threshold_results:
        print("\nComparison of validated performance across features:")
        print("-" * 60)
        
        feature_performance = []
        for score_col, results in robust_threshold_results.items():
            if 'test_roc_auc' in results:
                feature_performance.append({
                    'feature': score_col,
                    'test_roc_auc': results['test_roc_auc'],
                    'test_pr_auc': results['test_pr_auc'],
                    'cv_roc_auc': results['cv_results']['roc_auc']['mean'],
                    'cv_stability': results['cv_results']['roc_auc']['std']
                })
        
        # Sort by test ROC AUC
        feature_performance.sort(key=lambda x: x['test_roc_auc'], reverse=True)
        
        print(f"\nRanked by Test ROC AUC:")
        for i, perf in enumerate(feature_performance, 1):
            print(f"  {i}. {perf['feature']}:")
            print(f"     Test ROC AUC: {perf['test_roc_auc']:.4f}")
            print(f"     Test PR AUC: {perf['test_pr_auc']:.4f}")
            print(f"     CV Stability: {perf['cv_stability']:.4f}")
            
            if perf['test_roc_auc'] >= 0.9:
                performance_level = "Outstanding"
            elif perf['test_roc_auc'] >= 0.8:
                performance_level = "Excellent"
            elif perf['test_roc_auc'] >= 0.7:
                performance_level = "Good"
            else:
                performance_level = "Poor"
            
            print(f"     Performance: {performance_level}")
            print()
    
    return df, group_stats, test_results, basic_threshold_results, robust_threshold_results

if __name__ == "__main__":
    df, group_stats, test_results, basic_threshold_results, robust_threshold_results = main()