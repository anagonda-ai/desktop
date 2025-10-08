#!/usr/bin/env python3
"""
Statistical analysis and classification threshold optimization for E2P2 features.
Applies the same rigorous statistical methodology as foldseek analysis.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

def load_e2p2_features(feature_file):
    """Load E2P2 feature extraction results"""
    try:
        df = pd.read_csv(feature_file)
        print(f"Loaded {len(df)} clusters from {feature_file}")
        return df
    except Exception as e:
        print(f"Error loading features: {e}")
        return None

def create_binary_labels(df):
    """Create binary labels: 1 for BGC/MGC_CANDIDATE, 0 for RANDOM"""
    df['is_positive'] = df['classification_label'].astype(int)
    return df

def calculate_feature_statistics(df):
    """Calculate comprehensive feature statistics"""
    stats_dict = {}
    
    # Overall statistics
    stats_dict['total_clusters'] = len(df)
    stats_dict['positive_clusters'] = df['is_positive'].sum()
    stats_dict['negative_clusters'] = len(df) - df['is_positive'].sum()
    
    # Feature statistics
    for feature in ['num_distinct_enzyme_classes', 'num_distinct_enzyme_subfamilies', 'total_ec_numbers']:
        if feature in df.columns:
            stats_dict[f'{feature}_mean'] = df[feature].mean()
            stats_dict[f'{feature}_std'] = df[feature].std()
            stats_dict[f'{feature}_median'] = df[feature].median()
            stats_dict[f'{feature}_q25'] = df[feature].quantile(0.25)
            stats_dict[f'{feature}_q75'] = df[feature].quantile(0.75)
    
    return stats_dict

def group_comparison_analysis(df):
    """Perform group comparison analysis between positive and negative clusters"""
    positive_df = df[df['is_positive'] == 1]
    negative_df = df[df['is_positive'] == 0]
    
    print(f"\nGroup Comparison Analysis:")
    print(f"  Positive clusters (BGC/MGC): {len(positive_df)}")
    print(f"  Negative clusters (RANDOM): {len(negative_df)}")
    
    results = {}
    
    for feature in ['num_distinct_enzyme_classes', 'num_distinct_enzyme_subfamilies', 'total_ec_numbers']:
        if feature in df.columns:
            pos_values = positive_df[feature].values
            neg_values = negative_df[feature].values
            
            # Statistical tests
            t_stat, t_p = stats.ttest_ind(pos_values, neg_values)
            u_stat, u_p = stats.mannwhitneyu(pos_values, neg_values, alternative='two-sided')
            
            # Effect size (Cohen's d)
            pooled_std = np.sqrt(((len(pos_values) - 1) * np.var(pos_values) + 
                                 (len(neg_values) - 1) * np.var(neg_values)) / 
                                (len(pos_values) + len(neg_values) - 2))
            cohens_d = (np.mean(pos_values) - np.mean(neg_values)) / pooled_std
            
            results[feature] = {
                'positive_mean': np.mean(pos_values),
                'negative_mean': np.mean(neg_values),
                'positive_std': np.std(pos_values),
                'negative_std': np.std(neg_values),
                't_statistic': t_stat,
                't_p_value': t_p,
                'mannwhitney_u': u_stat,
                'mannwhitney_p': u_p,
                'cohens_d': cohens_d,
                'effect_size': 'large' if abs(cohens_d) > 0.8 else 'medium' if abs(cohens_d) > 0.5 else 'small'
            }
            
            print(f"\n{feature}:")
            print(f"  Positive: {np.mean(pos_values):.2f} ± {np.std(pos_values):.2f}")
            print(f"  Negative: {np.mean(neg_values):.2f} ± {np.std(neg_values):.2f}")
            print(f"  t-test: t={t_stat:.3f}, p={t_p:.3f}")
            print(f"  Mann-Whitney: U={u_stat:.3f}, p={u_p:.3f}")
            print(f"  Cohen's d: {cohens_d:.3f} ({results[feature]['effect_size']} effect)")
    
    return results

def find_optimal_threshold_single_feature(df, feature_name):
    """Find optimal threshold for a single feature using ROC analysis"""
    if feature_name not in df.columns:
        return None
    
    scores = df[feature_name].values
    labels = df['is_positive'].values
    
    # Remove NaN values
    valid_mask = ~np.isnan(scores)
    scores = scores[valid_mask]
    labels = labels[valid_mask]
    
    if len(np.unique(labels)) < 2:
        return None
    
    # Calculate ROC curve
    fpr, tpr, thresholds = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)
    
    # Calculate precision-recall curve
    precision, recall, pr_thresholds = precision_recall_curve(labels, scores)
    pr_auc = auc(recall, precision)
    
    # Find optimal thresholds
    results = {}
    
    # Youden's J statistic
    youden_j = tpr - fpr
    optimal_idx_youden = np.argmax(youden_j)
    optimal_threshold_youden = thresholds[optimal_idx_youden]
    
    # F1 score
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
    optimal_idx_f1 = np.argmax(f1_scores)
    optimal_threshold_f1 = pr_thresholds[optimal_idx_f1]
    
    # Balanced accuracy
    tnr = 1 - fpr
    balanced_acc = (tpr + tnr) / 2
    optimal_idx_balanced = np.argmax(balanced_acc)
    optimal_threshold_balanced = thresholds[optimal_idx_balanced]
    
    # Calculate metrics for each threshold
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
    
    # Store results
    methods = {
        'youden': optimal_threshold_youden,
        'f1': optimal_threshold_f1,
        'balanced': optimal_threshold_balanced
    }
    
    for method_name, threshold in methods.items():
        results[method_name] = calculate_metrics_at_threshold(threshold)
    
    results['roc_auc'] = roc_auc
    results['pr_auc'] = pr_auc
    results['feature_name'] = feature_name
    results['n_positives'] = np.sum(labels)
    results['n_negatives'] = len(labels) - np.sum(labels)
    
    return results

def robust_threshold_optimization_with_validation(df, feature_name, test_size=0.3, n_splits=5, random_state=42):
    """Perform robust threshold optimization with proper train/test splits"""
    if feature_name not in df.columns:
        return None
    
    scores = df[feature_name].values
    labels = df['is_positive'].values
    
    # Remove NaN values
    valid_mask = ~np.isnan(scores)
    scores = scores[valid_mask]
    labels = labels[valid_mask]
    
    if len(np.unique(labels)) < 2:
        return None
    
    # Train/Test Split
    X_train, X_test, y_train, y_test = train_test_split(
        scores.reshape(-1, 1), labels, 
        test_size=test_size, stratify=labels, 
        random_state=random_state
    )
    
    print(f"\nRobust threshold optimization for {feature_name}:")
    print(f"  Training set: {len(X_train)} samples ({np.sum(y_train)} positive, {len(y_train)-np.sum(y_train)} negative)")
    print(f"  Test set: {len(X_test)} samples ({np.sum(y_test)} positive, {len(y_test)-np.sum(y_test)} negative)")
    
    # Find optimal thresholds on training set
    fpr_train, tpr_train, thresholds_train = roc_curve(y_train, X_train.flatten())
    roc_auc_train = auc(fpr_train, tpr_train)
    
    precision_train, recall_train, pr_thresholds_train = precision_recall_curve(y_train, X_train.flatten())
    pr_auc_train = auc(recall_train, precision_train)
    
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
    
    # Evaluate on test set
    test_results = {}
    for method_name, threshold in [('youden', optimal_threshold_youden), 
                                  ('f1', optimal_threshold_f1), 
                                  ('balanced', optimal_threshold_balanced)]:
        
        y_pred_test = (X_test.flatten() >= threshold).astype(int)
        
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
        
        test_results[method_name] = {
            'threshold': threshold,
            'sensitivity': test_sensitivity,
            'specificity': test_specificity,
            'precision': test_precision,
            'recall': test_recall,
            'accuracy': test_accuracy,
            'f1_score': test_f1
        }
    
    # Cross-validation
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    lr_model = LogisticRegression(random_state=random_state, max_iter=1000)
    
    cv_accuracy = cross_val_score(lr_model, scores.reshape(-1, 1), labels, cv=skf, scoring='accuracy')
    cv_precision = cross_val_score(lr_model, scores.reshape(-1, 1), labels, cv=skf, scoring='precision')
    cv_recall = cross_val_score(lr_model, scores.reshape(-1, 1), labels, cv=skf, scoring='recall')
    cv_f1 = cross_val_score(lr_model, scores.reshape(-1, 1), labels, cv=skf, scoring='f1')
    cv_roc_auc = cross_val_score(lr_model, scores.reshape(-1, 1), labels, cv=skf, scoring='roc_auc')
    
    cv_results = {
        'accuracy': {'mean': cv_accuracy.mean(), 'std': cv_accuracy.std()},
        'precision': {'mean': cv_precision.mean(), 'std': cv_precision.std()},
        'recall': {'mean': cv_recall.mean(), 'std': cv_recall.std()},
        'f1': {'mean': cv_f1.mean(), 'std': cv_f1.std()},
        'roc_auc': {'mean': cv_roc_auc.mean(), 'std': cv_roc_auc.std()}
    }
    
    # Test set ROC analysis
    fpr_test, tpr_test, _ = roc_curve(y_test, X_test.flatten())
    roc_auc_test = auc(fpr_test, tpr_test)
    
    precision_test, recall_test, _ = precision_recall_curve(y_test, X_test.flatten())
    pr_auc_test = auc(recall_test, precision_test)
    
    return {
        'train_results': {
            'youden': {'threshold': optimal_threshold_youden},
            'f1': {'threshold': optimal_threshold_f1},
            'balanced': {'threshold': optimal_threshold_balanced},
            'roc_auc': roc_auc_train,
            'pr_auc': pr_auc_train
        },
        'test_results': test_results,
        'cv_results': cv_results,
        'test_roc_auc': roc_auc_test,
        'test_pr_auc': pr_auc_test,
        'train_size': len(X_train),
        'test_size': len(X_test),
        'feature_name': feature_name
    }

def evaluate_all_features(df):
    """Evaluate threshold optimization for all E2P2 features"""
    feature_columns = ['num_distinct_enzyme_classes', 'num_distinct_enzyme_subfamilies', 'total_ec_numbers']
    all_results = {}
    
    print("\n" + "="*80)
    print("E2P2 FEATURE THRESHOLD OPTIMIZATION")
    print("="*80)
    
    for feature in feature_columns:
        if feature not in df.columns:
            print(f"Skipping {feature} - not available in data")
            continue
            
        print(f"\n{'='*70}")
        print(f"ANALYZING {feature.upper()}")
        print("="*70)
        
        # Basic threshold optimization
        basic_results = find_optimal_threshold_single_feature(df, feature)
        if basic_results:
            print(f"\nBasic threshold optimization:")
            print(f"  ROC AUC: {basic_results['roc_auc']:.3f}")
            print(f"  PR AUC: {basic_results['pr_auc']:.3f}")
            print(f"  Data: {basic_results['n_positives']} positives, {basic_results['n_negatives']} negatives")
            
            print(f"\nOptimal thresholds by method:")
            for method in ['youden', 'f1', 'balanced']:
                if method in basic_results:
                    m = basic_results[method]
                    print(f"  {method.upper()}: threshold={m['threshold']:.4f}, "
                          f"F1={m['f1_score']:.3f}, MCC={m['mcc']:.3f}, "
                          f"Precision={m['precision']:.3f}, Recall={m['recall']:.3f}")
        
        # Robust validation
        robust_results = robust_threshold_optimization_with_validation(df, feature)
        if robust_results:
            print(f"\nRobust validation results:")
            print(f"  Test ROC AUC: {robust_results['test_roc_auc']:.4f}")
            print(f"  Test PR AUC: {robust_results['test_pr_auc']:.4f}")
            
            print(f"\nCross-validation stability:")
            cv_results = robust_results['cv_results']
            print(f"  CV ROC AUC: {cv_results['roc_auc']['mean']:.4f} ± {cv_results['roc_auc']['std']:.4f}")
            print(f"  CV F1 Score: {cv_results['f1']['mean']:.4f} ± {cv_results['f1']['std']:.4f}")
            
            print(f"\nTest set performance by threshold:")
            for method_name in ['youden', 'f1', 'balanced']:
                if method_name in robust_results['test_results']:
                    result = robust_results['test_results'][method_name]
                    print(f"\n{method_name.upper()} threshold ({result['threshold']:.4f}):")
                    print(f"  Test Accuracy: {result['accuracy']:.4f}")
                    print(f"  Test Precision: {result['precision']:.4f}")
                    print(f"  Test Recall: {result['recall']:.4f}")
                    print(f"  Test F1 Score: {result['f1_score']:.4f}")
            
            # Generalizability assessment
            test_roc_auc = robust_results['test_roc_auc']
            cv_roc_auc_mean = cv_results['roc_auc']['mean']
            cv_roc_auc_std = cv_results['roc_auc']['std']
            
            if cv_roc_auc_std < 0.05:
                stability = "Highly stable"
            elif cv_roc_auc_std < 0.1:
                stability = "Stable"
            else:
                stability = "Variable"
            
            if abs(test_roc_auc - cv_roc_auc_mean) < 0.05:
                generalization = "Excellent generalization"
            elif abs(test_roc_auc - cv_roc_auc_mean) < 0.1:
                generalization = "Good generalization"
            else:
                generalization = "Potential overfitting"
            
            print(f"\nModel Assessment:")
            print(f"  Stability: {stability}")
            print(f"  Generalization: {generalization}")
            
            # Production recommendations
            if test_roc_auc >= 0.8 and cv_roc_auc_std < 0.1:
                recommendation = "Highly recommended for production use"
            elif test_roc_auc >= 0.7 and cv_roc_auc_std < 0.15:
                recommendation = "Suitable for use with monitoring"
            else:
                recommendation = "Requires improvement before production use"
            
            print(f"  Recommendation: {recommendation}")
            
            if 'f1' in robust_results['test_results']:
                best_threshold = robust_results['test_results']['f1']['threshold']
                best_f1 = robust_results['test_results']['f1']['f1_score']
                print(f"\nRecommended Production Threshold: {best_threshold:.4f}")
                print(f"  Expected F1 Score: {best_f1:.4f}")
                print(f"  Expected Precision: {robust_results['test_results']['f1']['precision']:.4f}")
                print(f"  Expected Recall: {robust_results['test_results']['f1']['recall']:.4f}")
        
        all_results[feature] = {
            'basic': basic_results,
            'robust': robust_results
        }
    
    return all_results

def main():
    """Main analysis function"""
    # Load E2P2 features
    feature_file = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/feature_extraction_results/e2p2_feature_extraction_detailed.csv"
    
    print("E2P2 Statistical Analysis")
    print("=" * 50)
    
    df = load_e2p2_features(feature_file)
    if df is None:
        return
    
    # Create binary labels
    df = create_binary_labels(df)
    
    # Calculate feature statistics
    feature_stats = calculate_feature_statistics(df)
    
    print(f"\nDataset Overview:")
    print(f"  Total clusters: {feature_stats['total_clusters']}")
    print(f"  Positive clusters (BGC/MGC): {feature_stats['positive_clusters']}")
    print(f"  Negative clusters (RANDOM): {feature_stats['negative_clusters']}")
    
    # Group comparison analysis
    group_results = group_comparison_analysis(df)
    
    # Threshold optimization
    threshold_results = evaluate_all_features(df)
    
    # Final summary
    print(f"\n{'='*80}")
    print("FINAL SUMMARY - BEST PERFORMING FEATURES")
    print("="*80)
    
    if threshold_results:
        print(f"\nComparison of validated performance across features:")
        print("-" * 60)
        
        feature_performance = []
        for feature, results in threshold_results.items():
            if results['robust'] and 'test_roc_auc' in results['robust']:
                feature_performance.append({
                    'feature': feature,
                    'test_roc_auc': results['robust']['test_roc_auc'],
                    'test_pr_auc': results['robust']['test_pr_auc'],
                    'cv_roc_auc': results['robust']['cv_results']['roc_auc']['mean'],
                    'cv_stability': results['robust']['cv_results']['roc_auc']['std']
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
    
    return df, feature_stats, group_results, threshold_results

if __name__ == "__main__":
    df, feature_stats, group_results, threshold_results = main()
