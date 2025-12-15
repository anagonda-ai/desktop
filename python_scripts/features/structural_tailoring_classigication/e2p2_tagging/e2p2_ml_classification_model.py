#!/usr/bin/env python3
"""
E2P2 Feature ML Classification Model - Individual & Multi-feature
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, confusion_matrix, 
                            classification_report, matthews_corrcoef, f1_score)
from sklearn.model_selection import train_test_split, cross_validate
from sklearn.preprocessing import StandardScaler
import os
import warnings
warnings.filterwarnings('ignore')

E2P2_FEATURES = [
    'num_distinct_enzyme_classes',
    'num_distinct_enzyme_subclasses',
    'num_distinct_enzyme_families',
    'num_distinct_enzyme_subfamilies',
    'total_ec_numbers'
]

FEATURE_DESCRIPTIONS = {
    'num_distinct_enzyme_classes': 'Count of distinct first-level EC numbers (format: X.-.-.-). Represents diversity across the 6 main enzyme classes: 1=Oxidoreductases, 2=Transferases, 3=Hydrolases, 4=Lyases, 5=Isomerases, 6=Ligases. Low values suggest pathway-specific enzymatic function',
    'num_distinct_enzyme_subclasses': 'Count of distinct second-level EC numbers (format: X.Y.-.-). Represents diversity of enzyme subclasses within main categories. Measures functional breadth at intermediate specificity',
    'num_distinct_enzyme_families': 'Count of distinct third-level EC numbers (format: X.Y.Z.-). Represents diversity of enzyme families with specific reaction types. Higher granularity than subclasses, indicates mechanistic diversity',
    'num_distinct_enzyme_subfamilies': 'Count of distinct fourth-level EC numbers (format: X.Y.Z.W). Most specific enzyme classification level, representing exact substrate-specific enzymes. High diversity suggests metabolically promiscuous clusters; low diversity suggests pathway specialization',
    'total_ec_numbers': 'Total count of all E2P2 enzyme predictions including duplicates. High counts with low diversity indicate gene duplications and pathway amplification typical of MGCs. Positive predictor: more enzymes = higher MGC likelihood'
}

def train_model(df, features, feature_name=None, test_size=0.3, random_state=42):
    """Train classification model"""
    
    if feature_name:
        print(f"\n{'='*80}")
        print(f"INDIVIDUAL FEATURE MODEL: {feature_name}")
        print(f"Description: {FEATURE_DESCRIPTIONS[feature_name]}")
        print('='*80)
        features_to_use = [feature_name]
    else:
        print(f"\n{'='*80}")
        print(f"MULTI-FEATURE MODEL: {len(features)} features combined")
        print('='*80)
        features_to_use = features
    
    # Prepare data
    cluster_col = 'cluster_name'  # E2P2 uses 'cluster_name'
    data = df[[cluster_col] + features_to_use + ['label']].dropna()
    
    # Separate BGC clusters (validated MGCs) from training data
    bgc_mask = data[cluster_col].str.contains('BGC', case=False, na=False)
    bgc_data = data[bgc_mask].copy()
    train_data = data[~bgc_mask].copy()
    
    print(f"\nData: {len(data)} total samples")
    print(f"  Training set: {len(train_data)} samples ({train_data['label'].sum()} MGC, {(train_data['label']==0).sum()} Random)")
    print(f"  Validation set (BGC): {len(bgc_data)} samples ({bgc_data['label'].sum()} MGC, {(bgc_data['label']==0).sum()} Random)")
    
    # Prepare training/test split from non-BGC data
    X_train_full = train_data[features_to_use].values
    y_train_full = train_data['label'].values
    
    # Train/Test split (only on non-BGC data)
    X_train, X_test, y_train, y_test = train_test_split(
        X_train_full, y_train_full, test_size=test_size, stratify=y_train_full, random_state=random_state
    )
    
    # Prepare validation set (BGC clusters)
    X_val = bgc_data[features_to_use].values
    y_val = bgc_data['label'].values
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    X_val_scaled = scaler.transform(X_val)
    
    # Train model
    model = LogisticRegression(random_state=random_state, max_iter=1000)
    model.fit(X_train_scaled, y_train)
    
    # Get weights
    weights = model.coef_[0]
    intercept = model.intercept_[0]
    
    print(f"\n📊 LEARNED WEIGHTS:")
    if len(features_to_use) == 1:
        print(f"  Weight: {weights[0]:.6f}")
        print(f"  Intercept: {intercept:.6f}")
    else:
        print(f"  Intercept: {intercept:.6f}")
        for feat, weight in zip(features_to_use, weights):
            print(f"  {feat:45s}: {weight:+.6f}")
        
        # Feature importance
        feature_importance = sorted(zip(features_to_use, np.abs(weights)), key=lambda x: x[1], reverse=True)
        print(f"\n🏆 FEATURE IMPORTANCE:")
        for rank, (feat, importance) in enumerate(feature_importance, 1):
            print(f"  {rank}. {feat:45s}: {importance:.6f}")
    
    # Predictions on test set
    y_pred_proba_test = model.predict_proba(X_test_scaled)[:, 1]
    
    # ROC & PR on test set
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba_test)
    roc_auc = auc(fpr, tpr)
    precision, recall, pr_thresholds = precision_recall_curve(y_test, y_pred_proba_test)
    pr_auc = auc(recall, precision)
    
    # Find F1-optimal threshold on test set
    f1_scores = 2 * (precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-10)
    f1_optimal_idx = np.argmax(f1_scores)
    f1_optimal_threshold = pr_thresholds[f1_optimal_idx]
    
    print(f"\n🎯 TEST SET PERFORMANCE:")
    print(f"  PR AUC:  {pr_auc:.4f} (primary metric)")
    print(f"  ROC AUC: {roc_auc:.4f}")
    print(f"  F1-optimal threshold: {f1_optimal_threshold:.4f}")
    
    # Evaluate test set at F1-optimal threshold
    y_pred_test = (y_pred_proba_test >= f1_optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred_test).ravel()
    
    precision_val = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall_val = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = f1_score(y_test, y_pred_test)
    mcc = matthews_corrcoef(y_test, y_pred_test)
    
    print(f"\n📈 TEST SET METRICS (at F1-optimal threshold):")
    print(f"  Precision: {precision_val:.4f} | Recall: {recall_val:.4f} | F1: {f1:.4f} | MCC: {mcc:.4f}")
    print(f"  TP: {tp:4d} | FP: {fp:4d} | TN: {tn:5d} | FN: {fn:4d}")
    
    # Evaluate on validation set (BGC clusters)
    if len(X_val) > 0:
        y_pred_proba_val = model.predict_proba(X_val_scaled)[:, 1]
        
        # ROC & PR on validation set
        fpr_val, tpr_val, _ = roc_curve(y_val, y_pred_proba_val)
        roc_auc_val = auc(fpr_val, tpr_val)
        precision_val_curve, recall_val_curve, _ = precision_recall_curve(y_val, y_pred_proba_val)
        pr_auc_val = auc(recall_val_curve, precision_val_curve)
        
        # Evaluate validation set at F1-optimal threshold
        y_pred_val = (y_pred_proba_val >= f1_optimal_threshold).astype(int)
        tn_val, fp_val, fn_val, tp_val = confusion_matrix(y_val, y_pred_val).ravel()
        
        precision_val_metric = tp_val / (tp_val + fp_val) if (tp_val + fp_val) > 0 else 0
        recall_val_metric = tp_val / (tp_val + fn_val) if (tp_val + fn_val) > 0 else 0
        f1_val = f1_score(y_val, y_pred_val)
        mcc_val = matthews_corrcoef(y_val, y_pred_val)
        
        print(f"\n✅ VALIDATION SET (BGC) PERFORMANCE:")
        print(f"  PR AUC:  {pr_auc_val:.4f} (primary metric)")
        print(f"  ROC AUC: {roc_auc_val:.4f}")
        print(f"  📈 VALIDATION SET METRICS (at F1-optimal threshold):")
        print(f"  Precision: {precision_val_metric:.4f} | Recall: {recall_val_metric:.4f} | F1: {f1_val:.4f} | MCC: {mcc_val:.4f}")
        print(f"  TP: {tp_val:4d} | FP: {fp_val:4d} | TN: {tn_val:5d} | FN: {fn_val:4d}")
    else:
        roc_auc_val = np.nan
        pr_auc_val = np.nan
        precision_val_metric = np.nan
        recall_val_metric = np.nan
        f1_val = np.nan
        mcc_val = np.nan
    
    # Cross-validation with precision and recall (only on training data, excluding BGC)
    cv_scores = cross_validate(model, scaler.transform(X_train_full), y_train_full, cv=5, 
                               scoring=['roc_auc', 'average_precision', 'precision', 'recall', 'f1'], return_train_score=False)
    
    print(f"\n🔄 CROSS-VALIDATION (on training data):")
    print(f"  ROC AUC:   {cv_scores['test_roc_auc'].mean():.4f} ± {cv_scores['test_roc_auc'].std():.4f}")
    print(f"  PR AUC:    {pr_auc:.4f} (test set)")
    print(f"  Precision: {cv_scores['test_precision'].mean():.4f} ± {cv_scores['test_precision'].std():.4f}")
    print(f"  Recall:    {cv_scores['test_recall'].mean():.4f} ± {cv_scores['test_recall'].std():.4f}")
    print(f"  F1:        {cv_scores['test_f1'].mean():.4f} ± {cv_scores['test_f1'].std():.4f}")
    
    return {
        'features': features_to_use,
        'weights': weights,
        'intercept': intercept,
        'threshold': f1_optimal_threshold,
        'roc_auc': roc_auc,
        'pr_auc': pr_auc,
        'f1': f1,
        'mcc': mcc,
        'precision': precision_val,
        'recall': recall_val,
        'val_roc_auc': roc_auc_val if len(X_val) > 0 else np.nan,
        'val_pr_auc': pr_auc_val if len(X_val) > 0 else np.nan,
        'val_f1': f1_val if len(X_val) > 0 else np.nan,
        'val_mcc': mcc_val if len(X_val) > 0 else np.nan,
        'val_precision': precision_val_metric if len(X_val) > 0 else np.nan,
        'val_recall': recall_val_metric if len(X_val) > 0 else np.nan,
        'cv_roc_auc': cv_scores['test_roc_auc'].mean(),
        'cv_roc_auc_std': cv_scores['test_roc_auc'].std(),
        'cv_pr_auc': cv_scores['test_average_precision'].mean(),
        'cv_pr_auc_std': cv_scores['test_average_precision'].std(),
        'cv_precision': cv_scores['test_precision'].mean(),
        'cv_precision_std': cv_scores['test_precision'].std(),
        'cv_recall': cv_scores['test_recall'].mean(),
        'cv_recall_std': cv_scores['test_recall'].std(),
        'cv_f1': cv_scores['test_f1'].mean(),
        'cv_f1_std': cv_scores['test_f1'].std()
    }

def main():
    print("="*80)
    print("E2P2 ML CLASSIFICATION ANALYSIS")
    print("="*80)
    
    # Load data
    feature_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/e2p2_feature_extraction_detailed.csv"
    df = pd.read_csv(feature_file)
    
    # Create binary labels: 1 for MGC, 0 for RANDOM
    df['label'] = df['classification_label'].astype(int)
    
    print(f"\nLoaded {len(df)} clusters ({df['label'].sum()} MGC, {(df['label']==0).sum()} Random)")
    
    # Individual models
    print(f"\n{'='*80}")
    print("PART 1: INDIVIDUAL FEATURE MODELS")
    print("="*80)
    
    individual_results = []
    for feature in E2P2_FEATURES:
        if feature in df.columns:
            result = train_model(df, E2P2_FEATURES, feature_name=feature)
            individual_results.append(result)
        else:
            print(f"\n⚠️  Feature '{feature}' not found in data, skipping...")
    
    # Multi-feature model
    print(f"\n{'='*80}")
    print("PART 2: MULTI-FEATURE MODEL")
    print("="*80)
    
    # Use only available features
    available_features = [f for f in E2P2_FEATURES if f in df.columns]
    if len(available_features) >= 2:
        multi_result = train_model(df, available_features)
    else:
        print(f"\n⚠️  Not enough features available for multi-feature model")
        multi_result = None
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("="*80)
    
    if individual_results:
        print(f"\nIndividual Features (ranked by PR AUC - primary metric):")
        sorted_results = sorted(individual_results, key=lambda x: x.get('pr_auc', 0) if pd.notna(x.get('pr_auc')) else 0, reverse=True)
        for i, r in enumerate(sorted_results, 1):
            pr_auc_val = r.get('pr_auc', np.nan)
            precision_val = r.get('precision', np.nan)
            recall_val = r.get('recall', np.nan)
            pr_auc_str = f"PR AUC: {pr_auc_val:.4f}" if pd.notna(pr_auc_val) else "PR AUC: N/A"
            print(f"  {i}. {r['features'][0]:45s} {pr_auc_str} | Precision: {precision_val:.4f} | Recall: {recall_val:.4f} | F1: {r['f1']:.4f}")
    
    if multi_result:
        pr_auc_val = multi_result.get('pr_auc', np.nan)
        precision_val = multi_result.get('precision', np.nan)
        recall_val = multi_result.get('recall', np.nan)
        pr_auc_str = f"PR AUC: {pr_auc_val:.4f}" if pd.notna(pr_auc_val) else "PR AUC: N/A"
        print(f"\nMulti-Feature Model:")
        print(f"  {pr_auc_str} | Precision: {precision_val:.4f} | Recall: {recall_val:.4f} | F1: {multi_result['f1']:.4f} | MCC: {multi_result['mcc']:.4f}")
    
    # Save results
    if individual_results:
        summary = []
        for r in individual_results:
            summary.append({
                'model_type': 'individual',
                'feature': r['features'][0],
                'weight': r['weights'][0],
                'threshold': r['threshold'],
                'roc_auc': r['roc_auc'],
                'pr_auc': r.get('pr_auc', np.nan),
                'precision': r.get('precision', np.nan),
                'recall': r.get('recall', np.nan),
                'f1': r['f1'],
                'mcc': r['mcc'],
                'val_roc_auc': r.get('val_roc_auc', np.nan),
                'val_pr_auc': r.get('val_pr_auc', np.nan),
                'val_precision': r.get('val_precision', np.nan),
                'val_recall': r.get('val_recall', np.nan),
                'val_f1': r.get('val_f1', np.nan),
                'val_mcc': r.get('val_mcc', np.nan)
            })
        
        if multi_result:
            summary.append({
                'model_type': 'multi_feature',
                'feature': 'ALL_FEATURES_COMBINED',
                'weight': None,
                'threshold': multi_result['threshold'],
                'roc_auc': multi_result['roc_auc'],
                'pr_auc': multi_result.get('pr_auc', np.nan),
                'precision': multi_result.get('precision', np.nan),
                'recall': multi_result.get('recall', np.nan),
                'f1': multi_result['f1'],
                'mcc': multi_result['mcc'],
                'val_roc_auc': multi_result.get('val_roc_auc', np.nan),
                'val_pr_auc': multi_result.get('val_pr_auc', np.nan),
                'val_precision': multi_result.get('val_precision', np.nan),
                'val_recall': multi_result.get('val_recall', np.nan),
                'val_f1': multi_result.get('val_f1', np.nan),
                'val_mcc': multi_result.get('val_mcc', np.nan)
            })
        
        summary_df = pd.DataFrame(summary)
        output_dir = "/groups/itay_mayrose/alongonda/desktop/mibig_validate/e2p2"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, "e2p2_classification_models_summary.csv")
        summary_df.to_csv(output_file, index=False)
        print(f"\n💾 Saved: {output_file}")
        
        # Save multi-feature weights
        if multi_result:
            weights_df = pd.DataFrame({
                'feature': multi_result['features'],
                'weight': multi_result['weights'],
                'abs_importance': np.abs(multi_result['weights'])
            }).sort_values('abs_importance', ascending=False)
            
            weights_file = os.path.join(output_dir, "random_kegg_e2p2_multi_feature_weights.csv")
            weights_df.to_csv(weights_file, index=False)
            print(f"💾 Saved: {weights_file}")

if __name__ == "__main__":
    main()
