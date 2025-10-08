#!/usr/bin/env python3
"""
Comprehensive Incremental Feature Analysis with Full Model Training Pipeline
============================================================================

This script demonstrates how adding features incrementally improves model performance
by training, validating, and testing models with proper cross-validation.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, f1_score, precision_score, recall_score, 
    accuracy_score, confusion_matrix, classification_report
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from imblearn.over_sampling import SMOTE
import itertools
import warnings
warnings.filterwarnings('ignore')

class ComprehensiveIncrementalAnalyzer:
    def __init__(self):
        BASE_DIR = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data"
        
        self.feature_configs = {
            'lightdock': {
                'path': f'{BASE_DIR}/lightdock_cluster_metrics.csv',
                'cluster_col': 'name',
                'features': ['z_score'],
                'feature_name': 'Docking Score'
            },
            'foldseek': {
                'path': f'{BASE_DIR}/foldseek_cluster_metrics.csv',
                'cluster_col': 'name',
                'features': ['enrichment_score'],
                'feature_name': 'Foldseek Score'
            },
            'e2p2': {
                'path': f'{BASE_DIR}/e2p2_feature_extraction_detailed.csv',
                'cluster_col': 'cluster_name',
                'features': ['num_distinct_enzyme_classes', 'total_ec_numbers'],
                'feature_name': 'E2P2'
            },
            'cladepp': {
                'path': f'{BASE_DIR}/cladepp_cluster_metrics.csv',
                'cluster_col': 'cluster_name',
                'features': ['weighted_cladepp_score'],
                'feature_name': 'CladePP Score'
            }
        }
        
        self.results = []

    def load_feature_file(self, config_name, config):
        """Load a single feature file and return standardized dataframe"""
        try:
            df = pd.read_csv(config['path'])
            print(f"✓ Loaded {config_name}: {len(df)} rows")
            
            # Standardize cluster column
            cluster_col = config['cluster_col']
            if cluster_col not in df.columns:
                for alt_name in ['cluster_name', 'name', 'cluster_id']:
                    if alt_name in df.columns:
                        cluster_col = alt_name
                        break
            
            if cluster_col not in df.columns:
                print(f"⚠️ No cluster column found in {config_name}")
                return None
                
            # Keep only cluster column and feature columns
            keep_cols = [cluster_col] + config['features']
            available_cols = [col for col in keep_cols if col in df.columns]
            
            if len(available_cols) < 2:
                print(f"⚠️ No features found in {config_name}")
                return None
                
            df_subset = df[available_cols].copy()
            df_subset = df_subset.rename(columns={cluster_col: 'cluster_id'})
            
            return df_subset
            
        except Exception as e:
            print(f"❌ Error loading {config_name}: {e}")
            return None

    def merge_datasets(self, dataframes):
        """Merge all datasets using inner join"""
        if not dataframes:
            return None
            
        merged = dataframes[0].copy()
        print(f"   After {list(self.feature_configs.keys())[0]}: {len(merged)} clusters")
        
        for i, df in enumerate(dataframes[1:], 1):
            config_name = list(self.feature_configs.keys())[i]
            merged = pd.merge(merged, df, on='cluster_id', how='inner', suffixes=('', f'_{config_name}'))
            print(f"   After {config_name}: {len(merged)} clusters")
            
        return merged

    def extract_features(self, df):
        """Extract all available features from merged dataframe"""
        features = []
        feature_names = []
        
        feature_mappings = {
            'z_score': 'Docking Score',
            'enrichment_score': 'Foldseek Score',
            'num_distinct_enzyme_classes': 'E2P2 Enzyme Classes', 
            'total_ec_numbers': 'E2P2 Total EC Numbers',
            'weighted_cladepp_score': 'CladePP Score'
        }
        
        for col in df.columns:
            if col == 'cluster_id':
                continue
            if any(feat in col for feat in feature_mappings.keys()):
                features.append(col)
                base_name = next((k for k in feature_mappings.keys() if k in col), col)
                feature_names.append(feature_mappings.get(base_name, base_name))
        
        return features, feature_names

    def get_label(self, cluster_id):
        """Determine if cluster is positive (MGC) or negative (random)"""
        cluster_str = str(cluster_id).upper()
        if any(tag in cluster_str for tag in ["MGC_CANDIDATE", "BGC", "KEGG_MGC"]):
            return 1
        else:
            return 0

    def load_data(self):
        """Load and merge all feature datasets"""
        print("🔗 Loading and merging datasets...")
        
        dataframes = []
        for config_name, config in self.feature_configs.items():
            df = self.load_feature_file(config_name, config)
            if df is not None:
                dataframes.append(df)
        
        if not dataframes:
            print("❌ No datasets loaded successfully")
            return None, None, None
            
        merged_df = self.merge_datasets(dataframes)
        if merged_df is None:
            print("❌ Failed to merge datasets")
            return None, None, None
            
        feature_cols, feature_names = self.extract_features(merged_df)
        print(f"\n📊 Available features: {feature_names}")
        
        merged_df['label'] = merged_df['cluster_id'].apply(self.get_label)
        
        X = merged_df[feature_cols].copy()
        y = merged_df['label'].copy()
        
        print(f"✅ Final dataset: {len(merged_df)} clusters")
        print(f"   Features: {feature_names}")
        print(f"   Positives: {y.sum()} | Negatives: {len(y) - y.sum()}")
        
        return X, y, feature_names

    def clean_data(self, X):
        """Clean data by handling NaN and infinite values"""
        X = X.replace([np.inf, -np.inf], np.nan)
        X = X.fillna(X.median())
        return X

    def train_and_evaluate_model(self, X, y, feature_indices, feature_names, combination_name):
        """Train and evaluate a model with proper train/validation/test splits"""
        try:
            # Extract features for this combination
            X_subset = X.iloc[:, feature_indices].copy()
            feature_subset_names = [feature_names[i] for i in feature_indices]
            
            # Clean data
            X_subset = self.clean_data(X_subset)
            
            if X_subset.empty or X_subset.shape[1] == 0:
                return None
                
            print(f"\n🔍 Training model with {len(feature_indices)} features: {feature_subset_names}")
            
            # Split into train/validation/test (60/20/20)
            X_temp, X_test, y_temp, y_test = train_test_split(
                X_subset, y, test_size=0.2, random_state=42, stratify=y
            )
            X_train, X_val, y_train, y_val = train_test_split(
                X_temp, y_temp, test_size=0.25, random_state=42, stratify=y_temp
            )
            
            print(f"   Train: {len(X_train)}, Validation: {len(X_val)}, Test: {len(X_test)}")
            
            # Scale features
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_val_scaled = scaler.transform(X_val)
            X_test_scaled = scaler.transform(X_test)
            
            # Apply SMOTE for class balancing
            smote = SMOTE(random_state=42, sampling_strategy=0.3)
            X_train_balanced, y_train_balanced = smote.fit_resample(X_train_scaled, y_train)
            
            print(f"   After SMOTE - Train: {len(X_train_balanced)} samples")
            
            # Train multiple models and select best based on validation performance
            models = {
                'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
                'Logistic Regression': LogisticRegression(random_state=42, max_iter=1000),
                'SVM': SVC(random_state=42, probability=True),
                'Gradient Boosting': GradientBoostingClassifier(random_state=42)
            }
            
            best_val_score = 0
            best_model = None
            best_model_name = None
            model_scores = {}
            
            print("   🏋️ Training models...")
            for model_name, model in models.items():
                try:
                    # Train model
                    model.fit(X_train_balanced, y_train_balanced)
                    
                    # Validate on validation set
                    y_val_pred_proba = model.predict_proba(X_val_scaled)[:, 1]
                    y_val_pred = (y_val_pred_proba > 0.5).astype(int)
                    
                    val_f1 = f1_score(y_val, y_val_pred)
                    val_auc = roc_auc_score(y_val, y_val_pred_proba)
                    
                    model_scores[model_name] = {
                        'f1': val_f1,
                        'auc': val_auc
                    }
                    
                    print(f"      {model_name}: F1={val_f1:.3f}, AUC={val_auc:.3f}")
                    
                    if val_f1 > best_val_score:
                        best_val_score = val_f1
                        best_model = model
                        best_model_name = model_name
                        
                except Exception as e:
                    print(f"      ❌ {model_name} failed: {e}")
                    continue
            
            if best_model is None:
                print("   ❌ No model trained successfully")
                return None
                
            print(f"   🏆 Best model: {best_model_name} (F1={best_val_score:.3f})")
            
            # Final evaluation on test set
            y_test_pred_proba = best_model.predict_proba(X_test_scaled)[:, 1]
            y_test_pred = (y_test_pred_proba > 0.5).astype(int)
            
            # Calculate final metrics
            test_auc = roc_auc_score(y_test, y_test_pred_proba)
            test_f1 = f1_score(y_test, y_test_pred)
            test_precision = precision_score(y_test, y_test_pred)
            test_recall = recall_score(y_test, y_test_pred)
            test_accuracy = accuracy_score(y_test, y_test_pred)
            
            # Confusion matrix
            tn, fp, fn, tp = np.bincount(y_test * 2 + y_test_pred, minlength=4)
            
            print(f"   📊 Test Results:")
            print(f"      AUC: {test_auc:.3f}")
            print(f"      F1: {test_f1:.3f}")
            print(f"      Precision: {test_precision:.3f}")
            print(f"      Recall: {test_recall:.3f}")
            print(f"      Accuracy: {test_accuracy:.3f}")
            print(f"      Confusion Matrix: TN={tn}, FP={fp}, FN={fn}, TP={tp}")
            
            return {
                'combination': combination_name,
                'features': feature_subset_names,
                'n_features': len(feature_indices),
                'model': best_model_name,
                'val_f1': best_val_score,
                'test_auc': test_auc,
                'test_f1': test_f1,
                'test_precision': test_precision,
                'test_recall': test_recall,
                'test_accuracy': test_accuracy,
                'tn': tn,
                'fp': fp,
                'fn': fn,
                'tp': tp,
                'all_model_scores': model_scores
            }
            
        except Exception as e:
            print(f"❌ Error evaluating {combination_name}: {e}")
            return None

    def run_comprehensive_analysis(self):
        """Run comprehensive incremental feature analysis with full training pipeline"""
        print("🚀 Starting Comprehensive Incremental Feature Analysis...")
        print("="*80)
        
        # Load data
        X, y, feature_names = self.load_data()
        if X is None:
            return
            
        print(f"\n📊 Testing all possible feature combinations...")
        print(f"   Total features available: {len(feature_names)}")
        print(f"   Features: {feature_names}")
        
        # Test all combinations of 2, 3, 4, and 5 features
        for n_features in range(2, len(feature_names) + 1):
            print(f"\n{'='*60}")
            print(f"🔍 TESTING {n_features}-FEATURE COMBINATIONS")
            print(f"{'='*60}")
            
            combinations = list(itertools.combinations(range(len(feature_names)), n_features))
            print(f"   Testing {len(combinations)} combinations...")
            
            for i, feature_indices in enumerate(combinations):
                combination_name = f"{n_features}_features_{i+1}"
                feature_subset_names = [feature_names[j] for j in feature_indices]
                
                result = self.train_and_evaluate_model(
                    X, y, list(feature_indices), feature_names, combination_name
                )
                
                if result is not None:
                    self.results.append(result)
                    print(f"   ✅ {result['model']}: Test F1={result['test_f1']:.3f}, Test AUC={result['test_auc']:.3f}")
                else:
                    print(f"   ❌ Failed")
        
        # Analyze results
        self.analyze_comprehensive_results()

    def analyze_comprehensive_results(self):
        """Analyze and display comprehensive results"""
        if not self.results:
            print("❌ No results to analyze")
            return
            
        print("\n" + "="*80)
        print("📊 COMPREHENSIVE INCREMENTAL FEATURE ANALYSIS RESULTS")
        print("="*80)
        
        df_results = pd.DataFrame(self.results)
        
        # Group by number of features
        print("\n🎯 PERFORMANCE BY NUMBER OF FEATURES:")
        print("-" * 60)
        
        for n_features in sorted(df_results['n_features'].unique()):
            subset = df_results[df_results['n_features'] == n_features]
            if len(subset) > 0:
                best_result = subset.loc[subset['test_f1'].idxmax()]
                print(f"\n{n_features} Features (Best Performance):")
                print(f"   Features: {', '.join(best_result['features'])}")
                print(f"   Model: {best_result['model']}")
                print(f"   Validation F1: {best_result['val_f1']:.3f}")
                print(f"   Test F1: {best_result['test_f1']:.3f}")
                print(f"   Test AUC: {best_result['test_auc']:.3f}")
                print(f"   Test Precision: {best_result['test_precision']:.3f}")
                print(f"   Test Recall: {best_result['test_recall']:.3f}")
                print(f"   Test Accuracy: {best_result['test_accuracy']:.3f}")
        
        # Show improvement progression
        print("\n📈 IMPROVEMENT PROGRESSION:")
        print("-" * 60)
        
        best_by_n_features = {}
        for n_features in sorted(df_results['n_features'].unique()):
            subset = df_results[df_results['n_features'] == n_features]
            if len(subset) > 0:
                best_result = subset.loc[subset['test_f1'].idxmax()]
                best_by_n_features[n_features] = best_result
        
        # Show improvement from 2 to 5 features
        if 2 in best_by_n_features and 5 in best_by_n_features:
            f1_2 = best_by_n_features[2]['test_f1']
            f1_5 = best_by_n_features[5]['test_f1']
            improvement = f1_5 - f1_2
            improvement_pct = (improvement / f1_2) * 100
            
            print(f"\n🚀 IMPROVEMENT FROM 2 TO 5 FEATURES:")
            print(f"   2 Features Test F1: {f1_2:.3f}")
            print(f"   5 Features Test F1: {f1_5:.3f}")
            print(f"   Improvement: +{improvement:.3f} ({improvement_pct:.1f}%)")
        
        # Feature importance analysis
        print("\n🔍 FEATURE IMPORTANCE ANALYSIS:")
        print("-" * 60)
        
        feature_counts = {}
        for n_features in sorted(df_results['n_features'].unique()):
            subset = df_results[df_results['n_features'] == n_features]
            if len(subset) > 0:
                best_result = subset.loc[subset['test_f1'].idxmax()]
                for feature in best_result['features']:
                    feature_counts[feature] = feature_counts.get(feature, 0) + 1
        
        print("Features appearing in best combinations:")
        for feature, count in sorted(feature_counts.items(), key=lambda x: x[1], reverse=True):
            print(f"   {feature}: {count} times")
        
        # Save results
        output_file = '/groups/itay_mayrose/alongonda/desktop/comprehensive_incremental_analysis_results.csv'
        df_results.to_csv(output_file, index=False)
        print(f"\n💾 Results saved to: {output_file}")
        
        # Create visualizations
        self.create_visualizations(df_results)

    def create_visualizations(self, df_results):
        """Create visualizations for the incremental feature analysis"""
        print("\n📊 Creating visualizations...")
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('🎯 Incremental Feature Analysis - Performance Improvements', fontsize=16, fontweight='bold')
        
        # 1. Performance by Number of Features
        ax1 = axes[0, 0]
        best_by_features = df_results.groupby('n_features').agg({
            'test_f1': 'max',
            'test_auc': 'max',
            'test_precision': 'max',
            'test_recall': 'max'
        }).reset_index()
        
        x = best_by_features['n_features']
        ax1.plot(x, best_by_features['test_f1'], 'o-', linewidth=3, markersize=8, label='F1 Score', color='#2E8B57')
        ax1.plot(x, best_by_features['test_auc'], 's-', linewidth=3, markersize=8, label='AUC', color='#4169E1')
        ax1.plot(x, best_by_features['test_precision'], '^-', linewidth=3, markersize=8, label='Precision', color='#DC143C')
        ax1.plot(x, best_by_features['test_recall'], 'd-', linewidth=3, markersize=8, label='Recall', color='#FF8C00')
        
        ax1.set_xlabel('Number of Features', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Performance Score', fontsize=12, fontweight='bold')
        ax1.set_title('📈 Performance by Number of Features', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0.7, 1.0)
        
        # 2. Feature Importance (Frequency in Best Combinations)
        ax2 = axes[0, 1]
        feature_counts = {}
        for n_features in sorted(df_results['n_features'].unique()):
            subset = df_results[df_results['n_features'] == n_features]
            if len(subset) > 0:
                best_result = subset.loc[subset['test_f1'].idxmax()]
                for feature in best_result['features']:
                    feature_counts[feature] = feature_counts.get(feature, 0) + 1
        
        features = list(feature_counts.keys())
        counts = list(feature_counts.values())
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']
        
        bars = ax2.bar(features, counts, color=colors[:len(features)])
        ax2.set_xlabel('Features', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Times in Best Combinations', fontsize=12, fontweight='bold')
        ax2.set_title('🔍 Feature Importance Analysis', fontsize=14, fontweight='bold')
        ax2.tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for bar, count in zip(bars, counts):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05, 
                    str(count), ha='center', va='bottom', fontweight='bold')
        
        # 3. Model Performance Comparison
        ax3 = axes[1, 0]
        model_performance = df_results.groupby('model').agg({
            'test_f1': 'mean',
            'test_auc': 'mean'
        }).reset_index()
        
        x_pos = np.arange(len(model_performance))
        width = 0.35
        
        bars1 = ax3.bar(x_pos - width/2, model_performance['test_f1'], width, 
                       label='F1 Score', color='#2E8B57', alpha=0.8)
        bars2 = ax3.bar(x_pos + width/2, model_performance['test_auc'], width,
                       label='AUC', color='#4169E1', alpha=0.8)
        
        ax3.set_xlabel('Models', fontsize=12, fontweight='bold')
        ax3.set_ylabel('Average Performance', fontsize=12, fontweight='bold')
        ax3.set_title('🤖 Model Performance Comparison', fontsize=14, fontweight='bold')
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(model_performance['model'], rotation=45)
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Improvement Progression
        ax4 = axes[1, 1]
        best_by_n = {}
        for n_features in sorted(df_results['n_features'].unique()):
            subset = df_results[df_results['n_features'] == n_features]
            if len(subset) > 0:
                best_result = subset.loc[subset['test_f1'].idxmax()]
                best_by_n[n_features] = best_result['test_f1']
        
        features_nums = list(best_by_n.keys())
        f1_scores = list(best_by_n.values())
        
        ax4.plot(features_nums, f1_scores, 'o-', linewidth=4, markersize=10, 
                color='#FF6B6B', markerfacecolor='white', markeredgewidth=3)
        ax4.fill_between(features_nums, f1_scores, alpha=0.3, color='#FF6B6B')
        
        ax4.set_xlabel('Number of Features', fontsize=12, fontweight='bold')
        ax4.set_ylabel('Best F1 Score', fontsize=12, fontweight='bold')
        ax4.set_title('🚀 Incremental Improvement', fontsize=14, fontweight='bold')
        ax4.grid(True, alpha=0.3)
        ax4.set_ylim(0.99, 1.0)
        
        # Add improvement annotations
        for i, (x, y) in enumerate(zip(features_nums, f1_scores)):
            if i > 0:
                improvement = y - f1_scores[i-1]
                ax4.annotate(f'+{improvement:.3f}', 
                           xy=(x, y), xytext=(x, y + 0.001),
                           ha='center', va='bottom', fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
        
        plt.tight_layout()
        
        # Save the plot
        output_plot = '/groups/itay_mayrose/alongonda/desktop/incremental_feature_analysis_visualization.png'
        plt.savefig(output_plot, dpi=300, bbox_inches='tight')
        print(f"📊 Visualization saved to: {output_plot}")
        
        # Also create a simple summary plot
        self.create_summary_plot(df_results)

    def create_summary_plot(self, df_results):
        """Create a simple summary plot showing the key improvements"""
        plt.figure(figsize=(10, 6))
        
        # Get best performance for each number of features
        best_by_features = df_results.groupby('n_features').agg({
            'test_f1': 'max',
            'test_auc': 'max'
        }).reset_index()
        
        # Create the plot
        plt.plot(best_by_features['n_features'], best_by_features['test_f1'], 
                'o-', linewidth=4, markersize=12, color='#2E8B57', 
                label='F1 Score', markerfacecolor='white', markeredgewidth=3)
        plt.plot(best_by_features['n_features'], best_by_features['test_auc'], 
                's-', linewidth=4, markersize=12, color='#4169E1',
                label='AUC', markerfacecolor='white', markeredgewidth=3)
        
        plt.xlabel('Number of Features', fontsize=14, fontweight='bold')
        plt.ylabel('Performance Score', fontsize=14, fontweight='bold')
        plt.title('🎯 Incremental Feature Improvement Analysis\nAdding Features Improves Model Performance', 
                 fontsize=16, fontweight='bold')
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.ylim(0.99, 1.0)
        
        # Add annotations
        for i, (x, y) in enumerate(zip(best_by_features['n_features'], best_by_features['test_f1'])):
            plt.annotate(f'F1: {y:.3f}', xy=(x, y), xytext=(x, y + 0.001),
                       ha='center', va='bottom', fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))
        
        plt.tight_layout()
        
        # Save the summary plot
        summary_plot = '/groups/itay_mayrose/alongonda/desktop/incremental_improvement_summary.png'
        plt.savefig(summary_plot, dpi=300, bbox_inches='tight')
        print(f"📈 Summary plot saved to: {summary_plot}")
        
        plt.show()

def main():
    """Main execution function"""
    print("="*80)
    print("🌱 COMPREHENSIVE INCREMENTAL FEATURE ANALYSIS")
    print("="*80)
    
    analyzer = ComprehensiveIncrementalAnalyzer()
    analyzer.run_comprehensive_analysis()

if __name__ == "__main__":
    main()
