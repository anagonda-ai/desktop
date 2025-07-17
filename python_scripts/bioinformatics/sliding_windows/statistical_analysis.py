"""
Statistical Analysis Module for MGC Research

This module provides comprehensive statistical testing and analysis
for KEGG-MGC identification results and genome characteristics.
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import (
    pearsonr, spearmanr, kendalltau, mannwhitneyu, kruskal,
    shapiro, levene, chi2_contingency, fisher_exact
)
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import mcnemar
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import warnings
warnings.filterwarnings('ignore')


class MGCStatisticalAnalyzer:
    """Comprehensive statistical analysis for MGC research."""
    
    def __init__(self, genome_stats_file: str, mgc_results_file: str):
        """
        Initialize statistical analyzer.
        
        Args:
            genome_stats_file: Path to genome statistics CSV
            mgc_results_file: Path to MGC results CSV
        """
        self.genome_stats_file = genome_stats_file
        self.mgc_results_file = mgc_results_file
        self.output_dir = os.path.dirname(genome_stats_file)
        
        # Load data
        self.genome_stats = pd.read_csv(genome_stats_file)
        self.mgc_results = pd.read_csv(mgc_results_file)
        
        # Results storage
        self.test_results = {}
        self.correlation_results = {}
        self.regression_results = {}
        
    def test_normality(self, variables: list = None) -> dict:
        """
        Test normality of key variables using Shapiro-Wilk test.
        
        Args:
            variables: List of variables to test (default: all numeric)
            
        Returns:
            Dictionary with normality test results
        """
        if variables is None:
            variables = self.genome_stats.select_dtypes(include=[np.number]).columns
        
        normality_results = {}
        
        for var in variables:
            if var in self.genome_stats.columns:
                data = self.genome_stats[var].dropna()
                
                if len(data) >= 3:  # Minimum sample size for Shapiro-Wilk
                    statistic, p_value = shapiro(data)
                    
                    normality_results[var] = {
                        'statistic': statistic,
                        'p_value': p_value,
                        'is_normal': p_value > 0.05,
                        'n_samples': len(data)
                    }
        
        self.test_results['normality'] = normality_results
        return normality_results
    
    def comprehensive_correlation_analysis(self, method: str = 'all') -> dict:
        """
        Comprehensive correlation analysis with multiple methods.
        
        Args:
            method: Correlation method ('pearson', 'spearman', 'kendall', 'all')
            
        Returns:
            Dictionary with correlation results
        """
        # Select numeric variables
        numeric_vars = self.genome_stats.select_dtypes(include=[np.number]).columns
        target_var = 'mgc_count'
        
        correlation_results = {}
        
        methods = ['pearson', 'spearman', 'kendall'] if method == 'all' else [method]
        
        for corr_method in methods:
            method_results = {}
            
            for var in numeric_vars:
                if var != target_var and var in self.genome_stats.columns:
                    # Get clean data
                    clean_data = self.genome_stats[[var, target_var]].dropna()
                    
                    if len(clean_data) >= 3:
                        if corr_method == 'pearson':
                            r, p = pearsonr(clean_data[var], clean_data[target_var])
                        elif corr_method == 'spearman':
                            r, p = spearmanr(clean_data[var], clean_data[target_var])
                        elif corr_method == 'kendall':
                            r, p = kendalltau(clean_data[var], clean_data[target_var])
                        
                        method_results[var] = {
                            'correlation': r,
                            'p_value': p,
                            'n_samples': len(clean_data),
                            'effect_size': self._interpret_correlation(abs(r))
                        }
            
            correlation_results[corr_method] = method_results
        
        # Apply multiple testing correction
        if method == 'all':
            correlation_results = self._apply_multiple_testing_correction(correlation_results)
        
        self.correlation_results = correlation_results
        return correlation_results
    
    def _interpret_correlation(self, r: float) -> str:
        """Interpret correlation coefficient magnitude."""
        if r < 0.1:
            return 'negligible'
        elif r < 0.3:
            return 'small'
        elif r < 0.5:
            return 'medium'
        elif r < 0.7:
            return 'large'
        else:
            return 'very large'
    
    def _apply_multiple_testing_correction(self, correlation_results: dict) -> dict:
        """Apply Benjamini-Hochberg correction for multiple testing."""
        for method in correlation_results:
            if correlation_results[method]:
                p_values = [result['p_value'] for result in correlation_results[method].values()]
                var_names = list(correlation_results[method].keys())
                
                # Apply correction
                rejected, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')
                
                # Update results
                for i, var in enumerate(var_names):
                    correlation_results[method][var]['p_corrected'] = p_corrected[i]
                    correlation_results[method][var]['significant_corrected'] = rejected[i]
        
        return correlation_results
    
    def group_comparison_analysis(self) -> dict:
        """
        Compare MGC counts across different genome groups.
        
        Returns:
            Dictionary with group comparison results
        """
        comparison_results = {}
        
        # Create genome size groups
        if 'total_genes' in self.genome_stats.columns:
            # Tertile split
            tertiles = np.percentile(self.genome_stats['total_genes'].dropna(), [33.33, 66.67])
            
            def categorize_genome_size(genes):
                if genes <= tertiles[0]:
                    return 'Small'
                elif genes <= tertiles[1]:
                    return 'Medium'
                else:
                    return 'Large'
            
            self.genome_stats['genome_size_category'] = self.genome_stats['total_genes'].apply(categorize_genome_size)
            
            # Group MGC counts
            groups = {}
            for category in ['Small', 'Medium', 'Large']:
                mask = self.genome_stats['genome_size_category'] == category
                groups[category] = self.genome_stats[mask]['mgc_count'].dropna().values
            
            # Kruskal-Wallis test (non-parametric ANOVA)
            if all(len(group) >= 3 for group in groups.values()):
                h_stat, p_value = kruskal(*groups.values())
                
                comparison_results['genome_size_groups'] = {
                    'test': 'Kruskal-Wallis',
                    'statistic': h_stat,
                    'p_value': p_value,
                    'significant': p_value < 0.05,
                    'group_sizes': {cat: len(group) for cat, group in groups.items()},
                    'group_medians': {cat: np.median(group) for cat, group in groups.items()}
                }
                
                # Pairwise comparisons
                pairwise_results = {}
                categories = list(groups.keys())
                
                for i in range(len(categories)):
                    for j in range(i + 1, len(categories)):
                        cat1, cat2 = categories[i], categories[j]
                        u_stat, p_val = mannwhitneyu(groups[cat1], groups[cat2])
                        
                        pairwise_results[f'{cat1}_vs_{cat2}'] = {
                            'statistic': u_stat,
                            'p_value': p_val,
                            'significant': p_val < 0.05
                        }
                
                comparison_results['pairwise_comparisons'] = pairwise_results
        
        # KEGG annotation rate groups
        if 'kegg_annotation_rate' in self.genome_stats.columns:
            # Median split
            median_kegg = self.genome_stats['kegg_annotation_rate'].median()
            
            high_kegg = self.genome_stats[self.genome_stats['kegg_annotation_rate'] >= median_kegg]['mgc_count'].dropna()
            low_kegg = self.genome_stats[self.genome_stats['kegg_annotation_rate'] < median_kegg]['mgc_count'].dropna()
            
            if len(high_kegg) >= 3 and len(low_kegg) >= 3:
                u_stat, p_value = mannwhitneyu(high_kegg, low_kegg)
                
                comparison_results['kegg_annotation_groups'] = {
                    'test': 'Mann-Whitney U',
                    'statistic': u_stat,
                    'p_value': p_value,
                    'significant': p_value < 0.05,
                    'high_kegg_median': np.median(high_kegg),
                    'low_kegg_median': np.median(low_kegg),
                    'high_kegg_n': len(high_kegg),
                    'low_kegg_n': len(low_kegg)
                }
        
        self.test_results['group_comparisons'] = comparison_results
        return comparison_results
    
    def regression_analysis(self) -> dict:
        """
        Comprehensive regression analysis for MGC count prediction.
        
        Returns:
            Dictionary with regression results
        """
        regression_results = {}
        
        # Select predictor variables
        predictors = ['total_genes', 'gene_density', 'kegg_annotation_rate', 
                     'unique_kegg_pathways', 'total_kegg_annotations']
        
        available_predictors = [p for p in predictors if p in self.genome_stats.columns]
        
        if len(available_predictors) >= 2:
            # Prepare data
            data = self.genome_stats[available_predictors + ['mgc_count']].dropna()
            
            if len(data) >= 10:  # Minimum sample size for regression
                X = data[available_predictors]
                y = data['mgc_count']
                
                # Linear regression
                X_with_const = sm.add_constant(X)
                model = sm.OLS(y, X_with_const).fit()
                
                regression_results['linear_regression'] = {
                    'r_squared': model.rsquared,
                    'adj_r_squared': model.rsquared_adj,
                    'f_statistic': model.fvalue,
                    'f_p_value': model.f_pvalue,
                    'aic': model.aic,
                    'bic': model.bic,
                    'coefficients': model.params.to_dict(),
                    'p_values': model.pvalues.to_dict(),
                    'confidence_intervals': model.conf_int().to_dict(),
                    'n_observations': len(data)
                }
                
                # Random Forest regression for comparison
                rf_model = RandomForestRegressor(n_estimators=100, random_state=42)
                rf_scores = cross_val_score(rf_model, X, y, cv=5, scoring='r2')
                
                regression_results['random_forest'] = {
                    'mean_cv_r2': rf_scores.mean(),
                    'std_cv_r2': rf_scores.std(),
                    'cv_scores': rf_scores.tolist()
                }
                
                # Feature importance (Random Forest)
                rf_model.fit(X, y)
                feature_importance = dict(zip(available_predictors, rf_model.feature_importances_))
                regression_results['feature_importance'] = feature_importance
        
        self.regression_results = regression_results
        return regression_results
    
    def pathway_enrichment_analysis(self) -> dict:
        """
        Analyze pathway enrichment in high vs low MGC genomes.
        
        Returns:
            Dictionary with enrichment results
        """
        enrichment_results = {}
        
        # Define high vs low MGC groups
        mgc_median = self.genome_stats['mgc_count'].median()
        
        high_mgc_genomes = self.genome_stats[self.genome_stats['mgc_count'] >= mgc_median]['genome_name'].tolist()
        low_mgc_genomes = self.genome_stats[self.genome_stats['mgc_count'] < mgc_median]['genome_name'].tolist()
        
        # Clean genome names in MGC results for matching
        def clean_genome_name(source_file):
            name = os.path.basename(source_file).lower()
            suffixes = ['_updated_annotated.csv', '_annotated.csv', '.filtered.csv', '.filtered', '.csv']
            for suffix in suffixes:
                if name.endswith(suffix):
                    name = name[:-len(suffix)]
                    break
            return name.strip('._').replace('_', ' ').title()
        
        self.mgc_results['clean_genome_name'] = self.mgc_results['source_file'].apply(clean_genome_name)
        
        # Pathway counts by group
        high_mgc_pathways = self.mgc_results[self.mgc_results['clean_genome_name'].isin(high_mgc_genomes)]['pathway'].value_counts()
        low_mgc_pathways = self.mgc_results[self.mgc_results['clean_genome_name'].isin(low_mgc_genomes)]['pathway'].value_counts()
        
        # Get all pathways
        all_pathways = set(high_mgc_pathways.index) | set(low_mgc_pathways.index)
        
        # Fisher's exact test for each pathway
        pathway_tests = {}
        
        for pathway in all_pathways:
            high_with_pathway = high_mgc_pathways.get(pathway, 0)
            high_without_pathway = len(high_mgc_genomes) - high_with_pathway
            low_with_pathway = low_mgc_pathways.get(pathway, 0)
            low_without_pathway = len(low_mgc_genomes) - low_with_pathway
            
            # Create contingency table
            contingency_table = [
                [high_with_pathway, high_without_pathway],
                [low_with_pathway, low_without_pathway]
            ]
            
            # Fisher's exact test
            if sum(contingency_table[0]) > 0 and sum(contingency_table[1]) > 0:
                odds_ratio, p_value = fisher_exact(contingency_table)
                
                pathway_tests[pathway] = {
                    'odds_ratio': odds_ratio,
                    'p_value': p_value,
                    'high_mgc_count': high_with_pathway,
                    'low_mgc_count': low_with_pathway,
                    'enriched_in': 'high_mgc' if odds_ratio > 1 else 'low_mgc'
                }
        
        # Multiple testing correction
        if pathway_tests:
            p_values = [result['p_value'] for result in pathway_tests.values()]
            pathways = list(pathway_tests.keys())
            
            rejected, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')
            
            for i, pathway in enumerate(pathways):
                pathway_tests[pathway]['p_corrected'] = p_corrected[i]
                pathway_tests[pathway]['significant_corrected'] = rejected[i]
        
        enrichment_results['pathway_enrichment'] = pathway_tests
        enrichment_results['high_mgc_genomes'] = high_mgc_genomes
        enrichment_results['low_mgc_genomes'] = low_mgc_genomes
        
        return enrichment_results
    
    def save_statistical_results(self) -> str:
        """
        Save all statistical analysis results.
        
        Returns:
            Path to saved results file
        """
        # Combine all results
        all_results = {
            'normality_tests': self.test_results.get('normality', {}),
            'correlation_analysis': self.correlation_results,
            'group_comparisons': self.test_results.get('group_comparisons', {}),
            'regression_analysis': self.regression_results,
            'pathway_enrichment': self.test_results.get('pathway_enrichment', {})
        }
        
        # Save to file
        output_file = os.path.join(self.output_dir, 'statistical_analysis_results.txt')
        
        with open(output_file, 'w') as f:
            f.write("COMPREHENSIVE STATISTICAL ANALYSIS RESULTS\n")
            f.write("=" * 60 + "\n\n")
            
            # Normality tests
            f.write("NORMALITY TESTS (Shapiro-Wilk)\n")
            f.write("-" * 40 + "\n")
            for var, result in all_results['normality_tests'].items():
                f.write(f"{var}: W = {result['statistic']:.4f}, p = {result['p_value']:.4f}, ")
                f.write(f"Normal = {result['is_normal']}, n = {result['n_samples']}\n")
            f.write("\n")
            
            # Correlation analysis
            f.write("CORRELATION ANALYSIS\n")
            f.write("-" * 40 + "\n")
            for method, correlations in all_results['correlation_analysis'].items():
                f.write(f"\n{method.upper()} CORRELATIONS:\n")
                for var, result in correlations.items():
                    f.write(f"  {var}: r = {result['correlation']:.4f}, p = {result['p_value']:.4f}")
                    if 'p_corrected' in result:
                        f.write(f", p_corrected = {result['p_corrected']:.4f}")
                    f.write(f", effect_size = {result['effect_size']}, n = {result['n_samples']}\n")
            f.write("\n")
            
            # Group comparisons
            f.write("GROUP COMPARISONS\n")
            f.write("-" * 40 + "\n")
            for comparison, result in all_results['group_comparisons'].items():
                f.write(f"{comparison}: {result['test']}\n")
                f.write(f"  Statistic = {result['statistic']:.4f}, p = {result['p_value']:.4f}\n")
                f.write(f"  Significant = {result['significant']}\n")
                if 'group_medians' in result:
                    f.write(f"  Group medians: {result['group_medians']}\n")
                f.write("\n")
            
            # Regression analysis
            f.write("REGRESSION ANALYSIS\n")
            f.write("-" * 40 + "\n")
            if 'linear_regression' in all_results['regression_analysis']:
                lr = all_results['regression_analysis']['linear_regression']
                f.write(f"Linear Regression: R² = {lr['r_squared']:.4f}, adj R² = {lr['adj_r_squared']:.4f}\n")
                f.write(f"F-statistic = {lr['f_statistic']:.4f}, p = {lr['f_p_value']:.4f}\n")
                f.write(f"AIC = {lr['aic']:.2f}, BIC = {lr['bic']:.2f}\n")
            
            if 'random_forest' in all_results['regression_analysis']:
                rf = all_results['regression_analysis']['random_forest']
                f.write(f"Random Forest CV R² = {rf['mean_cv_r2']:.4f} ± {rf['std_cv_r2']:.4f}\n")
            f.write("\n")
        
        print(f"✅ Statistical analysis results saved to: {output_file}")
        return output_file
    
    def run_complete_analysis(self) -> dict:
        """
        Run complete statistical analysis pipeline.
        
        Returns:
            Dictionary with all analysis results
        """
        print("Running comprehensive statistical analysis...")
        
        # Run all analyses
        normality_results = self.test_normality()
        correlation_results = self.comprehensive_correlation_analysis()
        comparison_results = self.group_comparison_analysis()
        regression_results = self.regression_analysis()
        enrichment_results = self.pathway_enrichment_analysis()
        
        # Store enrichment results
        self.test_results['pathway_enrichment'] = enrichment_results
        
        # Save results
        results_file = self.save_statistical_results()
        
        print("✅ Complete statistical analysis finished.")
        
        return {
            'normality': normality_results,
            'correlations': correlation_results,
            'comparisons': comparison_results,
            'regression': regression_results,
            'enrichment': enrichment_results,
            'results_file': results_file
        }


def main():
    """Main function to run statistical analysis."""
    # Define paths
    base_dir = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3'
    genome_stats_file = os.path.join(base_dir, 'genome_statistics.csv')
    mgc_results_file = os.path.join(base_dir, 'potential_groups_w10_filtered.csv')
    
    # Check if files exist
    if not os.path.exists(genome_stats_file):
        print(f"Genome statistics file not found: {genome_stats_file}")
        print("Please run genome_statistics_extractor.py first.")
        return
    
    # Initialize analyzer
    analyzer = MGCStatisticalAnalyzer(genome_stats_file, mgc_results_file)
    
    # Run complete analysis
    results = analyzer.run_complete_analysis()
    
    return results


if __name__ == "__main__":
    main()