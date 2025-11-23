import os
import glob
import csv
import statistics
import concurrent.futures
import random

try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available. Statistical significance testing will be skipped.")

# Directory containing the merged CSV files
MERGED_DIR = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/aans_selected_1000/best_hits/merged"

def calculate_mean_std_for_file(csv_file, sample_size=100):
    """Calculate mean and standard deviation of bitscore column for a single CSV file.
    
    Returns statistics for both full dataset and a random sample.
    """
    bitscores = []
    filename = os.path.basename(csv_file)
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    bitscore = float(row['bitscore'])
                    bitscores.append(bitscore)
                except (ValueError, KeyError) as e:
                    # Skip rows with invalid bitscore values
                    continue
        
        if len(bitscores) == 0:
            return filename, None, None, 0, None, None, 0, None, None
        
        # Full dataset statistics
        full_mean = statistics.mean(bitscores)
        full_std = statistics.stdev(bitscores)  # Sample standard deviation
        full_count = len(bitscores)
        
        # Sample statistics (if we have enough data)
        sample_mean = None
        sample_std = None
        sample_count = 0
        sampled_bitscores = None
        t_statistic = None
        p_value = None
        
        if full_count >= sample_size:
            # Randomly sample 100 genes
            sampled_bitscores = random.sample(bitscores, sample_size)
            sample_mean = statistics.mean(sampled_bitscores)
            sample_std = statistics.stdev(sampled_bitscores)
            sample_count = sample_size
            
            # Perform statistical significance test (t-test)
            if SCIPY_AVAILABLE and len(bitscores) > 1 and len(sampled_bitscores) > 1:
                try:
                    # Use independent samples t-test
                    # Note: Since sample is a subset, we could also use one-sample t-test
                    # comparing sample mean to full mean, but independent t-test is more conservative
                    t_statistic, p_value = stats.ttest_ind(sampled_bitscores, bitscores, equal_var=False)
                except Exception as e:
                    # If test fails, continue without significance values
                    pass
        elif full_count > 0:
            # If we have less than sample_size but more than 0, use all available
            sample_mean = full_mean
            sample_std = full_std
            sample_count = full_count
        
        return filename, full_mean, full_std, full_count, sample_mean, sample_std, sample_count, t_statistic, p_value
    
    except Exception as e:
        print(f"Error processing {csv_file}: {e}")
        return filename, None, None, 0, None, None, 0, None, None

def main():
    # Set random seed for reproducibility
    random.seed(42)
    
    # Get all CSV files in the merged directory
    csv_files = glob.glob(os.path.join(MERGED_DIR, "*.csv"))
    csv_files.sort()
    
    print(f"Processing {len(csv_files)} CSV files with concurrent execution...")
    print("Calculating statistics for full dataset and 100-gene sample...")
    print()
    
    # Store results
    results = []
    
    # Process files concurrently using ThreadPoolExecutor (I/O bound operation)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit all tasks
        futures = {executor.submit(calculate_mean_std_for_file, csv_file): csv_file 
                   for csv_file in csv_files}
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(futures):
            filename, full_mean, full_std, full_count, sample_mean, sample_std, sample_count, t_statistic, p_value = future.result()
            
            if full_mean is not None:
                # Calculate differences
                mean_diff = None
                std_diff = None
                mean_diff_pct = None
                std_diff_pct = None
                is_significant = None
                significance_level = 0.05
                
                if sample_mean is not None:
                    mean_diff = sample_mean - full_mean
                    std_diff = sample_std - full_std
                    mean_diff_pct = (mean_diff / full_mean) * 100 if full_mean != 0 else 0
                    std_diff_pct = (std_diff / full_std) * 100 if full_std != 0 else 0
                    
                    # Determine if difference is statistically significant
                    if p_value is not None:
                        is_significant = p_value < significance_level
                
                results.append({
                    'filename': filename,
                    'full_mean': full_mean,
                    'full_std': full_std,
                    'full_count': full_count,
                    'sample_mean': sample_mean,
                    'sample_std': sample_std,
                    'sample_count': sample_count,
                    'mean_diff': mean_diff,
                    'std_diff': std_diff,
                    'mean_diff_pct': mean_diff_pct,
                    'std_diff_pct': std_diff_pct,
                    't_statistic': t_statistic,
                    'p_value': p_value,
                    'is_significant': is_significant
                })
                
                print(f"{filename}:")
                print(f"  Full dataset (n={full_count}):")
                print(f"    Mean bitscore: {full_mean:.2f}")
                print(f"    Std bitscore:  {full_std:.2f}")
                
                if sample_mean is not None:
                    print(f"  Sample (n={sample_count}):")
                    print(f"    Mean bitscore: {sample_mean:.2f}")
                    print(f"    Std bitscore:  {sample_std:.2f}")
                    print(f"  Comparison:")
                    print(f"    Mean difference: {mean_diff:+.2f} ({mean_diff_pct:+.2f}%)")
                    print(f"    Std difference:  {std_diff:+.2f} ({std_diff_pct:+.2f}%)")
                    
                    if p_value is not None:
                        significance_str = "SIGNIFICANT" if is_significant else "NOT SIGNIFICANT"
                        print(f"  Statistical test (t-test):")
                        print(f"    t-statistic: {t_statistic:.4f}")
                        print(f"    p-value: {p_value:.6f}")
                        print(f"    Significance (α={significance_level}): {significance_str}")
                    elif not SCIPY_AVAILABLE:
                        print(f"  Statistical test: scipy not available")
                else:
                    print(f"  Sample: Not enough data (need at least 100 genes)")
                print()
            else:
                print(f"{filename}: ERROR - Could not calculate statistics")
                print()
    
    # Sort results by filename for consistent output
    results.sort(key=lambda x: x['filename'])
    
    # Save results to summary files
    summary_file = os.path.join(MERGED_DIR, "bitscore_statistics_summary.csv")
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'filename', 'full_mean_bitscore', 'full_std_bitscore', 'full_num_rows',
            'sample_mean_bitscore', 'sample_std_bitscore', 'sample_num_rows',
            'mean_diff', 'std_diff', 'mean_diff_pct', 'std_diff_pct',
            't_statistic', 'p_value', 'is_significant'
        ])
        for res in results:
            writer.writerow([
                res['filename'],
                f"{res['full_mean']:.4f}" if res['full_mean'] else '',
                f"{res['full_std']:.4f}" if res['full_std'] else '',
                res['full_count'],
                f"{res['sample_mean']:.4f}" if res['sample_mean'] else '',
                f"{res['sample_std']:.4f}" if res['sample_std'] else '',
                res['sample_count'],
                f"{res['mean_diff']:.4f}" if res['mean_diff'] is not None else '',
                f"{res['std_diff']:.4f}" if res['std_diff'] is not None else '',
                f"{res['mean_diff_pct']:.4f}" if res['mean_diff_pct'] is not None else '',
                f"{res['std_diff_pct']:.4f}" if res['std_diff_pct'] is not None else '',
                f"{res['t_statistic']:.6f}" if res['t_statistic'] is not None else '',
                f"{res['p_value']:.6f}" if res['p_value'] is not None else '',
                'TRUE' if res['is_significant'] is True else ('FALSE' if res['is_significant'] is False else '')
            ])
    
    print(f"Summary saved to: {summary_file}")
    
    # Print overall comparison statistics
    print("\n" + "="*80)
    print("OVERALL COMPARISON SUMMARY")
    print("="*80)
    
    valid_comparisons = [r for r in results if r['sample_mean'] is not None]
    if valid_comparisons:
        mean_diffs = [r['mean_diff_pct'] for r in valid_comparisons]
        std_diffs = [r['std_diff_pct'] for r in valid_comparisons]
        
        # Statistical significance summary
        significant_tests = [r for r in valid_comparisons if r['is_significant'] is True]
        non_significant_tests = [r for r in valid_comparisons if r['is_significant'] is False]
        tests_with_pvalues = [r for r in valid_comparisons if r['p_value'] is not None]
        
        print(f"Files with valid comparisons: {len(valid_comparisons)}")
        print(f"\nMean bitscore difference (sample - full):")
        print(f"  Average: {statistics.mean(mean_diffs):.2f}%")
        print(f"  Std:     {statistics.stdev(mean_diffs):.2f}%")
        print(f"  Min:     {min(mean_diffs):.2f}%")
        print(f"  Max:     {max(mean_diffs):.2f}%")
        
        print(f"\nStd bitscore difference (sample - full):")
        print(f"  Average: {statistics.mean(std_diffs):.2f}%")
        print(f"  Std:     {statistics.stdev(std_diffs):.2f}%")
        print(f"  Min:     {min(std_diffs):.2f}%")
        print(f"  Max:     {max(std_diffs):.2f}%")
        
        if tests_with_pvalues:
            p_values = [r['p_value'] for r in tests_with_pvalues]
            print(f"\nStatistical significance (t-test, α=0.05):")
            print(f"  Tests performed: {len(tests_with_pvalues)}")
            print(f"  Significant differences: {len(significant_tests)} ({len(significant_tests)/len(tests_with_pvalues)*100:.1f}%)")
            print(f"  Non-significant differences: {len(non_significant_tests)} ({len(non_significant_tests)/len(tests_with_pvalues)*100:.1f}%)")
            print(f"  Average p-value: {statistics.mean(p_values):.6f}")
            print(f"  Median p-value: {statistics.median(p_values):.6f}")
            print(f"  Min p-value: {min(p_values):.6f}")
            print(f"  Max p-value: {max(p_values):.6f}")
    else:
        print("No valid comparisons available.")

if __name__ == "__main__":
    main()

