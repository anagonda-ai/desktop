#!/usr/bin/env python3
import os
import glob
import csv
import statistics
import concurrent.futures
from collections import defaultdict

# Base directory containing all dataset folders
BASE_DIR = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes"

def calculate_lpp_stats_for_file(txt_file):
    """Calculate mean, std, and count of LPP values from CSV-formatted TXT file."""
    lpp_values = []
    filename = os.path.basename(txt_file)
    
    try:
        with open(txt_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    lpp_str = row.get('LPP', '')
                    if lpp_str and lpp_str.lower() != 'none':
                        lpp = float(lpp_str)
                        lpp_values.append(lpp)
                except (ValueError, KeyError):
                    continue
        
        if len(lpp_values) == 0:
            return filename, None, None, 0
        
        mean_lpp = statistics.mean(lpp_values)
        std_lpp = statistics.stdev(lpp_values) if len(lpp_values) > 1 else 0.0
        count = len(lpp_values)
        
        return filename, mean_lpp, std_lpp, count
    
    except Exception as e:
        print(f"Error processing {txt_file}: {e}")
        return filename, None, None, 0

def process_dataset(dataset_dir):
    """Process all LPP CSV files in a dataset's best_hits_fixed directory."""
    dataset_name = os.path.basename(dataset_dir)
    best_hits_dir = os.path.join(dataset_dir, "best_hits_fixed")
    
    if not os.path.exists(best_hits_dir):
        return dataset_name, [], []
    
    # Find all LPP_matrix_*.txt files
    txt_files = glob.glob(os.path.join(best_hits_dir, "LPP_matrix_*.txt"))
    txt_files.sort()
    
    if not txt_files:
        return dataset_name, [], []
    
    print(f"\nProcessing dataset: {dataset_name}")
    print(f"Found {len(txt_files)} LPP matrix files")
    
    results = []
    organism_list = []
    
    # Process files concurrently
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(calculate_lpp_stats_for_file, txt_file): txt_file 
                   for txt_file in txt_files}
        
        for future in concurrent.futures.as_completed(futures):
            filename, mean_lpp, std_lpp, count = future.result()
            
            # Extract organism name from filename
            organism = filename.replace("LPP_matrix_", "").replace(".csv", "")
            organism_list.append(organism)
            
            if mean_lpp is not None:
                results.append({
                    'organism': organism,
                    'filename': filename,
                    'mean_lpp': mean_lpp,
                    'std_lpp': std_lpp,
                    'count': count
                })
                print(f"  ✓ {organism}: mean_LPP={mean_lpp:.4f}, std={std_lpp:.4f}, n={count}")
            else:
                print(f"  ✗ {organism}: No valid LPP values found")
    
    # Sort results by organism name
    results.sort(key=lambda x: x['organism'])
    
    # Save per-dataset summary
    summary_file = os.path.join(best_hits_dir, "LPP_statistics_summary.csv")
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['organism', 'filename', 'mean_lpp', 'std_lpp', 'num_records'])
        for res in results:
            writer.writerow([
                res['organism'],
                res['filename'],
                f"{res['mean_lpp']:.6f}",
                f"{res['std_lpp']:.6f}",
                res['count']
            ])
    
    print(f"  Summary saved to: {summary_file}")
    
    return dataset_name, results, organism_list

def main():
    # Get all dataset directories
    dataset_dirs = []
    for item in sorted(os.listdir(BASE_DIR)):
        item_path = os.path.join(BASE_DIR, item)
        if os.path.isdir(item_path):
            best_hits_dir = os.path.join(item_path, "best_hits_fixed")
            if os.path.exists(best_hits_dir):
                # Check for LPP matrix files
                lpp_files = glob.glob(os.path.join(best_hits_dir, "LPP_matrix_*.txt"))
                if lpp_files:
                    dataset_dirs.append(item_path)
    
    print("="*80)
    print("LPP STATISTICS CALCULATION FOR ALL DATASETS")
    print("="*80)
    print(f"Found {len(dataset_dirs)} datasets with LPP matrices\n")
    
    if len(dataset_dirs) == 0:
        print("ERROR: No datasets with LPP matrices found!")
        print("Make sure Step 2 (create LPP matrices) has completed successfully.")
        return
    
    # Store all results
    all_results = {}
    
    # Process each dataset
    for dataset_dir in dataset_dirs:
        dataset_name, results, organisms = process_dataset(dataset_dir)
        all_results[dataset_name] = results
    
    # Create comprehensive summary across all datasets
    print("\n" + "="*80)
    print("OVERALL SUMMARY ACROSS ALL DATASETS")
    print("="*80)
    
    # Collect statistics per dataset
    dataset_stats = []
    total_genes_all = 0
    total_organisms_all = 0
    
    for dataset_name in sorted(all_results.keys()):
        results = all_results[dataset_name]
        if results:
            all_means = [r['mean_lpp'] for r in results]
            all_stds = [r['std_lpp'] for r in results]
            all_counts = [r['count'] for r in results]
            
            total_genes = sum(all_counts)
            total_genes_all += total_genes
            total_organisms_all += len(results)
            
            dataset_stats.append({
                'dataset': os.path.basename(dataset_name),
                'num_organisms': len(results),
                'avg_mean_lpp': statistics.mean(all_means),
                'avg_std_lpp': statistics.mean(all_stds),
                'avg_records': statistics.mean(all_counts),
                'total_genes': total_genes,
                'min_mean_lpp': min(all_means),
                'max_mean_lpp': max(all_means)
            })
    
    if dataset_stats:
        # Sort by dataset name
        dataset_stats.sort(key=lambda x: x['dataset'])
        
        print(f"\n{'Dataset':<35} {'Org':<6} {'Mean LPP':<12} {'Std Dev':<12} {'Total Genes':<12}")
        print("-" * 80)
        for ds in dataset_stats:
            print(f"{ds['dataset']:<35} {ds['num_organisms']:<6} "
                  f"{ds['avg_mean_lpp']:<12.4f} {ds['avg_std_lpp']:<12.4f} {ds['total_genes']:<12}")
        
        # Grand totals
        print("-" * 80)
        grand_mean_lpp = statistics.mean([ds['avg_mean_lpp'] for ds in dataset_stats])
        grand_std_lpp = statistics.mean([ds['avg_std_lpp'] for ds in dataset_stats])
        print(f"{'TOTAL':<35} {total_organisms_all:<6} {grand_mean_lpp:<12.4f} "
              f"{grand_std_lpp:<12.4f} {total_genes_all:<12}")
        
        # Save comprehensive summary
        comprehensive_summary = os.path.join(BASE_DIR, "LPP_all_datasets_statistics_summary.csv")
        with open(comprehensive_summary, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['dataset', 'num_organisms', 'avg_mean_lpp', 'avg_std_lpp', 
                           'avg_records_per_organism', 'total_genes', 'min_mean_lpp', 'max_mean_lpp'])
            for ds in dataset_stats:
                writer.writerow([
                    ds['dataset'],
                    ds['num_organisms'],
                    f"{ds['avg_mean_lpp']:.6f}",
                    f"{ds['avg_std_lpp']:.6f}",
                    f"{ds['avg_records']:.2f}",
                    ds['total_genes'],
                    f"{ds['min_mean_lpp']:.6f}",
                    f"{ds['max_mean_lpp']:.6f}"
                ])
        
        print(f"\nComprehensive summary saved to: {comprehensive_summary}")
    else:
        print("No valid results to summarize.")
    
    print("\n" + "="*80)
    print("PROCESSING COMPLETE")
    print("="*80)

if __name__ == "__main__":
    main()