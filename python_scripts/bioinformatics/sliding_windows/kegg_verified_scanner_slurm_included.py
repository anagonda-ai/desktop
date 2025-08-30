import argparse
import os
import subprocess
import threading
import time
import numpy as np
import pandas as pd
import random
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# ------------------------------
# KEGG Database Path and Mapping
# ------------------------------

MAPPING_CSV = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping_with_fasta.csv"
mapping_df = pd.read_csv(MAPPING_CSV)
kegg_path_dict = dict(zip(mapping_df["filtered_path"], mapping_df["kegg_fasta"]))

# ------------------------------
# Argument Parsing
# ------------------------------

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='KEGG metabolic gene cluster analysis')
    
    parser.add_argument('--min_genes', 
                       type=int, 
                       default=3,
                       help='Minimum number of genes required in a group. Default: 3')
    
    parser.add_argument('--max_jobs', 
                       type=int, 
                       default=100,
                       help='Maximum number of SLURM jobs to run simultaneously. Default: 100')
    
    parser.add_argument('--window_sizes', 
                       type=int, 
                       nargs='+',
                       default=[10, 20],
                       help='Window sizes to process. Default: 10 20')
    
    # Analysis mode
    parser.add_argument('--run_shuffle_test',
                       action='store_true',
                       help='Run shuffle statistical test (200 iterations of shuffled data with sliding windows)')
    
    parser.add_argument('--num_shuffle_iterations',
                       type=int,
                       default=200,
                       help='Number of shuffle iterations. Default: 200')
    
    parser.add_argument('--shuffle_seed',
                       type=int,
                       help='Base seed for shuffle iterations (optional)')
    
    return parser.parse_args()

# ------------------------------
# Annotation Shuffling Functions
# ------------------------------

def shuffle_annotations_inplace(annotated_files, seed=None):
    """
    Shuffle the annotation data (kegg_ids, annotation, pathway) in-place across all files.
    Keeps the three annotation columns together as units.
    """
    if seed is not None:
        np.random.seed(seed)
    
    success_count = 0
    for annotated_file in annotated_files:
        try:
            # Read the file
            df = pd.read_csv(annotated_file)
            
            # Check if annotation columns exist
            annotation_columns = ['kegg_ids', 'annotation', 'pathway']
            missing_cols = [col for col in annotation_columns if col not in df.columns]
            if missing_cols:
                continue  # Skip files without proper annotation columns
            
            # Create shuffled indices for the annotation data
            shuffled_indices = np.random.permutation(len(df))
            
            # Apply shuffled annotation data
            original_annotations = df[annotation_columns].copy()
            df[annotation_columns] = original_annotations.iloc[shuffled_indices].reset_index(drop=True)
            
            # Save back to the same file
            df.to_csv(annotated_file, index=False)
            success_count += 1
            
        except Exception as e:
            print(f"⚠️  Error shuffling {annotated_file}: {str(e)}")
    
    return success_count

def backup_original_files(annotated_files, backup_dir):
    """Create backup copies of original annotated files"""
    os.makedirs(backup_dir, exist_ok=True)
    
    for file_path in annotated_files:
        filename = os.path.basename(file_path)
        backup_path = os.path.join(backup_dir, filename)
        shutil.copy2(file_path, backup_path)
    
    print(f"✅ Backed up {len(annotated_files)} files to: {backup_dir}")

def restore_original_files(annotated_files, backup_dir):
    """Restore original files from backup"""
    for file_path in annotated_files:
        filename = os.path.basename(file_path)
        backup_path = os.path.join(backup_dir, filename)
        if os.path.exists(backup_path):
            shutil.copy2(backup_path, file_path)

# ------------------------------
# String to Boolean Conversion
# ------------------------------

def str_to_bool(value):
    """Convert string representation to boolean"""
    if value.lower() in ['true', '1']:
        return True
    elif value.lower() in ['false', '0']:
        return False
    else:
        raise ValueError(f"Invalid boolean value: {value}")

# ------------------------------
# Count Running Jobs
# ------------------------------

def count_running_jobs(user="alongonda"):
    """Count number of running jobs for a user."""
    result = subprocess.run(
        ["squeue", "-u", user, "-h"],
        capture_output=True, text=True
    )
    jobs = result.stdout.strip().splitlines()
    return len(jobs)

# ------------------------------
# SLURM Submission
# ------------------------------

def submit_annotation_jobs(genome_files, temp_dir, annotated_dir, max_jobs=100):
    slurm_script = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/kegg_verified_annotation.sh"
    os.makedirs(annotated_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    for genome_file in genome_files:
        if genome_file not in kegg_path_dict:
            print(f"⚠️ No KEGG DB found for: {genome_file}")
            continue

        kegg_fasta = kegg_path_dict[genome_file]
        job_name = os.path.basename(genome_file).replace('.csv', '')
        output_file = os.path.join(annotated_dir, f"{job_name}_annotated.csv")
        if os.path.exists(output_file):
            print(f"✔️ Already exists: {output_file}")
            continue

        while True:
            current_jobs = count_running_jobs(user="alongonda")
            if current_jobs < max_jobs:
                break
            print(f"⏳ {current_jobs} jobs running (limit {max_jobs}). Waiting 30 seconds...")
            time.sleep(30)

        sbatch_cmd = [
            "sbatch",
            "--job-name", f"kegg_{job_name}",
            "--output", os.path.join(annotated_dir, f"{job_name}.out"),
            "--error", os.path.join(annotated_dir, f"{job_name}.err"),
            slurm_script,
            genome_file,
            annotated_dir,
            temp_dir,
            kegg_fasta
        ]
        result = subprocess.run(sbatch_cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✅ Submitted job for {job_name}: {result.stdout.strip()}")
        else:
            print(f"❌ Failed to submit job for {job_name}: {result.stderr.strip()}")

# ------------------------------
# Wait for SLURM Jobs
# ------------------------------

def wait_for_jobs(job_prefix="kegg_", user="alongonda", check_interval=60):
    print(f"⏳ Waiting for SLURM jobs starting with '{job_prefix}' by user '{user}'...")
    while True:
        result = subprocess.run(
            ["squeue", "-u", user, "-h", "-o", "%.18i %.9P %.100j %.8u %.2t %.10M %.6D %R"],
            capture_output=True, text=True
        )
        jobs = result.stdout.strip().splitlines()
        kegg_jobs = [j for j in jobs if job_prefix in j]

        if not kegg_jobs:
            print("✅ All KEGG annotation jobs are finished.")
            break
        else:
            print(f"🕒 {len(kegg_jobs)} jobs still running. Checking again in {check_interval} seconds...")
            time.sleep(check_interval)
            
# ------------------------------
# Sliding Window Clustering
# ------------------------------

def find_first_common_element(genes_and_pathways, min_genes):
    from itertools import combinations, product

    gene_sets = {gene: set([pathway for pathway in paths if pd.notna(pathway)])
                 for gene, paths in genes_and_pathways.items()}

    for group_size in range(len(gene_sets), min_genes - 1, -1):
        for selected_genes in combinations(gene_sets.keys(), group_size):
            selected_sets = [gene_sets[gene] for gene in selected_genes]
            for combination in product(*selected_sets):
                if len(set(combination)) == 1:
                    return combination[0], list(selected_genes)
    return None, []

# ------------------------------
# Count Num Of Sliding Window Functions
# ------------------------------

def count_sliding_windows(chr_data, window_size, min_genes):
    """
    Count how many windows would be generated by the sliding window algorithm
    """
    num_genes = len(chr_data)
    window_count = 0
    
    for i in range(num_genes):
        window = [chr_data.iloc[i]]
        start_index = chr_data.iloc[i]['index']
        for j in range(i + 1, num_genes):
            end_index = chr_data.iloc[j]['index']
            if (end_index - start_index <= window_size):
                window.append(chr_data.iloc[j])
            else:
                break
        
        if len(window) >= min_genes:
            window_count += 1
    
    return window_count

# ------------------------------
# Chromosome Random Processing Functions
# ------------------------------

def process_chromosome_random_groups(chr_data, window_size, min_genes, annotated_file, num_random_groups):
    """
    Process chromosome data using completely random groups of genes.
    """
    matches = []
    prev_matches = set()
    num_genes = len(chr_data)
    
    # If we have fewer genes than window size, skip
    if num_genes < window_size:
        return matches
    
    # Generate random groups
    for _ in range(num_random_groups):
        # Randomly select window_size genes from the chromosome
        random_indices = random.sample(range(num_genes), window_size)
        window = [chr_data.iloc[idx] for idx in random_indices]
        
        # Process the random group
        if len(window) >= min_genes:
            window_df = pd.DataFrame(window)
            genes_and_pathways = {row['id']: row['pathway'].split(",") for _, row in window_df.iterrows()}
            genes_and_annotations = {row['id']: row['annotation'].split(",") for _, row in window_df.iterrows()}
            
            pathway, metabolic_genes = find_first_common_element(genes_and_pathways, min_genes)
            
            # Only add if we found a pathway and haven't seen this exact gene combination
            if pathway and tuple(sorted(metabolic_genes)) not in prev_matches:
                prev_matches.add(tuple(sorted(metabolic_genes)))
                metabolic_annotations = [genes_and_annotations[gene][0] for gene in metabolic_genes]
                
                match = {
                    'pathway': pathway,
                    'genes': ','.join(window_df['id']),
                    'metabolic_genes': ','.join(metabolic_genes),
                    'metabolic_genes_annotations': ','.join(metabolic_annotations),
                    'start': window_df['start'].min(),
                    'end': window_df['end'].max(),
                    'source_file': annotated_file,
                    'group_type': 'random'
                }
                matches.append(match)
    
    return matches

# ------------------------------
# Chromosome Sliding Window Processing Functions
# ------------------------------

def process_chromosome_sliding(chr_data, window_size, min_genes, annotated_file):
    matches = []
    prev_matches = set()
    last_match = None
    num_genes = len(chr_data)
    i = 0
    while i < num_genes:
        window = [chr_data.iloc[i]]
        start_index = chr_data.iloc[i]['index']
        for j in range(i + 1, num_genes):
            end_index = chr_data.iloc[j]['index']
            if (end_index - start_index <= window_size):
                window.append(chr_data.iloc[j])
            else:
                break

        if len(window) >= min_genes:
            window_df = pd.DataFrame(window)
            genes_and_pathways = {row['id']: row['pathway'].split(",") for _, row in window_df.iterrows()}
            genes_and_annotations = {row['id']: row['annotation'].split(",") for _, row in window_df.iterrows()}
            pathway, metabolic_genes = find_first_common_element(genes_and_pathways, min_genes)
            if pathway and not tuple(metabolic_genes) in prev_matches:
                prev_matches.add(tuple(metabolic_genes))
                metabolic_annotations = [genes_and_annotations[gene][0] for gene in metabolic_genes]

                current_match = {
                    'pathway': pathway,
                    'genes': list(window_df['id']),
                    'metabolic_genes': metabolic_genes,
                    'metabolic_genes_annotations': metabolic_annotations,
                    'start': window_df['start'].min(),
                    'end': window_df['end'].max()
                }

                if (
                    last_match and
                    pathway == last_match['pathway'] and
                    metabolic_genes[:2] == last_match['metabolic_genes'][-2:]
                ):
                    new_genes = [g for g in current_match['genes'] if g not in last_match['genes']]
                    new_metabolic_genes = [g for g in current_match['metabolic_genes'] if g not in last_match['metabolic_genes']]
                    new_annotations = [
                        ann for gene, ann in zip(current_match['metabolic_genes'], current_match['metabolic_genes_annotations'])
                        if gene in new_metabolic_genes
                    ]
                    last_match['genes'] += new_genes
                    last_match['metabolic_genes'] += new_metabolic_genes
                    last_match['metabolic_genes_annotations'] += new_annotations
                    last_match['end'] = max(last_match['end'], current_match['end'])
                else:
                    if last_match:
                        matches.append({
                            'pathway': last_match['pathway'],
                            'genes': ','.join(last_match['genes']),
                            'metabolic_genes': ','.join(last_match['metabolic_genes']),
                            'metabolic_genes_annotations': ','.join(last_match['metabolic_genes_annotations']),
                            'start': last_match['start'],
                            'end': last_match['end'],
                            'source_file': annotated_file
                        })
                    last_match = current_match
        i += 1

    if last_match:
        matches.append({
            'pathway': last_match['pathway'],
            'genes': ','.join(last_match['genes']),
            'metabolic_genes': ','.join(last_match['metabolic_genes']),
            'metabolic_genes_annotations': ','.join(last_match['metabolic_genes_annotations']),
            'start': last_match['start'],
            'end': last_match['end'],
            'source_file': annotated_file
        })

    return matches

# ------------------------------
# Process Annotated Files
# ------------------------------

def process_annotated_file(annotated_file, output_file, window_size, min_genes, use_random_method=False, shuffle_annotations=False, shuffle_seed=None):
    file_lock = threading.Lock()
    df = pd.read_csv(annotated_file)
    
    # Shuffle annotations if requested
    if shuffle_annotations:
        annotation_columns = ['kegg_ids', 'annotation', 'pathway']
        
        # Check if all annotation columns exist
        missing_cols = [col for col in annotation_columns if col not in df.columns]
        if not missing_cols:  # Only shuffle if all columns exist
            if shuffle_seed is not None:
                np.random.seed(shuffle_seed)
            
            # Create shuffled indices for the annotation data
            shuffled_indices = np.random.permutation(len(df))
            
            # Apply shuffled annotation data
            original_annotations = df[annotation_columns].copy()
            df[annotation_columns] = original_annotations.iloc[shuffled_indices].reset_index(drop=True)
    
    df["index"] = df.index
    filtered_df = df[df["pathway"].notna()]
    chromosomes = filtered_df['chromosome'].unique()

    def wrapped(chr):
        chr_data = filtered_df[filtered_df['chromosome'] == chr]
        return process_chromosome_sliding(chr_data, window_size, min_genes, annotated_file)

    with ThreadPoolExecutor() as executor:
        results = executor.map(wrapped, chromosomes)

    all_matches = []
    for match_list in results:
        all_matches.extend(match_list)

    if all_matches and output_file:
        with file_lock:
            mode = 'a' if os.path.exists(output_file) else 'w'
            header = not os.path.exists(output_file)
            pd.DataFrame(all_matches).to_csv(output_file, mode=mode, header=header, index=False)

    return len(all_matches)

# ------------------------------
# Shuffle Statistical Test
# ------------------------------

def run_shuffle_statistical_test(annotated_files, window_size, min_genes, num_iterations, base_seed=None):
    """
    Run shuffle statistical test: 200 iterations of shuffled data analyzed with sliding windows.
    """
    shuffle_results = []
    
    print(f"🔄 Running {num_iterations} shuffle iterations with sliding windows...")
    
    # Use progress bar for iterations
    with tqdm(total=num_iterations, desc="Shuffle iterations", unit="iter") as pbar:
        for iteration in range(num_iterations):
            # Generate unique seed for this iteration
            iteration_seed = (base_seed or 42) + iteration
            
            # Run analysis with shuffling enabled (always sliding windows) - WITH CONCURRENCY
            total_matches = 0
            
            with ThreadPoolExecutor(max_workers=32) as executor:
                futures = [
                    executor.submit(
                        process_annotated_file,
                        annotated_file, 
                        None,  # No output file needed for shuffle test
                        window_size, 
                        min_genes, 
                        False,  # Always sliding windows for shuffle test
                        shuffle_annotations=True,  # Enable shuffling
                        shuffle_seed=iteration_seed + hash(annotated_file) % 1000  # Unique seed per file
                    )
                    for annotated_file in annotated_files
                ]
                
                for future in as_completed(futures):
                    total_matches += future.result()
            
            shuffle_results.append(total_matches)
            
            # Update progress bar with running statistics
            pbar.set_postfix({
                'Current': total_matches, 
                'Avg': f"{np.mean(shuffle_results):.1f}",
                'Std': f"{np.std(shuffle_results):.1f}" if len(shuffle_results) > 1 else "0.0"
            })
            pbar.update(1)
    
    return shuffle_results

def analyze_shuffle_results(shuffle_results, real_mgc_count=None):
    """
    Analyze the results from shuffle iterations and provide statistical summary.
    """
    shuffle_array = np.array(shuffle_results)
    
    stats = {
        'mean': np.mean(shuffle_array),
        'std': np.std(shuffle_array),
        'median': np.median(shuffle_array),
        'min': np.min(shuffle_array),
        'max': np.max(shuffle_array),
        'q25': np.percentile(shuffle_array, 25),
        'q75': np.percentile(shuffle_array, 75)
    }
    
    print("\n" + "="*60)
    print("SHUFFLE STATISTICAL TEST RESULTS")
    print("="*60)
    print(f"Number of iterations: {len(shuffle_results)}")
    print(f"Mean MGCs per iteration: {stats['mean']:.1f}")
    print(f"Standard deviation: {stats['std']:.1f}")
    print(f"Median: {stats['median']:.1f}")
    print(f"Range: {stats['min']} - {stats['max']}")
    print(f"25th percentile: {stats['q25']:.1f}")
    print(f"75th percentile: {stats['q75']:.1f}")
    
    if real_mgc_count is not None:
        print(f"\nCOMPARISON TO REAL DATA:")
        print(f"Real data MGCs: {real_mgc_count}")
        z_score = (real_mgc_count - stats['mean']) / stats['std'] if stats['std'] > 0 else 0
        print(f"Z-score: {z_score:.3f}")
        
        # Calculate empirical p-value
        if real_mgc_count >= stats['mean']:
            extreme_count = np.sum(shuffle_array >= real_mgc_count)
        else:
            extreme_count = np.sum(shuffle_array <= real_mgc_count)
        p_value = extreme_count / len(shuffle_results)
        print(f"Empirical p-value (one-tailed): {p_value:.4f}")
        
        # Two-tailed p-value
        extreme_count_two = np.sum(np.abs(shuffle_array - stats['mean']) >= np.abs(real_mgc_count - stats['mean']))
        p_value_two = extreme_count_two / len(shuffle_results)
        print(f"Empirical p-value (two-tailed): {p_value_two:.4f}")
    
    print("="*60)
    
    return stats

# ------------------------------
# Subset Filtering
# ------------------------------

def remove_subset_results(output_file):
    if not os.path.exists(output_file):
        return
        
    df = pd.read_csv(output_file)
    df['metabolic_genes_set'] = df['metabolic_genes'].apply(lambda x: set(x.split(',')))
    to_remove = set()
    for i, genes_i in enumerate(df['metabolic_genes_set']):
        for j, genes_j in enumerate(df['metabolic_genes_set']):
            if i != j and genes_i < genes_j:
                to_remove.add(i)
                break
    filtered_df = df.drop(list(to_remove)).drop(columns=['metabolic_genes_set'])
    filtered_output_file = output_file.replace(".csv", "_filtered.csv")
    filtered_df.to_csv(filtered_output_file, index=False)
    print(f"🧹 Removed {len(to_remove)} subset results: {filtered_output_file}")

# ------------------------------
# Main Function
# ------------------------------

def main():
    
    # Parse command line arguments
    args = parse_arguments()
    
    min_genes = args.min_genes
    max_jobs = args.max_jobs
    window_sizes = args.window_sizes
    
    # Print configuration
    print("="*50)
    print("CONFIGURATION:")
    print(f"  Analysis mode: {'Shuffle Test (200 iterations)' if args.run_shuffle_test else 'Normal Sliding Window'}")
    print(f"  Min genes: {min_genes}")
    print(f"  Window sizes: {window_sizes}")
    print(f"  Max SLURM jobs: {max_jobs}")
    if args.run_shuffle_test:
        print(f"  Shuffle iterations: {args.num_shuffle_iterations}")
    print("="*50)
    
    # Directories
    full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
    genome_dirs = [
        os.path.join(full_genome_dir, "final_dataset/filtered_no_mito")
    ]
    
    head_output_dir = f"/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_{min_genes}_overlap_merge"
    output_dir = os.path.join(head_output_dir, "kegg_scanner_min_genes_based_metabolic")
    temp_dir = os.path.join(head_output_dir, "blast_temp_annotated_metabolic")
    annotated_dir = os.path.join(head_output_dir, "annotated_genomes_metabolic")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(annotated_dir, exist_ok=True)

    # Genome files
    genome_files = []
    for genome_dir in genome_dirs:
        for file in os.listdir(genome_dir):
            full_path = os.path.join(genome_dir, file)
            if file.endswith('.csv') and full_path in kegg_path_dict and pd.notna(kegg_path_dict[full_path]):
                genome_files.append(full_path)

    print(f"🧬 {len(genome_files)} genome files with KEGG matches found.")

    # Submit SLURM jobs for annotation
    submit_annotation_jobs(genome_files=genome_files, temp_dir=temp_dir, annotated_dir=annotated_dir, max_jobs=max_jobs)

    # Wait until all jobs finish
    wait_for_jobs(job_prefix="kegg_", user="alongonda")

    # Get annotated files
    annotated_files = [os.path.join(annotated_dir, f) for f in os.listdir(annotated_dir) if f.endswith('_annotated.csv')]
    
    if not annotated_files:
        print("❌ No annotated files found!")
        return

    min_genes_subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    os.makedirs(min_genes_subdir, exist_ok=True)

    # Run analysis based on mode
    if args.run_shuffle_test:
        # Option 2: Shuffle statistical test
        print("\n" + "="*60)
        print("RUNNING SHUFFLE STATISTICAL TEST")
        print("="*60)
        
        for window_size in window_sizes:
            print(f"\n🧪 Running shuffle test for window size {window_size}")
            print(f"   Method: 200 iterations of shuffled data with sliding windows")
            
            # Run shuffle statistical test
            shuffle_results = run_shuffle_statistical_test(
                annotated_files=annotated_files,
                window_size=window_size,
                min_genes=min_genes,
                num_iterations=args.num_shuffle_iterations,
                base_seed=args.shuffle_seed
            )
            
            # Analyze and display results
            stats = analyze_shuffle_results(shuffle_results)
            
            # Save shuffle results
            shuffle_output_file = os.path.join(min_genes_subdir, f"shuffle_results_w{window_size}.csv")
            results_df = pd.DataFrame({
                'iteration': range(len(shuffle_results)),
                'mgc_count': shuffle_results
            })
            results_df.to_csv(shuffle_output_file, index=False)
            print(f"📊 Shuffle results saved to: {shuffle_output_file}")
        
        print(f"\n✅ Shuffle statistical test completed!")
        
    else:
        # Option 1: Normal sliding window analysis on real data
        print("\n" + "="*60)
        print("RUNNING NORMAL SLIDING WINDOW ANALYSIS")
        print("="*60)
        
        for window_size in window_sizes:
            output_file = os.path.join(min_genes_subdir, f"potential_groups_w{window_size}_real_data.csv")
            
            if os.path.exists(output_file):
                os.remove(output_file)

            total_matches = 0

            desc = f"Real data sliding window w{window_size}"
            with tqdm(total=len(annotated_files), desc=desc, unit='file') as pbar:
                with ThreadPoolExecutor() as executor:
                    futures = [
                        executor.submit(
                            process_annotated_file,
                            annotated_file,
                            output_file,
                            window_size,
                            min_genes,
                            False,  # Always sliding windows for real data
                            shuffle_annotations=False,  # No shuffling for real data
                            shuffle_seed=None
                        )
                        for annotated_file in annotated_files
                    ]
                    for future in as_completed(futures):
                        total_matches += future.result()
                        pbar.update(1)

            print(f"✔️ Total MGCs found for w{window_size}, min_genes={min_genes}: {total_matches}")
            remove_subset_results(output_file)
            
        print(f"\n✅ Normal analysis completed!")


if __name__ == "__main__":
    main()