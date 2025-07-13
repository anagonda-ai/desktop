import argparse
import os
import pandas as pd
import concurrent.futures
import re

def find_csv_files(root_dir):
    """Recursively find all CSV files in the given directory."""
    csv_files = {}
    for root, _, files in os.walk(root_dir):
        organism_name = root
        for file in files:
            if file.endswith(".csv"):
                if organism_name not in csv_files:
                    csv_files[organism_name] = []
                csv_files[organism_name].append(os.path.join(root, file))
    return csv_files

def detect_delimiter(file_path):
    """Detects the delimiter used in the CSV file."""
    with open(file_path, 'r', encoding='utf-8') as file:
        first_line = file.readline()
        if '\t' in first_line:
            return '\t'  # Tab-separated
        return ','  # Default to comma-separated

def calculate_gene_distances(df):
    """Calculate gene distances and return distance values, ensuring no zero distances."""
    df = df.sort_values(by='start')  # Ensure genes are ordered
    first_gene_start = df['start'].iloc[0]
    last_gene_end = df['end'].iloc[-1]
    total_distance = last_gene_end - first_gene_start
    return [total_distance]

def select_highest_score_genes(df):
    """Select the gene with the highest score if it appears in multiple chromosomes."""
    if 'bit_score' not in df.columns:
        print("Warning: Required column 'bit_score' not found. Using all genes as-is.")
        return df
    return df.loc[df.groupby('origin_file')['bit_score'].idxmax()]

def process_organism(organism, file_paths, output_dir, all_distances):
    print(f"Processing organism: {organism}")
    organism_output_dir = os.path.join(organism, output_dir)
    os.makedirs(organism_output_dir, exist_ok=True)
    
    # Dictionary to store data per chromosome
    chromosome_data = {}
    
    # Read each CSV and collect genes per chromosome
    for file_path in file_paths:
        if os.stat(file_path).st_size == 0:  # Check if file is empty
            print(f"Skipping empty file: {file_path}")
            continue
        
        delimiter = detect_delimiter(file_path)
        
        try:
            df = pd.read_csv(file_path, delimiter=delimiter, on_bad_lines="skip")
        except (pd.errors.EmptyDataError, pd.errors.ParserError):
            print(f"Skipping invalid or malformed file: {file_path}")
            continue
        
        # Normalize column names
        df.columns = df.columns.str.strip().str.lower()
        
        # Ensure required columns exist
        if 'chromosome' not in df.columns:
            print(f"Skipping file {file_path}, missing required columns.")
            continue
        
        # Add origin file column
        df['origin_file'] = os.path.basename(file_path).split(".")[0]
        
        # Group genes by chromosome
        for _, row in df.iterrows():
            chrom = row['chromosome']
            if chrom not in chromosome_data:
                chromosome_data[chrom] = []
            chromosome_data[chrom].append(row)
    
    # Create separate CSVs and collect distances for general statistics
    local_distances = []
    
    if chromosome_data:
        all_genes_df = pd.concat([pd.DataFrame(rows) for rows in chromosome_data.values()], ignore_index=True)
    else:
        all_genes_df = pd.DataFrame()  # Prevent error when concatenating empty list
    
    # Process chromosomes separately
    for chrom, rows in chromosome_data.items():
        if len(rows) < 2:
            print(f"Skipping chromosome {chrom} for {organism}, not enough data.")
            continue
        
        chrom_df = pd.DataFrame(rows)
        output_file = os.path.join(organism_output_dir, f"chromosome_{chrom}.csv")
        chrom_df.to_csv(output_file, index=False)
        print(f"Saved {output_file}")
        
        # Compute gene distances and store them
        distances = calculate_gene_distances(chrom_df)
        local_distances.extend(distances)
    
    # Save cross-chromosome clusters with highest scores
    if not all_genes_df.empty:
        highest_score_genes_df = select_highest_score_genes(all_genes_df)
        cross_chrom_output_file = os.path.join(organism_output_dir, "cross_chromosome_clusters.csv")
        highest_score_genes_df.to_csv(cross_chrom_output_file, index=False)
        print(f"Saved cross-chromosome clusters: {cross_chrom_output_file}")
    
    # Append to shared distance list
    all_distances.extend(local_distances)

def process_csv_files(root_dir):
    # Define input files
    best_hits_by_organism = os.path.join(root_dir, "blast_results_chromosome_separated/best_hits_by_organism")
    # Find all relevant CSV files per organism
    organism_files = find_csv_files(best_hits_by_organism)
    
    # Define output directory
    output_dir = "potential_clusters_by_chromosome"
    
    # List to store all gene distances for summary statistics
    all_distances = []
    
    # Use concurrency to process organisms in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for organism, file_paths in organism_files.items():
            futures.append(executor.submit(process_organism, organism, file_paths, output_dir, all_distances))
        
        # Wait for all futures to complete
        for future in concurrent.futures.as_completed(futures):
            future.result()
    
    # Save general summary statistics file
    if all_distances:
        stats = {
            'min_distance': min(all_distances),
            'max_distance': max(all_distances),
            'mean_distance': sum(all_distances) / len(all_distances),
            'median_distance': sorted(all_distances)[len(all_distances) // 2]
        }
        stats_df = pd.DataFrame([stats])
        stats_output_file = os.path.join(root_dir, "general_summary_statistics.csv")
        stats_df.to_csv(stats_output_file, index=False)
        print(f"Saved general summary statistics: {stats_output_file}")

def _extract_organism(origin_file):
    m = re.search(r'\.fasta_(.+)', origin_file)
    if m:
        return m.group(1)
    m = re.search(r'^[^_]+_(.+)', origin_file)
    if m:
        return m.group(1)
    return origin_file

def _parse_blast_file(file_path):
    hits = []
    origin_file = os.path.basename(file_path).split('.gene_transformed_filtered')[0]
    organism = _extract_organism(origin_file)
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
    except Exception as e:
        print(f"Skipping {file_path}: {e}")
        return hits
    print(f"Parsing {file_path} with {len(df)} rows")
    for _, row in df.iterrows():
        try:
            sseqid = str(row[1])
            sseqid_parts = sseqid.split('|')
            if len(sseqid_parts) >= 5:
                chrom = sseqid_parts[1]
                start = int(sseqid_parts[2])
                end = int(sseqid_parts[3])
                strand = sseqid_parts[4]
                row_index = int(sseqid_parts[5])
            else:
                chrom = 'Unknown'
                start = end = row_index = None
            hit = {
                'origin_file': origin_file,
                'organism': organism,
                'chromosome': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'row_index': row_index,
                'qseqid': row[0],
                'sseqid': row[1],
                'pident': row[2],
                'length': row[3],
                'mismatch': row[4],
                'gapopen': row[5],
                'qstart': row[6],
                'qend': row[7],
                'sstart': row[8],
                'send': row[9],
                'evalue': row[10],
                'bitscore': row[11],
            }
            hits.append(((organism, chrom), hit))
        except Exception:
            continue
    return hits

def _tightest_cluster_for_chromosome(args):
    (organism, chrom), hits, output_dir = args
    import pandas as pd
    import itertools
    from collections import defaultdict
    import os
    org_dir = os.path.join(output_dir, organism)
    os.makedirs(org_dir, exist_ok=True)
    print(f"[Organism {organism} | Chromosome {chrom}] Processing {len(hits)} hits...")
    origin_groups = defaultdict(list)
    for hit in hits:
        origin_groups[hit['origin_file']].append(hit)
    origins = list(origin_groups.keys())
    n = len(origins)
    if n < 2:
        print(f"[Organism {organism} | Chromosome {chrom}] Skipped: less than 2 origins.")
        return None
    best_cluster = None
    best_span = None
    best_size = 0
    for group_size in range(n, 1, -1):
        print(f"[Organism {organism} | Chromosome {chrom}] Trying group size {group_size}...")
        found = False
        for origin_subset in itertools.combinations(origins, group_size):
            groups = [origin_groups[origin] for origin in origin_subset]
            for combo in itertools.product(*groups):
                chrom_set = set(h['chromosome'] for h in combo)
                if len(chrom_set) != 1 or list(chrom_set)[0] != chrom:
                    continue
                starts = [h['start'] for h in combo if h['start'] is not None]
                ends = [h['end'] for h in combo if h['end'] is not None]
                if not starts or not ends:
                    continue
                span = max(ends) - min(starts)
                if (best_span is None) or (group_size > best_size) or (group_size == best_size and span < best_span):
                    best_span = span
                    best_cluster = combo
                    best_size = group_size
                    found = True
        if found:
            print(f"[Organism {organism} | Chromosome {chrom}] Found cluster with {best_size} genes, span {best_span}.")
            break
    if best_cluster is not None and len(best_cluster) >= 2:
        tightest_df = pd.DataFrame(best_cluster)
        tightest_file = os.path.join(org_dir, f"tightest_cluster_{organism}_{chrom}.csv")
        tightest_df.to_csv(tightest_file, index=False)
        print(f"[Organism {organism} | Chromosome {chrom}] Saved tightest cluster: {tightest_file}")
        # Return info for best-of-all selection
        return {'organism': organism, 'chromosome': chrom, 'file': tightest_file, 'span': best_span, 'size': best_size}
    print(f"[Organism {organism} | Chromosome {chrom}] No valid cluster found.")
    return None

def compute_tightest_clusters_from_raw_blast(directory, output_dir, max_workers=30):
    """
    State-of-the-art: Concurrently parse all *.csv_results.txt files and compute tightest clusters per chromosome.
    Uses ThreadPoolExecutor for I/O and ProcessPoolExecutor for computation.
    """
    import glob
    from collections import defaultdict
    import concurrent.futures
    import os
    import pandas as pd

    os.makedirs(output_dir, exist_ok=True)
    all_hits_chromosome_data = defaultdict(list)
    files = glob.glob(os.path.join(directory, '*.csv_results.txt'))

    print(f"[Main] Parsing {len(files)} BLAST result files with up to {max_workers} workers...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_parse_blast_file, files))
    for file_hits in results:
        for key, hit in file_hits:
            all_hits_chromosome_data[key].append(hit)

    print(f"[Main] Parsed hits for {len(all_hits_chromosome_data)} (organism, chromosome) pairs. Starting tightest cluster search...")
    cluster_args = [(key, hits, output_dir) for key, hits in all_hits_chromosome_data.items()]

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_tightest_cluster_for_chromosome, cluster_args))

    print(f"[Main] Finished tightest cluster search for all (organism, chromosome) pairs.")
    # Pick the best cluster for each organism (largest size, then smallest span)
    best_by_organism = {}
    for res in results:
        if res and 'organism' in res:
            org = res['organism']
            if org not in best_by_organism:
                best_by_organism[org] = res
            else:
                prev = best_by_organism[org]
                # Prefer larger size, then smaller span
                if (res['size'] > prev['size']) or (res['size'] == prev['size'] and res['span'] < prev['span']):
                    best_by_organism[org] = res
    # Save the best cluster for each organism
    for org, info in best_by_organism.items():
        org_dir = os.path.join(output_dir, org)
        best_file = os.path.join(org_dir, f"tightest_cluster_best_{org}.csv")
        df = pd.read_csv(info['file'])
        df.to_csv(best_file, index=False)
        print(f"[Main] Saved best tightest cluster for organism {org}: {best_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find homolog genes for Asaph cluster")
    parser.add_argument("--example_mgc", type=str, required=True, help="Path to example MGC directory")
    args = parser.parse_args()

    root_dir = args.example_mgc
    
    # Process files
    process_csv_files(root_dir)
    compute_tightest_clusters_from_raw_blast(directory=os.path.join(root_dir, "blast_results_chromosome_separated"),
                                              output_dir=os.path.join(root_dir, "tightest_clusters"),
                                              max_workers=30)
    print("Processing completed.")
