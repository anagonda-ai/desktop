import os
import pandas as pd
import concurrent.futures

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
    distances = df['start'].diff().dropna()  # Calculate distances and remove NaN
    distances = distances[distances > 0].tolist()  # Remove zero or negative distances
    return distances

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
        except pd.errors.EmptyDataError:
            print(f"Skipping invalid or empty file: {file_path}")
            continue
        except pd.errors.ParserError:
            print(f"Skipping malformed CSV file: {file_path}, possible inconsistent columns.")
            continue
        
        # Normalize column names (strip whitespace, lowercase)
        df.columns = df.columns.str.strip().str.lower()
        
        # Ensure 'chromosome' column exists
        if 'chromosome' not in df.columns:
            print(f"Skipping file {file_path}, missing 'chromosome' column.")
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
    for chrom, rows in chromosome_data.items():
        if not any(row['origin_file'].startswith(('adcs', 'cs')) for row in rows) or len(rows) < 2:
            print(f"Skipping chromosome {chrom} for {organism}, not enough data.")
            continue
        
        chrom_df = pd.DataFrame(rows)
        output_file = os.path.join(organism_output_dir, f"chromosome_{chrom}.csv")
        chrom_df.to_csv(output_file, index=False)
        print(f"Saved {output_file}")
        
        # Compute gene distances and store them
        distances = calculate_gene_distances(chrom_df)
        local_distances.extend(distances)
    
    # Append to shared distance list
    all_distances.extend(local_distances)

def process_csv_files(root_dir):
    # Find all relevant CSV files per organism
    organism_files = find_csv_files(root_dir)
    
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

if __name__ == "__main__":
    # Define input files
    root_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated/best_hits_by_organism"
    
    # Process files
    process_csv_files(root_dir)
