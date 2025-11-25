import os
import pandas as pd
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def collect_genes_from_organism(organism_dir):
    """
    Collect all genes from a single organism directory.
    
    Args:
        organism_dir: Path object for the organism directory
    
    Returns:
        Tuple: (organism_name, list of gene dictionaries)
    """
    organism_name = organism_dir.name
    genes = []
    
    print(f"Processing organism: {organism_name}")
    
    # Iterate through CSV files in the organism directory
    for csv_file in organism_dir.glob("*.csv"):
        try:
            df = pd.read_csv(csv_file)
            
            # Assuming each row represents a gene
            for idx, row in df.iterrows():
                gene_entry = row.to_dict()
                gene_entry['source_file'] = csv_file.name
                genes.append(gene_entry)
                
        except Exception as e:
            print(f"  Error reading {csv_file}: {e}")
    
    print(f"  Found {len(genes)} genes in {organism_name}")
    return organism_name, genes

def collect_genes_by_organism(base_dir, max_workers=None):
    """
    Collect genes organized by organism using concurrent processing.
    
    Args:
        base_dir: Path to the directory containing organism subdirectories
        max_workers: Maximum number of threads (None = use default)
    
    Returns:
        Dictionary: {organism_name: [gene_data_dicts]}
    """
    genes_by_organism = {}
    
    # Get all organism directories
    organism_dirs = [d for d in Path(base_dir).iterdir() if d.is_dir()]
    
    print(f"Found {len(organism_dirs)} organism directories")
    print("=" * 60)
    
    # Process organisms concurrently
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_organism = {
            executor.submit(collect_genes_from_organism, org_dir): org_dir
            for org_dir in organism_dirs
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_organism):
            try:
                organism_name, genes = future.result()
                genes_by_organism[organism_name] = genes
            except Exception as e:
                org_dir = future_to_organism[future]
                print(f"Error processing {org_dir.name}: {e}")
    
    return genes_by_organism

def process_single_organism(organism_name, genes, n_genes, output_path, base_dir):
    """
    Process a single organism: select random genes, save summary CSV, and create individual FASTA files.
    
    Args:
        organism_name: Name of the organism
        genes: List of gene dictionaries
        n_genes: Number of genes to select
        output_path: Path object for output directory
        base_dir: Base directory path (to locate original FASTA files)
    
    Returns:
        Tuple: (organism_name, number of genes selected)
    """
    print(f"\nProcessing {organism_name}:")
    print(f"  Total genes available: {len(genes)}")
    
    # Check if we have enough genes
    if len(genes) < n_genes:
        print(f"  Warning: Only {len(genes)} genes available. Selecting all.")
        selected_genes = genes
    else:
        # Randomly select n genes
        selected_genes = random.sample(genes, n_genes)
        print(f"  Randomly selected: {len(selected_genes)} genes")
    
    # Create organism-specific directory
    organism_output_dir = output_path / f"{organism_name}_selected_1000"
    organism_output_dir.mkdir(exist_ok=True)
    
    # Convert to DataFrame and save summary CSV
    df_output = pd.DataFrame(selected_genes)
    summary_file = organism_output_dir / f"{organism_name}_summary.csv"
    df_output.to_csv(summary_file, index=False)
    print(f"  ✓ Saved summary to: {summary_file}")
    
    # Create individual FASTA files for each gene
    print(f"  Creating individual FASTA files...")
    fasta_count = 0
    
    for idx, gene in enumerate(selected_genes):
        # Determine gene identifier (adjust column names as needed)
        # Common column names: 'gene_id', 'id', 'name', 'gene_name', 'sequence_id'
        gene_id = None
        for id_col in ['gene_id', 'id', 'name', 'gene_name', 'sequence_id', 'locus_tag']:
            if id_col in gene:
                gene_id = gene[id_col]
                break
        
        if gene_id is None:
            gene_id = f"gene_{idx + 1}"
        
        # Get sequence from CSV (check multiple possible column names)
        sequence = None
        for seq_col in ['Nucleotide_Sequence', 'sequence', 'seq', 'nucleotide', 'Sequence']:
            if seq_col in gene and pd.notna(gene[seq_col]):
                sequence = gene[seq_col]
                break
        
        if sequence and len(str(sequence)) > 0:
            # Sequence is in the CSV
            fasta_file = organism_output_dir / f"{gene_id}.fasta"
            with open(fasta_file, 'w') as f:
                f.write(f">{gene_id}\n")
                # Write sequence in 60 character lines (standard FASTA format)
                sequence_str = str(sequence)
                for i in range(0, len(sequence_str), 60):
                    f.write(sequence_str[i:i+60] + "\n")
            fasta_count += 1
    
    if fasta_count > 0:
        print(f"  ✓ Created {fasta_count} FASTA files")
    else:
        print(f"  ⚠ No sequences found in CSV. If sequences are in separate FASTA files,")
        print(f"    please specify the FASTA file location.")
    
    return organism_name, len(selected_genes)

def select_random_genes_per_organism(base_dir, n_genes=1000, output_dir="selected_genes", max_workers=None):
    """
    Randomly select n genes from EACH organism and save to separate CSV files using concurrent processing.
    
    Args:
        base_dir: Path to the directory containing organism subdirectories
        n_genes: Number of genes to randomly select per organism (default: 1000)
        output_dir: Directory to save output files
        max_workers: Maximum number of threads for concurrent processing (None = use default)
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Collect genes by organism (this step is already concurrent)
    print("Collecting genes from all organisms...")
    genes_by_organism = collect_genes_by_organism(base_dir, max_workers=max_workers)
    
    print(f"\nFound {len(genes_by_organism)} organisms")
    print("=" * 60)
    print("\nSelecting and saving genes...")
    
    # Process each organism concurrently
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks - now passing base_dir
        futures = {
            executor.submit(process_single_organism, org_name, genes, n_genes, output_path, base_dir): org_name
            for org_name, genes in genes_by_organism.items()
        }
        
        # Wait for all to complete
        results = []
        for future in as_completed(futures):
            try:
                org_name, num_selected = future.result()
                results.append((org_name, num_selected))
            except Exception as e:
                org_name = futures[future]
                print(f"Error processing {org_name}: {e}")
    
    print("\n" + "=" * 60)
    print(f"✓ All done! Files saved in '{output_dir}/' directory")
    print(f"✓ Processed {len(results)} organisms")

# Example usage
if __name__ == "__main__":
    # Set your base directory path here
    base_directory = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/origin"  # Change this to your actual path
    
    # Set random seed for reproducibility (optional)
    random.seed(42)
    
    # Select 1000 random genes from each organism (using all available CPU cores)
    select_random_genes_per_organism(
        base_directory, 
        n_genes=1000, 
        output_dir="/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes",
        max_workers=None  # None = use default (number of CPU cores)
    )