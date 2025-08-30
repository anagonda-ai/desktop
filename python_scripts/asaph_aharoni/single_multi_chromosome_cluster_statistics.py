import argparse
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def get_gene_order(df, base_genes):
    gene_order = {}
    for gene in base_genes:
        mask = df['origin_file'].astype(str).str.lower() == gene
        if mask.any() and 'start' in df.columns:
        # Get the minimum 'start' value for this gene
            gene_order[gene] = df.loc[mask, 'start'].min()
        else:
            gene_order[gene] = None
    # Assign order index (1-based) based on sorted 'start' values
    if 'start' in df.columns:
        # Filter genes that have a valid 'start' value
        valid_genes = {gene: val for gene, val in gene_order.items() if val is not None}
        # Sort genes by their 'start' value
        sorted_genes = sorted(valid_genes, key=lambda g: valid_genes[g])
        # Assign order index (1-based)
        for idx, gene in enumerate(sorted_genes, 1):
            gene_order[gene] = idx
        # Set None for genes not present
        for gene in gene_order:
            if gene not in sorted_genes:
                gene_order[gene] = None
    return gene_order

def get_gene_strands(df, base_genes):
    gene_strands = {}
    for gene in base_genes:
        mask = df['origin_file'].astype(str).str.lower() == gene
        if mask.any() and 'strand_value' in df.columns:
            # Get the most common strand for this gene
            gene_strands[gene] = df.loc[mask, 'strand_value'].iloc[0]
        else:
            gene_strands[gene] = None
    return gene_strands

def count_lines_in_csv(file_path, base_genes):
    """Counts the number of lines in a CSV file efficiently, excluding the header."""
    df = pd.read_csv(file_path)
    # Efficiently update gene_order based on 'origin_file' column values
    if 'origin_file' in df.columns:
        gene_order = get_gene_order(df, base_genes)
        gene_strands = get_gene_strands(df, base_genes)
    return len(df), gene_order, gene_strands

def find_largest_chromosome_file(dir_path, base_genes):
    """Finds the chromosome_*.csv file with the most lines in a given directory."""
    max_lines = 0
    largest_file = None
    best_gene_order = dict()
    best_gene_strands = dict()

    for file in os.listdir(dir_path):
        if file.startswith("chromosome_") and file.endswith(".csv"):
            file_path = os.path.join(dir_path, file)
            num_lines, current_gene_order, current_gene_strands = count_lines_in_csv(file_path, base_genes)
            if num_lines > max_lines:
                max_lines = num_lines
                largest_file = file
                best_gene_order = current_gene_order
                best_gene_strands = current_gene_strands

    return max_lines, largest_file, best_gene_order, best_gene_strands

def compare_csvs_in_each_dir(root_dir, output_file, base_genes):
    """Traverses 'x' to find all 'potential_clusters_by_chromosome' directories and compare CSVs."""
    results = []

    for dirpath, dirnames, filenames in os.walk(root_dir):
        if os.path.basename(dirpath) == "potential_clusters_by_chromosome":
            # Find the largest chromosome_*.csv file
            max_lines, largest_file, best_gene_order, best_gene_strands = find_largest_chromosome_file(dirpath, base_genes)
            
            # Find cross_chromosome_clusters.csv
            cross_chromosome_file = os.path.join(dirpath, "cross_chromosome_clusters.csv")
            cross_lines, cross_gene_order, cross_gene_strands = count_lines_in_csv(cross_chromosome_file, base_genes) if os.path.exists(cross_chromosome_file) else 0
            
            results.append([
                dirpath,
                largest_file if largest_file else "None",
                max_lines,
                cross_lines,
                cross_lines - max_lines,
                best_gene_order,
                best_gene_strands
                
            ])
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Returning results.")
    else:
        # Save results to a CSV file
        with open(output_file, mode="w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Directory", "Largest Chromosome File", "Largest Chromosome Lines", "Cross Chromosome Lines", "Difference", "MGC Genes Order", "MGC Genes Strands"])
            writer.writerows(results)

    print(f"Results saved to {output_file}")
    
    return results
    
def plot_results(results, output_dir, chunk_size=50):
    """Plots comparison between chromosome_*.csv and cross_chromosome_clusters.csv, 
       splitting into multiple plots if necessary for better visualization."""
    
    sub_dir_name = "comparison_plots"
    sub_dir = os.path.join(output_dir, sub_dir_name)
    # Ensure output directory exists
    os.makedirs(sub_dir, exist_ok=True)

    # Extract directories and values
    directories = [os.path.basename(res[0]) for res in results]
    chrom_lines = [res[2] for res in results]
    cross_lines = [res[3] for res in results]

    # Determine how many plots are needed
    total_plots = (len(directories) // chunk_size) + (1 if len(directories) % chunk_size != 0 else 0)

    for i in range(total_plots):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(directories))

        chunk_dirs = directories[start_idx:end_idx]
        chunk_chrom_lines = chrom_lines[start_idx:end_idx]
        chunk_cross_lines = cross_lines[start_idx:end_idx]

        # Create figure
        plt.figure(figsize=(16, 8))

        # Bar positions
        bar_width = 0.4
        x_indexes = np.arange(len(chunk_dirs))

        # Creating the bars
        plt.bar(x_indexes - bar_width / 2, chunk_chrom_lines, width=bar_width, 
                label="Largest Chromosome File", color="royalblue", alpha=0.8)
        plt.bar(x_indexes + bar_width / 2, chunk_cross_lines, width=bar_width, 
                label="Cross Chromosome File", color="coral", alpha=0.8)

        # Labels and title
        plt.xlabel("Directories", fontsize=12)
        plt.ylabel("Number of Lines", fontsize=12)
        plt.title(f"Comparison of Largest Chromosome File vs. Cross Chromosome File (Part {i+1})", fontsize=14)

        # Keep all x-tick labels but rotate them for better readability
        plt.xticks(x_indexes, chunk_dirs, rotation=90, ha="right", fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=12)
        plt.grid(axis="y", linestyle="--", alpha=0.7)

        # Adjust layout
        plt.tight_layout()

        # Save plot in the new directory
        plot_filename = os.path.join(sub_dir, f"comparison_plot_part_{i+1}.png")
        if not os.path.exists(plot_filename):
            print(f"Saving plot: {plot_filename}")
            plt.savefig(plot_filename)
            plt.close()  # Close the figure to save memory

    print(f"Plots saved in: {sub_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find homolog genes for Asaph cluster")
    parser.add_argument("--example_mgc", type=str, required=True, help="Path to example MGC directory")
    args = parser.parse_args()

    example_mgc = args.example_mgc
    base_genes = [os.path.basename(gene).replace(".fasta","").lower() for gene in os.listdir(example_mgc) if gene.endswith(".fasta")]

    best_hits_by_organism = os.path.join(example_mgc,"blast_results_chromosome_separated/best_hits_by_organism")  # Replace with the actual root directory
    output_filename = "comparison_results.csv"  # Output file name
    output_file = os.path.join(example_mgc, output_filename)
    
    
    results = compare_csvs_in_each_dir(best_hits_by_organism, output_file, base_genes)
    
    
    if results:
        plot_results(results, example_mgc)
