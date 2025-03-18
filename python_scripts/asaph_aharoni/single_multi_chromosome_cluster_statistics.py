import os
import csv
import matplotlib.pyplot as plt
import numpy as np

def count_lines_in_csv(file_path):
    """Counts the number of lines in a CSV file efficiently, excluding the header."""
    with open(file_path, 'r', encoding='utf-8') as f:
        return sum(1 for _ in f) - 1  # Subtract 1 to exclude the header

def find_largest_chromosome_file(dir_path):
    """Finds the chromosome_*.csv file with the most lines in a given directory."""
    max_lines = 0
    largest_file = None

    for file in os.listdir(dir_path):
        if file.startswith("chromosome_") and file.endswith(".csv"):
            file_path = os.path.join(dir_path, file)
            num_lines = count_lines_in_csv(file_path)
            if num_lines > max_lines:
                max_lines = num_lines
                largest_file = file

    return max_lines, largest_file

def compare_csvs_in_each_dir(root_dir, output_file):
    """Traverses 'x' to find all 'potential_clusters_by_chromosome' directories and compare CSVs."""
    results = []

    for dirpath, dirnames, filenames in os.walk(root_dir):
        if os.path.basename(dirpath) == "potential_clusters_by_chromosome":
            # Find the largest chromosome_*.csv file
            max_lines, largest_file = find_largest_chromosome_file(dirpath)
            
            # Find cross_chromosome_clusters.csv
            cross_chromosome_file = os.path.join(dirpath, "cross_chromosome_clusters.csv")
            cross_lines = count_lines_in_csv(cross_chromosome_file) if os.path.exists(cross_chromosome_file) else 0
            
            results.append([
                dirpath,
                largest_file if largest_file else "None",
                max_lines,
                cross_lines,
                cross_lines - max_lines
            ])

    # Save results to a CSV file
    with open(output_file, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Directory", "Largest Chromosome File", "Largest Chromosome Lines", "Cross Chromosome Lines", "Difference"])
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
        plt.savefig(plot_filename)
        plt.close()  # Close the figure to save memory

    print(f"Plots saved in: {sub_dir}")

if __name__ == "__main__":
    root_directory = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated/best_hits_by_organism"  # Replace with the actual root directory
    output_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output"  # Replace with the actual output directory
    output_filename = "comparison_results.csv"  # Output file name
    output_file = os.path.join(output_dir, output_filename)

    results = compare_csvs_in_each_dir(root_directory, output_file)
    
    
    if results:
        plot_results(results, output_dir)
