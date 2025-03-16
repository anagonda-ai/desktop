import os
import csv
import matplotlib.pyplot as plt

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
                max_lines - cross_lines
            ])

    # Save results to a CSV file
    with open(output_file, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Directory", "Largest Chromosome File", "Largest Chromosome Lines", "Cross Chromosome Lines", "Difference"])
        writer.writerows(results)

    print(f"Results saved to {output_file}")
    
    return results
    
def plot_results(results, output_dir):
    """Plots comparison between chromosome_*.csv and cross_chromosome_clusters.csv."""
    directories = [os.path.basename(res[0]) for res in results]
    chrom_lines = [res[2] for res in results]
    cross_lines = [res[3] for res in results]

    plt.figure(figsize=(12, 6))
    bar_width = 0.4
    x_indexes = range(len(directories))

    plt.bar(x_indexes, chrom_lines, width=bar_width, label="Largest Chromosome File", alpha=0.7)
    plt.bar([x + bar_width for x in x_indexes], cross_lines, width=bar_width, label="Cross Chromosome File", alpha=0.7)

    plt.xlabel("Directories")
    plt.ylabel("Number of Lines")
    plt.title("Comparison of Largest Chromosome File vs. Cross Chromosome File")
    plt.xticks([x + bar_width / 2 for x in x_indexes], directories, rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, "comparison_plot.png"))  # Save plot as an image file
    plt.show()  # Display the plot

if __name__ == "__main__":
    root_directory = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated/best_hits_by_organism"  # Replace with the actual root directory
    output_filename = "comparison_results.csv"  # Output file name
    output_file = os.path.join(root_directory, output_filename)

    results = compare_csvs_in_each_dir(root_directory, output_file)
    if results:
        plot_results(results, root_directory)
