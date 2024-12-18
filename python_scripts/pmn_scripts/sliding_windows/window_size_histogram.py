import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import re

directory_path = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs"

window_sizes = []
row_counts = []
num_of_genes = []

for file_name in os.listdir(directory_path):
    file_path = os.path.join(directory_path, file_name)
    if os.path.isfile(file_path):
        match = re.search(r'window_size_(\d+)', file_name)
        if match:
            window_size = int(match.group(1))
            try:
                df = pd.read_csv(file_path, sep=",")
                window_sizes.append(window_size)
                row_counts.append(len(df))
                
                num_genes = df['pathway_cluster_genes'].str.split(',').map(len).sum()
                num_of_genes.append(num_genes)
            except Exception as e:
                print(f"Could not process file {file_name}: {e}")



plt.figure(figsize=(12, 10))

plt.subplot(2, 1, 1)
bars = plt.bar(window_sizes, row_counts, color='skyblue')
plt.xticks(window_sizes)
plt.xlabel('Window Size')
plt.ylabel('Number of Unique Metabolic Gene Clusters Candidates')
plt.title('Distribution of Unique Metabolic Gene Candidates (From Single Pathway) Across Window Sizes')
plt.tight_layout()

# Adding specific values to each bar in the first subplot
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, int(yval), va='bottom')  # va='bottom' to place text above the bar

plt.subplot(2, 1, 2)
bars_genes = plt.bar(window_sizes, num_of_genes, color='lightgreen')
plt.xticks(window_sizes)
plt.xlabel('Window Size')
plt.ylabel('Number of Genes')
plt.title('Distribution of Genes (Sum From All The Candidates) Across Window Sizes')
plt.tight_layout()

# Adding specific values to each bar in the second subplot
for bar in bars_genes:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, int(yval), va='bottom')  # va='bottom' to place text above the bar


output_path_genes = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs/window_size_histogram_genes.png"
plt.savefig(output_path_genes, dpi=300)
plt.close()

print(f"Histogram of genes saved to {output_path_genes}")

print(f"Histogram saved to {output_path_genes}")
print(f"Current working directory: {os.getcwd()}")
print(f"Files in current directory: {os.listdir('.')}")

# פונקציה להצגת הקורלציה בין גודל החלון למספר המועמדים
def show_correlation(window_sizes, row_counts):
    correlation, _ = pearsonr(window_sizes, row_counts)
    print(f"Correlation between window size and number of candidates: {correlation:.2f}")

# הצגת הקורלציה
show_correlation(window_sizes, row_counts)