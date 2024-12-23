import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import re

directory_path = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs"

window_sizes = []
row_counts = []
num_of_genes = []
average_cluster_sizes = []

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
                
                avg_cluster_size = df['pathway_cluster_genes'].str.split(',').map(len).mean()
                average_cluster_sizes.append(avg_cluster_size)
                
            except Exception as e:
                print(f"Could not process file {file_name}: {e}")


plt.figure(figsize=(12, 20))

plt.subplot(4, 1, 1)
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

plt.subplot(4, 1, 2)
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

plt.subplot(4, 1, 3)
bars_avg_cluster = plt.bar(window_sizes, average_cluster_sizes, color='orange')
plt.xticks(window_sizes)
plt.xlabel('Window Size')
plt.ylabel('Average Cluster Size')
plt.title('Average Cluster Size Across Window Sizes')
plt.tight_layout()

# Adding specific values to each bar in the third subplot
for bar in bars_avg_cluster:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, round(yval, 2), va='bottom')  # va='bottom' to place text above the bar

# Calculate correlations
corr_candidates, _ = pearsonr(window_sizes, row_counts)
corr_genes, _ = pearsonr(window_sizes, num_of_genes)
corr_avg_cluster, _ = pearsonr(window_sizes, average_cluster_sizes)

# Add correlation information as a plot
plt.subplot(4, 1, 4)
plt.axis('off')
correlation_text = (
    f"Correlation between window size and number of candidates: {corr_candidates:.2f}\n"
    f"Correlation between window size and number of genes: {corr_genes:.2f}\n"
    f"Correlation between window size and average cluster size: {corr_avg_cluster:.2f}"
)
plt.text(0.5, 0.5, correlation_text, ha='center', va='center', fontsize=20, bbox={"facecolor": "white", "alpha": 0.5, "pad": 5})


output_path_genes = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs/window_size_histogram_genes.png"
plt.savefig(output_path_genes, dpi=300)
plt.close()