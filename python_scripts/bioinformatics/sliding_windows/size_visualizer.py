import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import re

directory_path = "/groups/itay_mayrose/alongonda/Plant_MGC/unique_basepair_based_only_metabolic_genes_input"

window_sizes = []
row_counts = []
num_of_genes = []
average_cluster_sizes = []

for file_name in os.listdir(directory_path):
    file_path = os.path.join(directory_path, file_name)
    if os.path.isfile(file_path):
        match = re.search(r'potential_groups_w(\d+)', file_name)
        if match:
            window_size = int(match.group(1))
            try:
                df = pd.read_csv(file_path, sep=",")
                window_sizes.append(window_size)
                row_counts.append(len(df))
                print(len(df))
                
                # num_genes = df['metabolic_genes'].str.split(',').map(len).sum()
                # Split the lists, flatten them into a single list, and get distinct values
                distinct_genes = set(gene for genes in df['metabolic_genes'].str.split(',') for gene in genes)

                # Count the total number of distinct genes
                num_distinct_genes = len(distinct_genes)
                num_of_genes.append(num_distinct_genes)
                
                avg_cluster_size = df['metabolic_genes'].str.split(',').map(len).mean()
                average_cluster_sizes.append(avg_cluster_size)
                
            except Exception as e:
                print(f"Could not process file {file_name}: {e}")


plt.figure(figsize=(12, 20))

plt.subplot(4, 1, 1)
bars = plt.bar(window_sizes, row_counts, color='skyblue')
# Fit a linear regression line to the data
coefficients = np.polyfit(window_sizes, row_counts, 1)
polynomial = np.poly1d(coefficients)
x_fit = np.linspace(min(window_sizes), max(window_sizes), 100)
y_fit = polynomial(x_fit)
plt.plot(x_fit, y_fit, color='red', linestyle='--')

# Add equation of the line to the plot
equation_text = f"y = {coefficients[0]:.2f}x + {coefficients[1]:.2f}"
plt.text(0.05, 0.95, equation_text, transform=plt.gca().transAxes, fontsize=12, verticalalignment='top', bbox={"facecolor": "white", "alpha": 0.5, "pad": 5})
plt.xticks(window_sizes)
plt.xlabel('Window Size (kbps)')
plt.ylabel('Number of Unique Metabolic Gene Clusters Candidates')
plt.title('Distribution of Unique Metabolic Gene Candidates (From Single Pathway) Across Window Sizes')
plt.tight_layout()

# Adding specific values to each bar in the first subplot
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, int(yval), va='bottom')  # va='bottom' to place text above the bar

plt.subplot(4, 1, 2)
bars_genes = plt.bar(window_sizes, num_of_genes, color='lightgreen')
# Fit a linear regression line to the data
coefficients_genes = np.polyfit(window_sizes, num_of_genes, 1)
polynomial_genes = np.poly1d(coefficients_genes)
x_fit_genes = np.linspace(min(window_sizes), max(window_sizes), 100)
y_fit_genes = polynomial_genes(x_fit_genes)
plt.plot(x_fit_genes, y_fit_genes, color='red', linestyle='--')

# Add equation of the line to the plot
equation_text_genes = f"y = {coefficients_genes[0]:.2f}x + {coefficients_genes[1]:.2f}"
plt.text(0.05, 0.95, equation_text_genes, transform=plt.gca().transAxes, fontsize=12, verticalalignment='top', bbox={"facecolor": "white", "alpha": 0.5, "pad": 5})
plt.xticks(window_sizes)
plt.xlabel('Window Size (kbps)')
plt.ylabel('Number of Genes')
plt.title('Distribution of Genes (Sum From All The Candidates) Across Window Sizes')
plt.tight_layout()

# Adding specific values to each bar in the second subplot
for bar in bars_genes:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, int(yval), va='bottom')  # va='bottom' to place text above the bar

plt.subplot(4, 1, 3)
bars_avg_cluster = plt.bar(window_sizes, average_cluster_sizes, color='orange')
# Fit a linear regression line to the data
coefficients_avg_cluster = np.polyfit(window_sizes, average_cluster_sizes, 1)
polynomial_avg_cluster = np.poly1d(coefficients_avg_cluster)
x_fit_avg_cluster = np.linspace(min(window_sizes), max(window_sizes), 100)
y_fit_avg_cluster = polynomial_avg_cluster(x_fit_avg_cluster)
plt.plot(x_fit_avg_cluster, y_fit_avg_cluster, color='red', linestyle='--')

# Add equation of the line to the plot
equation_text_avg_cluster = f"y = {coefficients_avg_cluster[0]:.2f}x + {coefficients_avg_cluster[1]:.2f}"
plt.text(0.05, 0.95, equation_text_avg_cluster, transform=plt.gca().transAxes, fontsize=12, verticalalignment='top', bbox={"facecolor": "white", "alpha": 0.5, "pad": 5})
plt.xticks(window_sizes)
plt.xlabel('Window Size (kbps)')
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


output_path_genes = os.path.join(directory_path, "window_size_histogram_genes.png")
plt.savefig(output_path_genes, dpi=300)
plt.close()