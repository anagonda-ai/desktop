import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram

def plot_dendrogram(linkage_matrix, labels):
    plt.figure(figsize=(8, 4))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=90)
    plt.tight_layout()
    plt.savefig("dendrogram.png")
    plt.close()

def plot_heatmap(npp_matrix):
    sns.clustermap(npp_matrix, cmap="vlag", metric="correlation", method="average", figsize=(10, 8))
    plt.savefig("heatmap.png")
    plt.close()
