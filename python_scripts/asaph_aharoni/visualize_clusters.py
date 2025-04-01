import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import gridspec, cm
from Bio import Phylo
import numpy as np
from matplotlib.colors import to_rgb

# Define file paths
MAPPING_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/dataset_organism_mapping.csv"
TREE_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/dataset_organism_mapping.nwk"
COMPARISON_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/comparison_results.csv"

# Load dataset organism mapping file
if os.path.exists(MAPPING_FILE):
    mapping_df = pd.read_csv(MAPPING_FILE)
    mapping_dict = dict(zip(mapping_df['Original Filename'], mapping_df['Organism']))
else:
    mapping_df = None
    mapping_dict = {}
    print("Warning: Dataset organism mapping file not found.")

# Load comparison results file
if os.path.exists(COMPARISON_FILE):
    comparison_df = pd.read_csv(COMPARISON_FILE)
    
    # Extract only the relevant part of the directory to match mapping_dict
    def clean_directory_name(directory):
        return directory.split("/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated/best_hits_by_organism/")[1].split(".gene_transformed_filtered")[0]
    
    comparison_df['Cleaned_Name'] = comparison_df['Directory'].apply(clean_directory_name)
    
    def calculate_largest_chromosome_length(row):
        if not pd.isna(row['Largest Chromosome File']):
            largest_chromosome_file = row['Directory'] + '/' + row['Largest Chromosome File']
            if os.path.exists(largest_chromosome_file):
                chrom_df = pd.read_csv(largest_chromosome_file)
                chrom_df = chrom_df.sort_values(by='start')
                length = chrom_df['end'].iloc[-1] - chrom_df['start'].iloc[0]
                index_diff = int(chrom_df['row_index'].iloc[-1]) - int(chrom_df['row_index'].iloc[0]) + 1
                return length, index_diff
        return None, None

    print(comparison_df)
    comparison_df[['Largest Chromosome Length', 'Index Difference']] = pd.DataFrame(
        comparison_df.apply(calculate_largest_chromosome_length, axis=1).tolist()
    )
    
    # Match Cleaned_Name with any substring in mapping_dict keys
    def find_matching_organism(cleaned_name):
        for key in mapping_dict:
            if cleaned_name in key:
                return mapping_dict[key]
        return 'Unknown'
    
    comparison_df['Organism'] = comparison_df['Cleaned_Name'].apply(find_matching_organism)
else:
    comparison_df = None
    print("Warning: Comparison results file not found.")

# Create a dictionary for organism to gene count mapping
if comparison_df is not None:
    gene_counts = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Cross Chromosome Lines']))
    chromosome_lines = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Largest Chromosome Lines']))
    chromosome_cluster_length = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Largest Chromosome Length']))
    
    max_value = 200  # Define the cap for Index Difference
    # Replace values in Index Difference column
    comparison_df['Index Difference'] = comparison_df['Index Difference'].apply(lambda x: min(x, max_value))
    
    # Create the dictionary after applying the cap
    chromosome_cluster_length_genes = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Index Difference']))
else:
    gene_counts = chromosome_lines = chromosome_cluster_length = chromosome_cluster_length_genes = {}

# Improved Color Function
def get_color(value, min_val, max_val, cmap_name="viridis"):
    """Returns a color from a matplotlib colormap based on normalized value."""
    norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
    cmap = cm.get_cmap(cmap_name)
    rgba = cmap(norm(value))
    return mcolors.to_hex(rgba)

def label_function_factory(data_dict,  cmap_name="viridis"):
    """Creates a function for labeling tree clades with colored values."""
    values = [v for v in data_dict.values() if isinstance(v, (int, float)) and not pd.isna(v)]
    min_val, max_val = (min(values), max(values)) if values else (0, 1)

    def label_function(clade):
        label = clade.name.replace("_", " ") if clade.name else None
        value = data_dict.get(label, min_val)
        color = get_color(value, min_val, max_val, cmap_name)
        return (f"{label}: {value}", color) if label else (None, color)

    return label_function, min_val, max_val

def normalize_branch_color(branch_color):
    return (
        branch_color.red / 255,
        branch_color.green / 255,
        branch_color.blue / 255
    )

# Tree Plotting with Colorbar
def plot_tree(tree, data_dict, filename, cmap_name="viridis"):
    """Plots a phylogenetic tree with colored labels and a color bar."""
    fig, ax = plt.subplots(figsize=(15, 15))

    # Generate label function
    label_func, min_val, max_val = label_function_factory(data_dict, cmap_name)

    # Assign colors and labels to clades
    for clade in tree.find_clades():
        label, color = label_func(clade)
        clade.color = color  # Assign color to clade
        clade.name = label if clade.is_terminal() else clade.name.split(":")[0] if clade.name else None  # Set formatted label

    # Clear internal node names (non-terminal)
    for clade in tree.get_nonterminals():
        clade.name = None
        clade.confidence = None  # This hides the number shown on the branch
        
    # Draw tree
    Phylo.draw(tree, axes=ax, do_show=False)

    # Create colorbar
    sm = cm.ScalarMappable(cmap=cm.get_cmap(cmap_name), norm=mcolors.Normalize(vmin=min_val, vmax=max_val))
    sm.set_array([])  # Empty array for color bar
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.04)
    cbar.set_label(f"Value")
    
    # Manually recolor the labels
    for clade in tree.get_terminals():
        for text in ax.texts:
            if text.get_text()[1:] == clade.name:
                try:
                    color = to_rgb(normalize_branch_color(clade.color))  # Convert to RGB tuple
                    text.set_color(color)
                except ValueError:
                        print(f"Invalid color: {clade.color}")

    # Save the figure
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    plt.tight_layout()
    
def get_color_discrete(value, discrete_values=[0, 1, 2, 3, 4], cmap_name="viridis"):
    """Returns a discrete color from a colormap for specific values."""
    cmap = plt.colormaps.get_cmap(cmap_name)

    # Define a colormap with exactly 5 colors
    listed_cmap = mcolors.ListedColormap([cmap(i / (len(discrete_values) - 1)) for i in range(len(discrete_values))])

    # Define boundaries for discrete bins
    boundaries = [-0.5] + [v + 0.5 for v in discrete_values]
    norm = mcolors.BoundaryNorm(boundaries, listed_cmap.N)

    # Assign color based on the value
    rgba = listed_cmap(norm(value))
    return mcolors.to_hex(rgba)

def label_function_factory_discrete(data_dict, cmap_name="viridis"):
    """Creates a function for labeling tree clades with discrete colored values."""
    values = [v for v in data_dict.values() if isinstance(v, (int, float)) and not pd.isna(v)]
    min_val, max_val = (min(values), max(values)) if values else (0, 1)

    def label_function(clade):
        label = clade.name.replace("_", " ") if clade.name else None
        value = data_dict.get(label, min_val)
        color = get_color_discrete(value, cmap_name=cmap_name)
        return (f"{label}: {value}", color) if label else (None, color)

    return label_function, min_val, max_val

# Tree Plotting with Colorbar
def plot_tree_discrete(tree, data_dict, filename, cmap_name="viridis"):
    """Plots a phylogenetic tree with colored labels and a color bar."""
    fig, ax = plt.subplots(figsize=(15, 15))

    # Generate label function
    label_func, min_val, max_val = label_function_factory_discrete(data_dict, cmap_name)

    # Assign colors and labels to clades
    for clade in tree.find_clades():
        label, color = label_func(clade)
        clade.color = color  # Assign color to clade
        clade.name = label if clade.is_terminal() else clade.name.split(":")[0] if clade.name else None  # Set formatted label

    # **Create Discrete Colorbar**
    discrete_values = [0, 1, 2, 3, 4]  # Defined distinct values

    cmap = plt.colormaps.get_cmap(cmap_name)
    listed_cmap = mcolors.ListedColormap([cmap(i / (len(discrete_values) - 1)) for i in range(len(discrete_values))])

    norm = mcolors.BoundaryNorm(boundaries=[-0.5] + [v + 0.5 for v in discrete_values], ncolors=listed_cmap.N)
    sm = cm.ScalarMappable(cmap=listed_cmap, norm=norm)
    sm.set_array([])

    # **Ensure discrete ticks (no gradient)**
    cbar = plt.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.04, ticks=discrete_values)
    cbar.set_label("Value")
    cbar.set_ticks(discrete_values)
    cbar.set_ticklabels([str(v) for v in discrete_values])  # Ensure only these labels appear
    
    # Clear internal node names (non-terminal)
    for clade in tree.get_nonterminals():
        clade.name = None
        clade.confidence = None  # This hides the number shown on the branch
    
    # Draw tree
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Manually recolor the labels
    for clade in tree.get_terminals():
        for text in ax.texts:
            if text.get_text()[1:] == clade.name:
                try:
                    color = to_rgb(normalize_branch_color(clade.color))  # Convert to RGB tuple
                    text.set_color(color)
                except ValueError:
                        print(f"Invalid color: {clade.color}")

    # Save the figure
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    plt.tight_layout()
    
    
def normalize_metric_colors(data_dict, cmap_name="viridis"):
    """Normalize values to colors using a colormap."""
    values = list(data_dict.values())
    norm = mcolors.Normalize(vmin=min(values), vmax=max(values))
    cmap = cm.get_cmap(cmap_name)
    return {k.replace("_", " "): cmap(norm(v)) for k, v in data_dict.items()}

def plot_metric_column(ax, metric_dict, colormap, title, y_positions, leaf_names):
    """Draw a column with colored bars and values, using a strikethrough line for NaN values."""
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_yticks(y_positions)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)

    for i, name in enumerate(leaf_names):
        val = metric_dict.get(name, None)
        color = colormap.get(name, (1, 1, 1, 0))
        ax.barh(y_positions[i], 1, color=color, edgecolor='none')

        if val is None or (isinstance(val, float) and np.isnan(val)):
            # Draw a horizontal line to indicate "deletion" (NaN)
            ax.plot([0.1, 0.9], [y_positions[i]] * 2, color="black", linewidth=1.5)
        else:
            display_text = f"≥{val}" if val == 200 else str(val)
            ax.text(0.5, y_positions[i], display_text, ha='center', va='center', fontsize=7, color='black')


def plot_combined_tree_with_metrics(tree_file, output_path, metrics, metric_titles, cmap_name="viridis"):
    """Creates a single tree plot with multiple value-annotated sidebars."""
    tree = Phylo.read(tree_file, "newick")
    leaf_names = [leaf.name.replace("_", " ") for leaf in tree.get_terminals()]
    num_leaves = len(leaf_names)

    scale_factor = max(10, len(tree.get_terminals()) * 0.25)  # auto-scale height
    fig = plt.figure(figsize=(18, scale_factor))
    gs = gridspec.GridSpec(1, len(metrics) + 1, width_ratios=[3] + [0.5] * len(metrics), wspace=0.05)

    # Clear internal node names (non-terminal)
    for clade in tree.get_nonterminals():
        clade.name = None
        clade.confidence = None  # This hides the number shown on the branch
    
    # Tree
    ax_tree = fig.add_subplot(gs[0, 0])
    Phylo.draw(tree, axes=ax_tree, do_show=False)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])

    # Metric columns
    for i, (metric_dict, title) in enumerate(zip(metrics, metric_titles)):
        ax = fig.add_subplot(gs[0, i + 1])
        color_map = normalize_metric_colors(metric_dict, cmap_name)
        plot_metric_column(ax, metric_dict, color_map, title, list(range(num_leaves)), leaf_names)

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    plt.tight_layout()
    print(f"✅ Combined tree saved: {output_path}")

# Save figures
output_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/summary_plots"
os.makedirs(output_dir, exist_ok=True)

# Load phylogenetic trees
MultiChromosomeGenesTree = Phylo.read(TREE_FILE, "newick")
MaxChromosomeGenesTree = Phylo.read(TREE_FILE, "newick")
ClusterLengthGenesTree = Phylo.read(TREE_FILE, "newick")
plot_combined_tree_with_metrics(
    tree_file=TREE_FILE,
    output_path=os.path.join(output_dir, "CombinedClusterMetrics.png"),
    metrics=[gene_counts, chromosome_lines, chromosome_cluster_length_genes],
    metric_titles=["Genes in the Genome", "Single-Chr Genes", "Cluster Length (Genes)"],
    cmap_name="viridis"
)
plot_tree_discrete(MultiChromosomeGenesTree, gene_counts, os.path.join(output_dir, "NumberOfHitsInGenome.png"))
plot_tree_discrete(MaxChromosomeGenesTree, chromosome_lines, os.path.join(output_dir, "SingleChromosomeGenes.png"))
plot_tree(ClusterLengthGenesTree, chromosome_cluster_length_genes, os.path.join(output_dir, "ClusterLengthGenes.png"))

print("✅ Phylogenetic trees saved successfully with separate metrics and scaled colors.")
