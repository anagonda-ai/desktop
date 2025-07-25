import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import gridspec, cm
from Bio import Phylo
import numpy as np
from matplotlib.colors import to_rgb
import ast
import argparse

parser = argparse.ArgumentParser(description="Find homolog genes for Asaph cluster")
parser.add_argument("--example_mgc", type=str, required=True, help="Path to example MGC directory")
args = parser.parse_args()

example_mgc = args.example_mgc

# Define file paths
MAPPING_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping.csv"
TREE_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/species.nwk"
COMPARISON_FILE = os.path.join(example_mgc, "comparison_results.csv")

def parse_dict(val):
    if isinstance(val, dict):
        return val
    if isinstance(val, str):
        try:
            return ast.literal_eval(val)
        except Exception:
            return {}
    return {}

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
        return directory.split(os.path.join(example_mgc, "blast_results_chromosome_separated/best_hits_by_organism/"))[1].split("/")[0]

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
    
    # Apply the function to calculate the largest chromosome length and index difference
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
    chromosome_cluster_length = dict(zip(comparison_df['Organism'].replace(" ", "_"), (comparison_df['Largest Chromosome Length'] / 1000).round())) # Convert to Kbp
    chromosome_cluster_gene_existance = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['MGC Genes Order'].apply(parse_dict)))
    chromosome_cluster_gene_strands = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['MGC Genes Strands'].apply(parse_dict)))

    max_value = 100  # Define the cap for Index Difference
    # Replace values in Index Difference column
    
    comparison_df['Index Difference'] = comparison_df['Index Difference'].apply(lambda x: min(x, max_value) if x is not None else None)
    
    # Create the dictionary after applying the cap
    chromosome_cluster_length_genes = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Index Difference']))
else:
    gene_counts = chromosome_lines = chromosome_cluster_length = chromosome_cluster_length_genes = chromosome_cluster_gene_existance = chromosome_cluster_gene_strands = {}

# Improved Color Function
def get_color(value, min_val, max_val, cmap_name="viridis"):
    """Returns a color from a matplotlib colormap based on normalized value."""
    norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
    cmap = plt.colormaps[cmap_name]
    
    if value is not None:
        rgba = cmap(norm(value))
    else:
        rgba = (1, 1, 1, 0)  # Default to transparent white if value is None
    return mcolors.to_hex(rgba)

def label_function_factory(data_dict,  cmap_name="viridis"):
    """Creates a function for labeling tree clades with colored values."""
    values = [v for v in data_dict.values() if isinstance(v, (int, float)) and not pd.isna(v)]
    min_val, max_val = (min(values), max(values)) if values else (0, 1)

    def label_function(clade):
        label = clade.name.replace("_", " ").lower() if clade.name else None
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
    sm = cm.ScalarMappable(cmap=plt.colormaps[cmap_name], norm=mcolors.Normalize(vmin=min_val, vmax=max_val))
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
        label = clade.name.replace("_", " ").lower() if clade.name else None
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
    
    filtered_values = [v for v in values if v is not None]
    
    if not filtered_values:
        norm = mcolors.Normalize(vmin=0, vmax=1)  # Default range if no valid values
    else:
        norm = mcolors.Normalize(vmin=min(filtered_values), vmax=max(filtered_values))
    cmap = plt.colormaps[cmap_name]
    
    return {k.replace("_", " "): cmap(norm(v)) for k, v in data_dict.items() if v is not None}

def plot_metric_column(ax, metric_dict, colormap, title, y_positions, leaf_names):
    """Draw a column with colored bars and values, using a strikethrough line for NaN values."""
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Set vertical limits to fully cover the tree height
    ax.set_ylim(min(y_positions) - 1, max(y_positions) + 1)

    ax.set_title(title, fontsize=10)

    for i, name in enumerate(leaf_names):
        y = y_positions[i]
        val = metric_dict.get(name, None)
        color = colormap.get(name, (1, 1, 1, 0))
        ax.barh(y, 1, color=color, edgecolor='none', height=1.0)

        if val is None or (isinstance(val, float) and np.isnan(val)):
            # Draw a horizontal line to indicate "deletion" (NaN)
            ax.plot([0.1, 0.9], [y] * 2, color="black", linewidth=1.5)
        else:
            display_text = f"≥{val}" if val >= 100 else str(val)
            ax.text(0.5, y, display_text, ha='center', va='center', fontsize=7, color='black')



def plot_combined_tree_with_metrics(tree_file, output_tree_path, output_data_path, metrics, metric_titles, chromosome_cluster_gene_existance, chromosome_cluster_gene_strands, cmap_name="viridis"):
    """Creates two plots: species tree and a value-annotated data matrix aligned with leaves."""
    tree = Phylo.read(tree_file, "newick")

    # Clear internal node labels
    for clade in tree.get_nonterminals():
        clade.name = None
        clade.confidence = None

    for clade in tree.get_terminals():
        clade.name = clade.name.lower()

    # Get leaf names in the drawn order
    all_gene_keys = set()
    for v in chromosome_cluster_gene_existance.values():
        if isinstance(v, dict):
            all_gene_keys.update(v.keys())
    all_gene_keys = sorted(all_gene_keys)

    if "haaap_transporters" in all_gene_keys:
        all_gene_keys.remove("haaap_transporters")

    total_columns = 1 + len(metrics) + len(all_gene_keys)
    n_leaves = len(tree.get_terminals())

    fig_height = max(20, n_leaves * 0.25)

    # --- Tree Figure ---
    fig_tree = plt.figure(figsize=(6, fig_height))
    ax_tree = fig_tree.add_subplot(111)
    label_positions = {}

    def custom_label_func(clade):
        if clade.is_terminal():
            label_positions[clade.name.replace("_", " ")] = len(label_positions)
            return clade.name.replace("_", " ")
        return None

    Phylo.draw(tree, do_show=False, axes=ax_tree, label_func=custom_label_func)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])
    fig_tree.savefig(output_tree_path, dpi=300, bbox_inches="tight")
    plt.close(fig_tree)

    # --- Data Matrix Figure ---
    leaf_names = sorted(label_positions, key=lambda name: label_positions[name], reverse=True)
    y_positions = np.linspace(0, len(leaf_names) - 1, len(leaf_names))

    fig_matrix = plt.figure(figsize=(max(5, 0.3 * total_columns), fig_height))
    gs = gridspec.GridSpec(1, total_columns, width_ratios=[0.4] * total_columns, wspace=0.05)

    for i, (metric_dict, title) in enumerate(zip(metrics, metric_titles)):
        ax = fig_matrix.add_subplot(gs[0, i])
        cmap = normalize_metric_colors(metric_dict, cmap_name)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=10)
        ax.set_ylim(len(leaf_names) - 0.5, -0.5)
        ax.invert_yaxis()

        for y_idx, name in enumerate(leaf_names):
            val = metric_dict.get(name, None)
            color = cmap.get(name, (1, 1, 1, 0))
            ax.barh(y_idx, 1, color=color, edgecolor='none')
            if val is None or (isinstance(val, float) and np.isnan(val)):
                ax.plot([0.1, 0.9], [y_idx] * 2, color="black", linewidth=1.5)
            else:
                display = f"≥{val}" if val == 200 else str(val)
                ax.text(0.5, y_idx, display, ha="center", va="center", fontsize=7)

    for gene_idx, gene_key in enumerate(all_gene_keys):
        ax = fig_matrix.add_subplot(gs[0, len(metrics) + gene_idx])
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(gene_key, fontsize=9)
        ax.set_ylim(len(leaf_names) - 0.5, -0.5)
        ax.invert_yaxis()
        ax.set_aspect('auto')

        for y_idx, name in enumerate(leaf_names):
            gene_dict = chromosome_cluster_gene_existance.get(name, {})
            gene_strand_dict = chromosome_cluster_gene_strands.get(name, {})
            val = gene_dict.get(gene_key, None) if isinstance(gene_dict, dict) else None
            strand = gene_strand_dict.get(gene_key, None) if isinstance(gene_strand_dict, dict) else None
            color_map = {1: "#4daf4a", 2: "#ffd92f", 3: "#ff7f00", 4: "#e41a1c"}
            color = color_map.get(val, (1, 1, 1, 0))
            ax.barh(y_idx, 1, color=color, edgecolor='none')
            if val is None:
                ax.plot([0.1, 0.9], [y_idx] * 2, color="black", linewidth=1.5)
            else:
                ax.text(0.5, y_idx, (val, strand) if val else "✗", ha="center", va="center", fontsize=10, color="black")

    fig_matrix.tight_layout()
    fig_matrix.savefig(output_data_path, dpi=300, bbox_inches="tight")
    plt.close(fig_matrix)
    print(f"✅ Species tree saved: {output_tree_path}")
    print(f"✅ Cluster metrics matrix saved: {output_data_path}")


# Save figures
output_dir = os.path.join(example_mgc, "summary_plots")
os.makedirs(output_dir, exist_ok=True)

# Load phylogenetic trees
MultiChromosomeGenesTree = Phylo.read(TREE_FILE, "newick")
MaxChromosomeGenesTree = Phylo.read(TREE_FILE, "newick")
ClusterLengthGenesTree = Phylo.read(TREE_FILE, "newick")

plot_combined_tree_with_metrics(
    tree_file=TREE_FILE,
    output_tree_path=os.path.join(output_dir, "SpeciesTree.png"),
    output_data_path=os.path.join(output_dir, "ClusterMetricsMatrix.png"),
    metrics=[gene_counts, chromosome_lines, chromosome_cluster_length_genes, chromosome_cluster_length],
    metric_titles=["Genome Hits", "Chr Hits", "Size (genes)", "Size (kbp)"],
    chromosome_cluster_gene_existance=chromosome_cluster_gene_existance,
    chromosome_cluster_gene_strands=chromosome_cluster_gene_strands,
    cmap_name="viridis"
)
plot_tree_discrete(MultiChromosomeGenesTree, gene_counts, os.path.join(output_dir, "NumberOfHitsInGenome.png"))
plot_tree_discrete(MaxChromosomeGenesTree, chromosome_lines, os.path.join(output_dir, "SingleChromosomeGenes.png"))
plot_tree(ClusterLengthGenesTree, chromosome_cluster_length_genes, os.path.join(output_dir, "ClusterLengthGenes.png"))

print("✅ Phylogenetic trees saved successfully with separate metrics and scaled colors.")