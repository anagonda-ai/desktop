from ete3 import NCBITaxa, Tree
import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from ete3.treeview import TreeStyle, NodeStyle
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

# **Force ETE3 to use offscreen rendering (for headless systems)**
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# **Initialize NCBI Taxonomy database**
ncbi = NCBITaxa()

# **Load species list from CSV**
df = pd.read_csv("/groups/itay_mayrose/alongonda/desktop/dataset_organism_mapping.csv")
species_list = df['Organism'].drop_duplicates().tolist()[:10] # Select top 10 species


# Step 1: Retrieve NCBI TaxIDs for all species
taxid_dict = {}
for species in species_list:
    taxid = ncbi.get_name_translator([species])
    if species in taxid:
        taxid_dict[str(taxid[species][0])] = species  # Store as string (TaxID -> Species)
        print(f"✔ Found NCBI TaxID for {species}: {taxid[species][0]}")
    else:
        print(f"❌ No TaxID found for {species}")

# Step 2: Get taxonomic tree topology
tree = ncbi.get_topology(list(map(int, taxid_dict.keys())))  # Convert TaxIDs to integers

# Step 3: Convert internal nodes (higher taxa) to readable names
def rename_internal_nodes(node):
    if node.is_leaf():
        return  # Skip terminal nodes (already labeled)
    
    try:
        lineage = ncbi.get_lineage(int(node.name))
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        
        # Prefer Family > Genus > Order for internal nodes
        for rank in ["family", "genus", "order"]:
            for taxid in lineage:
                if ranks[taxid] == rank:
                    node.name = names[taxid]
                    return  # Stop at the first match
        # Default to scientific name if no rank match
        node.name = names[lineage[-1]]
    
    except Exception as e:
        print(f"⚠️ Error retrieving taxonomic name for node {node.name}: {e}")

# Apply renaming to all internal nodes
for node in tree.traverse():
    rename_internal_nodes(node)

# Step 4: Save tree in Newick format
tree.write(outfile="ncbi_species_tree_named.nwk", format=1)
print("\n✅ NCBI Taxonomic Tree saved as 'ncbi_species_tree_named.nwk'")

# Step 5: Load and visualize tree
tree_file = "ncbi_species_tree_named.nwk"
tree = Phylo.read(tree_file, "newick")

# Step 6: Replace TaxIDs with species names in Biopython tree
for clade in tree.find_clades():
    if clade.name in taxid_dict:  # Check if it's a terminal node
        clade.name = taxid_dict[clade.name]  # Replace TaxID with species name

# Step 7: Create a high-quality circular phylogenetic tree
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=ax)

# Save visualization
plt.savefig("phylogenetic_tree_circular.png", dpi=300, bbox_inches="tight")
plt.show()

print("\n✅ Phylogenetic tree saved as 'phylogenetic_tree_circular.png'")
