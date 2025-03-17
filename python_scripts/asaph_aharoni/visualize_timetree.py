from ete3 import NCBITaxa, Tree
from ete3.treeview import TreeStyle, NodeStyle
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

# Load species list from CSV
df = pd.read_csv("/groups/itay_mayrose/alongonda/desktop/dataset_organism_mapping.csv")
species_list = df['Organism'].drop_duplicates().tolist()

# Initialize NCBI Taxonomy database
ncbi = NCBITaxa()

# Get NCBI taxon IDs for each species
taxid_dict = {}
for species in species_list:
    taxid = ncbi.get_name_translator([species])
    if species in taxid:
        taxid_dict[str(taxid[species][0])] = species  # Store as string (TaxID -> Species)
        print(f"✔ Found NCBI TaxID for {species}: {taxid[species][0]}")
    else:
        print(f"❌ No TaxID found for {species}")

# Get the taxonomic tree
tree = ncbi.get_topology(list(map(int, taxid_dict.keys())))  # Convert TaxIDs to integers

# Save the tree in Newick format
tree.write(outfile="ncbi_species_tree.nwk", format=1)
print("\n✅ NCBI Taxonomic Tree saved as ncbi_species_tree.nwk")

# Load the tree with Biopython
tree_file = "ncbi_species_tree.nwk"
tree = Phylo.read(tree_file, "newick")

# Replace TaxIDs with species names dynamically
for clade in tree.get_terminals():
    if clade.name in taxid_dict:  # Check if TaxID exists in mapping
        clade.name = taxid_dict[clade.name]  # Replace TaxID with species name

# Plot the tree
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=ax)

# Save the image
plt.savefig("phylogenetic_tree_named.png", dpi=300, bbox_inches="tight")
plt.show()

print("\n✅ Phylogenetic tree saved with species names as 'phylogenetic_tree_named.png'")
