from opentree import OT
import pandas as pd
import networkx as nx
from Bio import Phylo
import matplotlib.pyplot as plt

# Load species list from CSV file
df = pd.read_csv("/groups/itay_mayrose/alongonda/desktop/dataset_organism_mapping.csv")
species_list = df['Organism'].drop_duplicates().tolist()

# Step 1: Get OTT IDs using OpenTree's Taxonomy
ott_ids = []
for species in species_list:
    response = OT.tnrs_match(names=[species])  # Call the API
    response_dict = response.response_dict  # Extract response as a dictionary
    
    if "results" in response_dict and response_dict["results"]:
        match = response_dict["results"][0]["matches"]
        
        if match and isinstance(match, list):
            taxon_info = match[0].get("taxon", {})  # Extract "taxon" safely
            ott_id = taxon_info.get("ott_id")  # Get OTT ID safely
            
            if ott_id:
                ott_ids.append(str(ott_id))
                print(f"✔ Found OTT ID for {species}: {ott_id}")
            else:
                print(f"❌ No OTT ID found for {species}")
        else:
            print(f"❌ No match found for {species}")
    else:
        print(f"❌ No results for {species}")

print("\n✅ Valid OTT IDs:", ott_ids)

# Step 2: Get the Most Recent Common Ancestor (MRCA) for these species
if ott_ids:
    try:
        mrca = OT.taxon_mrca(ott_ids=ott_ids)  # Get MRCA taxon ID
        mrca_dict = mrca.response_dict

        if "mrca" in mrca_dict:
            mrca_ott_id = mrca_dict["mrca"]["ott_id"]
            print(f"\n✅ MRCA OTT ID: {mrca_ott_id}")

            # Step 3: Get the taxonomic subtree from MRCA
            tree = OT.taxon_subtree(ott_id=ott_ids)
            tree_dict = tree.response_dict

            # Save the tree in Newick format
            with open("taxonomic_species_tree.nwk", "w") as f:
                f.write(tree_dict.get("newick", ""))

            print("\n✅ Taxonomic Tree saved as taxonomic_species_tree.nwk")
        else:
            print("\n❌ Could not determine MRCA.")

    except Exception as e:
        print(f"\n❌ Error generating taxonomic tree: {e}")
else:
    print("\n❌ No valid OTT IDs found. Tree cannot be generated.")

# Load the Newick tree

# Load tree
tree_file = "taxonomic_species_tree.nwk"
tree = Phylo.read(tree_file, "newick")

# Convert tree to NetworkX graph
graph = Phylo.to_networkx(tree)

# Get node positions using a circular layout
pos = nx.circular_layout(graph)

# Create figure
plt.figure(figsize=(20, 20))

# Draw tree edges
nx.draw_networkx_edges(graph, pos, alpha=0.5, width=0.6)

# Draw labels for species (only terminal nodes)
labels = {n: n.name for n in tree.get_terminals() if n.name}
nx.draw_networkx_labels(graph, pos, labels, font_size=5, font_family="Arial", verticalalignment="center")

# Save and show
plt.savefig("phylogenetic_tree_circular.png", dpi=300, bbox_inches="tight")
plt.show()
