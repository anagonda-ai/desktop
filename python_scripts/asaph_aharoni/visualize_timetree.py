import os
import pandas as pd
import matplotlib.pyplot as plt
from ete3 import NCBITaxa, Tree
from Bio import Phylo
from multiprocessing import Pool, cpu_count
import requests

# **Force ETE3 to use offscreen rendering (for headless systems)**
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# **Initialize NCBI Taxonomy database**
ncbi = NCBITaxa()

# **Load species list from CSV**
df = pd.read_csv("/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/dataset_organism_mapping.csv")
species_list = df['Organism'].drop_duplicates().tolist()
print(species_list)

# **Step 1: Retrieve NCBI TaxIDs using Multiprocessing**
def get_taxid(species):
    """Retrieve NCBI TaxID for a given species."""
    local_ncbi = NCBITaxa()  # Create a local connection for each process
    try:
        taxid = local_ncbi.get_name_translator([species])
        if species in taxid:
            return str(taxid[species][0]), species  # Store as (TaxID, Species)
    except Exception as e:
        print(f"❌ Error fetching TaxID for {species}: {e}")
    return None

with Pool(processes=cpu_count()) as pool:
    results = pool.map(get_taxid, species_list)

# **Filter out failed lookups**
filtered_results = [result for result in results if result is not None]
# **Create a dictionary of TaxIDs and species names**
taxid_dict = {taxid: species for taxid, species in filtered_results if taxid}

print(f"\n✅ Retrieved {len(taxid_dict)} valid NCBI TaxIDs.")

# **Step 2: Get taxonomic tree topology**
tree = ncbi.get_topology(list(map(int, taxid_dict.keys())))

# **Step 3: Fetch taxonomic ranks using Multiprocessing**
def fetch_taxonomic_name(taxid):
    """Retrieve the highest available taxonomic rank for an internal node."""
    local_ncbi = NCBITaxa()  # Use separate instance
    try:
        lineage = local_ncbi.get_lineage(int(taxid))
        ranks = local_ncbi.get_rank(lineage)
        names = local_ncbi.get_taxid_translator(lineage)

        for rank in ["family", "genus", "order"]:  # Prefer higher ranks
            for lin_taxid in lineage:
                if ranks[lin_taxid] == rank:
                    return taxid, names[lin_taxid]
        return taxid, names[lineage[-1]]  # Default to last known rank
    except Exception:
        return taxid, None

internal_taxids = [node.name for node in tree.traverse() if not node.is_leaf()]
with Pool(processes=cpu_count()) as pool:
    results = pool.map(fetch_taxonomic_name, internal_taxids)

# **Create mapping for internal nodes**
internal_names = {taxid: name for taxid, name in results if name}

# **Step 4: Replace node names in tree**
for node in tree.traverse():
    if node.is_leaf():
        if node.name in taxid_dict:
            node.name = taxid_dict[node.name].replace(" ", "_")  # Replace TaxID with species name
    else:
        if node.name in internal_names:
            node.name = internal_names[node.name].replace(" ", "_")  # Replace internal nodes

# **Step 5: Retrieve Evolutionary Distances from NCBI**
def fetch_divergence_time(pair):
    """Retrieve evolutionary divergence time for a species pair."""
    species_1, species_2 = pair
    try:
        url = f"https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi?mode=All&taxid1={species_1}&taxid2={species_2}"
        response = requests.get(url).text
        
        # Extract time estimate
        if "Divergence Time: " in response:
            time_index = response.find("Divergence Time: ") + len("Divergence Time: ")
            divergence_time = response[time_index:].split(" ")[0]
            return (species_1, species_2, float(divergence_time))  # Return tuple
    except Exception:
        return None  # Ensure it returns None to avoid errors

# **Generate species pairs**
species_pairs = [(s1, s2) for i, s1 in enumerate(taxid_dict.keys()) for s2 in list(taxid_dict.keys())[i+1:]]

# **Retrieve divergence times in parallel**
with Pool(processes=cpu_count()) as pool:
    divergence_results = pool.map(fetch_divergence_time, species_pairs)

# **Filter out None results**
divergence_times = {(s1, s2): time for result in divergence_results if result is not None for s1, s2, time in [result]}

print(f"\n✅ Retrieved divergence times for {len(divergence_times)} species pairs.")

# **Step 6: Save tree in Newick format**
output_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/"
os.makedirs(output_dir, exist_ok=True)

tree_path = os.path.join(output_dir, "ncbi_species_tree_named.nwk")
tree.write(outfile=tree_path, format=1)
print("\n✅ NCBI Taxonomic Tree saved as 'ncbi_species_tree_named.nwk'")

# **Step 7: Load and visualize tree**
tree_path = os.path.join(output_dir, "species.nwk")
tree = Phylo.read(tree_path, "newick")

# **Step 9: Scale branch lengths by divergence times**
for clade in tree.get_nonterminals():
    for (s1, s2), time in divergence_times.items():
        if s1 in clade.name and s2 in clade.name:
            clade.branch_length = time  # Assign evolutionary distance

# **Step 10: Create a high-quality circular phylogenetic tree**
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=ax)

# **Add Evolutionary Distance Scale**
ax.set_xlabel("Evolutionary Distance (Millions of Years)")

# **Save visualization**
fih_path = os.path.join(output_dir, "phylogenetic_tree_circular.png")
plt.savefig(fih_path, dpi=300, bbox_inches="tight")
plt.show()

print("\n✅ Phylogenetic tree saved as 'phylogenetic_tree_circular.png'")
