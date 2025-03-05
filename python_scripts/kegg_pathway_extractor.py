import os
import requests
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

# Directory to store the data
root_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG"
os.makedirs(root_folder, exist_ok=True)

# Define the User-Agent header
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
}

# KEGG API limit: ~3 requests per second
BATCH_SIZE = 20  # Number of genes per batch request

# Fetch the list of all KEGG organisms
print("Fetching KEGG organism list...")
organisms_url = "https://rest.kegg.jp/list/organism"
organisms_data = requests.get(organisms_url, headers=headers).text.strip().split("\n")

# Filter only plant organisms (those with "Plants" in the category column)
plants = {}
for line in organisms_data:
    parts = line.split("\t")
    if len(parts) >= 4 and "Plants" in parts[3]:
        org_code, org_name = parts[1], parts[2]
        plants[org_code] = org_name

print(f"Found {len(plants)} plant organisms in KEGG.")

# Function to fetch nucleotide sequences for multiple genes at once
def fetch_nucleotide_sequences(gene_batch):
    try:
        gene_query = "+".join(gene_batch)  # Batch multiple genes in one request
        url = f"https://rest.kegg.jp/get/{gene_query}/ntseq"
        response = requests.get(url, headers=headers).text.strip()

        gene_sequences = {}
        if response.startswith(">"):
            genes_data = response.split("\n>")
            for entry in genes_data:
                lines = entry.split("\n")
                gene_id = lines[0].split()[0].replace(">", "").strip()
                sequence = "".join(lines[1:])
                gene_sequences[gene_id] = sequence

        return gene_sequences  # Dictionary of gene -> sequence
    except Exception as e:
        print(f"Error fetching nucleotide sequences: {e}")
        return {gene: "N/A" for gene in gene_batch}  # Return empty for failed requests

# Function to process a single pathway
def process_pathway(org_folder, org_code, pathway_line):
    if not pathway_line or "\t" not in pathway_line:
        return  # Skip invalid lines

    pathway_id, pathway_name = pathway_line.split("\t")

    # Extract only the pathway ID
    if ":" in pathway_id:
        pathway_id = pathway_id.split(":")[1]

    print(f"  Fetching genes for pathway: {pathway_id} - {pathway_name}")

    # Get genes for the pathway
    genes_url = f"https://rest.kegg.jp/link/{org_code}/{pathway_id}"
    genes_data = requests.get(genes_url, headers=headers).text.strip()
    if not genes_data:
        print(f"  ⚠️ No genes found for {pathway_id}. Skipping...")
        return

    genes = [line.split("\t")[1] for line in genes_data.split("\n") if "\t" in line]

    if not genes:
        return  # Skip empty pathways

    # Create a file for the pathway
    pathway_file = os.path.join(org_folder, f"{pathway_id}.txt")
    with open(pathway_file, "w") as f:
        f.write(f"# Pathway: {pathway_id} - {pathway_name}\n\n")
        f.write("# Gene_ID\tNucleotide_Sequence\n")

    # Process genes in batches
    with ThreadPoolExecutor() as executor:
        batch_futures = []
        for i in range(0, len(genes), BATCH_SIZE):
            gene_batch = genes[i : i + BATCH_SIZE]
            batch_futures.append(executor.submit(fetch_nucleotide_sequences, gene_batch))

        for future in as_completed(batch_futures):
            gene_sequences = future.result()
            with open(pathway_file, "a") as f:
                for gene, nt_seq in gene_sequences.items():
                    f.write(f"{gene}\t{nt_seq}\n")

            sleep(0.35)  # Respect KEGG rate limit (~3 requests/sec)

# Iterate over each plant organism
for org_code, org_name in plants.items():
    print(f"Processing plant: {org_code} - {org_name}")

    # Create a directory for the plant
    org_folder = os.path.join(root_folder, org_code)
    os.makedirs(org_folder, exist_ok=True)

    # Get all pathways for the plant
    pathways_url = f"https://rest.kegg.jp/list/pathway/{org_code}"
    pathways_data = requests.get(pathways_url, headers=headers).text.strip().split("\n")

    # Process each pathway
    for pathway_line in pathways_data:
        process_pathway(org_folder, org_code, pathway_line)

    print(f"  Done processing {org_code}")

print("✅ KEGG plant data collection complete!")
