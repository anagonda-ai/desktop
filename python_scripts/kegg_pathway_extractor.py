import os
import requests
import csv
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

# Directory to store the data
root_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_test2"
os.makedirs(root_folder, exist_ok=True)

# Define the User-Agent header
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
}

# KEGG API limit: ~3 requests per second
BATCH_SIZE = 10  # Number of genes per batch request

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
        
# Export the plants dictionary to a CSV file
plants_csv_file = os.path.join(root_folder, "plants_list.csv")
with open(plants_csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Organism_Code", "Organism_Name"])
    for org_code, org_name in plants.items():
        writer.writerow([org_code, org_name])

print(f"Found {len(plants)} plant organisms in KEGG.")

def parse_kegg_entry(entry_text):
    parsed = {}
    current_key = None
    for line in entry_text.splitlines():
        if not line.strip():
            continue
        if line[:12].strip():  # New field
            current_key = line[:12].strip()
            value = line[12:].strip() if not current_key == "AASEQ" else "" 
            parsed[current_key] = value
        else:  # Continuation of previous field
            if current_key:
                parsed[current_key] += line[12:].strip()

    return parsed

# Function to fetch nucleotide sequences for multiple genes at once
def fetch_nucleotide_sequences(gene_batch):
    try:
        gene_query = "+".join(gene_batch)
        url = f"https://rest.kegg.jp/get/{gene_query}"
        response = requests.get(url, headers=headers).text.strip()

        gene_info = {}
        if response.startswith("ENTRY"):
            entries = response.split("///")[:BATCH_SIZE]
            for entry in entries:
                parsed = parse_kegg_entry(entry)
                ORTHOLOGY = parsed.get("ORTHOLOGY", "")
                ENTRY = parsed.get("ENTRY", "").split(" ")[0]
                POSITION = parsed.get("POSITION", "")
                AASEQ = parsed.get("AASEQ", "")
                gene_info[ENTRY] = (ORTHOLOGY, AASEQ, POSITION)

        return gene_info  # Dictionary: gene_id -> (annotation, sequence)
    except Exception as e:
        print(f"Error fetching nucleotide sequences: {e}")
        return {gene: ("N/A", "N/A") for gene in gene_batch}


# Function to process a single pathway
def process_pathway(org_folder, org_code, pathway_line):
    if not pathway_line or "\t" not in pathway_line:
        return  # Skip invalid lines
    if "metabol" not in pathway_line.lower():
        print(f"  ⚠️ Skipping non-metabolic pathway: {pathway_line}")
        return  # Skip non-metabolic pathways
    
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
    pathway_file = os.path.join(org_folder, f"{pathway_id}.csv")
    # Create CSV file with header
    with open(pathway_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene_ID", "Annotation", "Nucleotide_Sequence", "Location"])

    # Process genes in batches
    with ThreadPoolExecutor() as executor:
        batch_futures = []
        for i in range(0, len(genes), BATCH_SIZE):
            gene_batch = genes[i : i + BATCH_SIZE]
            batch_futures.append(executor.submit(fetch_nucleotide_sequences, gene_batch))

        for future in as_completed(batch_futures):
            gene_info = future.result()
            with open(pathway_file, "a", newline="") as f:
                writer = csv.writer(f)
                for gene, (annotation, nt_seq, location) in gene_info.items():
                    writer.writerow([gene, annotation, nt_seq, location])


            sleep(0.35)  # Respect KEGG rate limit (~3 requests/sec)

# Iterate over each plant organism
for org_code, org_name in plants.items():
    org_folder = os.path.join(root_folder, org_code)
    
    # Skip if the organism folder already exists and contains pathway files
    if os.path.exists(org_folder) and any(fname.endswith(".txt") for fname in os.listdir(org_folder)):
        print(f"⏩ Skipping already processed plant: {org_code} - {org_name}")
        continue

    # Create a directory for the plant
    print(f"Processing plant: {org_code} - {org_name}")
    os.makedirs(org_folder, exist_ok=True)

    # Get all pathways for the plant
    pathways_url = f"https://rest.kegg.jp/list/pathway/{org_code}"
    pathways_data = requests.get(pathways_url, headers=headers).text.strip().split("\n")

    # Process each pathway
    for pathway_line in pathways_data:
        process_pathway(org_folder, org_code, pathway_line)

    print(f"  Done processing {org_code}")

print("✅ KEGG plant data collection complete!")
