import os
import requests
import csv
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup
root_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic_local"
os.makedirs(root_folder, exist_ok=True)

headers = {
    'User-Agent': 'Mozilla/5.0'
}

BATCH_SIZE = 10  # Genes per request
SLEEP_TIME = 0.35  # Respect KEGG API rate limit

# Get KEGG plant organisms
print("Fetching KEGG organism list...")
organisms_url = "https://rest.kegg.jp/list/organism"
organisms_data = requests.get(organisms_url, headers=headers).text.strip().split("\n")

plants = {}
for line in organisms_data:
    parts = line.split("\t")
    if len(parts) >= 4 and "Plants" in parts[3]:
        org_code, org_name = parts[1], parts[2]
        plants[org_code] = org_name

plants_csv = os.path.join(root_folder, "plants_list.csv")
with open(plants_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Organism_Code", "Organism_Name"])
    for code, name in plants.items():
        writer.writerow([code, name])

print(f"Found {len(plants)} plant organisms.")

# Parse KEGG flat entry
def parse_kegg_entry(entry_text):
    parsed = {}
    current_key = None
    for line in entry_text.splitlines():
        if not line.strip():
            continue
        if line[:12].strip():
            current_key = line[:12].strip()
            parsed[current_key] = line[12:].strip() if current_key != "AASEQ" else ""
        else:
            if current_key:
                parsed[current_key] += line[12:].strip()
    return parsed

# Check if a module is metabolic by querying its metadata
module_class_cache = {}

def is_metabolic_module(module_id):
    general_module_id = module_id.split("_")[1]
    if general_module_id in module_class_cache:
        return module_class_cache[general_module_id]
    try:
        url = f"https://rest.kegg.jp/get/{general_module_id}"
        response = requests.get(url, headers=headers).text
        is_metabolic = any("metabol" in line for line in response.splitlines() if line.startswith("CLASS"))
        module_class_cache[general_module_id] = is_metabolic
        sleep(SLEEP_TIME)  # throttle per KEGG rules
        return is_metabolic
    except Exception as e:
        print(f"⚠️ Failed to fetch module {module_id}: {e}")
        module_class_cache[module_id] = False
        return False

# Batch gene fetch
def fetch_gene_info(gene_batch, org_code):
    try:
        query = "+".join([f"{org_code}:{gene}" for gene in gene_batch])
        url = f"https://rest.kegg.jp/get/{query}"
        response = requests.get(url, headers=headers).text.strip()
        
        gene_info = {}
        for entry in response.split("///"):
            if entry.strip():
                parsed = parse_kegg_entry(entry)
                gene_id = parsed.get("ENTRY", "").split()[0]
                gene_info[gene_id] = (
                    parsed.get("ORTHOLOGY", ""),
                    parsed.get("AASEQ", ""),
                    parsed.get("POSITION", "")
                )
        return gene_info
    except Exception as e:
        print(f"Error fetching gene info: {e}")
        return {g: ("N/A", "N/A", "N/A") for g in gene_batch}

# Process one organism
def process_organism(org_code, org_name):
    org_folder = os.path.join(root_folder, org_code)
    if os.path.exists(org_folder) and any(fname.endswith(".csv") for fname in os.listdir(org_folder)):
        print(f"⏩ Skipping {org_code} (already done)")
        return

    print(f"Processing {org_code} - {org_name}")
    os.makedirs(org_folder, exist_ok=True)

    module_url = f"https://rest.kegg.jp/link/module/{org_code}"
    module_data = requests.get(module_url, headers=headers).text.strip().split("\n")

    # Build dict: module → list of genes (only metabolic)
    module_dict = {}
    for line in module_data:
        if "\t" not in line:
            continue
        gene, module = line.split("\t")
        gene_id = gene.split(":")[1]
        module_id = module.split(":")[1]

        if is_metabolic_module(module_id):
            module_dict.setdefault(module_id, []).append(gene_id)

    for module_id, genes in module_dict.items():
        module_file = os.path.join(org_folder, f"{module_id}.csv")
        with open(module_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Gene_ID", "Annotation", "Nucleotide_Sequence", "Location"])

        with ThreadPoolExecutor() as executor:
            futures = []
            for i in range(0, len(genes), BATCH_SIZE):
                batch = genes[i:i + BATCH_SIZE]
                futures.append(executor.submit(fetch_gene_info, batch, org_code))

            for future in as_completed(futures):
                result = future.result()
                with open(module_file, "a", newline="") as f:
                    writer = csv.writer(f)
                    for gene, (annot, aaseq, pos) in result.items():
                        writer.writerow([gene, annot, aaseq, pos])
                sleep(SLEEP_TIME)

    print(f"✅ Done: {org_code}")

# Run for all plants
for code, name in plants.items():
    process_organism(code, name)

print("✅ All metabolic module-based KEGG data collected.")
