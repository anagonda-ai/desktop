import os
import requests
import csv
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
import pandas as pd

headers = {
    'User-Agent': 'Mozilla/5.0'
}

module_class_cache = {}

BATCH_SIZE = 10  # Genes per request
SLEEP_TIME = 0.35  # Respect KEGG API rate limit

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

def is_metabolic_module(module_id_list, module_id):
    return module_id in module_id_list

# Process one organism
def process_organism(org_code, root_folder, metabolic_modules_path):
    # Load list of Module_IDs from CSV (["Module_ID", "Name", "Class"])
    module_id_list = []
    df = pd.read_csv(metabolic_modules_path)
    module_id_list = df["Module_ID"].tolist()            
        
    org_folder = os.path.join(root_folder, org_code)
    if os.path.exists(org_folder) and any(fname.endswith(".csv") for fname in os.listdir(org_folder)):
        print(f"⏩ Skipping {org_code} (already done)")
        return

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
        general_module_id = module.split("_")[1]

        if is_metabolic_module(module_id_list, general_module_id):
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
    
# Main function
if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Usage: python process_organism.py <org_code> <root_folder> <metabolic_modules_path>")
        print(sys.argv)
        sys.exit(1)

    org_code = sys.argv[1]
    root_folder = sys.argv[2]
    metabolic_modules_path = sys.argv[3]

    process_organism(org_code, root_folder, metabolic_modules_path)