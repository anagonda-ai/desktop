import os
import shutil
import requests
import csv
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup
root_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules"
metabolic_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic"
os.makedirs(metabolic_folder, exist_ok=True)

headers = {
    'User-Agent': 'Mozilla/5.0'
}

BATCH_SIZE = 10  # Genes per request
SLEEP_TIME = 0.35  # Respect KEGG API rate limit

metabolic_dict = {}


def is_metabolic_module(module_id):
    if module_id in metabolic_dict.keys():
        return metabolic_dict[module_id]
    try:
        url = f"https://rest.kegg.jp/get/{module_id}"
        response = requests.get(url, headers=headers).text
        is_metabolic = any("metabol" in line for line in response.splitlines() if line.startswith("CLASS"))
        metabolic_dict[module_id] = is_metabolic
        return is_metabolic
    except Exception as e:
        print(f"⚠️ Failed to fetch module {module_id}: {e}")
        metabolic_dict[module_id] = False
        return False

# Process one organism
def process_file(csv_path):
    sbudir = os.path.dirname(csv_path).replace("KEGG_annotations_modules", "KEGG_annotations_modules_metabolic")
    os.makedirs(sbudir, exist_ok=True)
    
    filename = os.path.basename(csv_path)
    module = filename.split("_")[1].replace(".csv", "")
    # Build set of metabolic modules
    if is_metabolic_module(module):
        shutil.copy(csv_path, csv_path.replace("KEGG_annotations_modules", "KEGG_annotations_modules_metabolic"))
        print(f"✅ Metabolic module found: {module}")


# Use ThreadPoolExecutor for concurrent processing
csv_files = []
for subdir, dirs, files in os.walk(root_folder):
    for file in files:
        if file.endswith(".csv") and not file.startswith("plants_list"):
            csv_files.append(os.path.join(subdir, file))

with ThreadPoolExecutor(max_workers=10) as executor:
    futures = {executor.submit(process_file, csv_path): csv_path for csv_path in csv_files}
    for future in as_completed(futures):
        try:
            future.result()
        except Exception as exc:
            print(f"⚠️ Exception occurred while processing {futures[future]}: {exc}")

print("✅ All metabolic module-based KEGG data collected.")
