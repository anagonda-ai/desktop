import os
import subprocess
import requests
import csv
import pandas as pd

BASE_URL = "https://rest.kegg.jp"

# Setup
root_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic_cluster"
os.makedirs(root_folder, exist_ok=True)

headers = {
    'User-Agent': 'Mozilla/5.0'
}

plants_csv = os.path.join(root_folder, "plants_list.csv")
if not os.path.exists(plants_csv):
    # Get KEGG plant organisms
    print("Fetching KEGG organism list...")
    organisms_url = f"{BASE_URL}/list/organism"
    organisms_data = requests.get(organisms_url, headers=headers).text.strip().split("\n")

    plants = {}
    for line in organisms_data:
        parts = line.split("\t")
        if len(parts) >= 4 and "Plants" in parts[3]:
            org_code, org_name = parts[1], parts[2]
            plants[org_code] = org_name

    with open(plants_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Organism_Code", "Organism_Name"])
        for code, name in plants.items():
            writer.writerow([code, name])

    print(f"Found {len(plants)} plant organisms.")
else:
    # Load existing plant organisms
    print("Loading existing KEGG organism list...")
    df = pd.read_csv(plants_csv)
    plants = dict(zip(df["Organism_Code"], df["Organism_Name"]))
    
    
metabolic_modules_path = os.path.join(root_folder, "metabolic_modules.csv")
if not os.path.exists(metabolic_modules_path):
    # Step 1: Get all modules
    print("Fetching all KEGG modules...")
    module_lines = requests.get(f"{BASE_URL}/list/module", headers=headers).text.strip().splitlines()
    module_ids = [line.split("\t")[0] for line in module_lines]

    print(f"Total modules found: {len(module_ids)}")

    # Step 2: Check which are metabolic
    metabolic_modules = []

    import concurrent.futures

    def check_metabolic(module_id):
        try:
            data = requests.get(f"{BASE_URL}/get/{module_id}", headers=headers).text
            lines = data.splitlines()
            module_class = ""
            for line in lines:
                if line.startswith("CLASS"):
                    module_class = line[12:].strip()
            if "metabolism" in module_class:
                return module_id
        except Exception as e:
            print(f"⚠️ Error with {module_id}: {e}")
        return None

    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(check_metabolic, module_ids))

    metabolic_modules = [mid for mid in results if mid]
    # Step 3: Save to CSV
    df = pd.DataFrame({"Module_ID": metabolic_modules})
    df.to_csv(metabolic_modules_path, index=False)

    print(f"\n✅ Saved {len(metabolic_modules)} metabolic modules to '{metabolic_modules_path}'")

# Run for all plants
for code, name in plants.items():
    cmd = f"sbatch /groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/kegg_organisms_download.sh {code} {root_folder} {metabolic_modules_path}"
    print(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {code}: {e}")
        continue

print("✅ All metabolic module-based KEGG data collected.")
