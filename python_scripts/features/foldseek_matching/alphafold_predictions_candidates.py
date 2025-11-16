import os
import requests
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import time
import threading
import shutil

# Configuration
API_SEMAPHORE = threading.Semaphore(3)
REQUEST_DELAY = 0.4
MAX_RETRIES = 3

LENGTH_THRESHOLD = 250  # adjustable

UNIPROT_SEARCH_API = "https://rest.uniprot.org/uniprotkb/search?query={}&fields=accession&format=json"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"

SOURCE_1 = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files"
SOURCE_2 = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/"

# ------------------------- Caching Functions -------------------------

def load_orf_cache(cache_file):
    cache = {}
    if cache_file.exists():
        with open(cache_file) as f:
            next(f)  # skip header
            for line in f:
                orf_name, fasta_stem, seq = line.strip().split("\t")
                cache[(orf_name, fasta_stem)] = seq
    return cache

def save_orf_cache(cache_file, cache_dict):
    with open(cache_file, "w") as f:
        f.write("orf_name\tfasta_stem\tsequence\n")
        for (orf, stem), seq in cache_dict.items():
            f.write(f"{orf}\t{stem}\t{seq}\n")

def update_orf_cache(fasta_file, cache_dict):
    source_path = str(fasta_file)
    index = 0
    if index is None:
        print(f"[!] Unknown source path: {fasta_file}")
        return {}

    new_entries = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        info = record.id if index == 0 else record.description
        parts = info.split(" | ")
        if len(parts) > index:
            orf_name = parts[index]
            key = (orf_name, fasta_file.stem)
            if key not in cache_dict:
                new_entries[key] = str(record.seq)
    return new_entries

def load_uniprot_cache(cache_file):
    cache = {}
    if cache_file.exists():
        with open(cache_file) as f:
            for line in f:
                orf, uid = line.strip().split("\t")
                cache[orf] = uid
    return cache

def save_uniprot_cache(cache_file, cache_dict):
    with open(cache_file, "w") as f:
        for orf, uid in cache_dict.items():
            f.write(f"{orf}\t{uid}\n")


# ------------------------- UniProt & AlphaFold -------------------------

def resolve_uniprot_id(orf_name):
    for attempt in range(1, MAX_RETRIES + 1):
        with API_SEMAPHORE:
            try:
                res = requests.get(UNIPROT_SEARCH_API.format(orf_name), timeout=10)
                res.raise_for_status()
                results = res.json().get("results", [])
                for entry in results:
                    xrefs = entry.get("uniProtKBCrossReferences", [])
                    if any(xref["database"] in ["Gramene", "EnsemblPlants", "Phytozome"] and orf_name in xref.get("id", "") for xref in xrefs):
                        print(f"[✔] Found {orf_name} in UniProt: {entry['primaryAccession']}")
                        return orf_name, entry["primaryAccession"]
                if results:
                    print(f"[!] No relevant cross-references found for {orf_name}, using first result: {results[0]['primaryAccession']}")
                    return orf_name, results[0]["primaryAccession"]
                break
            except Exception as e:
                print(f"[~] Retry {attempt} for {orf_name}: {e}")
                time.sleep(REQUEST_DELAY * 2)
            finally:
                time.sleep(REQUEST_DELAY)
    print(f"[!] Failed to resolve UniProt ID for {orf_name} after {MAX_RETRIES} attempts")
    return orf_name, None

def download_alphafold_pdb(uniprot_id, orf_name, subdir):
    if not uniprot_id:
        return f"[✖] Not found on AlphaFold: {orf_name}"
    try:
        subdir.mkdir(parents=True, exist_ok=True)
        out_path = subdir / f"{orf_name}.pdb"
        
        # Check if AlphaFold PDB already exists
        if out_path.exists():
            return f"[=] Exists: {out_path.name}"
        
        # First, get the PDB URL from the AlphaFold API
        api_url = ALPHAFOLD_API.format(uniprot_id)
        api_res = requests.get(api_url, timeout=10)
        if api_res.status_code != 200:
            return f"[✖] Not found on AlphaFold: {orf_name}"
        
        predictions = api_res.json()
        if not predictions:
            return f"[✖] Not found on AlphaFold: {orf_name}"
        
        # Get the PDB URL from the first prediction
        pdb_url = predictions[0].get("pdbUrl")
        if not pdb_url:
            return f"[✖] Not found on AlphaFold: {orf_name}"
        
        # Download the PDB file
        res = requests.get(pdb_url, timeout=15)
        if res.status_code == 200:
            out_path.write_text(res.text)
            return f"[✔] Downloaded: {out_path}"
        return f"[✖] Not found on AlphaFold: {orf_name}"
    except Exception as e:
        return f"[!] Error downloading {uniprot_id} ({orf_name}): {e}"
    
def run_on_gpu(fasta_path, pdb_out_dir):
    subprocess.run([
        "sbatch",
        "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/features_flow/colabfold_gpu.sh",
        str(fasta_path),
        str(pdb_out_dir)
    ], check=True)

    
def predict_structure_with_colabfold(orf_name, sequence, output_dir):
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = output_dir / f"{orf_name}.fasta"
        pdb_out_dir = output_dir / f"{orf_name}_alphafold"
        
        # Check if ColabFold prediction already exists
        existing_pdb = pdb_out_dir / f"{orf_name}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb"
        if existing_pdb.exists():
            return orf_name, sequence, output_dir  # Already exists

        # Remove stop codons (*) from sequence as ColabFold cannot process them
        cleaned_sequence = sequence.replace("*", "")
        fasta_path.write_text(f">{orf_name}\n{cleaned_sequence}\n")
        
        while True:
            result = subprocess.run(
                ["squeue", "-u", "alongonda"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            lines = result.stdout.strip().splitlines()
            num_jobs = len(lines) - 1  # subtract header
            if num_jobs < 100:
                break
            print(f"[~] {num_jobs} jobs running for 'alongonda'. Waiting...")
            time.sleep(60)

        run_on_gpu(fasta_path, pdb_out_dir)
        return orf_name, sequence, output_dir

    except Exception as e:
        print(f"[!] ColabFold error for {orf_name}: {e}")
        return None, None, None

# ------------------------- Main Pipeline -------------------------

def copy_to_all_appearances(orf, source_file, orf_to_fasta, output_root):
    """Copy a file to all directories where this ORF appears"""
    stems = orf_to_fasta[orf]
    if len(stems) <= 1:
        return  # Only one appearance, no copying needed
    
    copied_count = 0
    for stem in stems[1:]:  # Skip first one (original location)
        target_dir = output_root / stem
        target_dir.mkdir(parents=True, exist_ok=True)
        target_file = target_dir / source_file.name
        if not target_file.exists():
            shutil.copy2(source_file, target_file)
            copied_count += 1
    
    if copied_count > 0:
        print(f"[~] Copied {source_file.name} to {copied_count} additional directories for {orf}")

def main(fasta_list_path):
    output_dir = os.path.join(os.path.dirname(fasta_list_path), "kegg_random_mgc_pdb_files")
    output_root = Path(output_dir)
    output_root.mkdir(exist_ok=True)

    cache_file = output_root / "orf_cache.tsv"
    orf_cache = load_orf_cache(cache_file)

    # Update cache with new entries
    with open(fasta_list_path) as f:
        for line in f:
            fasta_path = Path(line.strip())
            if fasta_path.exists():
                new_entries = update_orf_cache(fasta_path, orf_cache)
                orf_cache.update(new_entries)
            else:
                print(f"[!] File not found: {fasta_path}")
    save_orf_cache(cache_file, orf_cache)

    # Reconstruct useful dictionaries - track ALL appearances of each ORF
    orf_to_fasta = {}
    sequences = {}
    for (orf, stem), seq in orf_cache.items():
        if orf not in orf_to_fasta:
            orf_to_fasta[orf] = []
            sequences[orf] = seq
        orf_to_fasta[orf].append(stem)

    print(f"[*] Total unique ORF names: {len(orf_to_fasta)}")
    print(f"[*] Total ORF appearances: {sum(len(stems) for stems in orf_to_fasta.values())}")

    uniprot_cache_file = output_root / "uniprot_cache.tsv"
    uniprot_cache = load_uniprot_cache(uniprot_cache_file)

    # Only resolve unknown ORFs
    to_resolve = [orf for orf in orf_to_fasta if orf not in uniprot_cache]

    print(f"[*] Resolving {len(to_resolve)} new ORFs from UniProt")

    uniprot_to_resolve = dict()
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(resolve_uniprot_id, orf): orf for orf in to_resolve}
        for future in as_completed(futures):
            orf, uid = future.result()
            uniprot_cache[orf] = uid
            uniprot_to_resolve[orf] = uid

    save_uniprot_cache(uniprot_cache_file, uniprot_cache)

    resolved_ids = [(orf, uid) for orf, uid in uniprot_to_resolve.items()]
    print(f"[+] Total UniProt IDs available: {len(resolved_ids)}")

    orf_for_colabfold = dict()
    
    # Separate ORFs with and without UniProt IDs
    orfs_with_uniprot = {orf: uid for orf, uid in uniprot_cache.items() if uid}
    orfs_without_uniprot = {orf: uid for orf, uid in uniprot_cache.items() if not uid}
    
    print(f"[*] ORFs with UniProt IDs: {len(orfs_with_uniprot)}")
    print(f"[*] ORFs without UniProt IDs: {len(orfs_without_uniprot)}")
    
    # Try AlphaFold download for ORFs with UniProt IDs
    with ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(
                download_alphafold_pdb,
                uid,
                orf,
                output_root / orf_to_fasta[orf][0]  # Use first appearance for download
            ): (orf, uid)
            for orf, uid in orfs_with_uniprot.items()
        }
        print(f"[*] Downloading AlphaFold PDBs for {len(futures)} ORFs")
        for future in as_completed(futures):
            result = future.result()
            print(result)
            if result.startswith("[✖] Not found on AlphaFold:"):
                orf = result.split(":")[1].strip()
                # Check ALL directories for this ORF
                pdb_exists = False
                for stem in orf_to_fasta[orf]:
                    subdir = output_root / stem
                    alphafold_pdb = subdir / f"{orf}.pdb"
                    if alphafold_pdb.exists():
                        pdb_exists = True
                        break
                
                if not pdb_exists:
                    # Use first directory for ColabFold prediction
                    subdir = output_root / orf_to_fasta[orf][0]
                    orf_for_colabfold[orf] = (orf, sequences[orf], subdir)
                    print(f"[!] {orf} not found on AlphaFold, will run ColabFold")
                else:
                    print(f"[~] {orf} has AlphaFold PDB, skipping ColabFold")
            elif result.startswith("[✔] Downloaded:"):
                # Extract ORF name and copy to all directories
                orf = result.split("/")[-1].replace(".pdb", "")
                pdb_file = Path(result.split(": ")[1])
                copy_to_all_appearances(orf, pdb_file, orf_to_fasta, output_root)
    
    # For ORFs without UniProt IDs, check if we need ColabFold
    for orf, uid in orfs_without_uniprot.items():
        # Check ALL directories for this ORF
        pdb_exists = False
        for stem in orf_to_fasta[orf]:
            subdir = output_root / stem
            alphafold_pdb = subdir / f"{orf}.pdb"
            if alphafold_pdb.exists():
                pdb_exists = True
                break
        
        if not pdb_exists:
            # Use first directory for ColabFold prediction
            subdir = output_root / orf_to_fasta[orf][0]
            orf_for_colabfold[orf] = (orf, sequences[orf], subdir)
            print(f"[!] {orf} has no UniProt ID, will run ColabFold")
        else:
            print(f"[~] {orf} has AlphaFold PDB, skipping ColabFold")
                
    with ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(
                predict_structure_with_colabfold,
                orf,
                seq,
                subdir
            ): (orf, seq, subdir)
            for orf, (orf, seq, subdir) in orf_for_colabfold.items()
        }
        print(f"[*] Running ColabFold for {len(futures)} ORFs")
        for future in as_completed(futures):
            orf, seq, subdir = future.result()
            if orf:
                print(f"[✔] ColabFold completed for {orf}")
                # Copy ColabFold results to all directories where this ORF appears
                colabfold_dir = subdir / f"{orf}_alphafold"
                if colabfold_dir.exists():
                    stems = orf_to_fasta[orf]
                    for stem in stems[1:]:  # Skip first one (already has it)
                        target_dir = output_root / stem / f"{orf}_alphafold"
                        if not target_dir.exists():
                            shutil.copytree(colabfold_dir, target_dir)
                            print(f"[~] Copied ColabFold results for {orf} to {stem}")
            else:
                print(f"[!] ColabFold failed for {orf}")

# ------------------------- Entry -------------------------

if __name__ == "__main__":
    fasta_list_file = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/kegg_random_mgc_candidates_fasta_files.txt"
    main(fasta_list_file)
