import requests
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import time
import threading

# Limit concurrent UniProt queries to 3 (UniProt recommends ~3/sec)
API_SEMAPHORE = threading.Semaphore(3)
REQUEST_DELAY = 0.4  # seconds between requests to avoid hammering the server
MAX_RETRIES = 3      # retry failed queries

UNIPROT_SEARCH_API = "https://rest.uniprot.org/uniprotkb/search?query={}&fields=accession&format=json"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb"

SOURCE_1 = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files"
SOURCE_2 = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_output/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"

def extract_orf_names_from_fasta(fasta_file):
    orf_data = []
    source_path = str(fasta_file)
    index = 0 if source_path.startswith(SOURCE_1) else 1 if source_path.startswith(SOURCE_2) else None
    if index is None:
        print(f"[!] Unknown source path: {fasta_file}")
        return orf_data

    for record in SeqIO.parse(fasta_file, "fasta"):
        info = record.id if index == 0 else record.description
        parts = info.split(" | ")
        if len(parts) > index:
            orf_name = parts[index]
            orf_data.append((orf_name, fasta_file.stem, str(record.seq)))
    return orf_data

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
                        return orf_name, entry["primaryAccession"]
                if results:
                    return orf_name, results[0]["primaryAccession"]
                break
            except Exception as e:
                print(f"[~] Retry {attempt} for {orf_name}: {e}")
                time.sleep(REQUEST_DELAY * 2)
            finally:
                time.sleep(REQUEST_DELAY)
    return orf_name, None

def download_alphafold_pdb(uniprot_id, orf_name, subdir):
    try:
        url = ALPHAFOLD_PDB_URL.format(uniprot_id)
        subdir.mkdir(parents=True, exist_ok=True)
        out_path = subdir / f"{orf_name}.pdb"
        if out_path.exists():
            return f"[=] Exists: {out_path.name}"
        res = requests.get(url, timeout=15)
        if res.status_code == 200:
            out_path.write_text(res.text)
            return f"[✔] Downloaded: {out_path}"
        return f"[✖] Not found on AlphaFold: {uniprot_id}"
    except Exception as e:
        return f"[!] Error downloading {uniprot_id} ({orf_name}): {e}"

def predict_structure_with_colabfold(orf_name, sequence, output_dir):
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = output_dir / f"{orf_name}.fasta"
        pdb_out_dir = output_dir / f"{orf_name}_colabfold"

        fasta_content = f">{orf_name}\n{sequence}\n"
        fasta_path.write_text(fasta_content)

        subprocess.run([
            "colabfold_batch", str(fasta_path), str(pdb_out_dir), "--msa-mode", "single_sequence"
        ], check=True)

        result_pdb = next(pdb_out_dir.glob("*.pdb"), None)
        if result_pdb:
            final_pdb = output_dir / f"{orf_name}_colabfold.pdb"
            result_pdb.rename(final_pdb)
            return f"[✔] Predicted with ColabFold: {final_pdb.name}"
        return f"[✖] ColabFold ran but no PDB found for {orf_name}"
    except Exception as e:
        return f"[!] ColabFold error for {orf_name}: {e}"

def main(fasta_list_path, output_root="alphafold_pdbs"):
    output_root = Path(output_root)
    output_root.mkdir(exist_ok=True)
    orf_to_fasta = {}
    sequences = {}
    with open(fasta_list_path) as f:
        for line in f:
            fasta_path = Path(line.strip())
            if fasta_path.exists():
                for orf_name, fasta_stem, seq in extract_orf_names_from_fasta(fasta_path):
                    orf_to_fasta[orf_name] = fasta_stem
                    sequences[orf_name] = seq
            else:
                print(f"[!] File not found: {fasta_path}")

    print(f"[*] Total unique ORF names: {len(orf_to_fasta)}")

    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(resolve_uniprot_id, orf): orf for orf in orf_to_fasta}
        resolved = [f.result() for f in as_completed(futures)]

    resolved_ids = [(orf, uid) for orf, uid in resolved if uid]
    print(f"[+] Resolved {len(resolved_ids)} UniProt IDs")

    with ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(
                download_alphafold_pdb,
                uid,
                orf,
                output_root / orf_to_fasta[orf]
            ): (orf, uid)
            for orf, uid in resolved_ids
        }
        for future in as_completed(futures):
            result = future.result()
            print(result)
            if result.startswith("[✖] Not found on AlphaFold:"):
                uid = result.split(":")[1].strip()
                orf = [o for o, u in resolved_ids if u == uid][0]
                subdir = output_root / orf_to_fasta[orf]
                colab_result = predict_structure_with_colabfold(orf, sequences[orf], subdir)
                print(colab_result)

if __name__ == "__main__":
    fasta_list_file = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_output/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/blast_all_vs_all/deduplicated_file_list.txt"
    main(fasta_list_file)
