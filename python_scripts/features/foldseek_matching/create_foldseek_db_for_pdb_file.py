import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
ROOT_DIR = Path(
    "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output/"
    "kegg_scanner_min_genes_based_metabolic/min_genes_3/"
    "mgc_candidates_fasta_files_without_e2p2_filtered_test/"
    "blast_all_vs_all/mgc_pdb_files"
)
OUTPUT_ROOT = ROOT_DIR.parent / "foldseek_predictions"
FOLDSEEK = "foldseek"
MAX_WORKERS = 16

def collect_pdbs(mgc_dir: Path):
    pdbs = list(mgc_dir.glob("*.pdb"))
    for sub in mgc_dir.glob("*_colabfold"):
        if sub.is_dir():
            top = list(sub.glob("*_unrelaxed_rank_001_*.pdb"))
            if top:
                pdbs.append(top[0])
    return pdbs

def merge_and_create_db(mgc_dir: Path):
    mgc_name = mgc_dir.name
    mgc_output_dir = OUTPUT_ROOT / mgc_name
    mgc_output_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = collect_pdbs(mgc_dir)
    if not pdb_files:
        return f"[!] No PDBs in {mgc_name}"

    merged_pdb = mgc_output_dir / f"{mgc_name}_merged_input.pdb"
    if not merged_pdb.exists():
        with open(merged_pdb, "w") as out_f:
            for pdb in pdb_files:
                with open(pdb) as in_f:
                    out_f.write(in_f.read())
        print(f"[+] Merged {len(pdb_files)} into {merged_pdb.name}")

    db_prefix = mgc_output_dir / f"{mgc_name}_merged"
    if (db_prefix.with_suffix(".dbtype")).exists():
        return f"[=] Already exists: {db_prefix.name}"

    try:
        subprocess.run(
            [FOLDSEEK, "createdb", merged_pdb.as_posix(), db_prefix.as_posix()],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        return f"[✔] Created DB for {mgc_name}"
    except subprocess.CalledProcessError as e:
        return f"[✖] DB failed for {mgc_name}: {e.stderr.decode()}"

def main():
    mgc_dirs = [d for d in ROOT_DIR.iterdir() if d.is_dir()]
    print(f"[*] Processing {len(mgc_dirs)} MGCs")

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(merge_and_create_db, d): d.name for d in mgc_dirs}
        for f in as_completed(futures):
            print(f.result())

if __name__ == "__main__":
    main()
