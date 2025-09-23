import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
FOLDSEEK_DB_ROOT = Path(
    "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mibig_kegg_foldseek_predictions/"
)
OUTPUT_DIR = FOLDSEEK_DB_ROOT / "foldseek_results"
FOLDSEEK = "foldseek"
MAX_WORKERS = 16  # Adjust based on your system

def run_foldseek_all_vs_all(mgc_dir: Path):
    """Run Foldseek all-vs-all comparison for a single MGC database"""
    mgc_name = mgc_dir.name
    db_prefix = mgc_dir / f"{mgc_name}_merged"
    
    # Check if database exists
    if not (db_prefix.with_suffix(".dbtype")).exists():
        return f"[!] Database not found: {db_prefix}"
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Output file
    output_file = OUTPUT_DIR / f"{mgc_name}_all_vs_all.tsv"
    tmp_dir = OUTPUT_DIR / f"tmp_{mgc_name}"
    
    # Skip if already processed
    if output_file.exists() and output_file.stat().st_size > 0:
        return f"[=] Already processed: {mgc_name}"
    
    try:
        # Run Foldseek easy-search (all-vs-all within the same database)
        cmd = [
            FOLDSEEK, "easy-search",
            db_prefix.as_posix(),  # Query database
            db_prefix.as_posix(),  # Target database (same as query for all-vs-all)
            output_file.as_posix(),
            tmp_dir.as_posix(),
            "-s", "7.5",  # Maximum sensitivity
            "--max-seqs", "300",  # Max hits per query
            "-e", "0.001",  # E-value threshold
            "--format-output", "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,lddt,prob"
        ]
        
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        # Clean up temp directory
        if tmp_dir.exists():
            subprocess.run(["rm", "-rf", tmp_dir.as_posix()], check=False)
        
        # Count results
        with open(output_file) as f:
            num_hits = sum(1 for _ in f)
        
        return f"[✔] {mgc_name}: {num_hits} hits found"
        
    except subprocess.CalledProcessError as e:
        # Clean up on failure
        if tmp_dir.exists():
            subprocess.run(["rm", "-rf", tmp_dir.as_posix()], check=False)
        return f"[✖] Failed {mgc_name}: {e.stderr}"
    except Exception as e:
        return f"[✖] Error {mgc_name}: {str(e)}"

def main():
    # Find all MGC directories with Foldseek databases
    mgc_dirs = []
    for d in FOLDSEEK_DB_ROOT.iterdir():
        if d.is_dir() and d.name.startswith("BGC"):
            # Check if the database files exist
            db_prefix = d / f"{d.name}_merged"
            if (db_prefix.with_suffix(".dbtype")).exists():
                mgc_dirs.append(d)
    
    print(f"[*] Found {len(mgc_dirs)} MGC databases to process")
    
    # Process with parallel execution
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(run_foldseek_all_vs_all, d): d.name for d in mgc_dirs}
        
        completed = 0
        for future in as_completed(futures):
            completed += 1
            mgc_name = futures[future]
            result = future.result()
            print(f"[{completed}/{len(mgc_dirs)}] {result}")
    
    print("[*] All Foldseek searches completed")

if __name__ == "__main__":
    main()