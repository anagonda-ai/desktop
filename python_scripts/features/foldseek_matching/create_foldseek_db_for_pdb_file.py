import os
import subprocess
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
ROOT_DIR = Path(
    "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_pdb_files/"
)
OLD_OUTPUT_ROOT = ROOT_DIR.parent / "mibig_kegg_foldseek_predictions"
FOLDSEEK = "foldseek"
MAX_WORKERS = 16

def collect_pdbs(mgc_dir: Path):
    """Collect all PDB files from the original MGC directory"""
    pdbs = list(mgc_dir.glob("*.pdb"))
    for sub in mgc_dir.glob("*_alphafold"):
        if sub.is_dir():
            top = list(sub.glob("*.pdb"))
            if top:
                pdbs.append(top[0])
    return pdbs

def fix_merged_pdb_and_recreate_db(mgc_name: str):
    """Fix the merged PDB file to properly separate proteins, then recreate the database"""
    
    # Get original PDB files from source directory
    original_mgc_dir = ROOT_DIR / mgc_name
    if not original_mgc_dir.exists():
        return f"[!] Original directory not found: {mgc_name}"
    
    pdb_files = collect_pdbs(original_mgc_dir)
    if not pdb_files:
        return f"[!] No PDB files found for {mgc_name}"
    
    # Work in the existing output directory
    mgc_output_dir = OLD_OUTPUT_ROOT / mgc_name
    mgc_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Backup the old merged file if it exists
    old_merged = mgc_output_dir / f"{mgc_name}_merged_input.pdb"
    if old_merged.exists():
        backup = mgc_output_dir / f"{mgc_name}_merged_input.pdb.backup"
        if not backup.exists():
            shutil.copy(old_merged, backup)
    
    # Create properly formatted merged PDB with clear protein separations
    fixed_merged_pdb = mgc_output_dir / f"{mgc_name}_merged_input.pdb"
    
    try:
        with open(fixed_merged_pdb, 'w') as out_f:
            for i, pdb_file in enumerate(pdb_files, 1):
                # Write a REMARK to identify each protein
                out_f.write(f"REMARK   1 Protein {i} from {pdb_file.name}\n")
                
                # Option 1: Use MODEL/ENDMDL tags (standard for multi-model PDB)
                out_f.write(f"MODEL     {i:4d}\n")
                
                with open(pdb_file) as in_f:
                    for line in in_f:
                        # Skip any existing END, MODEL, ENDMDL, or CONECT lines
                        if not line.startswith(('END', 'MODEL', 'ENDMDL', 'CONECT', 'MASTER')):
                            out_f.write(line)
                        elif line.startswith('TER'):
                            # Keep TER lines as they mark chain endings
                            out_f.write(line)
                
                out_f.write("ENDMDL\n")
                
                # Option 2: Alternative - use END between proteins (uncomment if MODEL doesn't work)
                # with open(pdb_file) as in_f:
                #     for line in in_f:
                #         if not line.startswith(('END', 'MODEL', 'ENDMDL')):
                #             out_f.write(line)
                # out_f.write("END\n")
            
            # Write final END for the entire file
            out_f.write("END\n")
        
        # Remove old database files
        db_prefix = mgc_output_dir / f"{mgc_name}_merged"
        for db_file in mgc_output_dir.glob(f"{mgc_name}_merged*"):
            if not db_file.suffix == '.pdb' and not str(db_file).endswith('.backup'):
                db_file.unlink()
        
        # Recreate the database with the fixed PDB file
        cmd = [FOLDSEEK, "createdb", fixed_merged_pdb.as_posix(), db_prefix.as_posix()]
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        # Verify the number of sequences in the new database
        verify_cmd = [
            FOLDSEEK, "createtsv",
            db_prefix.as_posix(),
            db_prefix.as_posix(),
            "-"  # Output to stdout
        ]
        verify_result = subprocess.run(verify_cmd, capture_output=True, text=True, check=False)
        
        protein_count = 0
        if verify_result.returncode == 0:
            protein_count = len(verify_result.stdout.strip().split('\n')) if verify_result.stdout else 0
        
        if protein_count > 1:
            return f"[✔] Fixed {mgc_name}: {protein_count} proteins recognized (from {len(pdb_files)} PDB files)"
        elif protein_count == 1:
            return f"[⚠] Warning {mgc_name}: Still shows 1 protein (expected {len(pdb_files)})"
        else:
            return f"[?] {mgc_name}: Could not verify protein count"
        
    except subprocess.CalledProcessError as e:
        return f"[✖] Failed to fix {mgc_name}: {e.stderr}"
    except Exception as e:
        return f"[✖] Error fixing {mgc_name}: {str(e)}"

def test_single_mgc(mgc_name: str):
    """Test fixing a single MGC to verify it works"""
    result = fix_merged_pdb_and_recreate_db(mgc_name)
    print(result)
    
    # Run a quick all-vs-all to check the results
    mgc_output_dir = OLD_OUTPUT_ROOT / mgc_name
    db_prefix = mgc_output_dir / f"{mgc_name}_merged"
    
    if (db_prefix.with_suffix(".dbtype")).exists():
        test_output = mgc_output_dir / "test_all_vs_all.tsv"
        tmp_dir = mgc_output_dir / "tmp_test"
        
        cmd = [
            FOLDSEEK, "easy-search",
            db_prefix.as_posix(),
            db_prefix.as_posix(),
            test_output.as_posix(),
            tmp_dir.as_posix(),
            "-s", "7.5",
            "--max-seqs", "300"
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            with open(test_output) as f:
                hits = sum(1 for _ in f)
            print(f"  Test search found {hits} hits")
            
            # Show first few hits
            with open(test_output) as f:
                print("  First few hits:")
                for i, line in enumerate(f):
                    if i < 5:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            print(f"    {parts[0]} vs {parts[1]}")
            
            # Cleanup
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)
            test_output.unlink()
                
        except Exception as e:
            print(f"  Test search failed: {e}")

def main():
    # Find all MGC directories that need fixing
    mgc_dirs_to_fix = []
    
    if OLD_OUTPUT_ROOT.exists():
        for d in OLD_OUTPUT_ROOT.iterdir():
            if d.is_dir() and d.name.startswith("BGC"):
                mgc_dirs_to_fix.append(d.name)
    
    if not mgc_dirs_to_fix:
        print("[!] No MGC databases found to fix")
        return
    
    print(f"[*] Found {len(mgc_dirs_to_fix)} MGC databases to fix")
    
    # Process all MGCs
    print(f"\n[*] Fixing all {len(mgc_dirs_to_fix)} MGC databases...")
    
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(fix_merged_pdb_and_recreate_db, mgc_name): mgc_name 
                  for mgc_name in mgc_dirs_to_fix}
        
        completed = 0
        successful = 0
        warnings = 0
        
        for future in as_completed(futures):
            completed += 1
            mgc_name = futures[future]
            result = future.result()
            
            if "[✔]" in result:
                successful += 1
            elif "[⚠]" in result:
                warnings += 1
                
            print(f"[{completed}/{len(mgc_dirs_to_fix)}] {result}")
    
    print(f"\n[*] Complete! Successfully fixed: {successful}, Warnings: {warnings}")
    print("[*] Databases have been fixed in place")
    print("[*] Old merged PDB files backed up as .backup")

if __name__ == "__main__":
    main()