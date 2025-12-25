#!/usr/bin/env python3
"""
Setup script for MGC Predictor Tool
Ensures all required data files and databases are ready for use.
"""

import os
import sys
import subprocess
from pathlib import Path

# Add parent directory to path for imports
FINAL_TOOL_DIR = Path(__file__).parent
sys.path.insert(0, str(FINAL_TOOL_DIR.parent))
sys.path.insert(0, str(FINAL_TOOL_DIR))

# Import config
import importlib.util
config_spec = importlib.util.spec_from_file_location("config", FINAL_TOOL_DIR / "config.py")
config = importlib.util.module_from_spec(config_spec)
config_spec.loader.exec_module(config)

def check_blast_installed():
    """Check if BLAST+ is installed and available"""
    try:
        result = subprocess.run(['blastp', '-version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0] if result.stdout else "BLAST+ installed"
            print(f"✓ BLAST+ is installed: {version_line}")
            return True
    except FileNotFoundError:
        pass
    except Exception as e:
        print(f"⚠ Warning: Could not check BLAST+ version: {e}")
    
    print("✗ BLAST+ is not installed or not in PATH")
    print("  Install via: conda install -c bioconda blast")
    return False

def check_kegg_fasta():
    """Check if KEGG FASTA file exists"""
    if os.path.exists(config.KEGG_DB):
        size_mb = os.path.getsize(config.KEGG_DB) / (1024 * 1024)
        print(f"✓ KEGG FASTA file exists: {config.KEGG_DB} ({size_mb:.1f} MB)")
        return True
    else:
        print(f"✗ KEGG FASTA file not found: {config.KEGG_DB}")
        return False

def check_blast_database():
    """Check if BLAST database files exist"""
    # BLAST database files: .phr, .pin, .psq (for protein databases)
    db_base = str(Path(config.KEGG_DB).with_suffix(''))
    db_files = [f"{db_base}.phr", f"{db_base}.pin", f"{db_base}.psq"]
    
    all_exist = all(os.path.exists(f) for f in db_files)
    
    if all_exist:
        print(f"✓ BLAST database files exist: {db_base}.*")
        return True
    else:
        missing = [f for f in db_files if not os.path.exists(f)]
        print(f"✗ BLAST database files missing: {', '.join([Path(f).name for f in missing])}")
        return False

def build_blast_database():
    """Build BLAST database from FASTA file"""
    if not check_blast_installed():
        print("  Cannot build database: BLAST+ not installed")
        return False
    
    if not check_kegg_fasta():
        print("  Cannot build database: KEGG FASTA file not found")
        return False
    
    db_base = str(Path(config.KEGG_DB).with_suffix(''))
    
    print(f"\nBuilding BLAST database from: {config.KEGG_DB}")
    print(f"Output database: {db_base}")
    print("This may take several minutes...")
    
    try:
        # Run makeblastdb
        cmd = [
            'makeblastdb',
            '-in', config.KEGG_DB,
            '-dbtype', 'prot',
            '-out', db_base,
            '-title', 'KEGG Metabolic Pathways'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        
        if result.returncode == 0:
            print("✓ BLAST database built successfully")
            return True
        else:
            print(f"✗ Error building BLAST database:")
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ Database build timed out (took longer than 1 hour)")
        return False
    except Exception as e:
        print(f"✗ Error building BLAST database: {e}")
        return False

def check_model_files():
    """Check if model files exist"""
    model_file = os.path.join(config.MODEL_WEIGHTS_DIR, config.COMBINED_MODEL)
    feature_order_file = os.path.join(config.MODEL_WEIGHTS_DIR, config.COMBINED_FEATURE_ORDER)
    weights_file = os.path.join(config.MODEL_WEIGHTS_DIR, config.COMBINED_WEIGHTS)
    
    all_exist = True
    
    if os.path.exists(model_file):
        size_kb = os.path.getsize(model_file) / 1024
        print(f"✓ Model file exists: {config.COMBINED_MODEL} ({size_kb:.1f} KB)")
    else:
        print(f"✗ Model file not found: {model_file}")
        all_exist = False
    
    if os.path.exists(feature_order_file):
        print(f"✓ Feature order file exists: {config.COMBINED_FEATURE_ORDER}")
    else:
        print(f"✗ Feature order file not found: {feature_order_file}")
        all_exist = False
    
    if os.path.exists(weights_file):
        print(f"✓ Weights file exists: {config.COMBINED_WEIGHTS}")
    else:
        print(f"⚠ Weights file not found: {weights_file} (optional)")
    
    return all_exist

def main():
    """Main setup function"""
    print("="*80)
    print("MGC Predictor Tool - Setup and Verification")
    print("="*80)
    print()
    
    # Check BLAST installation
    print("1. Checking BLAST+ installation...")
    blast_installed = check_blast_installed()
    print()
    
    # Check KEGG FASTA file
    print("2. Checking KEGG FASTA file...")
    kegg_exists = check_kegg_fasta()
    print()
    
    # Check BLAST database
    print("3. Checking BLAST database...")
    db_exists = check_blast_database()
    print()
    
    # Build database if needed (optional)
    if kegg_exists and not db_exists and blast_installed:
        print("4. BLAST database (optional)...")
        print("   Note: KEGG database is used as query FASTA, not as BLAST database")
        print("   Database build is optional and can be skipped")
        print()
    else:
        print("4. Skipping optional database build")
        print()
    
    # Check model files
    print("5. Checking model files...")
    models_exist = check_model_files()
    print()
    
    # Summary
    print("="*80)
    print("Setup Summary:")
    print("="*80)
    
    all_ready = True
    if not blast_installed:
        print("✗ BLAST+ is not installed")
        all_ready = False
    if not kegg_exists:
        print("✗ KEGG FASTA file is missing")
        all_ready = False
    # BLAST database is optional, don't fail if it doesn't exist
    # if not db_exists:
    #     print("✗ BLAST database is not built (optional)")
    if not models_exist:
        print("✗ Model files are missing")
        all_ready = False
    
    if all_ready:
        print("\n✓ All systems ready! You can run the MGC Predictor tool.")
        return 0
    else:
        print("\n✗ Some components are missing. Please install/configure them before running.")
        print("\nQuick fix commands:")
        if not blast_installed:
            print("  conda install -c bioconda blast")
        if not db_exists and kegg_exists and blast_installed:
            print(f"  python {__file__} --build-db")
        return 1

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Setup and verify MGC Predictor Tool')
    parser.add_argument('--build-db', action='store_true',
                       help='Build BLAST database from FASTA file')
    args = parser.parse_args()
    
    if args.build_db:
        if check_blast_installed() and check_kegg_fasta():
            success = build_blast_database()
            sys.exit(0 if success else 1)
        else:
            print("Cannot build database: prerequisites missing")
            sys.exit(1)
    else:
        sys.exit(main())

