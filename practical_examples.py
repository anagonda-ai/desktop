#!/usr/bin/env python3
"""
Practical examples showing how to run your bioinformatics tools 
using the new unified architecture.
"""

import sys
from pathlib import Path

# Add the src directory to the Python path
sys.path.insert(0, str(Path(__file__).parent / "src"))

def example_file_operations():
    """Show unified file operations."""
    print("📁 File Operations Examples:")
    print("=" * 40)
    
    from plant_mgc_analysis.utils.file_operations import file_manager
    from plant_mgc_analysis.config.settings import get_settings
    
    settings = get_settings()
    
    print("Instead of scattered CSV operations, use:")
    print("# Read any CSV file")
    print("data = file_manager.csv_handler.read_csv('path/to/file.csv')")
    print()
    
    print("# Read FASTA files")
    print("sequences = file_manager.fasta_handler.read_fasta('sequences.fasta')")
    print()
    
    print("# Convert FASTA to CSV")
    print("file_manager.convert_fasta_to_csv('input.fasta', 'output.csv')")
    print()
    
    print("# Merge multiple CSV files")
    print("file_manager.merge_csv_files(['file1.csv', 'file2.csv'], 'merged.csv')")
    print()
    
    print("# Use configured paths instead of hardcoded ones")
    print(f"arabidopsis_dir = settings.paths.arabidopsis_dir  # {settings.paths.arabidopsis_dir}")
    print(f"mgc_dir = settings.paths.mgc_dir  # {settings.paths.mgc_dir}")
    print()

def example_blast_analysis():
    """Show BLAST analysis using new system."""
    print("🔬 BLAST Analysis Example:")
    print("=" * 40)
    
    print("Instead of running blast_*.py scripts, use:")
    print()
    print("# Python API:")
    print("from plant_mgc_analysis.genomics.blast_analysis import BlastAnalyzer")
    print("from plant_mgc_analysis.core.types import AnalysisType")
    print()
    print("analyzer = BlastAnalyzer()")
    print("result = analyzer.run_analysis(genome_data, {")
    print("    'subject_file': 'database.fasta',")
    print("    'evalue': 1e-5,")
    print("    'max_targets': 100")
    print("})")
    print()
    
    print("# Command line (when CLI is working):")
    print("python -m plant_mgc_analysis blast sequences.fasta -d database.fasta --evalue 1e-5")
    print()

def example_kegg_analysis():
    """Show KEGG analysis using new system."""
    print("🧬 KEGG Analysis Example:")
    print("=" * 40)
    
    print("Instead of running kegg_*.py scripts, use:")
    print()
    print("# Python API:")
    print("from plant_mgc_analysis.utils.api_clients import KEGGAPIClient")
    print()
    print("kegg_client = KEGGAPIClient()")
    print("pathways = kegg_client.list_pathways('ath')  # Arabidopsis")
    print("organism_info = kegg_client.get_organism_info('ath')")
    print()
    
    print("# Command line (when CLI is working):")
    print("python -m plant_mgc_analysis kegg sequences.fasta --organisms ath")
    print()

def example_configuration():
    """Show configuration management."""
    print("⚙️  Configuration Management:")
    print("=" * 40)
    
    from plant_mgc_analysis.config.settings import get_settings
    
    settings = get_settings()
    
    print("Instead of hardcoded paths like:")
    print('"/groups/itay_mayrose/alongonda/datasets/"')
    print()
    print("Use configured paths:")
    print(f"settings.paths.datasets_dir = {settings.paths.datasets_dir}")
    print(f"settings.paths.desktop_dir = {settings.paths.desktop_dir}")
    print(f"settings.paths.mgc_dir = {settings.paths.mgc_dir}")
    print(f"settings.paths.arabidopsis_dir = {settings.paths.arabidopsis_dir}")
    print()
    
    print("Get organism-specific paths:")
    print("organism_config = settings.get_organism_config('arabidopsis')")
    print("# This gives you all paths for a specific organism")
    print()

def example_migration_workflow():
    """Show how to migrate legacy scripts."""
    print("🔄 Migration Workflow:")
    print("=" * 40)
    
    print("Step 1: Analyze your legacy scripts")
    print("from plant_mgc_analysis.utils.migration_tools import create_migration_plan")
    print("plan = create_migration_plan(Path('python_scripts/'))")
    print("print(plan.get_summary())")
    print()
    
    print("Step 2: Use the new unified functions")
    print("# Instead of your custom CSV reading:")
    print("from plant_mgc_analysis.utils.file_operations import file_manager")
    print("data = file_manager.csv_handler.read_csv('file.csv')")
    print()
    
    print("# Instead of hardcoded paths:")
    print("from plant_mgc_analysis.config.settings import get_settings")
    print("settings = get_settings()")
    print("data_dir = settings.paths.datasets_dir")
    print()
    
    print("Step 3: Use unified analysis")
    print("# Instead of separate analysis scripts:")
    print("from plant_mgc_analysis.genomics.blast_analysis import BlastAnalyzer")
    print("analyzer = BlastAnalyzer()")
    print("result = analyzer.run_analysis(data, parameters)")
    print()

def practical_replacement_guide():
    """Show specific script replacements."""
    print("🔄 Specific Script Replacements:")
    print("=" * 40)
    
    replacements = [
        ("blast_mibig_vs_all_genomes_parallel.py", "python -m plant_mgc_analysis blast *.fasta -d mibig.fasta --parallel"),
        ("kegg_pathway_extractor.py", "python -m plant_mgc_analysis kegg sequences.fasta --organisms ath"),
        ("sliding_window.py", "python -m plant_mgc_analysis sliding-window genome.fasta"),
        ("fasta_creator.py", "file_manager.fasta_handler.write_fasta(sequences, 'output.fasta')"),
        ("csv_converter.py", "file_manager.convert_fasta_to_csv('input.fasta', 'output.csv')"),
        ("file_merger.py", "file_manager.merge_csv_files(['file1.csv', 'file2.csv'], 'merged.csv')"),
    ]
    
    for old, new in replacements:
        print(f"❌ OLD: {old}")
        print(f"✅ NEW: {new}")
        print()

def main():
    """Run all examples."""
    print("🧬 Practical Examples - How to Use Your New Tools")
    print("=" * 60)
    print()
    
    example_file_operations()
    print()
    
    example_blast_analysis()
    print()
    
    example_kegg_analysis()
    print()
    
    example_configuration()
    print()
    
    example_migration_workflow()
    print()
    
    practical_replacement_guide()
    print()
    
    print("🎯 Summary - How to Run Your Tools:")
    print("=" * 40)
    print("1. ✅ Use the unified file operations instead of scattered CSV/FASTA code")
    print("2. ✅ Use configured paths instead of hardcoded ones")
    print("3. ✅ Use the unified analysis classes instead of separate scripts")
    print("4. ✅ Use the migration tools to convert remaining scripts")
    print("5. ✅ Use the CLI commands for common operations")
    print()
    
    print("📚 For complete documentation, see README.md")
    print("🔧 For migration help, run: python -m plant_mgc_analysis migrate plan python_scripts/")

if __name__ == "__main__":
    main()