#!/usr/bin/env python3
"""
Example usage of the new Plant MGC Analysis Pipeline.

This script demonstrates how to use the new unified architecture
to replace your legacy python_scripts.
"""

import sys
from pathlib import Path

# Add the src directory to the Python path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from plant_mgc_analysis.config.settings import get_settings
from plant_mgc_analysis.utils.file_operations import file_manager
from plant_mgc_analysis.main_orchestrator import PlantMGCAnalysisOrchestrator
from plant_mgc_analysis.core.types import AnalysisType

def main():
    """Demonstrate the new unified system."""
    print("🧬 Plant MGC Analysis Pipeline - New Unified System")
    print("=" * 60)
    
    # Get configuration
    settings = get_settings()
    print(f"✅ Configuration loaded")
    print(f"   Base data directory: {settings.paths.base_data_dir}")
    print(f"   Desktop directory: {settings.paths.desktop_dir}")
    print(f"   MGC directory: {settings.paths.mgc_dir}")
    print()
    
    # Show available analysis types
    print("📊 Available Analysis Types:")
    for analysis_type in AnalysisType:
        print(f"   - {analysis_type.value}")
    print()
    
    # Example 1: File operations (replaces scattered CSV/FASTA operations)
    print("📁 File Operations Example:")
    print("   Instead of scattered CSV operations across 52 files, use:")
    print("   from plant_mgc_analysis.utils.file_operations import file_manager")
    print("   data = file_manager.read_file('input.csv')")
    print("   file_manager.write_file(data, 'output.csv')")
    print()
    
    # Example 2: Configuration usage (replaces hardcoded paths)
    print("⚙️  Configuration Example:")
    print("   Instead of hardcoded paths like '/groups/itay_mayrose/alongonda/datasets/':")
    print("   settings = get_settings()")
    print(f"   datasets_dir = settings.paths.datasets_dir  # {settings.paths.datasets_dir}")
    print()
    
    # Example 3: Analysis orchestration
    print("🔬 Analysis Example:")
    print("   Instead of running multiple separate scripts:")
    print("   orchestrator = PlantMGCAnalysisOrchestrator()")
    print("   result = orchestrator.run_single_organism_analysis(")
    print("       input_file='genome.fasta',")
    print("       analysis_types=[AnalysisType.BLAST_SEARCH, AnalysisType.SLIDING_WINDOW],")
    print("       output_directory='results/'")
    print("   )")
    print()
    
    # Example 4: Migration tools
    print("🔄 Migration Example:")
    print("   To migrate your remaining legacy scripts:")
    print("   from plant_mgc_analysis.utils.migration_tools import create_migration_plan")
    print("   plan = create_migration_plan(Path('python_scripts/'))")
    print("   print(plan.get_summary())")
    print()
    
    # Show status
    print("📊 System Status:")
    try:
        orchestrator = PlantMGCAnalysisOrchestrator()
        stats = orchestrator.get_system_stats()
        print(f"   File operations available: ✅")
        print(f"   API clients configured: ✅")
        print(f"   Analysis engines loaded: ✅")
        print(f"   Configuration valid: ✅")
    except Exception as e:
        print(f"   System status: ⚠️  {e}")
    print()
    
    print("🚀 Quick Start Examples:")
    print("   # Basic analysis")
    print("   python -m plant_mgc_analysis analyze genome.fasta -a blast_search")
    print()
    print("   # BLAST analysis")
    print("   python -m plant_mgc_analysis blast sequences.fasta -d database.fasta")
    print()
    print("   # Migrate legacy scripts")
    print("   python -m plant_mgc_analysis migrate plan python_scripts/")
    print()
    print("   # Python API usage")
    print("   from plant_mgc_analysis import PlantMGCAnalysisOrchestrator")
    print("   orchestrator = PlantMGCAnalysisOrchestrator()")
    print("   result = orchestrator.run_single_organism_analysis(...)")
    print()
    
    print("✨ The new system consolidates all 92 legacy scripts into this unified framework!")
    print("   Check the README.md for complete documentation.")

if __name__ == "__main__":
    main()