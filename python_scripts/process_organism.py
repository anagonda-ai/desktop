#!/usr/bin/env python3
"""
State-of-the-art organism processing for KEGG metabolic modules.

This module processes organisms to extract metabolic module gene information
from KEGG database with professional error handling, rate limiting, and
comprehensive validation.

Usage:
    python process_organism.py <org_code> <root_folder> <metabolic_modules_path>

Example:
    python process_organism.py hsa /data/organisms modules.csv
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
import pandas as pd
import logging

# Import our new framework
sys.path.insert(0, str(Path(__file__).parent.parent))
from core import BaseProcessor, Config, get_config
from core.exceptions import ProcessingError, ValidationError, APIError
from core.types import ProcessingConfig, AnalysisResult, ProcessingStatus
from utils import KEGGClient, CSVHandler, BatchProcessor, ProgressTracker


@dataclass
class OrganismProcessingConfig(ProcessingConfig):
    """Enhanced configuration for organism processing."""
    batch_size: int = 10
    sleep_time: float = 0.35
    max_workers: int = 4
    skip_existing: bool = True
    validate_modules: bool = True


@dataclass
class GeneInfo:
    """Structured gene information."""
    gene_id: str
    annotation: str
    nucleotide_sequence: str
    location: str
    
    @classmethod
    def from_kegg_data(cls, gene_id: str, kegg_data: Tuple[str, str, str]) -> 'GeneInfo':
        """Create GeneInfo from KEGG API response tuple."""
        annotation, nucleotide_seq, location = kegg_data
        return cls(
            gene_id=gene_id,
            annotation=annotation or "N/A",
            nucleotide_sequence=nucleotide_seq or "N/A", 
            location=location or "N/A"
        )
    
    def to_csv_row(self) -> List[str]:
        """Convert to CSV row format."""
        return [self.gene_id, self.annotation, self.nucleotide_sequence, self.location]


@dataclass
class ModuleData:
    """Metabolic module data."""
    module_id: str
    genes: List[GeneInfo]
    organism_code: str
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to pandas DataFrame."""
        if not self.genes:
            return pd.DataFrame(columns=["Gene_ID", "Annotation", "Nucleotide_Sequence", "Location"])
        
        rows = [gene.to_csv_row() for gene in self.genes]
        return pd.DataFrame(rows, columns=["Gene_ID", "Annotation", "Nucleotide_Sequence", "Location"])


class KEGGParser:
    """
    Professional KEGG response parser.
    
    Handles parsing of KEGG flat file format with proper error handling
    and validation.
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def parse_kegg_entry(self, entry_text: str) -> Dict[str, str]:
        """
        Parse KEGG flat file entry format.
        
        Args:
            entry_text: Raw KEGG entry text
            
        Returns:
            Dictionary of parsed fields
        """
        parsed = {}
        current_key = None
        
        try:
            for line in entry_text.splitlines():
                if not line.strip():
                    continue
                
                # Check if line starts a new field (first 12 characters)
                field_name = line[:12].strip()
                if field_name:
                    current_key = field_name
                    # Special handling for AASEQ field
                    parsed[current_key] = line[12:].strip() if current_key != "AASEQ" else ""
                else:
                    # Continuation line
                    if current_key:
                        parsed[current_key] += line[12:].strip()
            
            return parsed
            
        except Exception as e:
            self.logger.error(f"Failed to parse KEGG entry: {e}")
            return {}
    
    def extract_gene_info(self, parsed_entry: Dict[str, str]) -> Tuple[str, str, str, str]:
        """
        Extract gene information from parsed entry.
        
        Args:
            parsed_entry: Parsed KEGG entry
            
        Returns:
            Tuple of (gene_id, annotation, sequence, location)
        """
        gene_id = parsed_entry.get("ENTRY", "").split()[0] if parsed_entry.get("ENTRY") else ""
        annotation = parsed_entry.get("ORTHOLOGY", "")
        sequence = parsed_entry.get("AASEQ", "")
        location = parsed_entry.get("POSITION", "")
        
        return gene_id, annotation, sequence, location


class OrganismProcessor(BaseProcessor):
    """
    State-of-the-art organism processor using professional architecture.
    
    Processes organisms to extract metabolic module information from KEGG
    with comprehensive error handling, progress tracking, and validation.
    """
    
    def __init__(self, config: Optional[OrganismProcessingConfig] = None):
        super().__init__(config)
        self.config = config or OrganismProcessingConfig()
        self.global_config = get_config()
        
        # Initialize components
        self.kegg_client = KEGGClient()
        self.csv_handler = CSVHandler()
        self.parser = KEGGParser()
        self.batch_processor = BatchProcessor(
            batch_size=self.config.batch_size,
            max_workers=self.config.max_workers
        )
        
        # Cache for metabolic modules
        self.metabolic_modules: Optional[List[str]] = None
    
    def validate_input(self, data: Dict[str, Any]) -> None:
        """
        Validate input parameters.
        
        Args:
            data: Input parameters dictionary
            
        Raises:
            ValidationError: If validation fails
        """
        required_fields = ["organism_code", "root_folder", "metabolic_modules_path"]
        for field in required_fields:
            if field not in data:
                raise ValidationError(f"Missing required field: {field}")
        
        # Validate organism code format
        org_code = data["organism_code"]
        if not org_code or not isinstance(org_code, str) or len(org_code) < 2:
            raise ValidationError(
                f"Invalid organism code: {org_code}",
                field_name="organism_code",
                invalid_value=org_code
            )
        
        # Validate paths
        root_folder = Path(data["root_folder"])
        if not root_folder.exists():
            if data.get("create_directories", True):
                root_folder.mkdir(parents=True, exist_ok=True)
                self.logger.info(f"Created root directory: {root_folder}")
            else:
                raise ValidationError(f"Root folder does not exist: {root_folder}")
        
        metabolic_modules_file = Path(data["metabolic_modules_path"])
        self.validate_file_exists(metabolic_modules_file, "metabolic modules file")
        
        # Validate file format
        if not metabolic_modules_file.suffix.lower() == '.csv':
            raise ValidationError(
                f"Metabolic modules file must be CSV format, got: {metabolic_modules_file.suffix}"
            )
    
    def load_metabolic_modules(self, modules_path: Path) -> List[str]:
        """
        Load metabolic module IDs from CSV file.
        
        Args:
            modules_path: Path to metabolic modules CSV file
            
        Returns:
            List of metabolic module IDs
            
        Raises:
            ValidationError: If file format is invalid
        """
        try:
            df = self.csv_handler.read_csv(modules_path)
            
            # Validate required columns
            if "Module_ID" not in df.columns:
                raise ValidationError(
                    "Metabolic modules CSV must contain 'Module_ID' column",
                    context={"columns": list(df.columns)}
                )
            
            module_ids = df["Module_ID"].astype(str).tolist()
            self.logger.info(f"Loaded {len(module_ids)} metabolic modules")
            
            return module_ids
            
        except Exception as e:
            if isinstance(e, ValidationError):
                raise
            raise ProcessingError(
                f"Failed to load metabolic modules: {e}",
                processor_name=self.__class__.__name__,
                stage="load_modules",
                original_error=e
            )
    
    def is_metabolic_module(self, module_id: str) -> bool:
        """Check if module ID is in the metabolic modules list."""
        if self.metabolic_modules is None:
            return False
        return module_id in self.metabolic_modules
    
    def check_existing_output(self, organism_folder: Path) -> bool:
        """
        Check if organism has already been processed.
        
        Args:
            organism_folder: Path to organism output folder
            
        Returns:
            True if organism already processed and skip_existing is True
        """
        if not self.config.skip_existing:
            return False
        
        if not organism_folder.exists():
            return False
        
        # Check for any CSV files in the folder
        csv_files = list(organism_folder.glob("*.csv"))
        if csv_files:
            self.logger.info(f"⏩ Skipping {organism_folder.name} (already processed)")
            return True
        
        return False
    
    def fetch_organism_modules(self, organism_code: str) -> Dict[str, List[str]]:
        """
        Fetch module-gene mappings for organism from KEGG.
        
        Args:
            organism_code: KEGG organism code
            
        Returns:
            Dictionary mapping module IDs to gene lists
            
        Raises:
            APIError: If KEGG API request fails
        """
        try:
            self.logger.info(f"Fetching modules for organism: {organism_code}")
            
            # Get module links for organism
            response = self.kegg_client.get(f"link/module/{organism_code}")
            
            if response.status_code != 200:
                raise APIError(
                    f"KEGG API returned status {response.status_code}",
                    api_name="KEGG",
                    endpoint=f"link/module/{organism_code}",
                    status_code=response.status_code
                )
            
            # Parse response
            module_dict = {}
            lines = response.data.strip().split("\n")
            
            for line in lines:
                if "\t" not in line:
                    continue
                
                try:
                    gene, module = line.split("\t")
                    gene_id = gene.split(":")[1]
                    module_id = module.split(":")[1]
                    
                    # Extract general module ID (remove organism prefix)
                    general_module_id = module_id.split("_")[1] if "_" in module_id else module_id
                    
                    # Only include metabolic modules
                    if self.is_metabolic_module(general_module_id):
                        module_dict.setdefault(module_id, []).append(gene_id)
                        
                except (IndexError, ValueError) as e:
                    self.logger.warning(f"Failed to parse line: {line} - {e}")
                    continue
            
            self.logger.info(f"Found {len(module_dict)} metabolic modules for {organism_code}")
            return module_dict
            
        except Exception as e:
            if isinstance(e, APIError):
                raise
            raise APIError(
                f"Failed to fetch organism modules: {e}",
                api_name="KEGG",
                endpoint=f"link/module/{organism_code}",
                original_error=e
            )
    
    def fetch_gene_batch_info(self, gene_batch: List[str], organism_code: str) -> Dict[str, GeneInfo]:
        """
        Fetch gene information for a batch of genes.
        
        Args:
            gene_batch: List of gene IDs
            organism_code: KEGG organism code
            
        Returns:
            Dictionary mapping gene IDs to GeneInfo objects
        """
        try:
            # Construct batch query
            query_genes = [f"{organism_code}:{gene}" for gene in gene_batch]
            query = "+".join(query_genes)
            
            self.logger.debug(f"Fetching info for {len(gene_batch)} genes")
            
            # Get gene information
            response = self.kegg_client.get(f"get/{query}")
            
            if response.status_code != 200:
                self.logger.error(f"KEGG API error for batch: status {response.status_code}")
                # Return N/A entries for failed batch
                return {gene: GeneInfo.from_kegg_data(gene, ("N/A", "N/A", "N/A")) 
                       for gene in gene_batch}
            
            # Parse response entries
            gene_info = {}
            entries = response.data.split("///")
            
            for entry in entries:
                if not entry.strip():
                    continue
                
                parsed = self.parser.parse_kegg_entry(entry)
                gene_id, annotation, sequence, location = self.parser.extract_gene_info(parsed)
                
                if gene_id:
                    gene_info[gene_id] = GeneInfo.from_kegg_data(
                        gene_id, (annotation, sequence, location)
                    )
            
            # Fill in missing genes with N/A
            for gene in gene_batch:
                if gene not in gene_info:
                    gene_info[gene] = GeneInfo.from_kegg_data(gene, ("N/A", "N/A", "N/A"))
            
            return gene_info
            
        except Exception as e:
            self.logger.error(f"Error fetching gene batch info: {e}")
            # Return N/A entries for all genes in failed batch
            return {gene: GeneInfo.from_kegg_data(gene, ("N/A", "N/A", "N/A")) 
                   for gene in gene_batch}
    
    def process_module_genes(self, module_id: str, genes: List[str], organism_code: str) -> ModuleData:
        """
        Process all genes for a single module.
        
        Args:
            module_id: Module ID
            genes: List of gene IDs in the module
            organism_code: KEGG organism code
            
        Returns:
            ModuleData object with processed gene information
        """
        self.logger.info(f"Processing module {module_id} with {len(genes)} genes")
        
        # Create batches for parallel processing
        def batch_function(gene_batch):
            return self.fetch_gene_batch_info(gene_batch, organism_code)
        
        # Process batches
        batch_results = self.batch_processor.process_parallel(
            genes,
            batch_function,
            description=f"Processing genes for module {module_id}"
        )
        
        # Collect all gene information
        all_gene_info = []
        for result in batch_results:
            if result.success:
                gene_info_dict = result.result
                all_gene_info.extend(gene_info_dict.values())
            else:
                self.logger.error(f"Batch failed: {result.error}")
        
        return ModuleData(
            module_id=module_id,
            genes=all_gene_info,
            organism_code=organism_code
        )
    
    def save_module_data(self, module_data: ModuleData, output_folder: Path) -> Path:
        """
        Save module data to CSV file.
        
        Args:
            module_data: Module data to save
            output_folder: Output folder path
            
        Returns:
            Path to saved file
        """
        output_file = output_folder / f"{module_data.module_id}.csv"
        
        # Convert to DataFrame and save
        df = module_data.to_dataframe()
        self.csv_handler.write_csv(df, output_file)
        
        self.logger.debug(f"Saved {len(module_data.genes)} genes to {output_file}")
        return output_file
    
    def process(self, data: Dict[str, Any], **kwargs) -> Dict[str, Any]:
        """
        Process organism to extract metabolic module information.
        
        Args:
            data: Input parameters
            **kwargs: Additional processing parameters
            
        Returns:
            Dictionary with processing results
        """
        organism_code = data["organism_code"]
        root_folder = Path(data["root_folder"])
        metabolic_modules_path = Path(data["metabolic_modules_path"])
        
        # Load metabolic modules
        self.metabolic_modules = self.load_metabolic_modules(metabolic_modules_path)
        
        # Setup organism output folder
        organism_folder = root_folder / organism_code
        
        # Check if already processed
        if self.check_existing_output(organism_folder):
            return {
                "organism_code": organism_code,
                "status": "skipped",
                "reason": "already_processed",
                "modules_processed": 0,
                "genes_processed": 0
            }
        
        # Create output folder
        organism_folder.mkdir(parents=True, exist_ok=True)
        
        # Fetch organism modules
        module_dict = self.fetch_organism_modules(organism_code)
        
        if not module_dict:
            self.logger.warning(f"No metabolic modules found for {organism_code}")
            return {
                "organism_code": organism_code,
                "status": "completed",
                "modules_processed": 0,
                "genes_processed": 0
            }
        
        # Process each module
        total_genes = 0
        processed_modules = 0
        
        tracker = ProgressTracker(len(module_dict), f"Processing modules for {organism_code}")
        
        for module_id, genes in module_dict.items():
            try:
                # Process module genes
                module_data = self.process_module_genes(module_id, genes, organism_code)
                
                # Save to file
                self.save_module_data(module_data, organism_folder)
                
                total_genes += len(genes)
                processed_modules += 1
                tracker.update(1)
                
            except Exception as e:
                self.logger.error(f"Failed to process module {module_id}: {e}")
                tracker.mark_failed(1)
        
        # Log completion
        self.logger.info(f"✅ Completed {organism_code}: {processed_modules} modules, {total_genes} genes")
        
        return {
            "organism_code": organism_code,
            "status": "completed",
            "modules_processed": processed_modules,
            "genes_processed": total_genes,
            "output_folder": str(organism_folder),
            "performance_stats": tracker.get_stats()
        }


def main():
    """
    Main entry point for organism processing script.
    
    Provides command-line interface with professional argument handling
    and comprehensive error reporting.
    """
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)
    
    # Validate command line arguments
    if len(sys.argv) < 4:
        print("Usage: python process_organism.py <org_code> <root_folder> <metabolic_modules_path>")
        print("\nArguments:")
        print("  org_code               : KEGG organism code (e.g., 'hsa' for human)")
        print("  root_folder           : Root directory for output files")
        print("  metabolic_modules_path: Path to CSV file with metabolic module IDs")
        print("\nExample:")
        print("  python process_organism.py hsa /data/organisms modules.csv")
        sys.exit(1)
    
    # Parse arguments
    organism_code = sys.argv[1]
    root_folder = sys.argv[2]
    metabolic_modules_path = sys.argv[3]
    
    # Create processor with configuration
    config = OrganismProcessingConfig(
        batch_size=10,
        max_workers=4,
        skip_existing=True
    )
    
    processor = OrganismProcessor(config)
    
    try:
        # Prepare input data
        input_data = {
            "organism_code": organism_code,
            "root_folder": root_folder,
            "metabolic_modules_path": metabolic_modules_path
        }
        
        # Run processing
        logger.info(f"Starting organism processing: {organism_code}")
        result = processor.run(input_data)
        
        if result.status == ProcessingStatus.COMPLETED:
            logger.info("Processing completed successfully")
            if result.results:
                stats = result.results
                print(f"✅ {stats['organism_code']}: {stats['modules_processed']} modules, "
                     f"{stats['genes_processed']} genes processed")
        else:
            logger.error(f"Processing failed: {result.error_message}")
            sys.exit(1)
    
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()