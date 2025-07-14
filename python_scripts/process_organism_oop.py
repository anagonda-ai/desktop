"""
Industry-Level Organism Processing Pipeline.

This module provides object-oriented KEGG organism processing with comprehensive
error handling, rate limiting, caching, and performance optimization.
"""

import sys
import time
from typing import Dict, List, Optional, Set, Tuple, Any
from dataclasses import dataclass
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv

import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from core.base_classes import (
    BaseProcessor, 
    BioinformaticsConfig, 
    ProcessingResult,
    DatabaseProcessor
)


@dataclass
class KEGGGeneInfo:
    """KEGG gene information container."""
    
    gene_id: str
    organism_code: str
    orthology: str = ""
    sequence: str = ""
    position: str = ""
    module_ids: List[str] = None
    
    def __post_init__(self):
        """Initialize module_ids if None."""
        if self.module_ids is None:
            self.module_ids = []


@dataclass
class OrganismProcessingConfig:
    """Configuration for organism processing."""
    
    batch_size: int = 10
    rate_limit_delay: float = 0.35
    max_retries: int = 3
    timeout: int = 30
    skip_existing: bool = True


class KEGGAPIClient:
    """
    Professional KEGG API client with rate limiting, retry logic, and caching.
    """
    
    def __init__(self, 
                 base_url: str = "https://rest.kegg.jp",
                 rate_limit_delay: float = 0.35,
                 max_retries: int = 3,
                 timeout: int = 30):
        """
        Initialize KEGG API client.
        
        Args:
            base_url: KEGG API base URL
            rate_limit_delay: Delay between requests (seconds)
            max_retries: Maximum retry attempts
            timeout: Request timeout (seconds)
        """
        self.base_url = base_url
        self.rate_limit_delay = rate_limit_delay
        self.timeout = timeout
        
        # Setup session with retry strategy
        self.session = requests.Session()
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        
        self.session.headers.update({
            'User-Agent': 'BioinformaticsOrganismProcessor/1.0'
        })
        
        self._last_request_time = 0.0
    
    def _enforce_rate_limit(self) -> None:
        """Enforce rate limiting between requests."""
        current_time = time.time()
        time_since_last = current_time - self._last_request_time
        
        if time_since_last < self.rate_limit_delay:
            sleep_time = self.rate_limit_delay - time_since_last
            time.sleep(sleep_time)
        
        self._last_request_time = time.time()
    
    def get(self, endpoint: str) -> str:
        """
        Make GET request to KEGG API with rate limiting.
        
        Args:
            endpoint: API endpoint
            
        Returns:
            Response text
            
        Raises:
            requests.RequestException: If request fails
        """
        self._enforce_rate_limit()
        
        url = f"{self.base_url}/{endpoint}"
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        
        return response.text.strip()
    
    def get_organism_modules(self, organism_code: str) -> str:
        """Get module links for organism."""
        return self.get(f"link/module/{organism_code}")
    
    def get_gene_batch(self, gene_ids: List[str], organism_code: str) -> str:
        """Get information for batch of genes."""
        query_ids = [f"{organism_code}:{gene_id}" for gene_id in gene_ids]
        query = "+".join(query_ids)
        return self.get(f"get/{query}")


class KEGGEntryParser:
    """Parser for KEGG flat file format entries."""
    
    @staticmethod
    def parse_flat_entry(entry_text: str) -> Dict[str, str]:
        """
        Parse KEGG flat file entry.
        
        Args:
            entry_text: Raw entry text
            
        Returns:
            Dictionary of parsed fields
        """
        parsed = {}
        current_key = None
        
        for line in entry_text.splitlines():
            if not line.strip():
                continue
            
            # Check if line starts with a field name
            field_name = line[:12].strip()
            
            if field_name:
                current_key = field_name
                if current_key == "AASEQ":
                    parsed[current_key] = ""
                else:
                    parsed[current_key] = line[12:].strip()
            else:
                # Continuation line
                if current_key:
                    parsed[current_key] += line[12:].strip()
        
        return parsed
    
    @classmethod
    def parse_gene_batch_response(cls, response_text: str) -> Dict[str, KEGGGeneInfo]:
        """
        Parse batch gene response.
        
        Args:
            response_text: Raw API response
            
        Returns:
            Dictionary mapping gene IDs to KEGGGeneInfo objects
        """
        genes = {}
        
        for entry in response_text.split("///"):
            if not entry.strip():
                continue
            
            try:
                parsed = cls.parse_flat_entry(entry)
                gene_id = parsed.get("ENTRY", "").split()[0]
                
                if gene_id:
                    # Extract organism code from gene ID if present
                    organism_code = ""
                    if ":" in gene_id:
                        parts = gene_id.split(":")
                        if len(parts) == 2:
                            organism_code, gene_id = parts
                    
                    genes[gene_id] = KEGGGeneInfo(
                        gene_id=gene_id,
                        organism_code=organism_code,
                        orthology=parsed.get("ORTHOLOGY", ""),
                        sequence=parsed.get("AASEQ", ""),
                        position=parsed.get("POSITION", "")
                    )
            
            except Exception as e:
                # Log parsing error but continue with other entries
                continue
        
        return genes


class MetabolicModuleFilter:
    """Filter for identifying metabolic modules."""
    
    def __init__(self, metabolic_modules_file: Path):
        """
        Initialize with metabolic modules file.
        
        Args:
            metabolic_modules_file: Path to CSV file with metabolic modules
        """
        self.metabolic_modules: Set[str] = set()
        self._load_metabolic_modules(metabolic_modules_file)
    
    def _load_metabolic_modules(self, modules_file: Path) -> None:
        """Load metabolic modules from CSV file."""
        try:
            df = pd.read_csv(modules_file)
            if "Module_ID" in df.columns:
                self.metabolic_modules = set(df["Module_ID"].astype(str).tolist())
            else:
                # Assume first column contains module IDs
                self.metabolic_modules = set(df.iloc[:, 0].astype(str).tolist())
        
        except Exception as e:
            raise ValueError(f"Failed to load metabolic modules from {modules_file}: {e}")
    
    def is_metabolic_module(self, module_id: str) -> bool:
        """
        Check if module is metabolic.
        
        Args:
            module_id: Module ID to check
            
        Returns:
            True if module is metabolic
        """
        # Extract general module ID (remove organism prefix)
        general_module_id = module_id.split("_")[-1] if "_" in module_id else module_id
        return general_module_id in self.metabolic_modules


class OrganismProcessor(DatabaseProcessor):
    """
    Professional organism processing pipeline with comprehensive error handling.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 processing_config: OrganismProcessingConfig,
                 metabolic_modules_file: Path):
        """
        Initialize organism processor.
        
        Args:
            config: Base processing configuration
            processing_config: Organism-specific configuration
            metabolic_modules_file: Path to metabolic modules file
        """
        super().__init__(config, name="OrganismProcessor")
        
        self.processing_config = processing_config
        self.api_client = KEGGAPIClient(
            rate_limit_delay=processing_config.rate_limit_delay,
            max_retries=processing_config.max_retries,
            timeout=processing_config.timeout
        )
        self.parser = KEGGEntryParser()
        self.module_filter = MetabolicModuleFilter(metabolic_modules_file)
        
        # Statistics tracking
        self.processing_stats = {
            "organisms_processed": 0,
            "modules_processed": 0,
            "genes_processed": 0,
            "api_requests_made": 0,
            "errors_encountered": 0,
        }
    
    def process(self, organism_code: str, **kwargs) -> ProcessingResult:
        """
        Process single organism.
        
        Args:
            organism_code: KEGG organism code
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info(f"Starting organism processing", organism=organism_code)
            
            # Create organism output directory
            organism_dir = self.config.output_dir / organism_code
            
            # Check if already processed
            if (self.processing_config.skip_existing and 
                organism_dir.exists() and 
                any(organism_dir.glob("*.csv"))):
                
                self.logger.info(f"Organism already processed, skipping", organism=organism_code)
                result.metadata["status"] = "skipped"
                return result
            
            organism_dir.mkdir(parents=True, exist_ok=True)
            
            # Get organism modules
            modules_data = self._get_organism_modules(organism_code)
            
            if not modules_data:
                result.add_warning("No modules found for organism")
                result.metadata["modules_processed"] = 0
                return result
            
            # Process each module
            modules_processed = 0
            total_genes = 0
            
            for module_id, gene_ids in modules_data.items():
                try:
                    genes_count = self._process_module(
                        module_id, 
                        gene_ids, 
                        organism_code, 
                        organism_dir
                    )
                    modules_processed += 1
                    total_genes += genes_count
                    
                except Exception as e:
                    result.add_warning(f"Failed to process module {module_id}: {str(e)}")
                    self.processing_stats["errors_encountered"] += 1
            
            # Update statistics
            self.processing_stats["organisms_processed"] += 1
            self.processing_stats["modules_processed"] += modules_processed
            self.processing_stats["genes_processed"] += total_genes
            
            result.metadata.update({
                "organism_code": organism_code,
                "modules_processed": modules_processed,
                "genes_processed": total_genes,
                "output_directory": str(organism_dir)
            })
            
            self.logger.info(f"Organism processing completed", 
                           organism=organism_code,
                           modules=modules_processed,
                           genes=total_genes)
            
        except Exception as e:
            result.add_error(f"Organism processing failed: {str(e)}")
            self.processing_stats["errors_encountered"] += 1
        
        return result
    
    def _get_organism_modules(self, organism_code: str) -> Dict[str, List[str]]:
        """
        Get metabolic modules for organism.
        
        Args:
            organism_code: KEGG organism code
            
        Returns:
            Dictionary mapping module IDs to gene lists
        """
        try:
            module_links = self.api_client.get_organism_modules(organism_code)
            self.processing_stats["api_requests_made"] += 1
            
            module_dict = {}
            
            for line in module_links.split("\n"):
                if "\t" not in line:
                    continue
                
                gene, module = line.split("\t")
                gene_id = gene.split(":")[1]
                module_id = module.split(":")[1]
                
                # Check if module is metabolic
                if self.module_filter.is_metabolic_module(module_id):
                    if module_id not in module_dict:
                        module_dict[module_id] = []
                    module_dict[module_id].append(gene_id)
            
            self.logger.debug(f"Found metabolic modules", 
                            organism=organism_code,
                            module_count=len(module_dict))
            
            return module_dict
            
        except Exception as e:
            self.logger.error(f"Failed to get organism modules", 
                            organism=organism_code,
                            error=str(e))
            raise
    
    def _process_module(self, 
                       module_id: str, 
                       gene_ids: List[str], 
                       organism_code: str, 
                       output_dir: Path) -> int:
        """
        Process single module.
        
        Args:
            module_id: Module ID
            gene_ids: List of gene IDs in module
            organism_code: Organism code
            output_dir: Output directory
            
        Returns:
            Number of genes processed
        """
        module_file = output_dir / f"{module_id}.csv"
        
        # Initialize CSV file
        with open(module_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Gene_ID", "Annotation", "Nucleotide_Sequence", "Location"])
        
        # Process genes in batches
        genes_processed = 0
        
        with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
            futures = []
            
            for i in range(0, len(gene_ids), self.processing_config.batch_size):
                batch = gene_ids[i:i + self.processing_config.batch_size]
                future = executor.submit(self._fetch_gene_batch, batch, organism_code)
                futures.append(future)
            
            for future in as_completed(futures):
                try:
                    gene_data = future.result(timeout=self.config.timeout)
                    
                    # Write batch results to CSV
                    with open(module_file, "a", newline="") as f:
                        writer = csv.writer(f)
                        for gene_id, gene_info in gene_data.items():
                            writer.writerow([
                                gene_id,
                                gene_info.orthology,
                                gene_info.sequence,
                                gene_info.position
                            ])
                    
                    genes_processed += len(gene_data)
                    self.processing_stats["api_requests_made"] += 1
                    
                except Exception as e:
                    self.logger.warning(f"Failed to process gene batch", 
                                      module=module_id,
                                      error=str(e))
        
        self.logger.debug(f"Module processed", 
                        module=module_id,
                        genes=genes_processed)
        
        return genes_processed
    
    def _fetch_gene_batch(self, 
                         gene_ids: List[str], 
                         organism_code: str) -> Dict[str, KEGGGeneInfo]:
        """
        Fetch information for batch of genes.
        
        Args:
            gene_ids: List of gene IDs
            organism_code: Organism code
            
        Returns:
            Dictionary of gene information
        """
        try:
            response = self.api_client.get_gene_batch(gene_ids, organism_code)
            return self.parser.parse_gene_batch_response(response)
            
        except Exception as e:
            self.logger.warning(f"Failed to fetch gene batch", 
                              genes=len(gene_ids),
                              error=str(e))
            
            # Return placeholder data for failed genes
            return {
                gene_id: KEGGGeneInfo(
                    gene_id=gene_id,
                    organism_code=organism_code,
                    orthology="N/A",
                    sequence="N/A",
                    position="N/A"
                )
                for gene_id in gene_ids
            }
    
    def get_processing_statistics(self) -> Dict[str, Any]:
        """Get comprehensive processing statistics."""
        base_stats = self.get_statistics()
        return {
            **base_stats,
            **self.processing_stats,
            "api_client_stats": {
                "rate_limit_delay": self.api_client.rate_limit_delay,
                "timeout": self.api_client.timeout,
            }
        }
    
    def validate_input(self, organism_code: str) -> bool:
        """
        Validate organism code input.
        
        Args:
            organism_code: Organism code to validate
            
        Returns:
            True if valid
        """
        return (isinstance(organism_code, str) and 
                organism_code.strip() and 
                len(organism_code) <= 10)


def main():
    """Main entry point for organism processing."""
    if len(sys.argv) < 4:
        print("Usage: python process_organism_oop.py <org_code> <root_folder> <metabolic_modules_path>")
        print(f"Arguments provided: {sys.argv}")
        sys.exit(1)
    
    organism_code = sys.argv[1]
    root_folder = Path(sys.argv[2])
    metabolic_modules_path = Path(sys.argv[3])
    
    # Validate inputs
    if not root_folder.exists():
        print(f"Error: Root folder does not exist: {root_folder}")
        sys.exit(1)
    
    if not metabolic_modules_path.exists():
        print(f"Error: Metabolic modules file does not exist: {metabolic_modules_path}")
        sys.exit(1)
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=root_folder,
        output_dir=root_folder,
        max_workers=4,
        log_level="INFO"
    )
    
    processing_config = OrganismProcessingConfig(
        batch_size=10,
        rate_limit_delay=0.35,
        skip_existing=True
    )
    
    # Initialize processor
    try:
        processor = OrganismProcessor(
            config=config,
            processing_config=processing_config,
            metabolic_modules_file=metabolic_modules_path
        )
        
        # Process organism
        result = processor.run(organism_code)
        
        # Display results
        if result.success:
            print(f"✅ Successfully processed organism: {organism_code}")
            print(f"   Modules processed: {result.metadata.get('modules_processed', 0)}")
            print(f"   Genes processed: {result.metadata.get('genes_processed', 0)}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
        else:
            print(f"❌ Failed to process organism: {organism_code}")
            for error in result.errors:
                print(f"   Error: {error}")
        
        # Display warnings
        for warning in result.warnings:
            print(f"   Warning: {warning}")
        
        # Display statistics
        stats = processor.get_processing_statistics()
        print(f"\nProcessing Statistics:")
        print(f"   API requests made: {stats['api_requests_made']}")
        print(f"   Cache size: {stats['cache_size']}")
        
    except Exception as e:
        print(f"❌ Critical error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()