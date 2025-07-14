"""
Industry-Level KEGG Metabolic Module Checker and Organism Processor.

This module provides comprehensive organism processing with metabolic module detection,
parallel processing, caching, and professional error handling.
"""

import sys
from typing import Dict, List, Optional, Set, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import time

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
class OrganismInfo:
    """Information about a KEGG organism."""
    
    code: str
    name: str
    division: str = ""
    lineage: str = ""
    modules_processed: int = 0
    genes_processed: int = 0
    processing_time: float = 0.0


@dataclass
class KEGGOrganismConfig:
    """Configuration for KEGG organism processing."""
    
    batch_size: int = 10
    rate_limit_delay: float = 0.35
    max_retries: int = 3
    timeout: int = 30
    skip_existing: bool = True
    organism_filter: str = "Plants"
    save_organism_list: bool = True


class KEGGOrganismAPIClient:
    """
    Professional KEGG API client for organism and module data retrieval.
    """
    
    def __init__(self, 
                 base_url: str = "https://rest.kegg.jp",
                 rate_limit_delay: float = 0.35,
                 max_retries: int = 3,
                 timeout: int = 30):
        """
        Initialize KEGG organism API client.
        
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
            'User-Agent': 'KEGGMetabolicModuleChecker/1.0'
        })
        
        self._last_request_time = 0.0
        self._module_cache: Dict[str, bool] = {}
        self._organism_cache: Dict[str, OrganismInfo] = {}
    
    def _enforce_rate_limit(self) -> None:
        """Enforce rate limiting between requests."""
        current_time = time.time()
        time_since_last = current_time - self._last_request_time
        
        if time_since_last < self.rate_limit_delay:
            sleep_time = self.rate_limit_delay - time_since_last
            time.sleep(sleep_time)
        
        self._last_request_time = time.time()
    
    def get_organisms(self, organism_filter: Optional[str] = None) -> Dict[str, OrganismInfo]:
        """
        Get KEGG organisms list with optional filtering.
        
        Args:
            organism_filter: Filter organisms by division (e.g., "Plants")
            
        Returns:
            Dictionary mapping organism codes to OrganismInfo objects
        """
        if self._organism_cache and organism_filter in ["Plants", None]:
            return self._organism_cache
        
        try:
            self._enforce_rate_limit()
            
            url = f"{self.base_url}/list/organism"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            
            organisms = {}
            
            for line in response.text.strip().split("\n"):
                parts = line.split("\t")
                if len(parts) >= 4:
                    org_code = parts[1]
                    org_name = parts[2]
                    division = parts[3]
                    
                    # Apply filter if specified
                    if organism_filter is None or organism_filter in division:
                        organisms[org_code] = OrganismInfo(
                            code=org_code,
                            name=org_name,
                            division=division
                        )
            
            # Cache results
            if organism_filter == "Plants":
                self._organism_cache = organisms
            
            return organisms
            
        except Exception as e:
            raise RuntimeError(f"Failed to fetch organisms: {e}")
    
    def is_metabolic_module(self, module_id: str) -> bool:
        """
        Check if module is metabolic using cached results.
        
        Args:
            module_id: KEGG module ID
            
        Returns:
            True if module is metabolic
        """
        # Extract general module ID (remove organism prefix)
        general_module_id = module_id.split("_")[-1] if "_" in module_id else module_id
        
        # Check cache first
        if general_module_id in self._module_cache:
            return self._module_cache[general_module_id]
        
        try:
            self._enforce_rate_limit()
            
            url = f"{self.base_url}/get/{general_module_id}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            
            # Check if module classification contains "metabol"
            is_metabolic = any(
                "metabol" in line.lower()
                for line in response.text.splitlines()
                if line.startswith("CLASS")
            )
            
            # Cache result
            self._module_cache[general_module_id] = is_metabolic
            return is_metabolic
            
        except Exception as e:
            # Cache as non-metabolic on failure
            self._module_cache[general_module_id] = False
            return False
    
    def get_organism_modules(self, organism_code: str) -> str:
        """Get module links for organism."""
        self._enforce_rate_limit()
        
        url = f"{self.base_url}/link/module/{organism_code}"
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        
        return response.text.strip()
    
    def get_gene_batch(self, gene_ids: List[str], organism_code: str) -> str:
        """Get information for batch of genes."""
        self._enforce_rate_limit()
        
        query_ids = [f"{organism_code}:{gene_id}" for gene_id in gene_ids]
        query = "+".join(query_ids)
        
        url = f"{self.base_url}/get/{query}"
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        
        return response.text.strip()
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        return {
            "module_cache_size": len(self._module_cache),
            "organism_cache_size": len(self._organism_cache),
            "metabolic_modules_cached": sum(self._module_cache.values()),
            "non_metabolic_modules_cached": len(self._module_cache) - sum(self._module_cache.values())
        }


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
    def parse_gene_batch_response(cls, response_text: str) -> Dict[str, Tuple[str, str, str]]:
        """
        Parse batch gene response.
        
        Args:
            response_text: Raw API response
            
        Returns:
            Dictionary mapping gene IDs to (orthology, sequence, position) tuples
        """
        genes = {}
        
        for entry in response_text.split("///"):
            if not entry.strip():
                continue
            
            try:
                parsed = cls.parse_flat_entry(entry)
                gene_id = parsed.get("ENTRY", "").split()[0]
                
                if gene_id:
                    genes[gene_id] = (
                        parsed.get("ORTHOLOGY", ""),
                        parsed.get("AASEQ", ""),
                        parsed.get("POSITION", "")
                    )
            
            except Exception:
                # Skip problematic entries
                continue
        
        return genes


class KEGGMetabolicModuleChecker(DatabaseProcessor):
    """
    Professional KEGG metabolic module checker and organism processor.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 organism_config: KEGGOrganismConfig):
        """
        Initialize KEGG metabolic module checker.
        
        Args:
            config: Base processing configuration
            organism_config: Organism-specific configuration
        """
        super().__init__(config, name="KEGGMetabolicModuleChecker")
        
        self.organism_config = organism_config
        self.api_client = KEGGOrganismAPIClient(
            rate_limit_delay=organism_config.rate_limit_delay,
            max_retries=organism_config.max_retries,
            timeout=organism_config.timeout
        )
        self.parser = KEGGEntryParser()
        
        # Statistics tracking
        self.processing_stats = {
            "organisms_found": 0,
            "organisms_processed": 0,
            "organisms_skipped": 0,
            "total_modules_processed": 0,
            "total_genes_processed": 0,
            "api_requests_made": 0,
            "errors_encountered": 0,
        }
    
    def process(self, input_data: Any = None, **kwargs) -> ProcessingResult:
        """
        Process all KEGG organisms and extract metabolic module data.
        
        Args:
            input_data: Optional input data (not used)
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info(f"Starting KEGG metabolic module checking",
                           organism_filter=self.organism_config.organism_filter)
            
            # Get organism list
            organisms = self._get_organisms()
            self.processing_stats["organisms_found"] = len(organisms)
            
            if not organisms:
                result.add_error("No organisms found matching filter criteria")
                return result
            
            self.logger.info(f"Found organisms to process", count=len(organisms))
            
            # Save organism list if requested
            if self.organism_config.save_organism_list:
                self._save_organism_list(organisms)
            
            # Process organisms
            organism_results = self._process_organisms_parallel(organisms)
            
            # Aggregate results
            successful_organisms = []
            total_modules = 0
            total_genes = 0
            
            for org_result in organism_results:
                if org_result.success:
                    successful_organisms.append(org_result.metadata.get("organism_info"))
                    total_modules += org_result.metadata.get("modules_processed", 0)
                    total_genes += org_result.metadata.get("genes_processed", 0)
                    
                    if org_result.metadata.get("status") == "skipped":
                        self.processing_stats["organisms_skipped"] += 1
                    else:
                        self.processing_stats["organisms_processed"] += 1
                else:
                    result.errors.extend(org_result.errors)
                    result.warnings.extend(org_result.warnings)
                    self.processing_stats["errors_encountered"] += 1
            
            # Update statistics
            self.processing_stats["total_modules_processed"] = total_modules
            self.processing_stats["total_genes_processed"] = total_genes
            
            result.metadata.update({
                "organisms_found": len(organisms),
                "organisms_processed": len(successful_organisms),
                "organisms_skipped": self.processing_stats["organisms_skipped"],
                "total_modules_processed": total_modules,
                "total_genes_processed": total_genes,
                "successful_organisms": successful_organisms,
                "processing_stats": self.processing_stats,
                "api_cache_stats": self.api_client.get_cache_stats()
            })
            
            self.logger.info(f"KEGG metabolic module checking completed",
                           processed=len(successful_organisms),
                           modules=total_modules,
                           genes=total_genes)
            
        except Exception as e:
            result.add_error(f"Processing failed: {str(e)}")
            self.processing_stats["errors_encountered"] += 1
        
        return result
    
    def process_single_organism(self, organism_code: str) -> ProcessingResult:
        """
        Process single organism.
        
        Args:
            organism_code: KEGG organism code
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info(f"Processing single organism", organism=organism_code)
            
            # Get organism info
            organisms = self.api_client.get_organisms(self.organism_config.organism_filter)
            
            if organism_code not in organisms:
                result.add_error(f"Organism not found: {organism_code}")
                return result
            
            organism_info = organisms[organism_code]
            
            # Process organism
            org_result = self._process_single_organism(organism_info)
            result = result.merge(org_result)
            
        except Exception as e:
            result.add_error(f"Single organism processing failed: {str(e)}")
        
        return result
    
    def _get_organisms(self) -> Dict[str, OrganismInfo]:
        """Get organisms matching filter criteria."""
        try:
            organisms = self.api_client.get_organisms(self.organism_config.organism_filter)
            self.processing_stats["api_requests_made"] += 1
            
            self.logger.debug(f"Retrieved organisms", 
                            count=len(organisms),
                            filter=self.organism_config.organism_filter)
            
            return organisms
            
        except Exception as e:
            self.logger.error(f"Failed to get organisms", error=str(e))
            raise
    
    def _save_organism_list(self, organisms: Dict[str, OrganismInfo]) -> None:
        """Save organism list to CSV file."""
        try:
            organism_file = self.config.output_dir / "plants_list.csv"
            
            with open(organism_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Organism_Code", "Organism_Name"])
                
                for organism_info in organisms.values():
                    writer.writerow([organism_info.code, organism_info.name])
            
            self.logger.info(f"Organism list saved", 
                           path=str(organism_file),
                           count=len(organisms))
            
        except Exception as e:
            self.logger.warning(f"Failed to save organism list", error=str(e))
    
    def _process_organisms_parallel(self, organisms: Dict[str, OrganismInfo]) -> List[ProcessingResult]:
        """
        Process organisms in parallel.
        
        Args:
            organisms: Dictionary of organisms to process
            
        Returns:
            List of processing results
        """
        organism_list = list(organisms.values())
        
        if self.config.max_workers == 1:
            # Sequential processing
            return [self._process_single_organism(org) for org in organism_list]
        
        # Parallel processing
        results = []
        
        with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
            future_to_organism = {
                executor.submit(self._process_single_organism, organism): organism
                for organism in organism_list
            }
            
            for future in as_completed(future_to_organism):
                organism = future_to_organism[future]
                try:
                    result = future.result(timeout=self.config.timeout)
                    results.append(result)
                except Exception as e:
                    error_result = ProcessingResult(success=False)
                    error_result.add_error(f"Processing failed for {organism.code}: {str(e)}")
                    results.append(error_result)
        
        return results
    
    def _process_single_organism(self, organism_info: OrganismInfo) -> ProcessingResult:
        """
        Process single organism.
        
        Args:
            organism_info: Organism information
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        start_time = time.time()
        
        try:
            organism_code = organism_info.code
            organism_dir = self.config.output_dir / organism_code
            
            # Check if already processed
            if (self.organism_config.skip_existing and
                organism_dir.exists() and
                any(organism_dir.glob("*.csv"))):
                
                self.logger.debug(f"Organism already processed, skipping", organism=organism_code)
                result.metadata.update({
                    "organism_info": organism_info,
                    "status": "skipped",
                    "modules_processed": 0,
                    "genes_processed": 0
                })
                return result
            
            organism_dir.mkdir(parents=True, exist_ok=True)
            
            # Get organism modules
            module_data = self.api_client.get_organism_modules(organism_code)
            self.processing_stats["api_requests_made"] += 1
            
            # Build metabolic module dictionary
            module_dict = {}
            
            for line in module_data.split("\n"):
                if "\t" not in line:
                    continue
                
                gene, module = line.split("\t")
                gene_id = gene.split(":")[1]
                module_id = module.split(":")[1]
                
                if self.api_client.is_metabolic_module(module_id):
                    if module_id not in module_dict:
                        module_dict[module_id] = []
                    module_dict[module_id].append(gene_id)
            
            # Process each metabolic module
            modules_processed = 0
            total_genes = 0
            
            for module_id, gene_ids in module_dict.items():
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
                    self.logger.warning(f"Failed to process module", 
                                      organism=organism_code,
                                      module=module_id,
                                      error=str(e))
            
            # Update organism info
            organism_info.modules_processed = modules_processed
            organism_info.genes_processed = total_genes
            organism_info.processing_time = time.time() - start_time
            
            result.metadata.update({
                "organism_info": organism_info,
                "status": "completed",
                "modules_processed": modules_processed,
                "genes_processed": total_genes
            })
            
            self.logger.debug(f"Organism processing completed",
                            organism=organism_code,
                            modules=modules_processed,
                            genes=total_genes)
            
        except Exception as e:
            result.add_error(f"Organism processing failed: {str(e)}")
        
        return result
    
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
            
            for i in range(0, len(gene_ids), self.organism_config.batch_size):
                batch = gene_ids[i:i + self.organism_config.batch_size]
                future = executor.submit(self._fetch_gene_batch, batch, organism_code)
                futures.append(future)
            
            for future in as_completed(futures):
                try:
                    gene_data = future.result(timeout=self.config.timeout)
                    
                    # Write batch results to CSV
                    with open(module_file, "a", newline="") as f:
                        writer = csv.writer(f)
                        for gene_id, (orthology, sequence, position) in gene_data.items():
                            writer.writerow([gene_id, orthology, sequence, position])
                    
                    genes_processed += len(gene_data)
                    self.processing_stats["api_requests_made"] += 1
                    
                except Exception as e:
                    self.logger.warning(f"Failed to process gene batch", 
                                      module=module_id,
                                      error=str(e))
        
        return genes_processed
    
    def _fetch_gene_batch(self, 
                         gene_ids: List[str], 
                         organism_code: str) -> Dict[str, Tuple[str, str, str]]:
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
            return {gene_id: ("N/A", "N/A", "N/A") for gene_id in gene_ids}
    
    def get_comprehensive_statistics(self) -> Dict[str, Any]:
        """Get comprehensive processing statistics."""
        base_stats = self.get_statistics()
        api_stats = self.api_client.get_cache_stats()
        
        return {
            **base_stats,
            **self.processing_stats,
            "api_cache_stats": api_stats,
            "processing_config": {
                "organism_filter": self.organism_config.organism_filter,
                "batch_size": self.organism_config.batch_size,
                "rate_limit_delay": self.organism_config.rate_limit_delay,
                "skip_existing": self.organism_config.skip_existing,
            }
        }
    
    def validate_input(self, input_data: Any = None) -> bool:
        """
        Validate input configuration.
        
        Args:
            input_data: Input data (not used)
            
        Returns:
            True if configuration is valid
        """
        return self.config.output_dir.exists()


def main():
    """Main entry point for KEGG metabolic module checking."""
    # Check for single organism processing
    if len(sys.argv) == 2:
        organism_code = sys.argv[1]
        single_organism = True
    else:
        organism_code = None
        single_organism = False
    
    # Configuration
    output_folder = Path("/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic_local")
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=output_folder,
        output_dir=output_folder,
        max_workers=4,
        log_level="INFO"
    )
    
    organism_config = KEGGOrganismConfig(
        batch_size=10,
        rate_limit_delay=0.35,
        organism_filter="Plants",
        skip_existing=True,
        save_organism_list=True
    )
    
    # Initialize checker
    try:
        checker = KEGGMetabolicModuleChecker(
            config=config,
            organism_config=organism_config
        )
        
        # Process organisms
        if single_organism:
            result = checker.process_single_organism(organism_code)
            print(f"Processing organism: {organism_code}")
        else:
            result = checker.run()
            print("Processing all organisms...")
        
        # Display results
        if result.success:
            if single_organism:
                print(f"✅ Successfully processed organism: {organism_code}")
                print(f"   Modules processed: {result.metadata.get('modules_processed', 0)}")
                print(f"   Genes processed: {result.metadata.get('genes_processed', 0)}")
            else:
                print("✅ All metabolic module-based KEGG data collected!")
                print(f"   Organisms found: {result.metadata.get('organisms_found', 0)}")
                print(f"   Organisms processed: {result.metadata.get('organisms_processed', 0)}")
                print(f"   Organisms skipped: {result.metadata.get('organisms_skipped', 0)}")
                print(f"   Total modules processed: {result.metadata.get('total_modules_processed', 0)}")
                print(f"   Total genes processed: {result.metadata.get('total_genes_processed', 0)}")
            
            print(f"   Processing time: {result.processing_time:.2f} seconds")
        else:
            print("❌ Processing failed!")
            for error in result.errors:
                print(f"   Error: {error}")
        
        # Display warnings
        for warning in result.warnings:
            print(f"   Warning: {warning}")
        
        # Display statistics
        stats = checker.get_comprehensive_statistics()
        print(f"\nProcessing Statistics:")
        print(f"   API requests made: {stats['api_requests_made']}")
        print(f"   Module cache size: {stats['api_cache_stats']['module_cache_size']}")
        print(f"   Metabolic modules cached: {stats['api_cache_stats']['metabolic_modules_cached']}")
        
    except Exception as e:
        print(f"❌ Critical error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())