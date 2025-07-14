"""
KEGG Database Integration for Plant MGC Analysis Pipeline.

This module provides object-oriented integration with the KEGG database for
metabolic pathway analysis and gene annotation retrieval.
"""

import asyncio
import time
from typing import Dict, List, Optional, Set, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import json

import requests
import aiohttp
import pandas as pd
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser

from ..core.base import DatabaseInterface, BioinformaticsProcessor
from ..core.types import DatabaseType, GeneInfo
from ..core.exceptions import DatabaseError, NetworkError, ValidationError
from ..utils.file_operations import CsvProcessor, SafeFileOperations


@dataclass
class KEGGModule:
    """Represents a KEGG metabolic module."""
    
    module_id: str
    name: str
    classification: str
    definition: str
    is_metabolic: bool = False
    genes: List[str] = field(default_factory=list)
    pathways: List[str] = field(default_factory=list)
    
    @property
    def general_id(self) -> str:
        """Get general module ID (without organism prefix)."""
        return self.module_id.split('_')[-1] if '_' in self.module_id else self.module_id


@dataclass
class KEGGGene:
    """Represents a KEGG gene entry."""
    
    gene_id: str
    organism: str
    orthology: str
    annotation: str
    sequence: str
    position: str
    pathways: List[str] = field(default_factory=list)
    modules: List[str] = field(default_factory=list)
    ec_numbers: List[str] = field(default_factory=list)
    
    @property
    def full_id(self) -> str:
        """Get full gene ID with organism prefix."""
        return f"{self.organism}:{self.gene_id}"


@dataclass
class KEGGPathway:
    """Represents a KEGG pathway."""
    
    pathway_id: str
    name: str
    classification: str
    genes: List[str] = field(default_factory=list)
    modules: List[str] = field(default_factory=list)
    compounds: List[str] = field(default_factory=list)
    is_metabolic: bool = False


class KEGGAPIClient:
    """Client for KEGG REST API with rate limiting and error handling."""
    
    def __init__(
        self,
        base_url: str = "https://rest.kegg.jp",
        rate_limit_delay: float = 0.35,
        max_retries: int = 3,
        timeout: int = 30,
    ):
        """
        Initialize KEGG API client.
        
        Args:
            base_url: KEGG API base URL
            rate_limit_delay: Delay between requests (seconds)
            max_retries: Maximum number of retries for failed requests
            timeout: Request timeout (seconds)
        """
        self.base_url = base_url
        self.rate_limit_delay = rate_limit_delay
        self.max_retries = max_retries
        self.timeout = timeout
        
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Plant-MGC-Analysis-Pipeline/1.0'
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
    
    def get(self, endpoint: str, **kwargs) -> str:
        """
        Make GET request to KEGG API.
        
        Args:
            endpoint: API endpoint
            **kwargs: Additional request parameters
            
        Returns:
            Response text
            
        Raises:
            NetworkError: If request fails
        """
        url = f"{self.base_url}/{endpoint}"
        
        for attempt in range(self.max_retries):
            try:
                self._enforce_rate_limit()
                
                response = self.session.get(
                    url,
                    timeout=self.timeout,
                    **kwargs
                )
                response.raise_for_status()
                
                return response.text.strip()
                
            except requests.exceptions.RequestException as e:
                if attempt == self.max_retries - 1:
                    raise NetworkError(
                        f"KEGG API request failed after {self.max_retries} attempts: {e}",
                        url=url,
                        status_code=getattr(e.response, 'status_code', None)
                    ) from e
                
                # Exponential backoff
                time.sleep(2 ** attempt)
        
        return ""
    
    def get_entry(self, entry_id: str) -> str:
        """Get KEGG entry by ID."""
        return self.get(f"get/{entry_id}")
    
    def list_entries(self, database: str) -> str:
        """List entries in KEGG database."""
        return self.get(f"list/{database}")
    
    def find_entries(self, database: str, query: str) -> str:
        """Find entries in KEGG database."""
        return self.get(f"find/{database}/{query}")
    
    def link_entries(self, target_db: str, source_db: str) -> str:
        """Get links between databases."""
        return self.get(f"link/{target_db}/{source_db}")
    
    def conv_entries(self, target_db: str, source_db: str) -> str:
        """Convert between database identifiers."""
        return self.get(f"conv/{target_db}/{source_db}")


class KEGGParser:
    """Parser for KEGG flat file format."""
    
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
                parsed[current_key] = line[12:].strip()
            else:
                # Continuation line
                if current_key:
                    parsed[current_key] += " " + line[12:].strip()
        
        return parsed
    
    @staticmethod
    def parse_gene_entry(entry_text: str, organism: str) -> KEGGGene:
        """
        Parse KEGG gene entry.
        
        Args:
            entry_text: Raw gene entry text
            organism: Organism code
            
        Returns:
            KEGGGene object
        """
        parsed = KEGGParser.parse_flat_entry(entry_text)
        
        gene_id = parsed.get("ENTRY", "").split()[0]
        
        # Parse pathways
        pathways = []
        if "PATHWAY" in parsed:
            pathway_text = parsed["PATHWAY"]
            for line in pathway_text.split():
                if line.startswith("ko"):
                    pathways.append(line)
        
        # Parse EC numbers
        ec_numbers = []
        if "ORTHOLOGY" in parsed:
            orthology = parsed["ORTHOLOGY"]
            # Extract EC numbers from orthology field
            import re
            ec_pattern = r'EC:(\d+\.\d+\.\d+\.\d+)'
            ec_numbers = re.findall(ec_pattern, orthology)
        
        return KEGGGene(
            gene_id=gene_id,
            organism=organism,
            orthology=parsed.get("ORTHOLOGY", ""),
            annotation=parsed.get("DEFINITION", ""),
            sequence=parsed.get("AASEQ", ""),
            position=parsed.get("POSITION", ""),
            pathways=pathways,
            ec_numbers=ec_numbers,
        )
    
    @staticmethod
    def parse_module_entry(entry_text: str) -> KEGGModule:
        """
        Parse KEGG module entry.
        
        Args:
            entry_text: Raw module entry text
            
        Returns:
            KEGGModule object
        """
        parsed = KEGGParser.parse_flat_entry(entry_text)
        
        module_id = parsed.get("ENTRY", "").split()[0]
        
        # Check if module is metabolic
        is_metabolic = False
        if "CLASS" in parsed:
            class_text = parsed["CLASS"].lower()
            is_metabolic = "metabol" in class_text
        
        return KEGGModule(
            module_id=module_id,
            name=parsed.get("NAME", ""),
            classification=parsed.get("CLASS", ""),
            definition=parsed.get("DEFINITION", ""),
            is_metabolic=is_metabolic,
        )


class KEGGDatabase(DatabaseInterface):
    """
    KEGG database interface for metabolic pathway analysis.
    
    This class provides comprehensive access to KEGG data including
    genes, pathways, modules, and organism-specific information.
    """
    
    def __init__(self, **kwargs):
        """Initialize KEGG database interface."""
        super().__init__(
            database_name="KEGG",
            connection_params={},
            **kwargs
        )
        
        self.client = KEGGAPIClient()
        self.parser = KEGGParser()
        
        # Caches
        self._metabolic_modules: Optional[Set[str]] = None
        self._organism_cache: Dict[str, Dict[str, Any]] = {}
        self._module_cache: Dict[str, KEGGModule] = {}
        self._pathway_cache: Dict[str, KEGGPathway] = {}
    
    def connect(self) -> None:
        """Establish connection (no-op for REST API)."""
        # Test connection with a simple request
        try:
            self.client.get("info/kegg")
            self.logger.info("Connected to KEGG database")
        except Exception as e:
            raise DatabaseError(
                f"Failed to connect to KEGG database: {e}",
                database_name="KEGG",
                operation="connect"
            ) from e
    
    def disconnect(self) -> None:
        """Close connection (no-op for REST API)."""
        if hasattr(self.client, 'session'):
            self.client.session.close()
    
    def query(self, query: str, parameters: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
        """
        Execute query against KEGG database.
        
        Args:
            query: Query string (KEGG API endpoint)
            parameters: Optional parameters
            
        Returns:
            Query results
        """
        try:
            response = self.client.get(query)
            return [{"response": response}]
        except Exception as e:
            raise DatabaseError(
                f"KEGG query failed: {e}",
                database_name="KEGG",
                operation="query"
            ) from e
    
    def load_metabolic_modules(self, modules_file: Optional[Path] = None) -> Set[str]:
        """
        Load metabolic modules from file or API.
        
        Args:
            modules_file: Optional CSV file with metabolic modules
            
        Returns:
            Set of metabolic module IDs
        """
        if self._metabolic_modules is not None:
            return self._metabolic_modules
        
        metabolic_modules = set()
        
        if modules_file and modules_file.exists():
            # Load from file
            processor = CsvProcessor()
            df = processor.read_file(modules_file)
            
            if "Module_ID" in df.columns:
                metabolic_modules = set(df["Module_ID"].tolist())
            else:
                metabolic_modules = set(df.iloc[:, 0].tolist())
            
            self.logger.info(f"Loaded {len(metabolic_modules)} metabolic modules from file")
        
        else:
            # Load from API
            try:
                modules_text = self.client.list_entries("module")
                
                for line in modules_text.split('\n'):
                    if '\t' in line:
                        module_id, module_name = line.split('\t', 1)
                        
                        # Check if module is metabolic
                        try:
                            module_entry = self.client.get_entry(module_id)
                            module = self.parser.parse_module_entry(module_entry)
                            
                            if module.is_metabolic:
                                metabolic_modules.add(module.general_id)
                        
                        except Exception as e:
                            self.logger.warning(f"Failed to check module {module_id}: {e}")
                
                self.logger.info(f"Loaded {len(metabolic_modules)} metabolic modules from API")
            
            except Exception as e:
                self.logger.error(f"Failed to load metabolic modules: {e}")
                metabolic_modules = set()
        
        self._metabolic_modules = metabolic_modules
        return metabolic_modules
    
    def get_organism_modules(self, organism_code: str) -> Dict[str, List[str]]:
        """
        Get metabolic modules for an organism.
        
        Args:
            organism_code: KEGG organism code
            
        Returns:
            Dictionary mapping module IDs to gene lists
        """
        if organism_code in self._organism_cache:
            return self._organism_cache[organism_code]
        
        try:
            # Get module-gene links for organism
            module_links = self.client.link_entries("module", organism_code)
            
            metabolic_modules = self.load_metabolic_modules()
            organism_modules = {}
            
            for line in module_links.split('\n'):
                if '\t' in line:
                    gene_id, module_id = line.split('\t')
                    
                    # Extract IDs
                    gene = gene_id.split(':')[1]
                    module = module_id.split(':')[1]
                    general_module = module.split('_')[-1] if '_' in module else module
                    
                    # Only include metabolic modules
                    if general_module in metabolic_modules:
                        if module not in organism_modules:
                            organism_modules[module] = []
                        organism_modules[module].append(gene)
            
            self._organism_cache[organism_code] = organism_modules
            self.logger.info(
                f"Found {len(organism_modules)} metabolic modules for {organism_code}"
            )
            
            return organism_modules
            
        except Exception as e:
            raise DatabaseError(
                f"Failed to get modules for organism {organism_code}: {e}",
                database_name="KEGG",
                operation="get_organism_modules"
            ) from e
    
    def get_gene_batch(
        self,
        gene_ids: List[str],
        organism_code: str,
        batch_size: int = 10
    ) -> Dict[str, KEGGGene]:
        """
        Get gene information in batches.
        
        Args:
            gene_ids: List of gene IDs
            organism_code: Organism code
            batch_size: Number of genes per batch
            
        Returns:
            Dictionary mapping gene IDs to KEGGGene objects
        """
        all_genes = {}
        
        # Process genes in batches
        for i in range(0, len(gene_ids), batch_size):
            batch = gene_ids[i:i + batch_size]
            
            try:
                # Build query for batch
                query_ids = [f"{organism_code}:{gene}" for gene in batch]
                query = "+".join(query_ids)
                
                # Get batch data
                response = self.client.get_entry(query)
                
                # Parse each entry
                for entry in response.split("///"):
                    if entry.strip():
                        try:
                            gene = self.parser.parse_gene_entry(entry, organism_code)
                            all_genes[gene.gene_id] = gene
                        except Exception as e:
                            self.logger.warning(f"Failed to parse gene entry: {e}")
                
                # Rate limiting
                time.sleep(self.client.rate_limit_delay)
                
            except Exception as e:
                self.logger.error(f"Failed to fetch batch {batch}: {e}")
                
                # Add placeholder entries for failed genes
                for gene_id in batch:
                    all_genes[gene_id] = KEGGGene(
                        gene_id=gene_id,
                        organism=organism_code,
                        orthology="N/A",
                        annotation="N/A",
                        sequence="N/A",
                        position="N/A",
                    )
        
        return all_genes
    
    def process_organism(
        self,
        organism_code: str,
        output_dir: Path,
        metabolic_modules_file: Optional[Path] = None,
        batch_size: int = 10,
        max_workers: int = 4,
    ) -> Dict[str, int]:
        """
        Process an organism and extract metabolic gene data.
        
        Args:
            organism_code: KEGG organism code
            output_dir: Output directory for results
            metabolic_modules_file: Optional file with metabolic modules
            batch_size: Genes per batch
            max_workers: Maximum worker threads
            
        Returns:
            Dictionary with processing statistics
        """
        self.logger.info(f"Processing organism: {organism_code}")
        
        # Create output directory
        org_output_dir = output_dir / organism_code
        org_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if already processed
        if any(org_output_dir.glob("*.csv")):
            self.logger.info(f"Organism {organism_code} already processed, skipping")
            return {"status": "skipped"}
        
        # Load metabolic modules
        if metabolic_modules_file:
            self.load_metabolic_modules(metabolic_modules_file)
        
        # Get organism modules
        organism_modules = self.get_organism_modules(organism_code)
        
        if not organism_modules:
            self.logger.warning(f"No metabolic modules found for {organism_code}")
            return {"status": "no_modules", "modules": 0}
        
        # Process each module
        total_genes = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            
            for module_id, gene_ids in organism_modules.items():
                future = executor.submit(
                    self._process_module,
                    module_id,
                    gene_ids,
                    organism_code,
                    org_output_dir,
                    batch_size
                )
                futures[future] = module_id
            
            for future in as_completed(futures):
                module_id = futures[future]
                try:
                    gene_count = future.result()
                    total_genes += gene_count
                    self.logger.debug(f"Processed module {module_id}: {gene_count} genes")
                except Exception as e:
                    self.logger.error(f"Failed to process module {module_id}: {e}")
        
        self.logger.info(
            f"Completed organism {organism_code}: "
            f"{len(organism_modules)} modules, {total_genes} genes"
        )
        
        return {
            "status": "completed",
            "modules": len(organism_modules),
            "genes": total_genes,
        }
    
    def _process_module(
        self,
        module_id: str,
        gene_ids: List[str],
        organism_code: str,
        output_dir: Path,
        batch_size: int,
    ) -> int:
        """Process a single module."""
        module_file = output_dir / f"{module_id}.csv"
        
        # Get gene information
        genes = self.get_gene_batch(gene_ids, organism_code, batch_size)
        
        # Write to CSV
        rows = []
        for gene_id, gene in genes.items():
            rows.append({
                "Gene_ID": gene.gene_id,
                "Annotation": gene.annotation,
                "Orthology": gene.orthology,
                "Sequence": gene.sequence,
                "Position": gene.position,
                "EC_Numbers": ";".join(gene.ec_numbers),
                "Pathways": ";".join(gene.pathways),
            })
        
        if rows:
            df = pd.DataFrame(rows)
            processor = CsvProcessor()
            processor.write_file(df, module_file)
        
        return len(rows)
    
    def validate_input(self, data: str) -> None:
        """Validate KEGG query input."""
        if not data or not isinstance(data, str):
            raise ValidationError(
                "KEGG query must be a non-empty string",
                field_name="query",
                field_value=data
            )
    
    def process(self, data: str, **kwargs) -> Dict[str, Any]:
        """Process KEGG query."""
        self.validate_input(data)
        
        try:
            response = self.client.get(data)
            return {"query": data, "response": response}
        except Exception as e:
            raise DatabaseError(
                f"Failed to process KEGG query: {e}",
                database_name="KEGG",
                operation="process"
            ) from e


class KEGGAnalyzer(BioinformaticsProcessor):
    """
    High-level analyzer for KEGG-based metabolic analysis.
    
    This class provides convenient methods for analyzing metabolic
    pathways and gene clusters using KEGG data.
    """
    
    def __init__(self, **kwargs):
        """Initialize KEGG analyzer."""
        super().__init__(**kwargs)
        self.database = KEGGDatabase(**kwargs)
        self._metabolic_annotations: Dict[str, Dict[str, str]] = {}
    
    def validate_input(self, data: Any) -> None:
        """Validate input data."""
        if not isinstance(data, (str, list, dict)):
            raise ValidationError(
                "Input must be organism code, gene list, or configuration dict",
                field_name="input_data",
                field_value=type(data).__name__
            )
    
    def load_metabolic_annotations(self, organisms: List[str], output_dir: Path) -> None:
        """
        Load metabolic annotations for multiple organisms.
        
        Args:
            organisms: List of organism codes
            output_dir: Output directory for organism data
        """
        self.logger.info(f"Loading metabolic annotations for {len(organisms)} organisms")
        
        for organism in organisms:
            try:
                stats = self.database.process_organism(organism, output_dir)
                self.logger.info(f"Processed {organism}: {stats}")
            except Exception as e:
                self.logger.error(f"Failed to process organism {organism}: {e}")
    
    def annotate_genes(self, gene_ids: List[str], organism_code: str) -> Dict[str, Dict[str, str]]:
        """
        Annotate genes with KEGG information.
        
        Args:
            gene_ids: List of gene IDs
            organism_code: Organism code
            
        Returns:
            Dictionary mapping gene IDs to annotation data
        """
        annotations = {}
        
        try:
            genes = self.database.get_gene_batch(gene_ids, organism_code)
            
            for gene_id, gene in genes.items():
                annotations[gene_id] = {
                    "annotation": gene.annotation,
                    "orthology": gene.orthology,
                    "ec_numbers": ";".join(gene.ec_numbers),
                    "pathways": ";".join(gene.pathways),
                    "is_metabolic": len(gene.ec_numbers) > 0 or len(gene.pathways) > 0,
                }
        
        except Exception as e:
            self.logger.error(f"Failed to annotate genes: {e}")
        
        return annotations
    
    def process(self, data: Any, **kwargs) -> Dict[str, Any]:
        """Process KEGG analysis request."""
        self.validate_input(data)
        
        if isinstance(data, str):
            # Single organism processing
            output_dir = kwargs.get("output_dir", self.temp_dir)
            return self.database.process_organism(data, Path(output_dir))
        
        elif isinstance(data, list):
            # Gene annotation
            organism = kwargs.get("organism", "")
            if not organism:
                raise ValidationError("Organism code required for gene annotation")
            
            return self.annotate_genes(data, organism)
        
        elif isinstance(data, dict):
            # Complex analysis configuration
            return self._process_complex_analysis(data, **kwargs)
        
        return {}