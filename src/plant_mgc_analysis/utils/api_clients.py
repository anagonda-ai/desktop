"""
Unified API clients for Plant MGC Analysis Pipeline.

This module consolidates all external API operations found across legacy scripts,
including KEGG, NCBI, Ensembl, and other bioinformatics databases.
"""

import time
import json
import urllib.parse
from typing import Any, Dict, List, Optional, Iterator, Union
from pathlib import Path
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from loguru import logger

from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..core.types import PathLike


@dataclass
class APIResponse:
    """Standardized API response wrapper."""
    
    status_code: int
    data: Any
    headers: Dict[str, str]
    url: str
    processing_time: float
    error_message: Optional[str] = None
    
    @property
    def success(self) -> bool:
        """Check if request was successful."""
        return 200 <= self.status_code < 300
    
    @property
    def is_rate_limited(self) -> bool:
        """Check if request was rate limited."""
        return self.status_code == 429


@dataclass
class APIClientConfig:
    """Configuration for API clients."""
    
    base_url: str
    timeout: int = 30
    max_retries: int = 3
    retry_delay: float = 1.0
    rate_limit_delay: float = 1.0
    max_rate_limit_retries: int = 5
    headers: Dict[str, str] = field(default_factory=dict)
    auth: Optional[Dict[str, str]] = None
    verify_ssl: bool = True
    
    def __post_init__(self):
        """Validate configuration."""
        if not self.base_url:
            raise ValidationError("Base URL is required")
        if self.timeout <= 0:
            raise ValidationError("Timeout must be positive")


class BaseAPIClient(BioinformaticsProcessor):
    """
    Base class for all API clients.
    
    Provides common functionality for HTTP requests, rate limiting,
    error handling, and response processing.
    """
    
    def __init__(self, config: APIClientConfig, **kwargs):
        """Initialize API client."""
        super().__init__(**kwargs)
        self.config = config
        self.session = self._create_session()
        self.request_count = 0
        self.last_request_time = 0.0
        
        # Statistics
        self.stats = {
            "total_requests": 0,
            "successful_requests": 0,
            "failed_requests": 0,
            "rate_limited_requests": 0,
            "total_processing_time": 0.0,
        }
    
    def _create_session(self) -> requests.Session:
        """Create configured requests session."""
        session = requests.Session()
        
        # Configure retries
        retry_strategy = Retry(
            total=self.config.max_retries,
            backoff_factor=self.config.retry_delay,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        
        # Set headers
        session.headers.update(self.config.headers)
        
        # Set auth if provided
        if self.config.auth:
            session.auth = (self.config.auth["username"], self.config.auth["password"])
        
        return session
    
    def _rate_limit_check(self) -> None:
        """Check and enforce rate limiting."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < self.config.rate_limit_delay:
            sleep_time = self.config.rate_limit_delay - time_since_last
            self.logger.debug(f"Rate limiting: sleeping for {sleep_time:.2f} seconds")
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
    
    def _make_request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
    ) -> APIResponse:
        """Make HTTP request with rate limiting and error handling."""
        self._rate_limit_check()
        
        url = urllib.parse.urljoin(self.config.base_url, endpoint)
        start_time = time.time()
        
        try:
            response = self.session.request(
                method=method,
                url=url,
                params=params,
                json=data,
                headers=headers,
                timeout=self.config.timeout,
                verify=self.config.verify_ssl,
            )
            
            processing_time = time.time() - start_time
            
            # Handle rate limiting
            if response.status_code == 429:
                self.stats["rate_limited_requests"] += 1
                self.logger.warning(f"Rate limited for endpoint: {endpoint}")
                
                # Extract retry-after header if available
                retry_after = response.headers.get('Retry-After')
                if retry_after:
                    sleep_time = int(retry_after)
                    self.logger.info(f"Sleeping for {sleep_time} seconds (Retry-After)")
                    time.sleep(sleep_time)
                    
                    # Retry the request
                    return self._make_request(method, endpoint, params, data, headers)
            
            # Parse response data
            try:
                response_data = response.json()
            except ValueError:
                response_data = response.text
            
            # Create API response
            api_response = APIResponse(
                status_code=response.status_code,
                data=response_data,
                headers=dict(response.headers),
                url=url,
                processing_time=processing_time,
                error_message=None if response.ok else response.text,
            )
            
            # Update statistics
            self.stats["total_requests"] += 1
            self.stats["total_processing_time"] += processing_time
            
            if response.ok:
                self.stats["successful_requests"] += 1
            else:
                self.stats["failed_requests"] += 1
                self.logger.error(f"API request failed: {response.status_code} - {response.text}")
            
            return api_response
            
        except requests.exceptions.Timeout:
            self.stats["failed_requests"] += 1
            raise ProcessingError(f"Request timeout for {url}")
        except requests.exceptions.ConnectionError:
            self.stats["failed_requests"] += 1
            raise ProcessingError(f"Connection error for {url}")
        except Exception as e:
            self.stats["failed_requests"] += 1
            raise ProcessingError(f"Request failed: {str(e)}")
    
    def get(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> APIResponse:
        """Make GET request."""
        return self._make_request("GET", endpoint, params=params)
    
    def post(self, endpoint: str, data: Optional[Dict[str, Any]] = None) -> APIResponse:
        """Make POST request."""
        return self._make_request("POST", endpoint, data=data)
    
    def get_stats(self) -> Dict[str, Any]:
        """Get API client statistics."""
        return self.stats.copy()
    
    def validate_input(self, data: Any) -> None:
        """Validate API input."""
        pass  # Override in subclasses
    
    def process(self, data: Any, **kwargs) -> Any:
        """Process API request."""
        pass  # Override in subclasses


class KEGGAPIClient(BaseAPIClient):
    """
    KEGG API client consolidating operations from 8+ legacy scripts.
    
    Replaces scattered KEGG operations found in:
    - process_organism.py
    - kegg_pathway_extractor_*.py
    - Various metabolic analysis scripts
    """
    
    def __init__(self, **kwargs):
        """Initialize KEGG API client."""
        config = APIClientConfig(
            base_url="https://rest.kegg.jp/",
            rate_limit_delay=1.0,  # KEGG requires 1 second between requests
            headers={"User-Agent": "Plant-MGC-Analysis/1.0"},
        )
        super().__init__(config, **kwargs)
    
    def get_organism_info(self, organism_code: str) -> Dict[str, Any]:
        """Get organism information from KEGG."""
        response = self.get(f"get/{organism_code}")
        
        if not response.success:
            raise ProcessingError(f"Failed to get organism info for {organism_code}")
        
        return self._parse_organism_info(response.data)
    
    def get_pathway_info(self, pathway_id: str) -> Dict[str, Any]:
        """Get pathway information from KEGG."""
        response = self.get(f"get/{pathway_id}")
        
        if not response.success:
            raise ProcessingError(f"Failed to get pathway info for {pathway_id}")
        
        return self._parse_pathway_info(response.data)
    
    def list_pathways(self, organism_code: str) -> List[Dict[str, str]]:
        """List all pathways for an organism."""
        response = self.get(f"list/pathway/{organism_code}")
        
        if not response.success:
            raise ProcessingError(f"Failed to list pathways for {organism_code}")
        
        return self._parse_pathway_list(response.data)
    
    def get_gene_info(self, gene_id: str) -> Dict[str, Any]:
        """Get gene information from KEGG."""
        response = self.get(f"get/{gene_id}")
        
        if not response.success:
            raise ProcessingError(f"Failed to get gene info for {gene_id}")
        
        return self._parse_gene_info(response.data)
    
    def find_genes_by_pathway(self, pathway_id: str) -> List[str]:
        """Find genes associated with a pathway."""
        response = self.get(f"link/gene/{pathway_id}")
        
        if not response.success:
            return []  # No genes found
        
        return self._parse_gene_list(response.data)
    
    def get_compound_info(self, compound_id: str) -> Dict[str, Any]:
        """Get compound information from KEGG."""
        response = self.get(f"get/{compound_id}")
        
        if not response.success:
            raise ProcessingError(f"Failed to get compound info for {compound_id}")
        
        return self._parse_compound_info(response.data)
    
    def search_compounds(self, query: str) -> List[Dict[str, str]]:
        """Search compounds by name or formula."""
        response = self.get(f"find/compound/{query}")
        
        if not response.success:
            return []  # No compounds found
        
        return self._parse_compound_list(response.data)
    
    def get_enzyme_info(self, enzyme_id: str) -> Dict[str, Any]:
        """Get enzyme information from KEGG."""
        response = self.get(f"get/{enzyme_id}")
        
        if not response.success:
            raise ProcessingError(f"Failed to get enzyme info for {enzyme_id}")
        
        return self._parse_enzyme_info(response.data)
    
    def batch_get_pathways(self, organism_codes: List[str]) -> Dict[str, List[Dict[str, str]]]:
        """Get pathways for multiple organisms."""
        results = {}
        
        for organism_code in organism_codes:
            try:
                pathways = self.list_pathways(organism_code)
                results[organism_code] = pathways
                self.logger.info(f"Retrieved {len(pathways)} pathways for {organism_code}")
            except Exception as e:
                self.logger.error(f"Failed to get pathways for {organism_code}: {e}")
                results[organism_code] = []
        
        return results
    
    def _parse_organism_info(self, data: str) -> Dict[str, Any]:
        """Parse organism information from KEGG response."""
        lines = data.split('\n')
        info = {}
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line.split('NAME')[1].strip()
            elif line.startswith('DEFINITION'):
                info['definition'] = line.split('DEFINITION')[1].strip()
            elif line.startswith('TAXONOMY'):
                info['taxonomy'] = line.split('TAXONOMY')[1].strip()
        
        return info
    
    def _parse_pathway_info(self, data: str) -> Dict[str, Any]:
        """Parse pathway information from KEGG response."""
        lines = data.split('\n')
        info = {}
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line.split('NAME')[1].strip()
            elif line.startswith('DESCRIPTION'):
                info['description'] = line.split('DESCRIPTION')[1].strip()
            elif line.startswith('CLASS'):
                info['class'] = line.split('CLASS')[1].strip()
        
        return info
    
    def _parse_pathway_list(self, data: str) -> List[Dict[str, str]]:
        """Parse pathway list from KEGG response."""
        pathways = []
        
        for line in data.split('\n'):
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    pathways.append({
                        'id': parts[0],
                        'name': parts[1]
                    })
        
        return pathways
    
    def _parse_gene_info(self, data: str) -> Dict[str, Any]:
        """Parse gene information from KEGG response."""
        lines = data.split('\n')
        info = {}
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line.split('NAME')[1].strip()
            elif line.startswith('DEFINITION'):
                info['definition'] = line.split('DEFINITION')[1].strip()
            elif line.startswith('PATHWAY'):
                info['pathway'] = line.split('PATHWAY')[1].strip()
        
        return info
    
    def _parse_gene_list(self, data: str) -> List[str]:
        """Parse gene list from KEGG response."""
        genes = []
        
        for line in data.split('\n'):
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 1:
                    genes.append(parts[0])
        
        return genes
    
    def _parse_compound_info(self, data: str) -> Dict[str, Any]:
        """Parse compound information from KEGG response."""
        lines = data.split('\n')
        info = {}
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line.split('NAME')[1].strip()
            elif line.startswith('FORMULA'):
                info['formula'] = line.split('FORMULA')[1].strip()
            elif line.startswith('EXACT_MASS'):
                info['exact_mass'] = line.split('EXACT_MASS')[1].strip()
        
        return info
    
    def _parse_compound_list(self, data: str) -> List[Dict[str, str]]:
        """Parse compound list from KEGG response."""
        compounds = []
        
        for line in data.split('\n'):
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    compounds.append({
                        'id': parts[0],
                        'name': parts[1]
                    })
        
        return compounds
    
    def _parse_enzyme_info(self, data: str) -> Dict[str, Any]:
        """Parse enzyme information from KEGG response."""
        lines = data.split('\n')
        info = {}
        
        for line in lines:
            if line.startswith('NAME'):
                info['name'] = line.split('NAME')[1].strip()
            elif line.startswith('CLASS'):
                info['class'] = line.split('CLASS')[1].strip()
            elif line.startswith('REACTION'):
                info['reaction'] = line.split('REACTION')[1].strip()
        
        return info


class NCBIAPIClient(BaseAPIClient):
    """
    NCBI API client for sequence and genome data.
    
    Consolidates NCBI operations found across multiple legacy scripts.
    """
    
    def __init__(self, email: str, api_key: Optional[str] = None, **kwargs):
        """Initialize NCBI API client."""
        config = APIClientConfig(
            base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
            rate_limit_delay=0.34,  # NCBI allows 3 requests per second
            headers={"User-Agent": "Plant-MGC-Analysis/1.0"},
        )
        super().__init__(config, **kwargs)
        self.email = email
        self.api_key = api_key
    
    def search_sequences(self, query: str, database: str = "protein") -> List[str]:
        """Search for sequences in NCBI database."""
        params = {
            "db": database,
            "term": query,
            "retmode": "json",
            "email": self.email,
        }
        
        if self.api_key:
            params["api_key"] = self.api_key
        
        response = self.get("esearch.fcgi", params=params)
        
        if not response.success:
            raise ProcessingError(f"Failed to search sequences: {query}")
        
        return response.data.get("esearchresult", {}).get("idlist", [])
    
    def fetch_sequences(self, ids: List[str], database: str = "protein") -> List[Dict[str, Any]]:
        """Fetch sequences by IDs."""
        params = {
            "db": database,
            "id": ",".join(ids),
            "retmode": "json",
            "email": self.email,
        }
        
        if self.api_key:
            params["api_key"] = self.api_key
        
        response = self.get("efetch.fcgi", params=params)
        
        if not response.success:
            raise ProcessingError(f"Failed to fetch sequences: {ids}")
        
        return self._parse_sequence_data(response.data)
    
    def _parse_sequence_data(self, data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse sequence data from NCBI response."""
        # Implementation depends on NCBI response format
        return []


class EnsemblAPIClient(BaseAPIClient):
    """
    Ensembl API client for genome annotations.
    
    Consolidates Ensembl operations found across legacy scripts.
    """
    
    def __init__(self, **kwargs):
        """Initialize Ensembl API client."""
        config = APIClientConfig(
            base_url="https://rest.ensembl.org/",
            headers={
                "Content-Type": "application/json",
                "User-Agent": "Plant-MGC-Analysis/1.0",
            },
        )
        super().__init__(config, **kwargs)
    
    def get_gene_info(self, gene_id: str, species: str = "arabidopsis_thaliana") -> Dict[str, Any]:
        """Get gene information from Ensembl."""
        response = self.get(f"lookup/id/{gene_id}", params={"species": species})
        
        if not response.success:
            raise ProcessingError(f"Failed to get gene info for {gene_id}")
        
        return response.data
    
    def get_sequence(self, gene_id: str, seq_type: str = "protein") -> str:
        """Get sequence for a gene."""
        response = self.get(f"sequence/id/{gene_id}", params={"type": seq_type})
        
        if not response.success:
            raise ProcessingError(f"Failed to get sequence for {gene_id}")
        
        return response.data.get("seq", "")


class UnifiedAPIManager:
    """
    Unified API manager that orchestrates all external API clients.
    
    This class provides a single interface for all API operations,
    automatically routing requests to appropriate clients.
    """
    
    def __init__(self, email: str, ncbi_api_key: Optional[str] = None):
        """Initialize unified API manager."""
        self.kegg_client = KEGGAPIClient()
        self.ncbi_client = NCBIAPIClient(email=email, api_key=ncbi_api_key)
        self.ensembl_client = EnsemblAPIClient()
        
        self.logger = logger
    
    def get_kegg_pathways(self, organism_code: str) -> List[Dict[str, str]]:
        """Get KEGG pathways for organism."""
        return self.kegg_client.list_pathways(organism_code)
    
    def get_ensembl_gene_info(self, gene_id: str, species: str = "arabidopsis_thaliana") -> Dict[str, Any]:
        """Get Ensembl gene information."""
        return self.ensembl_client.get_gene_info(gene_id, species)
    
    def search_ncbi_sequences(self, query: str, database: str = "protein") -> List[str]:
        """Search NCBI sequences."""
        return self.ncbi_client.search_sequences(query, database)
    
    def batch_process_organisms(self, organism_codes: List[str]) -> Dict[str, Dict[str, Any]]:
        """Process multiple organisms across all APIs."""
        results = {}
        
        for organism_code in organism_codes:
            try:
                organism_results = {
                    "kegg_pathways": self.kegg_client.list_pathways(organism_code),
                    "kegg_info": self.kegg_client.get_organism_info(organism_code),
                }
                results[organism_code] = organism_results
                self.logger.info(f"Processed organism: {organism_code}")
            except Exception as e:
                self.logger.error(f"Failed to process organism {organism_code}: {e}")
                results[organism_code] = {"error": str(e)}
        
        return results
    
    def get_combined_stats(self) -> Dict[str, Any]:
        """Get combined statistics from all API clients."""
        return {
            "kegg": self.kegg_client.get_stats(),
            "ncbi": self.ncbi_client.get_stats(),
            "ensembl": self.ensembl_client.get_stats(),
        }


# Global API manager instance
def create_api_manager(email: str, ncbi_api_key: Optional[str] = None) -> UnifiedAPIManager:
    """Create unified API manager instance."""
    return UnifiedAPIManager(email=email, ncbi_api_key=ncbi_api_key)