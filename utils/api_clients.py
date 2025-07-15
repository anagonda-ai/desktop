"""
Unified API clients for bioinformatics databases.

This module provides standardized, rate-limited, and cached access to
biological databases including KEGG, NCBI, and UniProt.
"""

import time
import requests
import json
from typing import Any, Dict, List, Optional, Union
from pathlib import Path
from dataclasses import dataclass
from enum import Enum
import logging
from Bio import Entrez
from urllib.parse import urljoin, quote
import asyncio
import aiohttp
from datetime import datetime, timedelta

from ..core.base import BatchAPIProcessor
from ..core.exceptions import APIError, ValidationError, ErrorCode
from ..core.types import DatabaseType, FilePath
from ..core.config import get_config


@dataclass
class APIResponse:
    """Standardized API response."""
    data: Any
    status_code: int
    headers: Dict[str, str]
    cached: bool = False
    timestamp: datetime = None
    
    def __post_init__(self):
        if self.timestamp is None:
            self.timestamp = datetime.now()


class CacheManager:
    """Simple file-based cache for API responses."""
    
    def __init__(self, cache_dir: Optional[Path] = None, ttl_hours: int = 24):
        self.cache_dir = cache_dir or Path("./cache/api")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.ttl = timedelta(hours=ttl_hours)
        self.logger = logging.getLogger(__name__)
    
    def _get_cache_path(self, key: str) -> Path:
        """Get cache file path for key."""
        # Simple hash to create filename
        import hashlib
        hash_key = hashlib.md5(key.encode()).hexdigest()
        return self.cache_dir / f"{hash_key}.json"
    
    def get(self, key: str) -> Optional[Any]:
        """Get cached data if valid."""
        cache_path = self._get_cache_path(key)
        
        if not cache_path.exists():
            return None
        
        try:
            with open(cache_path, 'r') as f:
                cache_data = json.load(f)
            
            # Check if cache is still valid
            cache_time = datetime.fromisoformat(cache_data['timestamp'])
            if datetime.now() - cache_time > self.ttl:
                cache_path.unlink()  # Remove expired cache
                return None
            
            self.logger.debug(f"Cache hit for key: {key[:50]}...")
            return cache_data['data']
            
        except Exception as e:
            self.logger.warning(f"Cache read error: {e}")
            return None
    
    def set(self, key: str, data: Any):
        """Cache data with timestamp."""
        cache_path = self._get_cache_path(key)
        
        try:
            cache_data = {
                'data': data,
                'timestamp': datetime.now().isoformat(),
                'key': key
            }
            
            with open(cache_path, 'w') as f:
                json.dump(cache_data, f)
            
            self.logger.debug(f"Cached data for key: {key[:50]}...")
            
        except Exception as e:
            self.logger.warning(f"Cache write error: {e}")


class KEGGClient(BatchAPIProcessor):
    """
    Professional KEGG API client with rate limiting and caching.
    
    Provides access to KEGG REST API with proper error handling,
    retry logic, and response caching.
    """
    
    def __init__(self, config=None):
        config = config or get_config()
        super().__init__(
            batch_size=10,
            rate_limit=config.database.kegg_rate_limit,
            max_workers=4
        )
        self.base_url = config.database.kegg_api_url
        self.cache = CacheManager(
            cache_dir=config.paths.cache_dir / "kegg",
            ttl_hours=config.compute.cache_ttl_hours
        )
        self.session = requests.Session()
        self.analysis_type = None  # For compatibility with base class
    
    def validate_input(self, data: Any) -> None:
        """Validate API request parameters."""
        if isinstance(data, str) and not data.strip():
            raise ValidationError("Empty query string provided")
    
    def process(self, query: str, **kwargs) -> APIResponse:
        """Process single API query."""
        return self.get(query, **kwargs)
    
    def process_batch(self, batch: List[str], **kwargs) -> List[APIResponse]:
        """Process batch of API queries."""
        results = []
        for query in batch:
            try:
                response = self.get(query, **kwargs)
                results.append(response)
            except Exception as e:
                self.logger.error(f"Failed to process query {query}: {e}")
                # Add error response
                results.append(APIResponse(
                    data=None,
                    status_code=500,
                    headers={},
                    cached=False
                ))
        return results
    
    def get(self, endpoint: str, use_cache: bool = True, **params) -> APIResponse:
        """
        Make GET request to KEGG API.
        
        Args:
            endpoint: API endpoint (e.g., 'list/organism', 'get/hsa:10458')
            use_cache: Whether to use cached responses
            **params: Additional query parameters
            
        Returns:
            APIResponse object
        """
        # Construct URL
        url = urljoin(self.base_url + "/", endpoint)
        if params:
            param_str = "&".join(f"{k}={quote(str(v))}" for k, v in params.items())
            url += f"?{param_str}"
        
        # Check cache first
        cache_key = f"kegg:{url}"
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data is not None:
                return APIResponse(
                    data=cached_data,
                    status_code=200,
                    headers={},
                    cached=True
                )
        
        # Apply rate limiting
        self._rate_limit_wait()
        
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            
            # Parse response
            content_type = response.headers.get('content-type', '')
            if 'json' in content_type:
                data = response.json()
            else:
                data = response.text
            
            # Cache successful response
            if use_cache and response.status_code == 200:
                self.cache.set(cache_key, data)
            
            return APIResponse(
                data=data,
                status_code=response.status_code,
                headers=dict(response.headers),
                cached=False
            )
            
        except requests.exceptions.RequestException as e:
            raise APIError(
                f"KEGG API request failed: {e}",
                api_name="KEGG",
                endpoint=endpoint,
                status_code=getattr(e.response, 'status_code', None),
                error_code=ErrorCode.API_CONNECTION_FAILED,
                original_error=e
            )
    
    def list_organisms(self) -> List[Dict[str, str]]:
        """Get list of all organisms in KEGG."""
        response = self.get("list/organism")
        
        organisms = []
        for line in response.data.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 3:
                organisms.append({
                    'code': parts[1],
                    'name': parts[2],
                    'full_name': parts[2] if len(parts) == 3 else parts[3]
                })
        
        return organisms
    
    def get_organism_info(self, organism_code: str) -> Dict[str, Any]:
        """Get detailed information about an organism."""
        response = self.get(f"get/{organism_code}")
        
        # Parse KEGG format
        info = {}
        current_key = None
        
        for line in response.data.split('\n'):
            if line.startswith('ENTRY'):
                info['entry'] = line.split()[1]
            elif line.startswith('NAME'):
                info['name'] = line[12:].strip()
            elif line.startswith('DEFINITION'):
                info['definition'] = line[12:].strip()
            elif line.startswith('TAXONOMY'):
                info['taxonomy'] = line[12:].strip()
        
        return info
    
    def get_pathway_info(self, pathway_id: str) -> Dict[str, Any]:
        """Get pathway information."""
        response = self.get(f"get/{pathway_id}")
        
        # Parse pathway data
        pathway_info = {
            'id': pathway_id,
            'genes': [],
            'compounds': [],
            'enzymes': []
        }
        
        lines = response.data.split('\n')
        current_section = None
        
        for line in lines:
            if line.startswith('NAME'):
                pathway_info['name'] = line[12:].strip()
            elif line.startswith('DESCRIPTION'):
                pathway_info['description'] = line[12:].strip()
            elif line.startswith('CLASS'):
                pathway_info['class'] = line[12:].strip()
            elif line.startswith('GENE'):
                current_section = 'genes'
                gene_line = line[12:].strip()
                if gene_line:
                    pathway_info['genes'].append(gene_line)
            elif line.startswith('COMPOUND'):
                current_section = 'compounds'
                compound_line = line[12:].strip()
                if compound_line:
                    pathway_info['compounds'].append(compound_line)
            elif line.startswith(' ') and current_section:
                # Continuation line
                continuation = line.strip()
                if continuation and current_section == 'genes':
                    pathway_info['genes'].append(continuation)
                elif continuation and current_section == 'compounds':
                    pathway_info['compounds'].append(continuation)
        
        return pathway_info
    
    def find_genes_by_function(self, organism_code: str, function_keywords: List[str]) -> List[Dict[str, str]]:
        """Find genes by functional annotation."""
        response = self.get(f"list/{organism_code}")
        
        genes = []
        for line in response.data.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 2:
                gene_id = parts[0]
                gene_desc = parts[1] if len(parts) > 1 else ""
                
                # Check if any keyword matches
                desc_lower = gene_desc.lower()
                if any(keyword.lower() in desc_lower for keyword in function_keywords):
                    genes.append({
                        'gene_id': gene_id,
                        'description': gene_desc
                    })
        
        return genes


class NCBIClient(BatchAPIProcessor):
    """
    Professional NCBI API client using Entrez utilities.
    
    Provides access to NCBI databases with proper authentication,
    rate limiting, and error handling.
    """
    
    def __init__(self, config=None):
        config = config or get_config()
        super().__init__(
            batch_size=100,  # NCBI allows larger batches
            rate_limit=config.database.ncbi_rate_limit,
            max_workers=2  # Be conservative with NCBI
        )
        
        # Configure Entrez
        Entrez.email = config.database.ncbi_email
        if config.database.ncbi_api_key:
            Entrez.api_key = config.database.ncbi_api_key
        
        self.cache = CacheManager(
            cache_dir=config.paths.cache_dir / "ncbi",
            ttl_hours=config.compute.cache_ttl_hours
        )
        self.analysis_type = None  # For compatibility
    
    def validate_input(self, data: Any) -> None:
        """Validate NCBI query parameters."""
        if not Entrez.email:
            raise ValidationError(
                "NCBI email is required. Set BIOINF_NCBI_EMAIL environment variable."
            )
    
    def process(self, query: str, database: str = "nucleotide", **kwargs) -> APIResponse:
        """Process single NCBI query."""
        return self.search(query, database, **kwargs)
    
    def process_batch(self, batch: List[str], database: str = "nucleotide", **kwargs) -> List[APIResponse]:
        """Process batch of NCBI queries."""
        # For NCBI, it's more efficient to join queries
        combined_query = " OR ".join(f"({query})" for query in batch)
        response = self.search(combined_query, database, **kwargs)
        
        # Return single response for the batch
        return [response]
    
    def search(self, query: str, database: str = "nucleotide", retmax: int = 100, use_cache: bool = True) -> APIResponse:
        """
        Search NCBI database.
        
        Args:
            query: Search query
            database: NCBI database name
            retmax: Maximum results to return
            use_cache: Whether to use cached responses
            
        Returns:
            APIResponse with search results
        """
        cache_key = f"ncbi:{database}:{query}:{retmax}"
        
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data is not None:
                return APIResponse(
                    data=cached_data,
                    status_code=200,
                    headers={},
                    cached=True
                )
        
        self._rate_limit_wait()
        
        try:
            # Search for IDs
            search_handle = Entrez.esearch(
                db=database,
                term=query,
                retmax=retmax,
                sort="relevance"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results["IdList"]
            
            if not id_list:
                data = {"records": [], "count": 0}
            else:
                # Fetch details
                fetch_handle = Entrez.efetch(
                    db=database,
                    id=id_list,
                    rettype="gb",
                    retmode="text"
                )
                records_text = fetch_handle.read()
                fetch_handle.close()
                
                data = {
                    "records": records_text,
                    "ids": id_list,
                    "count": len(id_list),
                    "query": query,
                    "database": database
                }
            
            if use_cache:
                self.cache.set(cache_key, data)
            
            return APIResponse(
                data=data,
                status_code=200,
                headers={},
                cached=False
            )
            
        except Exception as e:
            raise APIError(
                f"NCBI search failed: {e}",
                api_name="NCBI",
                endpoint=f"{database}/search",
                error_code=ErrorCode.API_CONNECTION_FAILED,
                original_error=e
            )
    
    def get_gene_coordinates(self, gene_ids: List[str], use_cache: bool = True) -> List[Dict[str, Any]]:
        """Get genomic coordinates for genes."""
        cache_key = f"ncbi:coordinates:{':'.join(gene_ids)}"
        
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data is not None:
                return cached_data
        
        coordinates = []
        
        # Process in batches
        for i in range(0, len(gene_ids), 100):
            batch = gene_ids[i:i+100]
            self._rate_limit_wait()
            
            try:
                # Get gene information
                handle = Entrez.efetch(
                    db="gene",
                    id=batch,
                    rettype="gene_table",
                    retmode="text"
                )
                result = handle.read()
                handle.close()
                
                # Parse coordinates (simplified)
                for line in result.split('\n'):
                    if line.strip() and '\t' in line:
                        parts = line.split('\t')
                        if len(parts) >= 5:
                            coordinates.append({
                                'gene_id': parts[0],
                                'chromosome': parts[1],
                                'start': int(parts[2]) if parts[2].isdigit() else None,
                                'end': int(parts[3]) if parts[3].isdigit() else None,
                                'strand': parts[4]
                            })
                
            except Exception as e:
                self.logger.error(f"Failed to get coordinates for batch: {e}")
        
        if use_cache:
            self.cache.set(cache_key, coordinates)
        
        return coordinates


class UniProtClient(BatchAPIProcessor):
    """
    Professional UniProt API client.
    
    Provides access to UniProt REST API with proper rate limiting
    and response caching.
    """
    
    def __init__(self, config=None):
        config = config or get_config()
        super().__init__(
            batch_size=25,
            rate_limit=config.database.uniprot_rate_limit,
            max_workers=3
        )
        self.base_url = config.database.uniprot_api_url
        self.cache = CacheManager(
            cache_dir=config.paths.cache_dir / "uniprot",
            ttl_hours=config.compute.cache_ttl_hours
        )
        self.session = requests.Session()
        self.analysis_type = None  # For compatibility
    
    def validate_input(self, data: Any) -> None:
        """Validate UniProt query parameters."""
        if isinstance(data, str) and not data.strip():
            raise ValidationError("Empty query string provided")
    
    def process(self, query: str, **kwargs) -> APIResponse:
        """Process single UniProt query."""
        return self.search(query, **kwargs)
    
    def process_batch(self, batch: List[str], **kwargs) -> List[APIResponse]:
        """Process batch of UniProt queries."""
        results = []
        for query in batch:
            try:
                response = self.search(query, **kwargs)
                results.append(response)
            except Exception as e:
                self.logger.error(f"Failed to process query {query}: {e}")
                results.append(APIResponse(
                    data=None,
                    status_code=500,
                    headers={},
                    cached=False
                ))
        return results
    
    def search(self, query: str, format: str = "json", size: int = 100, use_cache: bool = True) -> APIResponse:
        """
        Search UniProt database.
        
        Args:
            query: Search query
            format: Response format (json, tsv, fasta)
            size: Maximum results to return
            use_cache: Whether to use cached responses
            
        Returns:
            APIResponse with search results
        """
        cache_key = f"uniprot:{query}:{format}:{size}"
        
        if use_cache:
            cached_data = self.cache.get(cache_key)
            if cached_data is not None:
                return APIResponse(
                    data=cached_data,
                    status_code=200,
                    headers={},
                    cached=True
                )
        
        self._rate_limit_wait()
        
        try:
            url = f"{self.base_url}/uniprotkb/search"
            params = {
                'query': query,
                'format': format,
                'size': size
            }
            
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse response
            if format == 'json':
                data = response.json()
            else:
                data = response.text
            
            if use_cache:
                self.cache.set(cache_key, data)
            
            return APIResponse(
                data=data,
                status_code=response.status_code,
                headers=dict(response.headers),
                cached=False
            )
            
        except requests.exceptions.RequestException as e:
            raise APIError(
                f"UniProt API request failed: {e}",
                api_name="UniProt",
                endpoint="search",
                status_code=getattr(e.response, 'status_code', None),
                error_code=ErrorCode.API_CONNECTION_FAILED,
                original_error=e
            )


class DatabaseClientFactory:
    """Factory for creating database API clients."""
    
    @staticmethod
    def create_client(database_type: DatabaseType, config=None):
        """Create appropriate API client for database type."""
        clients = {
            DatabaseType.KEGG: KEGGClient,
            DatabaseType.NCBI: NCBIClient,
            DatabaseType.UNIPROT: UniProtClient
        }
        
        client_class = clients.get(database_type)
        if not client_class:
            raise ValidationError(f"Unsupported database type: {database_type}")
        
        return client_class(config)
    
    @staticmethod
    def create_all_clients(config=None) -> Dict[str, Any]:
        """Create all available API clients."""
        return {
            'kegg': KEGGClient(config),
            'ncbi': NCBIClient(config),
            'uniprot': UniProtClient(config)
        }