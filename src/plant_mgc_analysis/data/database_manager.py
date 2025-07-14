"""
Database Management Layer for Plant MGC Analysis Pipeline.

This module provides unified access to all biological databases (KEGG, MiBIG, 
PlantCyc, etc.) with caching, connection pooling, and comprehensive error handling.
"""

import asyncio
import sqlite3
import time
from typing import Dict, List, Optional, Any, Type, Union, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import pickle
from collections import defaultdict

import pandas as pd
import aiofiles
import aiosqlite

from ..core.base import DatabaseInterface, BioinformaticsProcessor, LoggerMixin
from ..core.types import DatabaseType, GeneInfo, MGCCandidate
from ..core.exceptions import DatabaseError, ValidationError, ConfigurationError
from ..metabolic.kegg_integration import KEGGDatabase, KEGGGene, KEGGModule
from ..utils.file_operations import SafeFileOperations, CsvProcessor
from ..config.settings import get_settings


@dataclass
class DatabaseConfig:
    """Configuration for database connections."""
    
    database_type: DatabaseType
    connection_params: Dict[str, Any] = field(default_factory=dict)
    cache_enabled: bool = True
    cache_ttl: int = 3600  # Cache TTL in seconds
    max_connections: int = 10
    timeout: int = 30
    retry_attempts: int = 3
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        if self.cache_ttl < 0:
            raise ConfigurationError("Cache TTL must be non-negative")
        if self.max_connections <= 0:
            raise ConfigurationError("Max connections must be positive")


@dataclass
class CacheEntry:
    """Cache entry with expiration tracking."""
    
    data: Any
    timestamp: float
    ttl: int
    hit_count: int = 0
    
    @property
    def is_expired(self) -> bool:
        """Check if cache entry has expired."""
        return time.time() - self.timestamp > self.ttl
    
    def access(self) -> Any:
        """Access cached data and increment hit count."""
        self.hit_count += 1
        return self.data


class DatabaseCache:
    """High-performance caching layer for database operations."""
    
    def __init__(self, max_size: int = 10000, default_ttl: int = 3600):
        """
        Initialize database cache.
        
        Args:
            max_size: Maximum number of cached entries
            default_ttl: Default TTL for cache entries
        """
        self.max_size = max_size
        self.default_ttl = default_ttl
        self._cache: Dict[str, CacheEntry] = {}
        self._access_order: List[str] = []
        
    def get(self, key: str) -> Optional[Any]:
        """Get cached value by key."""
        if key not in self._cache:
            return None
        
        entry = self._cache[key]
        if entry.is_expired:
            self.delete(key)
            return None
        
        # Update access order (LRU)
        if key in self._access_order:
            self._access_order.remove(key)
        self._access_order.append(key)
        
        return entry.access()
    
    def set(self, key: str, value: Any, ttl: Optional[int] = None) -> None:
        """Set cached value with optional TTL."""
        if len(self._cache) >= self.max_size:
            self._evict_lru()
        
        entry = CacheEntry(
            data=value,
            timestamp=time.time(),
            ttl=ttl or self.default_ttl
        )
        
        self._cache[key] = entry
        if key in self._access_order:
            self._access_order.remove(key)
        self._access_order.append(key)
    
    def delete(self, key: str) -> None:
        """Delete cached entry."""
        if key in self._cache:
            del self._cache[key]
        if key in self._access_order:
            self._access_order.remove(key)
    
    def clear(self) -> None:
        """Clear all cached entries."""
        self._cache.clear()
        self._access_order.clear()
    
    def _evict_lru(self) -> None:
        """Evict least recently used entry."""
        if self._access_order:
            lru_key = self._access_order[0]
            self.delete(lru_key)
    
    def stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        total_hits = sum(entry.hit_count for entry in self._cache.values())
        return {
            "size": len(self._cache),
            "max_size": self.max_size,
            "total_hits": total_hits,
            "expired_entries": sum(1 for entry in self._cache.values() if entry.is_expired),
        }


class DatabaseConnectionPool:
    """Connection pool for database management."""
    
    def __init__(self, config: DatabaseConfig):
        """Initialize connection pool."""
        self.config = config
        self._connections: List[Any] = []
        self._available: List[Any] = []
        self._lock = asyncio.Lock()
    
    async def get_connection(self) -> Any:
        """Get available connection from pool."""
        async with self._lock:
            if self._available:
                return self._available.pop()
            
            if len(self._connections) < self.config.max_connections:
                conn = await self._create_connection()
                self._connections.append(conn)
                return conn
            
            # Wait for connection to become available
            while not self._available:
                await asyncio.sleep(0.1)
            
            return self._available.pop()
    
    async def return_connection(self, connection: Any) -> None:
        """Return connection to pool."""
        async with self._lock:
            self._available.append(connection)
    
    async def _create_connection(self) -> Any:
        """Create new database connection."""
        # This is a placeholder - specific implementations would create
        # actual database connections based on database type
        return None
    
    async def close_all(self) -> None:
        """Close all connections in pool."""
        async with self._lock:
            for conn in self._connections:
                if hasattr(conn, 'close'):
                    if asyncio.iscoroutinefunction(conn.close):
                        await conn.close()
                    else:
                        conn.close()
            
            self._connections.clear()
            self._available.clear()


class UnifiedDatabaseManager(BioinformaticsProcessor):
    """
    Unified database manager providing access to all biological databases.
    
    This class manages connections to multiple databases (KEGG, MiBIG, etc.)
    and provides a unified interface for data access with caching and 
    performance optimization.
    """
    
    def __init__(self, **kwargs):
        """Initialize unified database manager."""
        super().__init__(**kwargs)
        
        self.settings = get_settings()
        self.cache = DatabaseCache(
            max_size=self.settings.cache_max_size,
            default_ttl=self.settings.cache_default_ttl
        )
        
        # Database instances
        self._databases: Dict[DatabaseType, DatabaseInterface] = {}
        self._connection_pools: Dict[DatabaseType, DatabaseConnectionPool] = {}
        
        # Statistics
        self._query_stats = defaultdict(int)
        self._error_stats = defaultdict(int)
        
        # Initialize databases
        self._initialize_databases()
    
    def _initialize_databases(self) -> None:
        """Initialize all database connections."""
        # KEGG Database
        kegg_config = DatabaseConfig(
            database_type=DatabaseType.KEGG,
            connection_params={"rate_limit_delay": 0.35}
        )
        self._databases[DatabaseType.KEGG] = KEGGDatabase()
        
        # Add other databases as needed
        # self._databases[DatabaseType.MIBIG] = MiBIGDatabase()
        # self._databases[DatabaseType.PLANTCYC] = PlantCycDatabase()
    
    def get_database(self, db_type: DatabaseType) -> DatabaseInterface:
        """
        Get database instance by type.
        
        Args:
            db_type: Database type
            
        Returns:
            Database interface instance
            
        Raises:
            ValidationError: If database type not supported
        """
        if db_type not in self._databases:
            raise ValidationError(
                f"Database type {db_type} not supported",
                field_name="database_type",
                field_value=db_type
            )
        
        return self._databases[db_type]
    
    async def query_database(
        self,
        db_type: DatabaseType,
        query: str,
        parameters: Optional[Dict[str, Any]] = None,
        use_cache: bool = True,
    ) -> List[Dict[str, Any]]:
        """
        Execute query against specified database.
        
        Args:
            db_type: Database type
            query: Query string
            parameters: Optional query parameters
            use_cache: Whether to use caching
            
        Returns:
            Query results
        """
        # Generate cache key
        cache_key = f"{db_type}:{query}:{hash(str(parameters))}"
        
        # Check cache first
        if use_cache:
            cached_result = self.cache.get(cache_key)
            if cached_result is not None:
                self._query_stats[f"{db_type}_cache_hit"] += 1
                return cached_result
        
        try:
            # Execute query
            database = self.get_database(db_type)
            result = database.query(query, parameters)
            
            # Cache result
            if use_cache:
                self.cache.set(cache_key, result)
            
            self._query_stats[f"{db_type}_query"] += 1
            return result
            
        except Exception as e:
            self._error_stats[f"{db_type}_error"] += 1
            raise DatabaseError(
                f"Query failed for {db_type}: {e}",
                database_name=str(db_type),
                operation="query"
            ) from e
    
    def annotate_genes(
        self,
        genes: List[GeneInfo],
        organism_code: str,
        databases: Optional[List[DatabaseType]] = None,
    ) -> Dict[str, Dict[str, Any]]:
        """
        Annotate genes using multiple databases.
        
        Args:
            genes: List of genes to annotate
            organism_code: Organism code
            databases: List of databases to use (default: all)
            
        Returns:
            Dictionary mapping gene IDs to annotation data
        """
        if databases is None:
            databases = list(self._databases.keys())
        
        annotations = {}
        gene_ids = [gene.gene_id for gene in genes]
        
        for db_type in databases:
            try:
                if db_type == DatabaseType.KEGG:
                    kegg_db = self.get_database(db_type)
                    kegg_annotations = kegg_db.get_gene_batch(gene_ids, organism_code)
                    
                    for gene_id, kegg_gene in kegg_annotations.items():
                        if gene_id not in annotations:
                            annotations[gene_id] = {}
                        
                        annotations[gene_id][f"{db_type}_annotation"] = kegg_gene.annotation
                        annotations[gene_id][f"{db_type}_pathways"] = kegg_gene.pathways
                        annotations[gene_id][f"{db_type}_ec_numbers"] = kegg_gene.ec_numbers
                
                # Add other database annotations here
                
            except Exception as e:
                self.logger.error(f"Failed to annotate with {db_type}: {e}")
        
        return annotations
    
    def find_metabolic_genes(
        self,
        genes: List[GeneInfo],
        organism_code: str,
        confidence_threshold: float = 0.7,
    ) -> List[GeneInfo]:
        """
        Identify metabolic genes using database annotations.
        
        Args:
            genes: List of genes to analyze
            organism_code: Organism code
            confidence_threshold: Confidence threshold for classification
            
        Returns:
            List of metabolic genes
        """
        annotations = self.annotate_genes(genes, organism_code)
        metabolic_genes = []
        
        for gene in genes:
            gene_annotation = annotations.get(gene.gene_id, {})
            
            # Check for metabolic indicators
            is_metabolic = False
            confidence = 0.0
            
            # KEGG-based classification
            if f"{DatabaseType.KEGG}_pathways" in gene_annotation:
                pathways = gene_annotation[f"{DatabaseType.KEGG}_pathways"]
                if pathways:
                    is_metabolic = True
                    confidence += 0.5
            
            if f"{DatabaseType.KEGG}_ec_numbers" in gene_annotation:
                ec_numbers = gene_annotation[f"{DatabaseType.KEGG}_ec_numbers"]
                if ec_numbers:
                    is_metabolic = True
                    confidence += 0.5
            
            # Add annotation data to gene
            if hasattr(gene, 'annotations'):
                gene.annotations.update(gene_annotation)
            else:
                gene.annotations = gene_annotation
            
            if is_metabolic and confidence >= confidence_threshold:
                metabolic_genes.append(gene)
        
        return metabolic_genes
    
    def build_metabolic_network(
        self,
        genes: List[GeneInfo],
        organism_code: str,
    ) -> Dict[str, List[str]]:
        """
        Build metabolic network from gene annotations.
        
        Args:
            genes: List of genes
            organism_code: Organism code
            
        Returns:
            Dictionary representing metabolic network
        """
        annotations = self.annotate_genes(genes, organism_code)
        network = defaultdict(list)
        
        # Group genes by pathway
        pathway_genes = defaultdict(list)
        
        for gene in genes:
            gene_annotation = annotations.get(gene.gene_id, {})
            pathways = gene_annotation.get(f"{DatabaseType.KEGG}_pathways", [])
            
            for pathway in pathways:
                pathway_genes[pathway].append(gene.gene_id)
        
        # Build network connections
        for pathway, gene_list in pathway_genes.items():
            for gene in gene_list:
                network[gene].extend([g for g in gene_list if g != gene])
        
        return dict(network)
    
    def export_annotations(
        self,
        annotations: Dict[str, Dict[str, Any]],
        output_file: Path,
        format: str = "csv",
    ) -> None:
        """
        Export annotations to file.
        
        Args:
            annotations: Annotation data
            output_file: Output file path
            format: Export format (csv, json, excel)
        """
        if format == "csv":
            # Flatten nested annotations
            rows = []
            for gene_id, gene_annotations in annotations.items():
                row = {"gene_id": gene_id}
                
                for key, value in gene_annotations.items():
                    if isinstance(value, list):
                        row[key] = ";".join(str(v) for v in value)
                    else:
                        row[key] = str(value)
                
                rows.append(row)
            
            df = pd.DataFrame(rows)
            processor = CsvProcessor()
            processor.write_file(df, output_file)
        
        elif format == "json":
            with open(output_file, 'w') as f:
                json.dump(annotations, f, indent=2, default=str)
        
        elif format == "excel":
            # Flatten and write to Excel
            rows = []
            for gene_id, gene_annotations in annotations.items():
                row = {"gene_id": gene_id}
                
                for key, value in gene_annotations.items():
                    if isinstance(value, list):
                        row[key] = ";".join(str(v) for v in value)
                    else:
                        row[key] = str(value)
                
                rows.append(row)
            
            df = pd.DataFrame(rows)
            df.to_excel(output_file, index=False)
        
        else:
            raise ValidationError(f"Unsupported export format: {format}")
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database manager statistics."""
        return {
            "query_stats": dict(self._query_stats),
            "error_stats": dict(self._error_stats),
            "cache_stats": self.cache.stats(),
            "available_databases": list(self._databases.keys()),
        }
    
    def clear_cache(self) -> None:
        """Clear all cached data."""
        self.cache.clear()
        self.logger.info("Database cache cleared")
    
    def validate_input(self, data: Any) -> None:
        """Validate input data."""
        if not isinstance(data, (list, dict, str)):
            raise ValidationError(
                "Input must be gene list, query dict, or database name",
                field_name="input_data",
                field_value=type(data).__name__
            )
    
    def process(self, data: Any, **kwargs) -> Dict[str, Any]:
        """
        Process database management request.
        
        Args:
            data: Input data (genes, query, etc.)
            **kwargs: Additional parameters
            
        Returns:
            Processing results
        """
        self.validate_input(data)
        
        if isinstance(data, list):
            # Gene annotation request
            organism_code = kwargs.get("organism_code", "")
            if not organism_code:
                raise ValidationError("Organism code required for gene annotation")
            
            annotations = self.annotate_genes(data, organism_code)
            return {"annotations": annotations}
        
        elif isinstance(data, dict):
            # Complex query request
            db_type = data.get("database")
            query = data.get("query")
            parameters = data.get("parameters")
            
            if not db_type or not query:
                raise ValidationError("Database type and query required")
            
            result = asyncio.run(
                self.query_database(
                    DatabaseType(db_type),
                    query,
                    parameters
                )
            )
            return {"results": result}
        
        elif isinstance(data, str):
            # Statistics or info request
            if data == "stats":
                return self.get_statistics()
            elif data == "clear_cache":
                self.clear_cache()
                return {"status": "cache_cleared"}
        
        return {}
    
    def __enter__(self):
        """Context manager entry."""
        for db in self._databases.values():
            if hasattr(db, 'connect'):
                db.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        for db in self._databases.values():
            if hasattr(db, 'disconnect'):
                db.disconnect()


class DatabaseFactory:
    """Factory for creating database instances."""
    
    _database_types = {
        DatabaseType.KEGG: KEGGDatabase,
        # Add other database types here
    }
    
    @classmethod
    def create_database(
        self,
        db_type: DatabaseType,
        config: Optional[DatabaseConfig] = None,
        **kwargs
    ) -> DatabaseInterface:
        """
        Create database instance.
        
        Args:
            db_type: Database type
            config: Optional database configuration
            **kwargs: Additional parameters
            
        Returns:
            Database interface instance
        """
        if db_type not in self._database_types:
            raise ValidationError(f"Unsupported database type: {db_type}")
        
        database_class = self._database_types[db_type]
        return database_class(**kwargs)
    
    @classmethod
    def register_database(
        self,
        db_type: DatabaseType,
        database_class: Type[DatabaseInterface]
    ) -> None:
        """Register new database type."""
        self._database_types[db_type] = database_class