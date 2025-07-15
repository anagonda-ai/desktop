"""
BLAST Analysis Workflow orchestrating database creation, sequence searches,
and result processing into a unified pipeline.
"""

import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional
import pandas as pd
from Bio.Blast import NCBIXML

from .base_workflow import BaseWorkflow
from core.exceptions import ProcessingError, ValidationError
from core.types import ProcessingStatus
from utils import BatchProcessor, ProgressTracker, FileManager, ExecutorType


class BLASTAnalysisWorkflow(BaseWorkflow):
    """
    Unified BLAST analysis workflow orchestrating database creation,
    sequence searches, and comprehensive result processing.
    
    This workflow coordinates:
    1. Database preparation and indexing
    2. Query sequence processing
    3. BLAST search execution
    4. Result parsing and analysis
    5. Hit filtering and reporting
    """
    
    def __init__(self):
        super().__init__("BLAST")
        
        # Workflow components
        self.file_manager = FileManager()
        self.batch_processor = BatchProcessor(
            executor_type=ExecutorType.PROCESS,
            max_workers=4
        )
        
        # Workflow state
        self.database_prepared = False
        self.searches_completed = 0
        self.total_hits_found = 0
    
    def get_required_parameters(self) -> List[str]:
        """Get required parameters for BLAST workflow."""
        return ["workspace", "query_sequences", "target_databases"]
    
    def get_required_subdirectories(self) -> List[str]:
        """Get required subdirectories for BLAST workflow."""
        return [
            "databases", 
            "queries", 
            "results",
            "filtered_results",
            "reports",
            "temp"
        ]
    
    def validate_parameters(self, params: Dict[str, Any]) -> None:
        """Validate BLAST workflow parameters."""
        super().validate_parameters(params)
        
        # Validate query sequences
        query_sequences = params.get('query_sequences')
        if isinstance(query_sequences, (str, Path)):
            query_path = Path(query_sequences)
            if not query_path.exists():
                raise ValidationError(f"Query sequences file not found: {query_path}")
        elif not isinstance(query_sequences, list):
            raise ValidationError("query_sequences must be a file path or list of sequences")
        
        # Validate target databases
        target_databases = params.get('target_databases')
        if isinstance(target_databases, (str, Path)):
            db_path = Path(target_databases)
            if not db_path.exists():
                raise ValidationError(f"Target database not found: {db_path}")
        elif not isinstance(target_databases, list):
            raise ValidationError("target_databases must be a file path or list of database paths")
        
        # Validate BLAST parameters
        blast_type = params.get('blast_type', 'blastp')
        valid_blast_types = ['blastp', 'blastn', 'blastx', 'tblastn', 'tblastx']
        if blast_type not in valid_blast_types:
            raise ValidationError(f"Invalid blast_type: {blast_type}. Must be one of: {valid_blast_types}")
    
    def _prepare_databases(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare and index BLAST databases."""
        self.logger.info("Starting database preparation phase")
        
        output_paths = self.get_output_paths()
        databases_dir = output_paths['databases']
        
        target_databases = params['target_databases']
        if isinstance(target_databases, (str, Path)):
            target_databases = [target_databases]
        
        prepared_databases = []
        
        for db_source in target_databases:
            db_source = Path(db_source)
            
            # Determine database name and type
            db_name = db_source.stem
            db_type = params.get('db_type', 'prot')  # 'prot' or 'nucl'
            
            # Target database path
            target_db_path = databases_dir / f"{db_name}_blastdb"
            
            # Check if database already exists
            if self._database_exists(target_db_path, db_type):
                self.logger.info(f"✅ Database already exists: {target_db_path}")
                prepared_databases.append(str(target_db_path))
                continue
            
            try:
                # Build BLAST database
                self.logger.info(f"Building BLAST database from {db_source}")
                
                cmd = [
                    "makeblastdb",
                    "-in", str(db_source),
                    "-dbtype", db_type,
                    "-out", str(target_db_path),
                    "-parse_seqids"
                ]
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                self.logger.info(f"✅ Created database: {target_db_path}")
                prepared_databases.append(str(target_db_path))
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Failed to create database {db_name}: {e.stderr}")
                raise ProcessingError(f"Database creation failed: {e.stderr}")
        
        self.database_prepared = True
        
        return {
            'prepared_databases': prepared_databases,
            'total_databases': len(prepared_databases)
        }
    
    def _database_exists(self, db_path: Path, db_type: str) -> bool:
        """Check if BLAST database files exist."""
        if db_type == 'prot':
            extensions = ['.pin', '.phr', '.psq']
        else:
            extensions = ['.nin', '.nhr', '.nsq']
        
        # Check for regular or multi-volume database files
        for ext in extensions:
            if (db_path.parent / f"{db_path.name}{ext}").exists():
                return True
            if (db_path.parent / f"{db_path.name}.00{ext}").exists():
                return True
        
        return False
    
    def _prepare_queries(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare and validate query sequences."""
        self.logger.info("Preparing query sequences")
        
        output_paths = self.get_output_paths()
        queries_dir = output_paths['queries']
        
        query_sequences = params['query_sequences']
        
        if isinstance(query_sequences, (str, Path)):
            # Copy query file to workspace
            query_file = Path(query_sequences)
            target_query = queries_dir / query_file.name
            
            import shutil
            shutil.copy2(query_file, target_query)
            
            # Validate FASTA format
            query_count = self._validate_fasta(target_query)
            
            return {
                'query_file': str(target_query),
                'query_count': query_count
            }
        else:
            # Create FASTA file from sequence list
            target_query = queries_dir / "queries.fasta"
            
            with open(target_query, 'w') as f:
                for i, seq in enumerate(query_sequences):
                    f.write(f">query_{i+1}\n{seq}\n")
            
            return {
                'query_file': str(target_query),
                'query_count': len(query_sequences)
            }
    
    def _validate_fasta(self, fasta_file: Path) -> int:
        """Validate FASTA file and count sequences."""
        sequence_count = 0
        
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        sequence_count += 1
            
            if sequence_count == 0:
                raise ValidationError(f"No sequences found in FASTA file: {fasta_file}")
            
            self.logger.info(f"Validated {sequence_count} sequences in {fasta_file}")
            return sequence_count
            
        except Exception as e:
            raise ValidationError(f"Invalid FASTA file {fasta_file}: {e}")
    
    def _execute_blast_searches(self, params: Dict[str, Any], 
                               query_info: Dict[str, Any], 
                               db_info: Dict[str, Any]) -> Dict[str, Any]:
        """Execute BLAST searches against all databases."""
        self.logger.info("Starting BLAST search execution")
        
        output_paths = self.get_output_paths()
        results_dir = output_paths['results']
        
        query_file = query_info['query_file']
        databases = db_info['prepared_databases']
        
        # BLAST parameters
        blast_type = params.get('blast_type', 'blastp')
        evalue = params.get('evalue', 1e-3)
        num_threads = params.get('num_threads', 8)
        max_target_seqs = params.get('max_target_seqs', 500)
        
        search_results = []
        
        tracker = ProgressTracker(len(databases), "BLAST searches")
        
        for database in databases:
            db_name = Path(database).name
            result_file = results_dir / f"{db_name}_results.xml"
            
            try:
                self.logger.info(f"Running {blast_type} against {db_name}")
                
                cmd = [
                    blast_type,
                    "-query", query_file,
                    "-db", database,
                    "-out", str(result_file),
                    "-outfmt", "5",  # XML format
                    "-evalue", str(evalue),
                    "-num_threads", str(num_threads),
                    "-max_target_seqs", str(max_target_seqs)
                ]
                
                # Add additional parameters
                if 'word_size' in params:
                    cmd.extend(["-word_size", str(params['word_size'])])
                if 'gapopen' in params:
                    cmd.extend(["-gapopen", str(params['gapopen'])])
                if 'gapextend' in params:
                    cmd.extend(["-gapextend", str(params['gapextend'])])
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Parse result summary
                hit_count = self._count_blast_hits(result_file)
                
                search_results.append({
                    'database': db_name,
                    'result_file': str(result_file),
                    'hits_found': hit_count,
                    'status': 'success'
                })
                
                self.total_hits_found += hit_count
                self.searches_completed += 1
                
                self.logger.info(f"✅ {blast_type} vs {db_name}: {hit_count} hits")
                tracker.update(1)
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"BLAST search failed for {db_name}: {e.stderr}")
                search_results.append({
                    'database': db_name,
                    'result_file': None,
                    'hits_found': 0,
                    'status': 'failed',
                    'error': e.stderr
                })
                tracker.mark_failed(1)
        
        return {
            'search_results': search_results,
            'total_searches': len(databases),
            'successful_searches': self.searches_completed,
            'total_hits': self.total_hits_found
        }
    
    def _count_blast_hits(self, xml_file: Path) -> int:
        """Count total hits in BLAST XML file."""
        try:
            hit_count = 0
            with open(xml_file) as handle:
                for record in NCBIXML.parse(handle):
                    hit_count += len(record.alignments)
            return hit_count
        except Exception as e:
            self.logger.warning(f"Could not count hits in {xml_file}: {e}")
            return 0
    
    def _process_and_filter_results(self, params: Dict[str, Any], 
                                   search_results: Dict[str, Any]) -> Dict[str, Any]:
        """Process and filter BLAST results."""
        self.logger.info("Processing and filtering BLAST results")
        
        output_paths = self.get_output_paths()
        filtered_dir = output_paths['filtered_results']
        
        # Filtering criteria
        min_identity = params.get('min_identity', 30.0)
        min_coverage = params.get('min_coverage', 50.0)
        max_evalue = params.get('max_evalue', 1e-5)
        
        all_filtered_hits = []
        processing_summary = []
        
        for search_result in search_results['search_results']:
            if search_result['status'] != 'success':
                continue
            
            result_file = Path(search_result['result_file'])
            db_name = search_result['database']
            
            try:
                # Process XML results
                filtered_hits = self._filter_blast_xml(
                    result_file, 
                    min_identity, 
                    min_coverage, 
                    max_evalue
                )
                
                # Save filtered results
                if filtered_hits:
                    filtered_file = filtered_dir / f"{db_name}_filtered.csv"
                    df = pd.DataFrame(filtered_hits)
                    df.to_csv(filtered_file, index=False)
                    
                    all_filtered_hits.extend(filtered_hits)
                    
                    processing_summary.append({
                        'database': db_name,
                        'original_hits': search_result['hits_found'],
                        'filtered_hits': len(filtered_hits),
                        'filtered_file': str(filtered_file)
                    })
                    
                    self.logger.info(
                        f"Filtered {db_name}: "
                        f"{search_result['hits_found']} -> {len(filtered_hits)} hits"
                    )
                
            except Exception as e:
                self.logger.error(f"Failed to process results for {db_name}: {e}")
        
        # Save combined filtered results
        if all_filtered_hits:
            combined_file = filtered_dir / "all_filtered_hits.csv"
            df = pd.DataFrame(all_filtered_hits)
            df.to_csv(combined_file, index=False)
            
            self.logger.info(f"Saved {len(all_filtered_hits)} total filtered hits to {combined_file}")
        
        return {
            'processing_summary': processing_summary,
            'total_filtered_hits': len(all_filtered_hits),
            'filtering_criteria': {
                'min_identity': min_identity,
                'min_coverage': min_coverage,
                'max_evalue': max_evalue
            }
        }
    
    def _filter_blast_xml(self, xml_file: Path, min_identity: float, 
                         min_coverage: float, max_evalue: float) -> List[Dict[str, Any]]:
        """Filter BLAST XML results based on criteria."""
        filtered_hits = []
        
        with open(xml_file) as handle:
            for record in NCBIXML.parse(handle):
                query_id = record.query
                query_length = record.query_length
                
                for alignment in record.alignments:
                    hit_id = alignment.hit_def
                    hit_length = alignment.length
                    
                    for hsp in alignment.hsps:
                        # Calculate metrics
                        identity = (hsp.identities / hsp.align_length) * 100
                        coverage = (hsp.align_length / query_length) * 100
                        evalue = hsp.expect
                        
                        # Apply filters
                        if (identity >= min_identity and 
                            coverage >= min_coverage and 
                            evalue <= max_evalue):
                            
                            filtered_hits.append({
                                'query_id': query_id,
                                'hit_id': hit_id,
                                'identity_percent': identity,
                                'coverage_percent': coverage,
                                'evalue': evalue,
                                'bit_score': hsp.bits,
                                'alignment_length': hsp.align_length,
                                'query_start': hsp.query_start,
                                'query_end': hsp.query_end,
                                'hit_start': hsp.sbjct_start,
                                'hit_end': hsp.sbjct_end
                            })
        
        return filtered_hits
    
    def _execute_workflow(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Execute the BLAST analysis workflow."""
        
        # Phase 1: Prepare databases
        db_results = self._prepare_databases(params)
        
        # Phase 2: Prepare queries
        query_results = self._prepare_queries(params)
        
        # Phase 3: Execute BLAST searches
        search_results = self._execute_blast_searches(params, query_results, db_results)
        
        # Phase 4: Process and filter results
        filter_results = self._process_and_filter_results(params, search_results)
        
        # Integrate all results
        integrated_results = {
            'workflow_summary': {
                'databases_prepared': db_results['total_databases'],
                'queries_processed': query_results['query_count'],
                'searches_completed': self.searches_completed,
                'total_hits_found': self.total_hits_found,
                'filtered_hits': filter_results['total_filtered_hits']
            },
            'database_preparation': db_results,
            'query_preparation': query_results,
            'blast_searches': search_results,
            'result_filtering': filter_results,
            'output_paths': {k: str(v) for k, v in self.get_output_paths().items()}
        }
        
        return integrated_results
    
    def post_execution_hook(self, result) -> None:
        """Generate comprehensive BLAST workflow report."""
        self._generate_blast_report(result)
    
    def _generate_blast_report(self, result) -> None:
        """Generate detailed BLAST analysis report."""
        output_paths = self.get_output_paths()
        report_file = output_paths['reports'] / "blast_workflow_summary.md"
        
        with open(report_file, 'w') as f:
            f.write("# BLAST Analysis Workflow Summary\n\n")
            f.write(f"**Execution Time:** {result.execution_time:.2f} seconds\n")
            f.write(f"**Status:** {result.status.value}\n\n")
            
            if result.results:
                summary = result.results.get('workflow_summary', {})
                f.write("## Analysis Summary\n\n")
                f.write(f"- **Databases Prepared:** {summary.get('databases_prepared', 0)}\n")
                f.write(f"- **Queries Processed:** {summary.get('queries_processed', 0)}\n")
                f.write(f"- **Searches Completed:** {summary.get('searches_completed', 0)}\n")
                f.write(f"- **Total Hits Found:** {summary.get('total_hits_found', 0)}\n")
                f.write(f"- **Filtered Hits:** {summary.get('filtered_hits', 0)}\n\n")
                
                # Add filtering criteria
                filter_info = result.results.get('result_filtering', {})
                if 'filtering_criteria' in filter_info:
                    criteria = filter_info['filtering_criteria']
                    f.write("## Filtering Criteria\n\n")
                    f.write(f"- **Minimum Identity:** {criteria.get('min_identity', 'N/A')}%\n")
                    f.write(f"- **Minimum Coverage:** {criteria.get('min_coverage', 'N/A')}%\n")
                    f.write(f"- **Maximum E-value:** {criteria.get('max_evalue', 'N/A')}\n\n")
            
            # Add output locations
            f.write("## Output Locations\n\n")
            for path_name, path in self.get_output_paths().items():
                f.write(f"- **{path_name.replace('_', ' ').title()}:** `{path}`\n")
        
        self.logger.info(f"Generated BLAST report: {report_file}")