#!/usr/bin/env python3
"""
Main application entry point for the bioinformatics analysis system.

This unified system orchestrates KEGG pathway analysis, BLAST searches, 
metabolic gene cluster analysis, and phylogenetic studies through a 
cohesive architecture.
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from core import Config, get_config
from core.exceptions import BioinformaticsError
from workflows import (
    KEGGAnalysisWorkflow, 
    BLASTAnalysisWorkflow, 
    MGCAnalysisWorkflow,
    PhylogeneticWorkflow
)
from utils import FileManager, setup_logging


class BioinformaticsApplication:
    """
    Main application orchestrator for bioinformatics analysis workflows.
    
    Provides unified entry point for all analysis pipelines with shared
    configuration, logging, and error handling.
    """
    
    def __init__(self):
        self.config = get_config()
        self.file_manager = FileManager()
        self.logger = logging.getLogger(__name__)
        
        # Initialize available workflows
        self.workflows = {
            'kegg': KEGGAnalysisWorkflow(),
            'blast': BLASTAnalysisWorkflow(),
            'mgc': MGCAnalysisWorkflow(),
            'phylo': PhylogeneticWorkflow()
        }
    
    def setup_workspace(self, workspace_path: Path) -> None:
        """Setup unified workspace for all analyses."""
        self.logger.info(f"Setting up workspace: {workspace_path}")
        
        # Create standard directory structure
        directories = [
            workspace_path / "data" / "raw",
            workspace_path / "data" / "processed", 
            workspace_path / "results" / "kegg",
            workspace_path / "results" / "blast",
            workspace_path / "results" / "mgc",
            workspace_path / "results" / "phylo",
            workspace_path / "logs",
            workspace_path / "cache",
            workspace_path / "config"
        ]
        
        for directory in directories:
            self.file_manager.ensure_directory(directory)
    
    def run_workflow(self, workflow_name: str, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run a specific workflow with unified error handling."""
        if workflow_name not in self.workflows:
            raise BioinformaticsError(f"Unknown workflow: {workflow_name}")
        
        self.logger.info(f"Starting {workflow_name} workflow")
        
        try:
            workflow = self.workflows[workflow_name]
            result = workflow.execute(params)
            
            self.logger.info(f"Completed {workflow_name} workflow successfully")
            return result
            
        except Exception as e:
            self.logger.error(f"Workflow {workflow_name} failed: {e}")
            raise
    
    def run_integrated_pipeline(self, pipeline_config: Dict[str, Any]) -> Dict[str, Any]:
        """Run integrated pipeline with multiple workflows."""
        self.logger.info("Starting integrated bioinformatics pipeline")
        
        results = {}
        
        # Execute workflows in dependency order
        workflow_order = pipeline_config.get('workflow_order', [])
        
        for workflow_name in workflow_order:
            if workflow_name in pipeline_config:
                params = pipeline_config[workflow_name]
                
                # Inject results from previous workflows
                params['previous_results'] = results
                
                # Run workflow
                workflow_result = self.run_workflow(workflow_name, params)
                results[workflow_name] = workflow_result
        
        self.logger.info("Integrated pipeline completed successfully")
        return results


def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Unified bioinformatics analysis system",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        'workflow',
        choices=['kegg', 'blast', 'mgc', 'phylo', 'pipeline'],
        help='Workflow to execute'
    )
    
    parser.add_argument(
        '--config',
        type=Path,
        help='Configuration file path'
    )
    
    parser.add_argument(
        '--workspace',
        type=Path,
        default=Path.cwd() / "workspace",
        help='Workspace directory (default: ./workspace)'
    )
    
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Logging level'
    )
    
    # Workflow-specific arguments
    parser.add_argument('--organism', help='Organism code for KEGG analysis')
    parser.add_argument('--query-file', type=Path, help='Query file for BLAST')
    parser.add_argument('--database', help='Database path for BLAST')
    parser.add_argument('--mgc-candidates', type=Path, help='MGC candidates file')
    parser.add_argument('--tree-file', type=Path, help='Tree file for phylogenetic analysis')
    
    return parser


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(level=args.log_level, log_dir=args.workspace / "logs")
    
    try:
        # Initialize application
        app = BioinformaticsApplication()
        
        # Setup workspace
        app.setup_workspace(args.workspace)
        
        # Load configuration if provided
        if args.config:
            app.config.load_from_file(args.config)
        
        # Prepare workflow parameters
        params = {
            'workspace': args.workspace,
            'organism': args.organism,
            'query_file': args.query_file,
            'database': args.database,
            'mgc_candidates': args.mgc_candidates,
            'tree_file': args.tree_file
        }
        
        # Filter None values
        params = {k: v for k, v in params.items() if v is not None}
        
        # Run workflow
        if args.workflow == 'pipeline':
            if not args.config:
                raise BioinformaticsError("Pipeline workflow requires --config file")
            
            pipeline_config = app.config.load_pipeline_config(args.config)
            result = app.run_integrated_pipeline(pipeline_config)
        else:
            result = app.run_workflow(args.workflow, params)
        
        print(f"✅ {args.workflow} workflow completed successfully")
        if result:
            print(f"Results: {result}")
    
    except BioinformaticsError as e:
        print(f"❌ Error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n⚠️  Operation cancelled by user")
        sys.exit(130)
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()