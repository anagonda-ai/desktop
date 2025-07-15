"""
Command-line interface for Plant MGC Analysis Pipeline.

This module provides a unified CLI that replaces all 92 legacy scripts
with a single, cohesive interface.
"""

import sys
import time
from pathlib import Path
from typing import List, Optional, Dict, Any

import click
from loguru import logger

from .main_orchestrator import (
    PlantMGCAnalysisOrchestrator,
    AnalysisJobConfig,
    AnalysisMode,
    run_single_organism_analysis,
    run_batch_analysis,
)
from .core.types import AnalysisType
from .config.settings import get_settings, create_environment_config
from .utils.migration_tools import create_migration_plan, migrate_legacy_scripts


@click.group()
@click.option('--config', type=click.Path(exists=True), help='Configuration file path')
@click.option('--env', type=click.Choice(['development', 'production', 'cluster']), default='production', help='Environment')
@click.option('--debug', is_flag=True, help='Enable debug logging')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.pass_context
def cli(ctx, config, env, debug, verbose):
    """
    Plant MGC Analysis Pipeline - Unified CLI
    
    This tool consolidates all 92 legacy scripts into a single,
    cohesive command-line interface for plant metabolic gene cluster analysis.
    """
    ctx.ensure_object(dict)
    
    # Setup configuration
    if config:
        settings = get_settings().load_config(Path(config))
    else:
        settings = create_environment_config(env)
    
    if debug:
        settings.debug = True
        settings.logging.level = "DEBUG"
    
    ctx.obj['settings'] = settings
    ctx.obj['verbose'] = verbose
    
    # Setup logging
    logger.configure(
        handlers=[
            {
                "sink": sys.stderr,
                "level": settings.logging.level,
                "format": settings.logging.format,
            }
        ]
    )
    
    if verbose:
        click.echo(f"Using environment: {env}")
        click.echo(f"Debug mode: {debug}")
        click.echo(f"Configuration: {config or 'default'}")


@cli.command()
@click.argument('input_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--analysis-types', '-a', multiple=True, 
              type=click.Choice([at.value for at in AnalysisType]), 
              default=['blast_search'], help='Analysis types to run')
@click.option('--output-dir', '-o', type=click.Path(), default='output', help='Output directory')
@click.option('--organism', help='Organism name')
@click.option('--mode', type=click.Choice([m.value for m in AnalysisMode]), 
              default='single_organism', help='Analysis mode')
@click.option('--parallel/--no-parallel', default=True, help='Run analyses in parallel')
@click.option('--max-workers', type=int, default=4, help='Maximum number of parallel workers')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def analyze(ctx, input_files, analysis_types, output_dir, organism, mode, parallel, max_workers, dry_run):
    """
    Run comprehensive analysis on genome files.
    
    This command replaces multiple legacy scripts by providing a unified
    interface for all analysis types.
    """
    settings = ctx.obj['settings']
    verbose = ctx.obj['verbose']
    
    # Convert analysis types
    analysis_type_enums = [AnalysisType(at) for at in analysis_types]
    
    # Create job configuration
    job_config = AnalysisJobConfig(
        job_id=f"cli_analysis_{int(time.time())}",
        analysis_types=analysis_type_enums,
        input_files=list(input_files),
        output_directory=Path(output_dir),
        parameters={
            'organism_name': organism,
            'verbose': verbose,
        },
        mode=AnalysisMode(mode),
        parallel_execution=parallel,
        max_workers=max_workers,
    )
    
    if dry_run:
        click.echo("Dry run mode - showing what would be executed:")
        click.echo(f"Job ID: {job_config.job_id}")
        click.echo(f"Analysis types: {[at.value for at in analysis_type_enums]}")
        click.echo(f"Input files: {list(input_files)}")
        click.echo(f"Output directory: {output_dir}")
        click.echo(f"Mode: {mode}")
        click.echo(f"Parallel: {parallel}")
        return
    
    # Run analysis
    orchestrator = PlantMGCAnalysisOrchestrator(config=settings.model_dump())
    
    with click.progressbar(length=100, label='Running analysis') as bar:
        result = orchestrator.run_analysis_job(job_config)
        bar.update(100)
    
    # Display results
    if result.success:
        click.echo(click.style("✓ Analysis completed successfully!", fg='green'))
        click.echo(f"Results: {result.total_results}")
        click.echo(f"Processing time: {result.processing_time:.2f}s")
        click.echo(f"Output files: {len(result.output_files)}")
        
        if verbose:
            click.echo("\nOutput files:")
            for file_path in result.output_files:
                click.echo(f"  - {file_path}")
    else:
        click.echo(click.style("✗ Analysis failed!", fg='red'))
        for error in result.errors:
            click.echo(f"Error: {error}")


@cli.command()
@click.argument('input_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--database', '-d', help='BLAST database path')
@click.option('--program', type=click.Choice(['blastp', 'blastn', 'blastx', 'tblastn', 'tblastx']), 
              default='blastp', help='BLAST program')
@click.option('--evalue', type=float, default=1e-5, help='E-value threshold')
@click.option('--max-targets', type=int, default=100, help='Maximum target sequences')
@click.option('--output-dir', '-o', type=click.Path(), default='blast_output', help='Output directory')
@click.option('--threads', type=int, default=4, help='Number of threads')
@click.pass_context
def blast(ctx, input_files, database, program, evalue, max_targets, output_dir, threads):
    """
    Run BLAST analysis on sequences.
    
    This command consolidates all BLAST-related legacy scripts into a single interface.
    """
    settings = ctx.obj['settings']
    
    # Prepare parameters
    parameters = {
        'blast_search': {
            'program': program,
            'evalue': evalue,
            'max_targets': max_targets,
            'num_threads': threads,
            'subject_file': database,
        }
    }
    
    # Run analysis
    if len(input_files) == 1:
        result = run_single_organism_analysis(
            input_files[0],
            ['blast_search'],
            output_dir,
            parameters,
            config=settings.model_dump(),
        )
    else:
        result = run_batch_analysis(
            input_files,
            ['blast_search'],
            output_dir,
            parameters,
            config=settings.model_dump(),
        )
    
    # Display results
    if result.success:
        click.echo(click.style("✓ BLAST analysis completed!", fg='green'))
        click.echo(f"Processing time: {result.processing_time:.2f}s")
    else:
        click.echo(click.style("✗ BLAST analysis failed!", fg='red'))
        for error in result.errors:
            click.echo(f"Error: {error}")


@cli.command()
@click.argument('input_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--organisms', '-org', multiple=True, help='Organism codes for KEGG')
@click.option('--pathways', '-p', multiple=True, help='Specific pathway IDs')
@click.option('--output-dir', '-o', type=click.Path(), default='kegg_output', help='Output directory')
@click.option('--batch-size', type=int, default=10, help='Batch size for API requests')
@click.pass_context
def kegg(ctx, input_files, organisms, pathways, output_dir, batch_size):
    """
    Analyze sequences using KEGG pathways.
    
    This command consolidates all KEGG-related legacy scripts.
    """
    settings = ctx.obj['settings']
    
    # Prepare parameters
    parameters = {
        'kegg_analysis': {
            'organisms': list(organisms),
            'pathways': list(pathways),
            'batch_size': batch_size,
        }
    }
    
    # Run analysis
    result = run_batch_analysis(
        input_files,
        ['kegg_analysis'],
        output_dir,
        parameters,
        config=settings.model_dump(),
    )
    
    # Display results
    if result.success:
        click.echo(click.style("✓ KEGG analysis completed!", fg='green'))
        click.echo(f"Processing time: {result.processing_time:.2f}s")
    else:
        click.echo(click.style("✗ KEGG analysis failed!", fg='red'))
        for error in result.errors:
            click.echo(f"Error: {error}")


@cli.command()
@click.argument('input_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--window-size', type=int, default=50000, help='Window size in base pairs')
@click.option('--step-size', type=int, default=10000, help='Step size in base pairs')
@click.option('--min-genes', type=int, default=3, help='Minimum genes per window')
@click.option('--output-dir', '-o', type=click.Path(), default='sliding_window_output', help='Output directory')
@click.pass_context
def sliding_window(ctx, input_files, window_size, step_size, min_genes, output_dir):
    """
    Perform sliding window analysis on genomes.
    
    This command consolidates all sliding window analysis scripts.
    """
    settings = ctx.obj['settings']
    
    # Prepare parameters
    parameters = {
        'sliding_window': {
            'window_size': window_size,
            'step_size': step_size,
            'min_genes': min_genes,
        }
    }
    
    # Run analysis
    result = run_batch_analysis(
        input_files,
        ['sliding_window'],
        output_dir,
        parameters,
        config=settings.model_dump(),
    )
    
    # Display results
    if result.success:
        click.echo(click.style("✓ Sliding window analysis completed!", fg='green'))
        click.echo(f"Processing time: {result.processing_time:.2f}s")
    else:
        click.echo(click.style("✗ Sliding window analysis failed!", fg='red'))
        for error in result.errors:
            click.echo(f"Error: {error}")


@cli.group()
def migrate():
    """Migration tools for legacy scripts."""
    pass


@migrate.command()
@click.argument('directory', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), default='migration_plan.md', help='Output file for migration plan')
@click.pass_context
def plan(ctx, directory, output):
    """
    Create migration plan for legacy scripts.
    
    This analyzes all Python scripts in the directory and creates a comprehensive
    migration plan showing priorities, duplications, and effort estimates.
    """
    directory_path = Path(directory)
    
    click.echo(f"Analyzing scripts in: {directory_path}")
    
    with click.progressbar(length=100, label='Analyzing scripts') as bar:
        migration_plan = create_migration_plan(directory_path)
        bar.update(100)
    
    # Save migration plan
    from .utils.migration_tools import ScriptMigrator
    migrator = ScriptMigrator()
    report = migrator.create_migration_report(migration_plan)
    
    output_path = Path(output)
    with open(output_path, 'w') as f:
        f.write(report)
    
    # Display summary
    summary = migration_plan.get_summary()
    click.echo(f"Migration plan created: {output_path}")
    click.echo(f"Total scripts: {summary['total_scripts']}")
    click.echo(f"High priority: {summary['high_priority']}")
    click.echo(f"Medium priority: {summary['medium_priority']}")
    click.echo(f"Low priority: {summary['low_priority']}")
    click.echo(f"Estimated effort: {summary['estimated_effort_hours']} hours")


@migrate.command()
@click.argument('directory', type=click.Path(exists=True))
@click.option('--priority', type=click.Choice(['high', 'medium', 'low', 'all']), 
              default='high', help='Migration priority')
@click.option('--dry-run', is_flag=True, help='Show what would be migrated without executing')
@click.pass_context
def scripts(ctx, directory, priority, dry_run):
    """
    Migrate legacy scripts to new architecture.
    
    This command migrates legacy scripts based on priority level.
    """
    directory_path = Path(directory)
    
    click.echo(f"Migrating {priority} priority scripts in: {directory_path}")
    
    if dry_run:
        click.echo("Dry run mode - showing what would be migrated:")
    
    with click.progressbar(length=100, label='Migrating scripts') as bar:
        results = migrate_legacy_scripts(directory_path, priority, dry_run)
        bar.update(100)
    
    # Display results
    successful = len([r for r in results if r['transformed']])
    failed = len(results) - successful
    
    click.echo(f"Migration completed:")
    click.echo(f"  Successful: {successful}")
    click.echo(f"  Failed: {failed}")
    
    if ctx.obj['verbose']:
        for result in results:
            status = "✓" if result['transformed'] else "✗"
            click.echo(f"  {status} {result['original_path']}")


@cli.command()
@click.pass_context
def status(ctx):
    """Show system status and statistics."""
    settings = ctx.obj['settings']
    
    # Create orchestrator to get stats
    orchestrator = PlantMGCAnalysisOrchestrator(config=settings.model_dump())
    stats = orchestrator.get_system_stats()
    
    # Display status
    click.echo("Plant MGC Analysis Pipeline Status")
    click.echo("=" * 40)
    
    # Jobs
    jobs = stats.get('jobs', {})
    click.echo(f"Active jobs: {jobs.get('active_count', 0)}")
    click.echo(f"Completed jobs: {jobs.get('completed_count', 0)}")
    click.echo(f"Total jobs: {jobs.get('total_count', 0)}")
    
    # File operations
    file_ops = stats.get('file_operations', {})
    click.echo(f"Files processed: {file_ops.get('files_processed', 0)}")
    click.echo(f"Processing errors: {file_ops.get('error_count', 0)}")
    
    # Configuration
    click.echo(f"\nConfiguration:")
    click.echo(f"  Environment: {'cluster' if settings.is_cluster else 'local'}")
    click.echo(f"  Debug mode: {settings.debug}")
    click.echo(f"  Max workers: {settings.compute.max_workers}")
    click.echo(f"  Cache directory: {settings.cache_dir}")


@cli.command()
@click.pass_context
def info(ctx):
    """Show information about the analysis pipeline."""
    settings = ctx.obj['settings']
    
    click.echo(f"Plant MGC Analysis Pipeline")
    click.echo(f"Version: {settings.app_version}")
    click.echo(f"")
    click.echo(f"This unified pipeline consolidates 92 legacy scripts into a single,")
    click.echo(f"cohesive command-line interface for plant metabolic gene cluster analysis.")
    click.echo(f"")
    click.echo(f"Available analysis types:")
    for analysis_type in AnalysisType:
        click.echo(f"  - {analysis_type.value}")
    click.echo(f"")
    click.echo(f"Configuration:")
    click.echo(f"  Data directory: {settings.data_dir}")
    click.echo(f"  Output directory: {settings.output_dir}")
    click.echo(f"  Cache directory: {settings.cache_dir}")
    click.echo(f"  Environment: {'cluster' if settings.is_cluster else 'local'}")


if __name__ == '__main__':
    cli()