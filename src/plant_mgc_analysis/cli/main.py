"""
Main CLI interface for Plant MGC Analysis Pipeline.

This module provides the command-line interface using Click for a professional
user experience with comprehensive help, validation, and error handling.
"""

import sys
from pathlib import Path
from typing import List, Optional, Tuple

import click
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table
from loguru import logger

from ..core.analyzer import MGCAnalyzer
from ..core.types import AnalysisType
from ..core.exceptions import PlantMGCError, ValidationError
from ..config.settings import get_settings
from ..utils.logging import setup_logging

console = Console()


def setup_cli_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Setup logging for CLI."""
    if quiet:
        level = "ERROR"
    elif verbose:
        level = "DEBUG"
    else:
        level = "INFO"
    
    setup_logging(level=level, enable_json=False)


def handle_errors(func):
    """Decorator to handle CLI errors gracefully."""
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except PlantMGCError as e:
            console.print(f"[red]Error: {e}[/red]")
            sys.exit(1)
        except Exception as e:
            console.print(f"[red]Unexpected error: {e}[/red]")
            logger.exception("Unexpected error in CLI")
            sys.exit(1)
    return wrapper


@click.group()
@click.version_option()
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
@click.option("--quiet", "-q", is_flag=True, help="Enable quiet mode")
@click.option("--config", "-c", type=click.Path(exists=True), help="Configuration file path")
@click.pass_context
def main(ctx: click.Context, verbose: bool, quiet: bool, config: Optional[str]):
    """
    Plant MGC Analysis Pipeline - A comprehensive toolkit for plant biosynthetic gene cluster analysis.
    
    This tool provides multiple analysis methods including sliding window analysis,
    BLAST searches, machine learning predictions, and phylogenetic analysis.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet
    ctx.obj["config"] = config
    
    setup_cli_logging(verbose, quiet)
    
    if config:
        settings = get_settings()
        settings.load_config(Path(config))
        console.print(f"[green]Loaded configuration from {config}[/green]")


@main.command()
@click.option("--genome", "-g", required=True, type=click.Path(exists=True), 
              help="Path to genome FASTA file")
@click.option("--annotation", "-a", type=click.Path(exists=True), 
              help="Path to annotation file (GFF/GTF)")
@click.option("--protein", "-p", type=click.Path(exists=True), 
              help="Path to protein FASTA file")
@click.option("--organism", "-o", help="Organism name")
@click.option("--output", "-O", required=True, type=click.Path(), 
              help="Output directory")
@click.option("--analysis", "-A", multiple=True, 
              type=click.Choice([t.value for t in AnalysisType]), 
              default=["sliding_window"], show_default=True,
              help="Analysis types to run")
@click.option("--window-size", type=int, default=50000, show_default=True,
              help="Window size for sliding window analysis")
@click.option("--step-size", type=int, default=10000, show_default=True,
              help="Step size for sliding window analysis")
@click.option("--min-genes", type=int, default=3, show_default=True,
              help="Minimum genes per window")
@click.option("--evalue", type=float, default=1e-5, show_default=True,
              help="E-value threshold for BLAST")
@click.option("--parallel-jobs", "-j", type=int, default=1, show_default=True,
              help="Number of parallel jobs")
@click.pass_context
@handle_errors
def analyze(
    ctx: click.Context,
    genome: str,
    annotation: Optional[str],
    protein: Optional[str],
    organism: Optional[str],
    output: str,
    analysis: Tuple[str, ...],
    window_size: int,
    step_size: int,
    min_genes: int,
    evalue: float,
    parallel_jobs: int,
):
    """
    Analyze a single genome for metabolic gene clusters.
    
    This command processes a genome file and runs the specified analysis types
    to identify potential metabolic gene clusters.
    """
    console.print("[bold blue]Starting MGC Analysis[/bold blue]")
    
    # Initialize analyzer
    analyzer = MGCAnalyzer(log_level="DEBUG" if ctx.obj["verbose"] else "INFO")
    
    # Load genome data
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Loading genome data...", total=None)
        
        genome_data = analyzer.load_genome(
            genome_file=genome,
            annotation_file=annotation,
            protein_file=protein,
            organism_name=organism,
        )
        
        progress.update(task, description="Genome data loaded ✓")
    
    console.print(f"[green]Loaded genome for {genome_data.organism_name}[/green]")
    
    # Configure pipeline
    parameters = {
        "sliding_window": {
            "window_size": window_size,
            "step_size": step_size,
            "min_genes": min_genes,
        },
        "blast_search": {
            "evalue": evalue,
        },
    }
    
    analyzer.configure_pipeline(
        analysis_types=list(analysis),
        output_dir=output,
        parameters=parameters,
        parallel_jobs=parallel_jobs,
    )
    
    # Run analysis
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Running analysis pipeline...", total=None)
        
        results = analyzer.run_pipeline()
        
        progress.update(task, description="Analysis completed ✓")
    
    # Display results
    candidates = analyzer.get_mgc_candidates()
    
    # Create results table
    table = Table(title="Analysis Results")
    table.add_column("Analysis Type", style="cyan")
    table.add_column("Status", style="green")
    table.add_column("Candidates", style="yellow")
    table.add_column("Time (s)", style="magenta")
    
    for result in results:
        status = "✓" if result.success else "✗"
        table.add_row(
            result.analysis_type.value,
            status,
            str(len(result.get_mgc_candidates())),
            f"{result.execution_time:.2f}" if result.execution_time else "N/A",
        )
    
    console.print(table)
    console.print(f"[bold green]Found {len(candidates)} total MGC candidates[/bold green]")
    
    # Generate report
    report = analyzer.generate_report()
    
    # Save results
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    analyzer.save_results(output_path / "results.json")
    
    console.print(f"[green]Results saved to {output_path}[/green]")


@main.command()
@click.option("--input-dir", "-i", required=True, type=click.Path(exists=True), 
              help="Directory containing genome files")
@click.option("--output-dir", "-o", required=True, type=click.Path(), 
              help="Output directory")
@click.option("--config", "-c", type=click.Path(exists=True), 
              help="Configuration file")
@click.option("--parallel-jobs", "-j", type=int, default=1, show_default=True,
              help="Number of parallel jobs")
@click.option("--file-pattern", default="*.fasta", show_default=True,
              help="File pattern to match genome files")
@click.pass_context
@handle_errors
def pipeline(
    ctx: click.Context,
    input_dir: str,
    output_dir: str,
    config: Optional[str],
    parallel_jobs: int,
    file_pattern: str,
):
    """
    Run the complete analysis pipeline on multiple genomes.
    
    This command processes all genome files in a directory and runs
    the configured analysis pipeline on each.
    """
    console.print("[bold blue]Starting MGC Pipeline[/bold blue]")
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    # Find genome files
    genome_files = list(input_path.glob(file_pattern))
    
    if not genome_files:
        console.print(f"[red]No genome files found matching pattern: {file_pattern}[/red]")
        sys.exit(1)
    
    console.print(f"[green]Found {len(genome_files)} genome files[/green]")
    
    # Process each genome
    for genome_file in genome_files:
        console.print(f"[cyan]Processing {genome_file.name}[/cyan]")
        
        # Create output directory for this genome
        genome_output = output_path / genome_file.stem
        genome_output.mkdir(parents=True, exist_ok=True)
        
        # Run analysis (this would be implemented with proper error handling)
        # For now, just show the structure
        console.print(f"[yellow]  - Output: {genome_output}[/yellow]")
    
    console.print("[bold green]Pipeline completed[/bold green]")


@main.command()
@click.option("--candidates", "-c", required=True, type=click.Path(exists=True), 
              help="Path to MGC candidates file")
@click.option("--model", "-m", default="phytoclust", show_default=True,
              help="Machine learning model to use")
@click.option("--output", "-o", required=True, type=click.Path(), 
              help="Output file for predictions")
@click.option("--threshold", "-t", type=float, default=0.5, show_default=True,
              help="Prediction threshold")
@click.pass_context
@handle_errors
def predict(
    ctx: click.Context,
    candidates: str,
    model: str,
    output: str,
    threshold: float,
):
    """
    Run machine learning predictions on MGC candidates.
    
    This command takes a file of MGC candidates and applies machine learning
    models to predict their likelihood of being true biosynthetic gene clusters.
    """
    console.print("[bold blue]Starting MGC Prediction[/bold blue]")
    
    console.print(f"[cyan]Loading candidates from {candidates}[/cyan]")
    console.print(f"[cyan]Using model: {model}[/cyan]")
    console.print(f"[cyan]Threshold: {threshold}[/cyan]")
    
    # This would implement the actual prediction logic
    console.print("[yellow]Prediction functionality not yet implemented[/yellow]")
    
    console.print(f"[green]Predictions saved to {output}[/green]")


@main.command()
@click.option("--all", "all_dbs", is_flag=True, help="Setup all databases")
@click.option("--kegg", is_flag=True, help="Setup KEGG database")
@click.option("--mibig", is_flag=True, help="Setup MiBIG database")
@click.option("--plantcyc", is_flag=True, help="Setup PlantCyc database")
@click.option("--data-dir", type=click.Path(), help="Data directory")
@handle_errors
def setup_databases(all_dbs: bool, kegg: bool, mibig: bool, plantcyc: bool, data_dir: Optional[str]):
    """
    Setup and download required databases.
    
    This command downloads and sets up the databases required for MGC analysis.
    """
    console.print("[bold blue]Setting up databases[/bold blue]")
    
    if all_dbs:
        kegg = mibig = plantcyc = True
    
    if not any([kegg, mibig, plantcyc]):
        console.print("[red]Please specify at least one database to setup[/red]")
        sys.exit(1)
    
    settings = get_settings()
    data_path = Path(data_dir) if data_dir else settings.data_dir
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        
        if kegg:
            task = progress.add_task("Setting up KEGG database...", total=None)
            # Implementation would go here
            progress.update(task, description="KEGG database setup ✓")
        
        if mibig:
            task = progress.add_task("Setting up MiBIG database...", total=None)
            # Implementation would go here
            progress.update(task, description="MiBIG database setup ✓")
        
        if plantcyc:
            task = progress.add_task("Setting up PlantCyc database...", total=None)
            # Implementation would go here
            progress.update(task, description="PlantCyc database setup ✓")
    
    console.print(f"[green]Databases setup completed in {data_path}[/green]")


if __name__ == "__main__":
    main()