"""
Main entry point for Plant MGC Analysis Pipeline.

This allows the package to be run as a module:
python -m plant_mgc_analysis
"""

from .cli_main import cli

if __name__ == '__main__':
    cli()