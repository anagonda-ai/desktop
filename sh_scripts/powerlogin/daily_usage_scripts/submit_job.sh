#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q itaym
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N statistical_scanner
#PBS -e /groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#PBS -o /groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#PBS -l select=ncpus=32:mem=32gb

# Run your script
python "/groups/itay_mayrose/alongonda/desktop/python_scripts/move_files.py"