#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q itaym
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N statistical_scanner
#PBS -e /groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#PBS -o /groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#PBS -l select=ncpus=32:mem=32gb

conda deactivate
conda activate plantismash

# Run your script
python /groups/itay_mayrose/alongonda/tools/plant_mgc_existing_tools/plantismash/run_antismash.py /groups/itay_mayrose/alongonda/datasets/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_rm.chromosome.1.fa