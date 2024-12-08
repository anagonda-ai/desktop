#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q itaym
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N test
#PBS -e /groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#PBS -o /groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#PBS -l select=ncpus=8:mem=8gb
/groups/itay_mayrose/alongonda/desktop/alphafold/scripts/download_bfd.sh /groups/itay_mayrose_nosnap/alongonda/alphafold_database_script/