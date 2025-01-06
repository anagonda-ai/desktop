#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q itaym
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N filter_by_similarity_plaza
#PBS -e /groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#PBS -o /groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
#PBS -l select=ncpus=32:mem=32gb
python /groups/itay_mayrose/alongonda/desktop/python_scripts/pmn_scripts/mgc_candidate_scripts/filter_by_similarity.py --targets_dir /groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/organisms --filter_by /groups/itay_mayrose_nosnap/alongonda/full_genomes/chloroplast_prot_genes --output /groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza_without_chloroplast