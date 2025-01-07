#!/bin/bash

#PBS -S /bin/bash
#PBS -r y
#PBS -q itaym
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N test
#PBS -e /groups/itay_mayrose_nosnap/alongonda/desktop/example_jobs/error.ER
#PBS -o /groups/itay_mayrose_nosnap/alongonda/desktop/example_jobs/out.OU
#PBS -l select=ncpus=1:mem=1gb
cd /groups/itay_mayrose_nosnap/alongonda/desktop/example_jobs/
touch test.txt
