#!/bin/bash
#PBS -S /bin/bash
#PBS -q itaym
#PBS -l select=ncpus=32:mem=32gb
#PBS -e /groups/itay_mayrose/alongonda/desktop/example_jobs/error.ER
#PBS -o /groups/itay_mayrose/alongonda/desktop/example_jobs/out.OU
cd /groups/itay_mayrose/alongonda/desktop
source deactivate
source activate AlphaPulldown
python AlphaPulldown/alphapulldown/scripts/create_individual_features.py --fasta_paths="/groups/itay_mayrose/alongonda/desktop/Alphapulldown_TestRun/sequences.fasta" --data_dir=/groups/itay_mayrose/alongonda/alphafold_database_script/ --output_dir="/groups/itay_mayrose/alongonda/desktop/Alphapulldown_TestRun/output_dir/" --max_template_date=2050-01-01 --save_msa_files=False