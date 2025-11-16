#!/bin/bash
#SBATCH --job-name=actual_mgc_cladepp_analysis_%j
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/out_actual_mgc_cladepp_analysis_%j.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/error_actual_mgc_cladepp_analysis_%j.ER
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=7-00:00:00

tree_path=$1
comparison_csv=$2
mgc_output_dir=$3
mapping_file=$4
compute_gain_loss_coevolution=$5
max_workers=$6
use_processes=$7

python /groups/itay_mayrose/alongonda/desktop/python_scripts/features/cladepp_phylo_profiling/cladepp_phylo_profiling_helpfuncs/analyze_tree_clades_dynamic.py --tree_path "$tree_path" --comparison_csv "$comparison_csv" --mgc_output_dir "$mgc_output_dir" --mapping_file "$mapping_file" --compute_gain_loss_coevolution "$compute_gain_loss_coevolution" --max_workers "$max_workers" --use_processes "$use_processes"