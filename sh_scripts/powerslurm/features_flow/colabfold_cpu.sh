#!/bin/bash
#SBATCH --job-name=cpu_colabfold_${SLURM_ARRAY_TASK_ID}
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/example_jobs/cpu_out.OU
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/example_jobs/cpu_error.ER
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --time=24:00:00

# Get the input FASTA file from the command line argument
fasta_path=$1

# Define the output directory
pdb_out_dir=$2

echo "Starting my SLURM job"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Allocated CPUs: $SLURM_JOB_CPUS_PER_NODE"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
module purge
module load alphafold/alphafold_non_docker_2.3.1
echo "ALPHAFOLD_SCRIPT_PATH: $ALPHAFOLD_SCRIPT_PATH"
echo "ALPHAFOLD_DB_PATH: $ALPHAFOLD_DB_PATH"

FILE=$pdb_out_dir/ranked_0.pdb
if test -f "$FILE"; then
    echo "$FILE exists."
else 
    bash $ALPHAFOLD_SCRIPT_PATH/run_alphafold.sh -d $ALPHAFOLD_DB_PATH -o $pdb_out_dir -f $fasta_path -t $(date +%Y-%m-%d) -g false
fi