#!/bin/bash
#SBATCH --job-name=alon_interactive
#SBATCH --partition=itaym
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=24:00:00  # Adjust as needed

srun --partition=itaym --ntasks=1 --cpus-per-task=32 --mem=32G --time=04:00:00 --pty bash

