#!/bin/bash
#SBATCH -J motif_array
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --array=1-50
#SBATCH -o logfiles/motif_all-default-%x-%A-%a.out
#SBATCH -e logfiles/motif_all-default-%x-%A-%a.err

# Args
# Dependencies
ml R/3.6.1
ml gcc/5.5.0
ml atlas
# Command
Rscript --vanilla motif_break_array.R "$SLURM_ARRAY_TASK_ID" variant_in_all.rds 50 motif_all_JASPAR default 
