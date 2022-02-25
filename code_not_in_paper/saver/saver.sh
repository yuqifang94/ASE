#!/bin/bash -l
#SBATCH --time=10:0:0
#SBATCH --mem=80G
#SBATCH --partition=lrgmem
ml R
Rscript saver.R $1
