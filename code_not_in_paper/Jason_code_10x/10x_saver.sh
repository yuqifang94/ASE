#!/bin/bash -l
#SBATCH --time=60:0:0
#SBATCH --mem=200G
#SBATCH --partition=lrgmem
ml R/4.0.2
Rscript 10x_saver.R $1
