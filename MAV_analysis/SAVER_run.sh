#!/bin/bash -l
#SBATCH --time=60:0:0
#SBATCH --mem=200G
#SBATCH --partition=lrgmem
#SBATCH -o logfile/SAVER-%A.out
#SBATCH -e logfile/SAVER-%A.err
ml R/4.0.2
echo $1
Rscript 10x_saver.R $1

#for i in `seq 1 7`
# do
# qsub 10x_saver.sh $i
# done