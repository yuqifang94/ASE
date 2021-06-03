#!/bin/bash
#SBATCH -J mbias
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH -o mbias_logfile/mbias-%A.out
#SBATCH -e mbias_logfile/mbias-%A.err
bam_in=$1
echo processing: $bam_in
bismark_methylation_extractor -o /ibox/afeinbe2/yfang/allele_agnostic_mouse_all/mbias --parallel 12 --genome_folder  /ibox/afeinbe2/yfang/allele_agnostic_mouse_all/CPEL/CPEL_agnostic/cpelasm/fasta --mbias_only -s $bam_in

