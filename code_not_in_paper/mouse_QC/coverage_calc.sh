#!/bin/bash
#SBATCH -J coverage_calc
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH -o ../../downstream/logfiles/coverage-ChIP-%A.out
#SBATCH -e ../../downstream/logfiles/coverage-ChIP-%A.err
#bedtools v2.29.2
bed_in=$1
bam_in=$2
cov_out=$3

#bedtools coverage -sorted -counts -a ../../downstream/input/mouse_analysis/enhancer_selection/bin_enhancer.bed -b ../../downstream/data/mouse_ChIP/bam_files/E16_5_midbrain_H3K27ac_2.bam  > ../../downstream/data/mouse_ChIP/bam_files/E16_5_midbrain_H3K27ac_2.bam.cov
bedtools coverage -sorted -counts -a ${bed_in} -b ${bam_in} > ${cov_out}
