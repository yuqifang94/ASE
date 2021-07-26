#!/bin/bash
#SBATCH -J coverage_calc
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH -o ../../downstream/output/mouse_analysis/data_collection_log/coverage-%A.out
#SBATCH -e ../../downstream/output/mouse_analysis/data_collection_log/coverage-%A.err
#bedtools v2.29.2
bed_in=$1
bam_in=$2
cov_out=$3


bedtools coverage -f 1 -sorted -counts -a ${bed_in} -b ${bam_in} > ${cov_out}
