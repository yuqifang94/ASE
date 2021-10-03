#!/bin/bash
#$ -l cee,mf=3G,h_vmem=4G
#$ -pe local 1
#$ -o "../logfiles/bedGrpah-bw-$JOB_ID.err"
#$ -j y
#$ -m ebas
#$ -M "yfang27@jhmi.edu"
#$ -N bedGraph-bw
#$ -cwd
# Args
bedGraph_in="$1"
output_dir="$2"
intermediate_dir="$3"
chrom_size_file="$4"
fn=${bedGraph_in/*\//}
bw_out=${fn/\.bedGraph/.bw}

# Cut bedfile for first row (data)
echo cut -f1,2,3,4 $bedGraph_in > $intermediate_dir$fn.cut
echo bedGraphToBigWig $intermediate_dir$fn.cut $chrom_size_file $output_dir$bw_out

#Example submission
#qsub data_IO/bedGraph_to_bw.sh ../data_submission/compliment_MML_NME_model_mouse/mm10_heart_day10_5_all_allele_agnostic_mml.bedGraph ../data_submission/compliment_MML_NME_model_mouse/ ../data_submission/bedGraph_to_bw_intermediate/ ../data_submission/mm10.chrom.size