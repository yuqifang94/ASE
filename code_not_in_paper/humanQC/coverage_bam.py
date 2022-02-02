#!/usr/bin/env python3
import sys
import os
import argparse
import pathlib
import subprocess 
import shlex
def main():
	rf = "ls ~/work/shared/CpelAsm/data/{149,150,112,HuFGM02,H1}/bam/*{genome1,genome2}.bam"
	gff="het"
	bam_in = read_file(rf)
	for fl in bam_in:
		fn=fl.split("/")[-1]
		sp=fn.split("_")[0]
		po=fn.split(".")[-2]
		fn_red=fn.rstrip(".bam")
		cov="sub_cov/"+fn_red+gff+".sh"
		sub=open(cov,"w")
		sub.write(heads(fn_red,"/scratch/users/yfang27@jhu.edu/yfang/allele/alignment/R/gff/","10:00:00","24","shared,parallel,skylake","Coverage"))
		sub.write("module load intel/18.0\n")
		sub.write("samtools view -@24 -bhL gff/%s_%s.bed -o temp/%s_sub.bam %s \n" %(sp,gff,fn_red,fl))
		sub.write("Rscript --vanilla ~/data/yfang/bin/allele-specific/bin/bam_coverage_exe.R -s %s -f temp/%s_sub.bam -o out_%s -g %s" %(sp,fn_red,gff,gff))
		sub.close()
def read_file(cmd):

	check = subprocess.call(cmd,shell=True)
	if check==0:
		gzf = subprocess.check_output(cmd,shell=True)
	else:
		sys.exit("No samples found in fastq directory")
	#Convert files into the list
	out = gzf.decode("utf-8").split("\n")
	out=list(filter(None,out))
	return(out)
def heads(name,sdir,time,ncore,par,jn): #header of the submission script
	out = "#!/bin/sh\n"
	out = out+"#SBATCH --job-name=%s-%s\n" %(jn,name)
	out = out+"#SBATCH --cpus-per-task=%s\n" %ncore
	out = out+"#SBATCH --time=%s\n" %time
	out = out+"#SBATCH -p %s\n" %par
	out = out+"#SBATCH -o %ssub_cov/%s_%s_%%A.out\n" %(sdir,jn,name)
	out = out+"#SBATCH -e %ssub_cov/%s_%s_%%A.err\n" %(sdir,jn,name)
	out = out+"#SBATCH --chdir %s\n"%sdir
	#out=out+ "echo $PWD\n"
	return(out)
if __name__=="__main__":main()
