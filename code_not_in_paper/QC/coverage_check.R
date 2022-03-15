library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(exomeCopy)
library(Mus.musculus)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#mouse use all the CpG sites
chrs_mm10 <- names(Mmusculus)[1:19]#2276796
cgs_mm10 <- lapply(chrs_mm10, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr_mm10 <- do.call(c, lapply(1:19, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs_mm10[[x]], width = 2)))) #use both location
export.bed(cpgr_mm10,'../downstream/output/mouse_analysis/QC/mm10_all_CpG.bed')
#In bam file folder: 

#for fn in mm10*all.sort.dup.bam; do sbatch coverage_calc.sh mm10_all_CpG.sort.bed $fn coverage_file_all_CG/$fn.agnostic.cov; done


# Human coverage ----------------------------------------------------------

chrs <- names(Hsapiens)[1:25]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:25, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2)))) #use first location
seqlevels(cpgr)=gsub('chr','',seqlevels(cpgr))
export.bed(sort(cpgr),'../downstream/output/human_analysis/QC/hg19_all_CpG.bed')
#In coverage folder: /ibox/afeinbe2/yfang/human_analysis_bam_all/coverage
#for fn in ../bam_all/*phased.all.sorted.bam; do sbatch coverage_calc.sh hg19_all_CpG.bed $fn "${fn/\.\.\/bam_all\//}".agnostic.cov; done
#Need to rerun human coverage analysis
bed_out=GRanges()
cov_dir='../downstream/data/coverage_human/'
for (fn in dir(cov_dir,pattern='all.sorted.bam.agnostic.cov')){
  bed_in=import.bedGraph(fn)
  sample=gsub('_phased.all.sorted.bam.agnostic.cov','',fn)
  bed_in$coverage=bed_in$NA.2
  bed_in$Sample=sample
  bed_in=bed_in[seqlevels(bed_in)%in%1:22]
  mcols(bed_in)=mcols(bed_in)[,c('coverage','Sample')]
  bed_out=c(bed_out,bed_in)
  
}
#Mouse coverage
cov_dir='../downstream/data/coverage_mouse/'
bed_out=GRanges()
for (fn in dir(cov_dir,pattern='*.cov')){
  bed_in=import.bedGraph(paste0(cov_dir,fn))
  sample=gsub('_all.sort.dup.bam.agnostic.cov','',fn)
  bed_in$coverage=bed_in$NA.2
  bed_in$Sample=sample
  bed_in=bed_in[seqlevels(bed_in)%in%paste0('chr',1:19)]
  mcols(bed_in)=mcols(bed_in)[,c('coverage','Sample')]
  bed_out=c(bed_out,bed_in)
  
}
saveRDS(bed_out,'../downstream/output/mouse_analysis/mouse_coverage.rds')
#ASM coverage
#for fn in ../bam_asm/*.bam; do sbatch coverage_calc.sh hg19_all_CpG.bed $fn "${fn/\.\.\/bam_asm\//}".agnostic.cov; donebed_out=GRanges()
cov_dir='../downstream/data/coverage_human/'
bed_out=GRanges()
for (fn in dir(cov_dir,pattern='sort.genome')){
  bed_in=import.bedGraph(paste0(cov_dir,fn))
  sample=gsub('_phased.sort*.bam.agnostic.cov','',fn)
  cat("Processing: ",fn,"\n")
  bed_in$coverage=bed_in$NA.2
  bed_in$Sample=sample
  bed_in=bed_in[seqlevels(bed_in)%in%1:22]
  mcols(bed_in)=mcols(bed_in)[,c('coverage','Sample')]
  bed_out=c(bed_out,bed_in)
  
}
#Mbias
#for fn in ../bam_all/{149,150,HuFGM02,H1}*.bam; do echo sbatch mbias.sh $fn $PWD ~/data/yfang/referenceGenome/hg19_Arioc/; done