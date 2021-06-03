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



# Human coverage ----------------------------------------------------------

chrs <- names(Hsapiens)[1:24]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2)))) #use first location
export.bed(cpgr,'../downstream/output/human_analysis/QC/hg19_all_CpG.bed')
#Need to rerun human coverage analysis
bed_out=GRanges()
for (fn in dir(pattern='*.cov')){
  bed_in=import.bedGraph(fn)
  sample=gsub('_phased.all.sorted.bam.agnostic.cov','',fn)
  bed_in$coverage=bed_in$NA.2
  bed_in$Sample=sample
  bed_in=bed_in[seqlevels(bed_in)%in%1:22]
  mcols(bed_in)=mcols(bed_in)[,c('coverage','Sample')]
  bed_out=c(bed_out,bed_in)
  
  
}