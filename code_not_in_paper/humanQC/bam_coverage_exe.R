#!/usr/bin/env Rscript
library(bsseq)
library(rtracklayer)
library(Rsamtools)
library(rlist)
library("optparse")
bam_gr<-function(aln){
  aln_sub=aln[c('qname','rname','strand','pos')]
  aln_sub=as.data.frame(aln_sub)
  aln_sub$start=aln_sub$pos
  aln_sub$end=aln[['qwidth']]+aln_sub$start
  aln_sub$strand='*'
  colnames(aln_sub)=c('qname','seqnames','strand', 'pos','start','end')
  makeGRangesFromDataFrame(aln_sub)
}
cov_bam<-function(CpG_gff,bam){countOverlaps(CpG_gff,bam)}

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="Patient name", metavar="character"),
  make_option(c("-f", "--input"), type="character", default=NULL, 
              help="input file name", metavar="character"),
  make_option(c("-o","--output"),type="character",default=NULL,
              metavar="character",help="Output directory"),
  make_option(c("-g","--het"),type="character",default=NULL,
              metavar="character",help="gff file type")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
fn=unlist(strsplit(opt$input,"/"))
fn=fn[length(fn)]
cat(paste(fn,'\n'))
tt1=proc.time()[[3]]
p_tag=ScanBamParam(what=c('qname',"rname", "strand", "pos", "qwidth"))
bam_in=scanBam(BamFile(opt$input),param=p_tag)
bam_in=bam_gr(bam_in[[1]])
tt2=proc.time()[[3]]
cat(paste('Finish reading in bam file',opt$file,'in',tt2-tt1,'\n'))
bam_all=c()
bam_apply=list()
cat('generating all coverage list\n')
t1=proc.time()[[3]]
gff_out=readRDS(paste('gff/',opt$sample,"_",opt$het,"_gff_all.rds",sep=""))
bam_all=countOverlaps(gff_out,bam_in)
saveRDS(bam_all,paste(opt$output,gsub('.bam',paste(opt$het,'all.rds',sep='.'),fn),sep='/'))
t2=proc.time()[[3]]
cat(paste('time for generating all coverage list =',t2-t1,'\n'))
#for (chr in 1:22){
#  cat(paste('reading in gff file',opt$sample,'\n'))
#  gff_out=readRDS(paste('gff/chr/',opt$sample,'_het_gff','_chr',chr,'.rds',sep=''))
#  cat(paste('generating haplotype coverage list for chr',chr,'\n'))
#  bam_apply=c(bam_apply,mclapply(gff_out,cov_bam,bam=bam_in,mc.cores=8))
  
  
#}
#cat(paste('saving result',paste(opt$output,gsub('.bam','.cov.rds',fn),sep='/')))
#saveRDS(bam_apply,paste(opt$output,gsub('.bam','.cov.rds',fn),sep='/'))


##read in rds file
# ASM_size=list()
# for (fn in dir()){
#   cov=readRDS(fn)
#   fn_sp=strsplit(fn,'_')[[1]]
#   fn_spg=strsplit(fn_sp[length(fn_sp)],'\\.')[[1]][3]
#   fnc=paste(paste(fn_sp[-c(length(fn_sp),length(fn_sp)-1)],collapse = '_'),fn_spg,sep="_")
#   cat(paste(fn,':\n'),sep="")
#   cat(paste(mean(cov[cov!=0]),'\n',sep=''))
#   ASM_size[[fnc]]=c(fnc,sum(cov>=6))
#   
# }
# do.call(rbind,ASM_size)
