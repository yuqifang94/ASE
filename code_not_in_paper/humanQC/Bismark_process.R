library(bsseq)
library(rtracklayer)
library(Rsamtools)
library(rlist)
setwd("CpG_report")
#you need to separate positive and negative strain then calculate coverage
coverage_report<-function(file,strandCollapse=TRUE,nThread=1,subverbose=TRUE){
  dt <- bsseq:::.readBismarkAsDT(
    file = file,
    col_spec = "BSseq",
    check = TRUE,
    nThread = nThread,
    verbose = subverbose)
  
  if (strandCollapse && !is.null(dt[["strand"]])) {
    # Shift loci on negative strand by 1 to the left and then remove
    # strand since no longer valid.
    #select minus strain, those start -1 , make strand =NULL
    dt[strand == "-", start := start - 1L][, strand := NULL]
    # Aggregate counts at loci with the same 'seqnames' and 'start'. data table aggregate
    dt <- dt[, list(M = sum(M), U = sum(U)), by = c("seqnames", "start")]
  }
  #Note this does not remove ambigous CpG locations: forward strain != reverse strain -1
  dt$cov=dt$M+dt$U
  
  
  print(file)
  print(mean(dt$cov[dt$cov!=0]))
  
  return(dt[dt$cov!=0,])
  
}
fn=dir()#file_name
sample_name=NULL
tissue=NULL
read_type=NULL
coverage=list()
read_gff=function(gff_file){
  gff_table=readGFF(gff_file)
  #construct gff for each single CpG
  tt1=proc.time()
  print('separating CpGs')
  CpG_df=CpG_df(gff_table)
  tt2=proc.time()
  print(paste('CpG separation done in',tt2-tt1))
  
  tt1=proc.time()
  print('make Granges object')
  mclapply(CpG_df,makeGRangesFromDataFrame,mc.cores=24)
  tt2=proc.time()
  print(paste('make Granges object done in',tt2-tt1))
}
list.append2<-function(list_in,list_append){
  l=length(list_in)
  k=length(list_append)
  list_in[(l+1):(l+k)]=list_append
  list_in
}
CpG_df<-function(gff_df){
  #out=list()
  #for (i in 1:nrow(gff_df)){out=list.append2(out,CpG_split(gff_df[i,]))}
  #out
  out=mclapply(1:nrow(gff_df), function(x){CpG_split(gff_df[x,])},mc.cores=24)
  do.call('c',out)
}
CpG_split<-function(gff_single){
  chr=gff_single[1]
  CpG=gff_single$CpGs
  CpG=gsub('\\[','',CpG)
  CpG=as.numeric(unlist(gsub('\\]','',CpG)))
  out=list()
  out=lapply(CpG,function(x){data.frame(seqlevels=chr,start=x ,end=x+1)})
  return(out)
}

bam_gr<-function(bamfile){
  bamFile=BamFile(bamfile)
  aln=scanBam(bamFile)[[1]]
  aln_sub=aln[c('qname','rname','strand','pos')]
  aln_sub=as.data.frame(aln_sub)
  aln_sub$start=aln_sub$pos
  aln_sub$end=aln[['qwidth']]+aln_sub$start
  aln_sub$strand='*'
  colnames(aln_sub)=c('qname','seqnames','strand', 'pos','start','end')
  makeGRangesFromDataFrame(aln_sub)
}
cov_bam<-function(CpG_gff,bam){length(subsetByOverlaps(bam,CpG_gff))}
gff_out=read_gff('STL001_het.juliasm.gff')
bam_in=bam_gr('')
bam_apply=mclappy(gff_out,cov_bam,bam=bam_in,mc.cores=24)
#change spleen file name. stem 27 not ready
for (samples in fn){
  sn=gsub(".CpG_report.txt","",samples)
  sn_split=strsplit(sn,"_")[[1]]
  type_loc=grep("single|paired",sn_split)
  read_type=c(read_type,sn_split[type_loc])
  sample_name=c(sample_name,sn_split[1])
  tissue=c(tissue,gsub(", ","_",toString(sn_split[2:(type_loc-1)])))
  coverage[[sn]]=coverage_report(samples)
  gff_in=read_gff(gff_file)
  coverage[[sn]]$end=coverage[[sn]]$start+1
  coverage_df=makeGRangesFromDataFrame(coverage[[sn]],keep.extra.columns = T)
  coverage_het[[sn]]=subsetByOverlaps(coverage_df,gff_table)
}
# df_col=data.frame(sample_name=sample_name,tissue=tissue,read_type=read_type)
# rownames(df_col)=paste(sample_name,tissue,sep="_")
# data_input = read.bismark(file=file,
#                           colData =df_col,
#                           BPPARAM =MulticoreParam(workers=6,progressbar = TRUE),
#                           rmZeroCov = FALSE,
#                           verbose = TRUE,
#                           strandCollapse=TRUE,
#                           BACKEND="HDF5Array",
#                           dir="../bis_raw_allele",
#                           replace=TRUE)

bamfile='~/data/yfang/allele_specific/alignment/STL001_Bladder_single_finish/bam_output/STL001_Bladder_single_phased.sort.genome1.bam'

