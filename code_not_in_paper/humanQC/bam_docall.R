library(bsseq)
library(rtracklayer)
library(Rsamtools)
library(rlist)
library("optparse")
read_gff=function(gff_file,ft){
  gff_table=read.delim2(gff_file,header=F)[,c(1,4,5,6,9)]
  colnames(gff_table)=c('chr','start','end','number of CpG','CpGs')
   gff_table$CpGs=gsub('CpGs=','',gff_table$CpGs)
   if(ft=='het'){gff_table$CpGs=do.call('c',mclapply(strsplit(gff_table$CpGs,";"),function(x){x[2]},mc.cores=12))}
  #construct gff for each single CpG
  tt1=proc.time()[[3]]
  cat('separating CpGs\n')
  CpG_out=CpG_df(gff_table)
  CpG_out=lapply(CpG_out,makeGRangesFromDataFrame)
  CpG_out=do.call('c',CpG_out)
  cat('Generating bed files')
  #export(CpG_out,paste('gff/',gsub('.cpel.gff','.bed',gff_file),sep=''),format="BED")
  tt2=proc.time()[[3]]
  cat(paste('CpG separation done in',tt2-tt1,'\n'))
  gc()
  tt1=proc.time()[[3]]
  cat('make Granges object')
  #print(CpG_out)
  #out=makeGRangesFromDataFrame(CpG_out,keep.extra.columns=TRUE)
  print(paste('gff/',gsub('.cpel.gff','.bed',gff_file),sep=''))
  export(CpG_out,paste('gff/',gsub('.cpel.gff','.bed',gff_file),sep=''),format="BED")
  tt2=proc.time()[[3]]
  tt2=proc.time()[[3]]
  cat(paste('make Granges object done in',tt2-tt1,'\n'))
  return(CpG_out)
  #saveRDS(out,'out.rds')
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
  out=mclapply(1:nrow(gff_df), function(x){CpG_split(gff_df[x,])},mc.cores=20)
  out=do.call('c',out)
  out

}
CpG_split<-function(gff_single){
  chr=gff_single[1]
  CpG=gff_single$CpGs
  CpG=gsub('\\[','',CpG)
  CpG=as.numeric(strsplit((unlist(gsub('\\]','',CpG))),',')[[1]])
  out=list()
  out=lapply(CpG,function(x){data.frame(seqlevels=chr,start=x ,end=x+1)})
  return(out)
}
option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="Patient name", metavar="character"),
  make_option(c("-t","--hetgff"),type="character",default=NULL,help="gff type",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#print(opt)
cat('reading gff file\n')
gff_out=read_gff(opt$sample,ft=opt$hetgff)
#gff_out=readRDS(opt$sample)
cat('start combine regions \n')
t1=proc.time()[3]
gc()
print(gsub('.cpel.gff','_gff_all.rds',opt$sample))
saveRDS(gff_out,gsub('.cpel.gff','_gff_all.rds',opt$sample))
cat(paste('combining done in ',proc.time()[3]-t1,'\n',sep=''))
# chr_sum=0
# for(chr in 1:22){
#   t1=proc.time()[3]
#   chr_loc=mclapply(gff_out,function(x,chr){all(seqnames(x)==chr)},chr=chr,mc.cores=5)
#   cat(paste('chromosome sub finish in ',proc.time()[3]-t1,'\n',sep=''))
#   chr_out=gff_out[which(do.call('c',chr_loc))]
#   chr_sum=length(chr_out)+chr_sum
#   chr_out_name=gsub('het_gff',paste('het_gff_chr',chr,sep=''),gsub('gff/','',opt$sample))
#   cat(paste('saving output to ','gff/chr/',chr_out_name,'\n',sep=''))
#   saveRDS(chr_out,paste('gff/chr/',chr_out_name))
# }
# print(chr_sum)
# print(length(gff_out))
