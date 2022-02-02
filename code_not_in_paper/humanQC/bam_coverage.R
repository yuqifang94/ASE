library(bsseq)
library(rtracklayer)
library(Rsamtools)
library(rlist)
read_gff=function(gff_file){
  gff_table=readGFF(gff_file)
  #construct gff for each single CpG
  print('separating CpG')
  CpG_out=CpG_df(gff_table)
  print(CpG_out[1:10])
  print('make Granges object')
  return(mclapply(CpG_out,makeGRangesFromDataFrame,mc.cores=24))
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
  #do.call('c',out)
  return(out)
}
CpG_split<-function(gff_single){
  chr=gff_single[1]
  CpG=gff_single$CpGs
  CpG=gsub('\\[','',CpG)
  CpG=as.numeric(unlist(gsub('\\]','',CpG)))
  out=list()
  out=lapply(CpG,function(x){data.frame(seqlevels=chr,start=x ,end=x+1)})
  return(do.call(rbind,out))
}

gff_list=c('STL001','STL002','STL003','STL011','HUES64','H9','skin03','HuFGM02')
for (gff_file in gff_list){
  #gff_out=read_gff(paste(gff_file,'_het.juliasm.gff',sep=''))
  gff_out=readRDS(paste('gff/',gff_file,'_het_gff.rds',sep=''))
  chr_sum=0
  for(chr in 1:22){
    t1=proc.time()[3]
    chr_loc=mclapply(gff_out,function(x,chr){all(seqnames(x)==chr)},chr=chr,mc.cores=8)
    cat(paste('chromosome sub finish at',proc.time()[3]-t1,sep=''))
    chr_out=gff_out[which(do.call('c',chr_loc))]
    chr_sum=length(chr_out)+chr_sum
    saveRDS(chr_out,paste('gff/',gff_file,'_chr',chr,'_het_gff.rds',sep=''))
  }
}

