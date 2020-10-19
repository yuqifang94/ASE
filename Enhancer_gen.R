rm(list=ls())
library(data.table)
library(rtracklayer)
library(Gmisc)
#Read in FeDMR 
dir_in='../downstream/input/FeDMR/'
FeDMR_fn=dir(dir_in,pattern='tsv')
FeDMR_in=fastDoCall('c',lapply(dir(dir_in,pattern='tsv'),function(fn){
  sp=gsub('.tsv','',fn)
  sp=gsub("feDMR_","",sp)
  tissue=gsub(".*_","",sp)
  stage=gsub(paste0("_",tissue),"",sp)
  stage=gsub("_5",".5",stage)
  tissue=gsub("craniofacial","EFP",tissue)
  tissue=gsub("tube","NT",tissue)
  fn_in=fread(paste0(dir_in,fn))
  fn_in=makeGRangesFromDataFrame(fn_in,seqnames.field = "chrom",keep.extra.columns = T)
  fn_in$tissue=tissue
  fn_in$stage=stage
  fn_in$sample=paste(tissue,stage,sep='-')
  return(fn_in)
  
})
)
FeDMR_in_mcols=as.data.table(mcols(FeDMR_in))
FeDMR_in_mcols$tissue_exist=TRUE
FeDMR_in_gr=unique(FeDMR_in)
FeDMR_in_gr=FeDMR_in_gr[order(FeDMR_in_gr$dmr_id,decreasing=T)]
mcols(FeDMR_in_gr)=mcols(FeDMR_in_gr)[,c("dmr_id","score","tissue_specificity","stage_specificity")]
FeDMR_in_gr$dmr_id_original=FeDMR_in_gr$dmr_id
FeDMR_in_gr$dmr_id=NULL
FeDMR_in_mcols_dc=dcast.data.table(FeDMR_in_mcols,dmr_id~sample,value.var = "tissue_exist",fill=FALSE)
FeDMR_in_mcols_dc=FeDMR_in_mcols_dc[order(FeDMR_in_mcols_dc$dmr_id,decreasing = T)]
mcols(FeDMR_in_gr)=cbind(mcols(FeDMR_in_gr),FeDMR_in_mcols_dc)
print(which(FeDMR_in_gr$dmr_id_original!=FeDMR_in_gr$dmr_id))
FeDMR_in_gr$dmr_id_original=NULL
saveRDS(FeDMR_in_gr,'../downstream/output/FeDMR.rds')
# percent_cov=NULL
# UC_in=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')
# for(tissue in unique(FeDMR_in_mcols$tissue)){
# percent_cov=rbind(percent_cov,data.frame(tissue=tissue,
#   percent=length(subsetByOverlaps(FeDMR_in_gr[apply(mcols(FeDMR_in_gr)[,gsub('-.*','',colnames(mcols(FeDMR_in_gr)))==tissue],1,any)],UC_in[[tissue]]))/
#     length(FeDMR_in_gr[apply(mcols(FeDMR_in_gr)[,gsub('-.*','',colnames(mcols(FeDMR_in_gr)))==tissue],1,any)])
#   ))
# }

#Read in chromHMM

#ChromHMM
read_chromHMM_bed<-function(bed_dir,rep){
  bed_out=GRanges()
  for(fn in dir(bed_dir,pattern='.bed.gz')){
    #get sample name etc
    fn_split=strsplit(fn,'_')[[1]]
    stage=gsub('e','E',fn_split[1])
    tissue=gsub('facial-prominence','EFP',fn_split[2])
    tissue=gsub('neural-tube','NT',tissue)
    bed_in=read.table(paste(bed_dir,fn,sep=''))
    colnames(bed_in)=c('chr','start','end','chrom_num','chrom_state')
    bed_in=makeGRangesFromDataFrame(bed_in,keep.extra.columns = T)
    bed_in$stage=stage
    bed_in$tissue=tissue
    bed_in$rep=rep
    bed_in$Sample=paste(stage,tissue,sep='-')
    bed_out=c(bed_out,bed_in)
  }
  return(bed_out)
}
chromHMM_pooled= read_chromHMM_bed('../downstream/input/chromHMM_mm10/pooled/','pooled')
saveRDS(chromHMM_pooled,'../downstream/output/chromHMM.rds')
chromHMM_enhancer=chromHMM_pooled[sub('-.*','',as.character(chromHMM_pooled$chrom_state))=="En"]
saveRDS(chromHMM_enhancer,'../downstream/output/chromHMM_enhancer.rds')
