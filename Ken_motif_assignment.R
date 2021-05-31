rm(list=ls())
library(data.table)
library(dplyr)
library(glue)
library(GenomicRanges)
library(parallel)
convert_GR<-function(x,direction="GR"){
  if(direction=="GR"){
    strand=unlist(lapply(strsplit(x,','),function(x)x[2]))
    strand=ifelse(is.na(strand),"*",strand)
    x=ifelse(is.na(strand),x,sub(paste0(',\\','-'),'',sub(paste0(',\\','+'),'',x)))
    
    gr=GRanges(seqnames=sub(':.*','',x),
               IRanges(start=as.numeric(sub('-.*','',sub('.*:','',x))),
                       end=as.numeric(sub('.*-','',x))),strand=strand)
    return(gr)}else
      if(direction=="DT"){
        
        x_dt=as.data.table(mcols(x))
        x_dt$region=paste0(seqnames(x),':',start(x),'-',end(x))
        return(x_dt)
      }else
        if(direction=="matrix"){
          x_mt=as.matrix(mcols(x))
          rownames(x_mt)=paste0(seqnames(x),':',start(x),'-',end(x))
          return(x_mt)
        }
  
}

# motif_lib_mm10 <- readRDS("../downstream/input/mouse_analysis/motif_analysis/motif_JASPAR_mm10.rds")
# 
# cluster_gr_anno <- readRDS(file="../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_all_region.rds")
motif_lib_mm10 <- readRDS("../downstream/input/mouse_analysis/motif_analysis/motif_JASPAR_mm10.rds")

cluster_gr_anno <- readRDS(file="../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_all_region.rds")

cluster_gr_anno=lapply(cluster_gr_anno,function(x) {
  
  
  x_gr=convert_GR(x$region)
  mcols(x_gr)=x[,list(tissue,region_type,cluster,enhancer)]
  return(x_gr)
  })

temp <- lapply(c(1:length(cluster_gr_anno)), function(i) {
  
  select_motif_site <- mclapply(c(1:length(motif_lib_mm10)), function(m){
    cluster_data_select <- cluster_gr_anno[[i]]
    overlap_region <- findOverlaps(motif_lib_mm10[[m]], cluster_data_select)
    overlap_gene <- motif_lib_mm10[[m]][queryHits(overlap_region)]
    temp_region <- as.data.frame(cluster_data_select[subjectHits(overlap_region)])
    overlap_gene$region <- glue("{temp_region[,1]}:{temp_region[,2]}-{temp_region[,3]}")
    return(overlap_gene)
  },mc.cores=20)
  names(select_motif_site) <- names(motif_lib_mm10)
  return(select_motif_site)
  
})
names(temp) <- names(cluster_gr_anno)
saveRDS(temp, file="../downstream/output/mouse_analysis/motif_analysis/tissue_region_motif_all_regions.rds")

