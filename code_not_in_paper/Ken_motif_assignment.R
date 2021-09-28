rm(list=ls())
library(data.table)
library(dplyr)
library(glue)
library(GenomicRanges)
library(parallel)
source('mainFunctions_sub.R')
motif_lib_mm10 <- readRDS(JASPAR_motif_mm10_file)

cluster_gr_anno <- readRDS(file=tissue_out_filtered_fn)

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
saveRDS(temp, file=tissue_region_motif_all_regions_fn)

