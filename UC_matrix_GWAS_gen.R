rm(list=ls())
source('mainFunctions_sub.R')
#Get all regions we analyzed UC in at least one samples
all_regions=readRDS(UC_in_matrix_ls_file)
regions_analyzed=lapply(all_regions,granges)
names(regions_analyzed)=NULL
regions_analyzed=unique(do.call('c',regions_analyzed))
cluster_01_in=readRDS(paste0(dir_cluster_in_01,'uc_0.1_1.rds'))
for(ts in names(cluster_01_in)){
  mcols(regions_analyzed)[,ts]=0
  mcols(regions_analyzed)[convert_GR(regions_analyzed,direction='DT')$region %in% names(cluster_01_in[[ts]]),ts]=1
  
  
}
saveRDS(regions_analyzed,'../downstream/output/mouse_analysis/GWAS_prep/region_analyzed_01_select_all.rds')