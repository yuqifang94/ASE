rm(list=ls())
source('mainFunctions_sub.R')
#Get all regions we analyzed UC in at least one samples
all_regions=readRDS(UC_in_matrix_ls_file)
regions_analyzed=lapply(all_regions,granges)
names(regions_analyzed)=NULL
regions_analyzed=unique(do.call('c',regions_analyzed))
regions_analyzed_all=regions_analyzed
cluster_01_in=readRDS(paste0(dir_cluster_in_01,'uc_0.1_1.rds'))
for(ts in names(cluster_01_in)){
  mcols(regions_analyzed_all)[,ts]=0
  mcols(regions_analyzed_all)[convert_GR(regions_analyzed_all,direction='DT')$region %in% names(cluster_01_in[[ts]]),ts]=1
  
  
}
saveRDS(regions_analyzed_all,paste0(GWAS_prep_dir,'region_analyzed_01_select_all.rds'))

#prepare for NME driven ones
region_catogrization=readRDS(tissue_out_filtered_fn)
regions_analyzed_NME=regions_analyzed
mcols(regions_analyzed_NME)=NULL
for(ts in names(region_catogrization)){
  mcols(regions_analyzed_NME)[,ts]=0
  mcols(regions_analyzed_NME)[convert_GR(regions_analyzed_NME,direction='DT')$region %in% 
                              region_catogrization[[ts]][region_type%in% c('NME only','Both')]$region,ts]=1
  
  
}
saveRDS(regions_analyzed_NME,paste0(GWAS_prep_dir,'region_analyzed_01_select_NME.rds'))

#prepare for MML driven ones

regions_analyzed_MML=regions_analyzed
mcols(regions_analyzed_MML)=NULL
for(ts in names(region_catogrization)){
  mcols(regions_analyzed_MML)[,ts]=0
  mcols(regions_analyzed_MML)[convert_GR(regions_analyzed_MML,direction='DT')$region %in% 
                              region_catogrization[[ts]][region_type%in% c('MML only','Both')]$region,ts]=1
  
  
}
saveRDS(regions_analyzed_MML,paste0(GWAS_prep_dir,'region_analyzed_01_select_MML.rds'))