rm(list=ls())
source('mainFunctions_sub.R')
# Correlation matrix forming ----------------------------------------------
#Clustering assignment
dir_in_cor_Jason='../downstream/input/mouse_analysis/correlation_analysis/all_regions/'
dMML_cor=readRDS(dmml_cor_file)
dNME_cor=readRDS(dnme_cor_file)
dmml_perm=readRDS(paste0(dir_in_cor_Jason,'fullpermudmmlcor.rds'))
dnme_perm=readRDS(paste0(dir_in_cor_Jason,'fullpermudnmecor.rds'))
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.text.x=element_text(size=16),
                                 axis.text.y=element_text(size=16))


cluster_assignment(dir_cluster_in_01,dir_out_cluster01,cutoffs=0.1,cluster_region_out_fn=cluster_01_region_out_fn,figure_path=figure_path)
#Merge into data table
#Filtered result
cor_dt_pre_process_fn=paste0(dir_out_rds_correlation,'correlation_dt_N17_kmeans_10run_filtered_all_regions.rds')
if(!file.exists(cor_dt_pre_process_fn)){
  cor_dt_filtered=lapply(names(dNME_cor),cor_dt_preprocessing,dMML_cor=dMML_cor,
                         dNME_cor=dNME_cor,dmml_perm=dmml_perm,dnme_perm=dnme_perm,
                         filtered=TRUE,folder_input=dir_out_cluster01)
  names(cor_dt_filtered)=names(dNME_cor)
  saveRDS(cor_dt_filtered,cor_dt_pre_process_fn)
}
cor_dt_filtered=readRDS(cor_dt_pre_process_fn)
#Plot the density for each one
tissue_out_filtered=lapply(names(cor_dt_filtered),correlation_processing,cor_dt=cor_dt_filtered,filtered=T,
                           dir_figure=paste0(dir_out_rds_correlation,'/correlation_figure/'))
names(tissue_out_filtered)=names(cor_dt_filtered)
tissue_out_filtered_fn=paste0(dir_out_rds_correlation,'tissue_out_N17_kmeans_10run_filtered_all_region.rds')
saveRDS(tissue_out_filtered,tissue_out_filtered_fn)
lapply(tissue_out_filtered,nrow)
cluster_out=readRDS('../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_1/cluster_all_region_assignment_filtered_0_1.rds')
lapply(cluster_out,function(x) nrow(x[!grepl('chrX|chrY',regions)]))
tissue_out_filtered_enhancer=lapply(tissue_out_filtered,function(x) x[queryHits(findOverlaps(convert_GR(x$region),bin_enhancer))])
# Plot the distribution of each category ----------------------------------
#Promoter
plot_correlation(tissue_out_filtered,pdf_fn=paste0(dir_out_rds_correlation,'correlation_main.pdf'))

mm10_TSS=get_mm10_tss()
tissue_out_filtered_promoter=lapply(tissue_out_filtered,function(tissue_out_filtered_ts){
  olap=findOverlaps(convert_GR(tissue_out_filtered_ts$region,direction = 'GR'),mm10_TSS,maxgap = 2000)
  tissue_out_filtered_ts$promoter=FALSE
  tissue_out_filtered_ts$promoter[queryHits(olap)]=TRUE
  return(tissue_out_filtered_ts)
}
)
tissue_out_filtered_promoter=do.call(rbind,tissue_out_filtered_promoter)
plot_correlation(tissue_out_filtered_promoter,pdf_fn=paste0(dir_out_rds_correlation,'correlation_main_promoter.pdf'))



# Assign region type to each region in csv_folder -------------------------

#Assign DNase region
DNAase=readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp.rds')


assign_regions(tissue_out_filtered,dir_out_cluster,DNAase)
