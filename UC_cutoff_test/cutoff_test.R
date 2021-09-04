# nolint start
source('mainFunctions_sub.R')
#We have already done the clustering for each cutoff, we need to run cluster assignment and GO and correlation analysis# nolint
#Input folder 
cutoff_char=as.character(gsub('\\.','',cutoff))
cluster_dir_in=paste0('../downstream/input/mouse_analysis/clustering/tissue_specific/uc_',cutoff_char,'/')
#creating output folder
output_dir=paste0('../downstream/output/mouse_analysis/cutoff_testing/',cutoff_char,'/')
if(!dir.exists(output_dir)){dir.create(output_dir)}
cluster_assigned_dir=paste0(output_dir,'cluster_assigned/')
if(!dir.exists(cluster_assigned_dir)){dir.create(cluster_assigned_dir)}
  UC_merge=readRDS(UC_merge_max_loc_file)
#Clustering

figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_all',cutoff_char,'.tiff')
cluster_assignment(cluster_dir_in,cluster_assigned_dir,cutoffs=cutoff,
            cluster_region_out_fn=paste0(cluster_assigned_dir,'cluster_assginment_filtered_',cutoff_char,'.rds'),figure_name=figure_name,UC_merge)

#Reading in all correlation
dMML_cor=readRDS(dmml_cor_file)
dNME_cor=readRDS(dnme_cor_file)
dmml_perm=readRDS(paste0(dir_in_cor_Jason,'fullpermudmmlcor.rds'))
dnme_perm=readRDS(paste0(dir_in_cor_Jason,'fullpermudnmecor.rds'))
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"), # nolint
                                 axis.text.x=element_text(size=16),
                                 axis.text.y=element_text(size=16))

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
                           dir_figure=paste0(dir_out_rds_correlation,'correlation_figure/'))
names(tissue_out_filtered)=names(cor_dt_filtered)

saveRDS(tissue_out_filtered,tissue_out_filtered_fn)



# nolint send