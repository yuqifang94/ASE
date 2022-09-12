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
NME_only_name = "Predominantly\nNME"
MML_only_name = "Predominantly\nMML"
#588188 used for clustering
#
figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_all',gsub('.','',cutoffs),'.tiff')
cluster_assignment(dir_cluster_in_01,dir_out_cluster01,cutoffs=0.1,cluster_region_out_fn=cluster_01_region_out_fn,
                    figure_path=figure_path,figure_width=1800,figure_height=2000,res=200)
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
tissue_out_filtered=lapply(names(cor_dt_filtered),correlation_processing,cor_dt=cor_dt_filtered,filtered=T,NME_only_name=NME_only_name,MML_only_name=MML_only_name,
                           dir_figure=paste0(dir_out_rds_correlation,'correlation_figure/'))
names(tissue_out_filtered)=names(cor_dt_filtered)

saveRDS(tissue_out_filtered,tissue_out_filtered_fn)

# Plot the distribution of each category ----------------------------------
#all
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
region_type_count=lapply(tissue_out_filtered,function(x) table(x$region_type))
region_type_count=do.call(rbind,region_type_count)
apply(apply(region_type_count,1,function(x) x/sum(x)),1,range)
t.test(region_type_count[,'NME only'],region_type_count[,'MML only'],alternative='greater')
plot_correlation(tissue_out_filtered,pdf_fn=paste0(dir_out_rds_correlation,'correlation_main.pdf'))
#enhancer
bin_enhancer=readRDS(bin_enhancer_rds)
tissue_out_filtered_enhancer=lapply(tissue_out_filtered,function(x) x[queryHits(findOverlaps(convert_GR(x$region),bin_enhancer))])
plot_correlation(tissue_out_filtered_enhancer,pdf_fn=paste0(dir_out_rds_correlation,'correlation_main_enhancer.pdf'))
#Promoter
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
assign_regions(tissue_out_filtered,dir_out_cluster01,DNAase)

#Plotting non-tissue-specific heatmap
figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_non_ts_01.jpeg')
cluster_assignment(dir_cluster_in_01_non_ts,dir_out_cluster01_non_ts,cutoffs=0.1,
        cluster_region_out_fn=cluster_01_region_out_non_ts_fn,figure_name=figure_name)

# Plot example regions for each catogries -------------------------
UC_merge=readRDS(UC_merge_file)
UC_merge=UC_merge$heart
UC_merge=UC_merge[,grepl("-E1",colnames(UC_merge))]
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
tissue_out_filtered=tissue_out_filtered$heart
tissue_out_filtered=tissue_out_filtered[order(abs(cor_diff),decreasing=T)][region%in%rownames(UC_merge)]
 preProcessingRegion<-function(datIn){
    datOut=data.table(value=as.numeric(datIn),metaDat=names(datIn))
    datOut$statType=gsub("-.*","",datOut$metaDat)
    datOut$stage=gsub("dNME-|dMML-|UC-heart-|-all","",datOut$metaDat)
    datOutUCdNME=datOut[statType=="UC"]
    colnames(datOutUCdNME)[1]="UC"
    datOutUCdNME$value=datOut[statType=="dNME"][match(datOutUCdNME$stage,stage)]$value
    datOutUCdNME$statType="dNME"
    datOutUCdMML=datOut[statType=="UC"]
    colnames(datOutUCdMML)[1]="UC"
    datOutUCdMML$value=datOut[statType=="dMML"][match(datOutUCdMML$stage,stage)]$value
    datOutUCdMML$statType="dMML"
    return(rbind(datOutUCdNME,datOutUCdMML))
 }
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.text.x=element_text(size=10),
                                 axis.text.y=element_text(size=10),legend.title=element_blank())


 pdf("../downstream/output/mouse_analysis/QC/region_cat_example.pdf")
  NME_only_region=preProcessingRegion(UC_merge[tissue_out_filtered[region_type==NME_only_name][order(dNME_cor,decreasing=T)][3]$region,])
  NME_only_fig=ggplot(NME_only_region,aes(x=UC,y=value,color=statType))+geom_point(size=0.5)+geom_smooth(se=FALSE,method="lm")+ggtitle(NME_only_name)+theme_glob
  MML_only_region=preProcessingRegion(UC_merge[tissue_out_filtered[region_type==MML_only_name][order(round(dMML_cor,digits=3),decreasing=T)][3]$region,])
  MML_only_fig=ggplot(MML_only_region,aes(x=UC,y=value,color=statType))+geom_point(size=0.5)+geom_smooth(se=FALSE,method="lm")+ggtitle(MML_only_name)+theme_glob
  both_region=preProcessingRegion(UC_merge[tissue_out_filtered[region_type=="Both"][order(round(dMML_cor,digits=3),round(dNME_cor,digits=3),decreasing=T)][12]$region,])
  both_fig=ggplot(both_region,aes(x=UC,y=value,color=statType))+geom_point(size=0.5)+geom_smooth(se=FALSE,method="lm")+ggtitle("Both")+theme_glob
  neither_region=preProcessingRegion(UC_merge[tissue_out_filtered[region_type=="Neither"][order(round(dMML_cor,digits=3),round(dNME_cor,digits=3),decreasing=F)][3]$region,])
  neither_fig=ggplot(neither_region,aes(x=UC,y=value,color=statType))+geom_point(size=0.5)+geom_smooth(se=FALSE,method="lm")+ggtitle("Neither")+theme_glob
  print(ggarrange(NME_only_fig,MML_only_fig,both_fig,neither_fig,nrow=2,ncol=2,common.legend = TRUE, legend="bottom"))
 
 dev.off()
