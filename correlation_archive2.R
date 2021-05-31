rm(list=ls())
source('mainFunctions_sub.R')
# Correlation matrix forming ----------------------------------------------
#Getting fitted data and correlation outside the fitted data
dir_in='../downstream/input/mouse_analysis/clustering/tissue_specific/currently_in_use/ts_cluster_0_1/'
dMML_cor=readRDS('../downstream/input/mouse_analysis/correlation_analysis/all_regions/fulldmmlcor.rds')
dNME_cor=readRDS('../downstream/input/mouse_analysis/correlation_analysis/all_regions/fulldnmecor.rds')
dmml_perm=readRDS('../downstream/input/mouse_analysis/correlation_analysis/all_regions/fullpermudmmlcor.rds')
dnme_perm=readRDS('../downstream/input/mouse_analysis/correlation_analysis/all_regions/fullpermudnmecor.rds')
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.text.x=element_text(size=16),
                                 axis.text.y=element_text(size=16))
#Merge into data table
#Filtered result
cor_dt_filtered=lapply(names(dNME_cor),cor_dt_preprocessing,dMML_cor=dMML_cor,
                       dNME_cor=dNME_cor,dmml_perm=dmml_perm,dnme_perm=dnme_perm,
                       filtered=TRUE,folder_input=dir_in)
names(cor_dt_filtered)=names(dNME_cor)
saveRDS(cor_dt_filtered,'../downstream/output/mouse_analysis/correlation/correlation_dt_N17_kmeans_10run_filtered_all_regions.rds')
#Plot the density for each one
cor_dt_filtered=readRDS('../downstream/output/mouse_analysis/correlation/correlation_dt_N17_kmeans_10run_filtered_all_regions.rds')
tissue_out_filtered=lapply(names(cor_dt_filtered),correlation_processing,cor_dt=cor_dt_filtered,filtered=T,
                           dir_in="../downstream/output/mouse_analysis/correlation/kmeans_10_run_filtered_all_regions/")
names(tissue_out_filtered)=names(cor_dt_filtered)
tissue_out_filtered_fn='../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_all_region_V2.rds'
saveRDS(tissue_out_filtered,tissue_out_filtered_fn)
tissue_out_filtered_comp=lapply(names(tissue_out_filtered_old),function(x) {
  tissue_out_filtered[[x]]$region_type_old=tissue_out_filtered_old[[x]][match(tissue_out_filtered[[x]]$region,region)]$region_type
  return(tissue_out_filtered[[x]])
  })
lapply(tissue_out_filtered_comp,function(x) sum(x$region_type==x$region_type_old)/nrow(x))

lapply(tissue_out_filtered,function(x) {
  olap=findOverlaps(convert_GR(x$region),enhancer)
  return(unique(x[queryHits(olap)]$enhancer))
  
})

#Plot the density for each one remove repeat overlap
dir_in='../downstream/input/mouse_analysis/clustering/tissue_specific/currently_in_use/ts_cluster_0_1_non_repeats/'
cor_dt_filtered=lapply(names(dNME_cor),cor_dt_preprocessing,dMML_cor=dMML_cor,
                       dNME_cor=dNME_cor,dmml_perm=dmml_perm,dnme_perm=dnme_perm,
                       filtered=TRUE,folder_input=dir_in)
names(cor_dt_filtered)=names(dNME_cor)
saveRDS(cor_dt_filtered,'../downstream/output/mouse_analysis/correlation/correlation_dt_N17_kmeans_10run_filtered_all_regions_non_repeats.rds')
cor_dt_filtered=readRDS('../downstream/output/mouse_analysis/correlation/correlation_dt_N17_kmeans_10run_filtered_all_regions_non_repeats.rds')

tissue_out_filtered=lapply(names(cor_dt_filtered),correlation_processing,cor_dt=cor_dt_filtered,filtered=T,
                           dir_in="../downstream/output/mouse_analysis/correlation/kmeans_10_run_filtered_all_regions_non_repeats/")
names(tissue_out_filtered)=names(cor_dt_filtered)
tissue_out_filtered_fn='../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_all_region_non_repeats.rds'
saveRDS(tissue_out_filtered,tissue_out_filtered_fn)


# Plot the distribution of each category ----------------------------------
#Overall
tissue_out_filtered_fn='../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_all_region.rds'
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
mm10_TSS=get_mm10_tss()
tissue_out_filtered=lapply(tissue_out_filtered,function(tissue_out_filtered_ts){
  olap=findOverlaps(convert_GR(tissue_out_filtered_ts$region,direction = 'GR'),mm10_TSS,maxgap = 2000)
  tissue_out_filtered_ts$promoter=FALSE
  tissue_out_filtered_ts$promoter[queryHits(olap)]=TRUE
  return(tissue_out_filtered_ts)
}
)
tissue_out_filtered=do.call(rbind,tissue_out_filtered)
#Change name
tissue_out_filtered[region_type=="dMML_only"]$region_type="dMML only"
tissue_out_filtered[region_type=="dNME_only"]$region_type="dNME only"
tissue_out_filtered[region_type=="Both"]$region_type="Both"
tissue_out_filtered[region_type=="Neither"]$region_type="Neither"

tissue_out_filtered$region_type=factor(tissue_out_filtered$region_type,
                                                 levels=c("dMML only","dNME only","Both","Neither"))
#All regions
tissue_out_filtered_frequency=tissue_out_filtered[tissue!="NT",list(count=length(region)),by=list(tissue,region_type)]
tissue_out_filtered_frequency_p=tissue_out_filtered_frequency[,list(percentage=round(count/sum(count)*100,digits=0),region_type=region_type),by=list(tissue)]
tissue_out_filtered_frequency_p[,list(minp=min(percentage),maxp=max(percentage)),by=list(region_type)]
t.test(tissue_out_filtered_frequency_p[region_type=="dMML only"]$percentage,
       tissue_out_filtered_frequency_p[region_type=="dNME only"]$percentage,alternative = "less")


#Enhancer only
tissue_out_filtered_frequency_enc=tissue_out_filtered[tissue!="NT"&enhancer==T,list(count=length(region)),by=list(tissue,region_type)]
tissue_out_filtered_frequency_enc_p=tissue_out_filtered_frequency_enc[,list(percentage=round(count/sum(count)*100,digits=0),region_type=region_type),by=list(tissue)]
tissue_out_filtered_frequency_enc_p[,list(minp=min(percentage),maxp=max(percentage)),by=list(region_type)]
t.test(tissue_out_filtered_frequency_enc_p[region_type=="dMML only"]$percentage,
       tissue_out_filtered_frequency_enc_p[region_type=="dNME only"]$percentage,alternative = "less")
#Promoter only
tissue_out_filtered_frequency_prom=tissue_out_filtered[tissue!="NT"&promoter==T,list(count=length(region)),by=list(tissue,region_type)]
tissue_out_filtered_frequency_prom_p=tissue_out_filtered_frequency_prom[,list(percentage=round(count/sum(count)*100,digits=0),region_type=region_type),by=list(tissue)]
tissue_out_filtered_frequency_prom_p[,list(minp=min(percentage),maxp=max(percentage)),by=list(region_type)]
t.test(tissue_out_filtered_frequency_prom_p[region_type=="dMML only"]$percentage,
       tissue_out_filtered_frequency_prom_p[region_type=="dNME only"]$percentage,alternative = "less")

pdf('../downstream/output/graphs/Figure5/Figure5C_dMML_dNME_catogrization_all_non_repeats.pdf',width=7,height=16)
print(ggplot(tissue_out_filtered_frequency, aes(y=count, x=tissue,fill=region_type)) + 
        geom_bar( stat="identity",position="fill")+ylab("")+xlab("")+
        
        theme_classic()+theme(axis.text.x=element_text(size=32,angle=90),
                              axis.text.y=element_text(size=32))+
        scale_fill_manual(values=brewer.pal(4,'Set1'))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        theme(legend.position = "bottom",legend.title = element_blank(),legend.text =element_text(size=28)))
dev.off()


# Correlation using all regions -------------------------------------------

#Filtered result
cor_dt_unfiltered=lapply(names(dNME_cor),cor_dt_preprocessing,dMML_cor=dMML_cor,
                       dNME_cor=dNME_cor,dmml_perm=dmml_perm,dnme_perm=dnme_perm,
                       filtered=FALSE,folder_input=dir_in)
names(cor_dt_unfiltered)=names(dNME_cor)
saveRDS(cor_dt_unfiltered,'../downstream/output/mouse_analysis/correlation/correlation_dt_unfiltered_N17_all_regions.rds')
#Plot the density for each one
cor_dt_unfiltered=readRDS('../downstream/output/mouse_analysis/correlation/correlation_dt_unfiltered_N17_all_regions.rds')
tissue_out_unfiltered=mclapply(names(cor_dt_unfiltered),correlation_processing,cor_dt=cor_dt_unfiltered,filtered=F,subsmple_plot=1,
                           dir_in="../downstream/output/mouse_analysis/correlation/kmeans_10_run_unfiltered_all_regions/",mc.cores=7)
names(tissue_out_unfiltered)=names(cor_dt_unfiltered)
saveRDS(tissue_out_unfiltered,'tissue_out_unfiltered_N17_all_regions_clustering.rds')

#All regions 
tissue_out_all=readRDS('../downstream/output/correlation/tissue_out_unfiltered.rds')
tissue_out_all_frequency=do.call(rbind,tissue_out_all)
tissue_out_all_frequency=tissue_out_all_frequency[tissue!="NT",list(count=length(region)),by=list(tissue,region_type)]
#Change name
tissue_out_all_frequency[region_type=="dMML_only"]$region_type="UC_M+N-"
tissue_out_all_frequency[region_type=="dNME_only"]$region_type="UC_M-N+"
tissue_out_all_frequency[region_type=="Both"]$region_type="UC_M+N+"
tissue_out_all_frequency[region_type=="Neither"]$region_type="UC_M-N-"
tissue_out_all_frequency$region_type=factor(tissue_out_all_frequency$region_type,levels=c("UC_M+N-","UC_M-N+","UC_M+N+","UC_M-N-"))

pdf('../downstream/output/graphs/Figure5/Figure5CS_dMML_dNME_catogrization.pdf',width=10,height=7)
print(ggplot(tissue_out_all_frequency, aes(y=count, x=tissue,fill=region_type)) + 
        geom_bar( stat="identity",position="fill")+ylab("")+xlab("")+theme_glob+
        theme(legend.position = "bottom",legend.title = element_blank(),legend.text =element_text(size=16))+
        scale_fill_manual(values=brewer.pal(4,'Set1')))
dev.off()

tissue_out_all_frequency_p=tissue_out_all_frequency[,list(percentage=round(count/sum(count)*100,digits=0),region_type=region_type),by=list(tissue)]
tissue_out_all_frequency_p[,list(minp=min(percentage),maxp=max(percentage)),by=list(region_type)]
t.test(tissue_out_all_frequency_p[region_type=="UC_M+N-"]$percentage,
       tissue_out_all_frequency_p[region_type=="UC_M-N+"]$percentage,alternative = "less")
# Assign region type to each region in csv_folder -------------------------
tissue_out_filtered=readRDS('../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_all_region.rds')
#Assign DNase region
DNAase=readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp.rds')
DNAase=convert_GR(DNAase,direction="DT")

folder_in="../downstream/input/mouse_analysis/clustering/tissue_specific/currently_in_use/ts_cluster_0_1"
#For repeats

folder_in="../downstream/input/mouse_analysis/clustering/tissue_specific/currently_in_use/ts_cluster_0_1_non_repeats/"
folder_out=paste0(folder_in,'region_assigned/')
ifelse(!dir.exists(file.path(folder_out)), dir.create(file.path(folder_out)), FALSE)
#assign regions 
lapply(tissue_out_filtered,function(cor_dt_in){
  tissue_in=unique(cor_dt_in$tissue)
  csv_in=fread(paste0(folder_in,tissue_in,'.csv'))
  csv_in$region_type=cor_dt_in[match(csv_in$regions,cor_dt_in$region)]$region_type
  csv_in$DNAase=FALSE
  DNAase_region=which(csv_in$regions%in% DNAase$region)
  csv_in[DNAase_region]$DNAase=TRUE
  tss=get_mm10_tss()
  
 
  dt_nearest=GenomicRanges::distanceToNearest(convert_GR(csv_in$regions),tss)
  csv_in$distance=mcols(dt_nearest)$distance
  csv_in$gene=names(tss)[subjectHits(dt_nearest)]
  write.csv(csv_in,paste0(folder_out,tissue_in,'.csv'),row.names = F)
  return(NULL)
})

tissue_out_filtered=readRDS('../downstream/output/correlation/tissue_out_filtered.rds')
tissue_out_filtered_frequency=do.call(rbind,tissue_out_filtered)
tissue_out_filtered_frequency_promoter=tissue_out_filtered_frequency[abs(distance)<=2000,list(count=length(region)),by=list(tissue,region_type)]
ggplot(tissue_out_filtered_frequency_promoter, aes(y=count, x=tissue,fill=region_type)) + 
  geom_bar( stat="identity",position="fill")+ylab("")+xlab("")




tissue_out_filtered_frequency_enhancer=tissue_out_filtered_frequency[enhancer==TRUE,list(count=length(region)),by=list(tissue,region_type)]
ggplot(tissue_out_filtered_frequency_enhancer, aes(y=count, x=tissue,fill=region_type)) + 
  geom_bar( stat="identity",position="fill")+ylab("")+xlab("")

tissue_out_filtered_frequency_promoter=tissue_out_filtered_frequency[abs(distance)<=2000,list(count=length(region)),by=list(region_type)]
tissue_out_filtered_frequency_promoter$state="promoter"
tissue_out_filtered_frequency_enhancer=tissue_out_filtered_frequency[enhancer==TRUE,list(count=length(region)),by=list(region_type)]
tissue_out_filtered_frequency_enhancer$state="enhancer"
tissue_out_filtered_frequency_all=rbind(tissue_out_filtered_frequency_promoter,tissue_out_filtered_frequency_enhancer)
ggplot(tissue_out_filtered_frequency_all, aes(y=count, x=state,fill=region_type)) + 
  geom_bar( stat="identity",position="fill")+ylab("")+xlab("")

