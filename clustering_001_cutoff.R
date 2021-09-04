rm(list=ls())
source('mainFunctions_sub.R')
# clustering K-means --------------------------------------------------------------
dir_in='../downstream/input/mouse_analysis/clustering/tissue_specific/uc_001/'
total_run=10
#Convert into df with major
cutoffs=0.01
cluster_out=list()
for(fn in c(paste0('uc_',cutoffs,'_',1:10,'.rds'))){
  cluster_in=readRDS(paste0(dir_in,fn))
  for(ts in names(cluster_in)){
    if(fn==paste0('uc_',cutoffs,'_',1,'.rds')){
      cluster_out[[ts]]=data.table(regions=names(cluster_in[[ts]]),cluster_1=cluster_in[[ts]])
    }else{
      cluster_out[[ts]][[paste0("cluster_",gsub(paste0('uc_|.rds|',cutoffs,'_'),'',fn))]]=cluster_in[[ts]][ cluster_out[[ts]]$region]
      
      
    }
  }
  
  
}
#Find major cluster use cluster_1 as reference
cluster_out=lapply(cluster_out,function(x){
  for(i in 1:10){
    
    x[[paste0("major_cluster_",i)]]=as.numeric(NA)
    x[[paste0("major_cluster_in_",i)]]=as.numeric(NA)
    #For each cluster, find major cluster
    for(j in 1:10){
      x[cluster_1==j][[paste0("major_cluster_",i)]]= as.numeric(names(which.max(table(x[cluster_1==j][[paste0("cluster_",i)]]))))  
      x[cluster_1==j][[paste0("major_cluster_in_",i)]]=x[cluster_1==j][[paste0("major_cluster_",i)]]==x[cluster_1==j][[paste0("cluster_",i)]]#If in major cluster
    }
  }
  x$percent_cluster_in=rowSums(x[,grepl("major_cluster_in",colnames(x)),with=FALSE])/(total_run)
  return(x)
  
})




pdf('../downstream/output/mouse_analysis/clustering/proportion_run_kmeans_10_all_regions_001.pdf',width=3,height=3)
for(ts in names(cluster_out)){
  hist(cluster_out[[ts]]$percent_cluster_in,xlab="Proportion of runs in major cluster",main=ts)
  
  
}
dev.off()
#cluster$EFP['chr1:93970429-93970678',]
cluster=readRDS('../downstream/input/mouse_analysis/clustering/tissue_specific/uc_001/uc_0.01_1.rds')
UC_merge_max_loc=readRDS('../downstream/input/mouse_analysis/CPEL_output/UC_merge_max_loc_all_regions.rds')
UC_merge_max_loc_sub=lapply(names(UC_merge_max_loc),function(x) {
  print(x)
  return(UC_merge_max_loc[[x]][names(cluster[[x]]),])
  
})
names(UC_merge_max_loc_sub)=names(UC_merge_max_loc)
saveRDS(UC_merge_max_loc_sub,'../downstream/input/mouse_analysis/clustering/UC_merge_max_loc_cluster001.rds')
UC_merge=readRDS('../downstream/input/mouse_analysis/clustering/UC_merge_max_loc_cluster001.rds')
dir_out='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_01/'
cluster_region_out=list()
for(ts in names(cluster_out)){
  cluster_out_ts=cluster_out[[ts]]
  UC_ts=UC_merge[[ts]][,grepl('UC-',colnames(UC_merge[[ts]]))]
  
  #Define core clusters
  core_cluster=cluster_out_ts[percent_cluster_in==1]
  core_cluster=core_cluster[,list(regions,cluster_1,percent_cluster_in)]
  core_cluster=cbind(core_cluster,UC_ts[core_cluster$regions,])
  cols=colnames(core_cluster)[grepl(".5-E",colnames(core_cluster))]
  #find patterns of core clusters
  core_cluster_pattern=core_cluster[,lapply(.SD,mean),.SDcols=cols,by=cluster_1]
  core_cluster_pattern=core_cluster_pattern[order(cluster_1)]
  #Find cluster to assign
  cluster_to_assign=cluster_out_ts[percent_cluster_in>=0.5&percent_cluster_in<1]
  
  cluster_to_assign_UC=UC_ts[cluster_to_assign$regions,]
  #Each row is a region, each column is a cluster
  core_cluster_pattern_mt=as.matrix(core_cluster_pattern[,-1])
  rownames(core_cluster_pattern_mt)=core_cluster_pattern$cluster_1
  cor_cluster_out=cor(t(cluster_to_assign_UC),t(core_cluster_pattern_mt))
  #prepare to assign,make sure rows are consistent
  cor_cluster_out=cor_cluster_out[cluster_to_assign$regions,]
  cluster_to_assign$correlation=rowMax(cor_cluster_out)
  cluster_to_assign$cluster=colnames(cor_cluster_out)[apply(cor_cluster_out,1,which.max)]
  #Summary
  cluster_to_assign=cluster_to_assign[,list(regions,cluster,correlation)]
  cluster_to_assign$region_type="noncore_cluster"
  core_cluster$cluster=core_cluster$cluster_1
  core_cluster$correlation=1
  core_cluster$region_type="core_cluster"
  region_out=rbind(core_cluster[,list(regions,cluster,correlation,region_type)],cluster_to_assign[,list(regions,cluster,correlation,region_type)])
  region_out$tissue=ts
  region_out=region_out
  UC_max_ts=UC_merge[[ts]][,grepl('max',colnames(UC_merge[[ts]]))]
  
  UC_max_ts$UC_max_time_adj =gsub(paste0(ts,'-|-all'),'',UC_max_ts$UC_max_time_adj )
  UC_max_ts$UC_max_time  =gsub(paste0(ts,'-|-all'),'',UC_max_ts$UC_max_time  )
  region_out=cbind(region_out,UC_max_ts[region_out$regions,])
  cluster_region_out[[ts]]=region_out
  
  cat("Percent left for:",ts,nrow(region_out)/nrow(cluster_out_ts),'\n')
  write.csv(region_out,paste0(dir_out,'cluster_assigned/',ts,'.csv'))
}
cluster_region_out_fn='../downstream/output/mouse_analysis/clustering/cluster_all_region_assignment_filtered_001.rds'
saveRDS(cluster_region_out,cluster_region_out_fn)

# Plot heatmap ------------------------------------------------------------

library(RColorBrewer)
library(pheatmap)
library(gplots)
UC_merge=readRDS('../downstream/input/mouse_analysis/UC_only_all_regions.rds')
UC_merge=lapply(UC_merge,function(x) x[,!grepl('max',colnames(x))])
d <- lapply(UC_merge,function(x) x[,grepl('UC-',colnames(x))])
names(d)=names(UC_merge)
#dmml <- sapply(dmml,abs)
#dnme <- sapply(dnme,abs)
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
#clu <- readRDS('../downstream/input/Jason_UC_cluster/uc_0.1.rds')
clu=readRDS(cluster_region_out_fn)
clu=lapply(clu,function(x){
  out=as.numeric(x$cluster)
  names(out)=x$regions
  return(out)
  
})
d=d[tissue_all]
d=lapply(d,function(x) {
  
  colnames(x)=gsub(paste0('UC-|-all|',paste(tissue_all,'-',sep='',collapse = '|')),'',colnames(x))
  return(x)
})
d <- sapply(d,function(i) {
  
  i <- i[rowSums(i) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  i <- i[,order(match(colnames(i),timeorder))]
  
  #i <- scalematrix(i)
  i <- i[complete.cases(i),]
})

mat_out=matrix(ncol=39,nrow=0)
rowann_out=data.frame()
row_gap=c(0)
for (n in names(d)) {
  cl <- clu[[n]]
  cl <- sort(cl)
  mat <- do.call(cbind,sapply(tissue_all,function(i) {
    tmp <- matrix(NA,nrow=length(cl),ncol=ncol(d[[i]]),dimnames = list(names(cl),colnames(d[[i]])))
    
    rn <- intersect(names(cl),rownames(d[[i]]))
    tmp[rn,] <- as.matrix(d[[i]][rn,])
    
    colnames(tmp) <- paste0(i,':',colnames(tmp))
    
    tmp
  }))
  mat_out=rbind(mat_out,mat)
  print(head(mat_out))
  rowann <- data.frame(tissue_r=n,cluster=sub(':.*','',cl),
                       #dMMLJSDcor=dmml[[n]][rownames(mat)],
                       #dNMEJSDcor=dnme[[n]][rownames(mat)],
                       stringsAsFactors = F)
  rownames(rowann) <- rownames(mat)
  rowann <- rowann[,ncol(rowann):1]
  rowann_out=rbind(rowann_out,rowann)
  #row_gap=c(row_gap,row_gap[length(row_gap)]+cumsum(rle(sub(':.*','',cl))$lengths))
  row_gap=c(row_gap,row_gap[length(row_gap)]+nrow(mat))
}
#Refine plotting parameters
colann <- data.frame(time=sub('.*:','',colnames(mat_out)),tissue=sub(':.*','',colnames(mat_out)),stringsAsFactors = F)
rownames(colann) <- colnames(mat_out)


c1 <- mouse_color()
c2 <- brewer.pal(10,'Set3')
names(c2) <- 1:10
c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
names(c4) <- sort(unique(colann[,1]))
# tiff(paste0('../downstream/output/heatmap_acrosstissue/',n,'.tiff'),width=3000,height=3000,res=300)
# #png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
# pheatmap(mat,cluster_rows = F,annotation_row = rowann,cluster_cols = F,
#          annotation_col = colann,show_colnames = F,show_rownames = F,
#          gaps_row = row_gap,gaps_col = cumsum(rle(colann[,1])$lengths),
#          annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4,dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)))
# dev.off()

tiff('../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_01/all_sc_N17_ft_kmeans_10run_filtered_all_001.tiff',width=5000,height=5000,res=300)
#png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
pheatmap(scalematrix(mat_out),cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
         annotation_col = colann,show_colnames = F,show_rownames = F,
         gaps_row = row_gap,gaps_col = cumsum(rle(colann[,2])$lengths),
         annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4
                                  #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)
         ))
dev.off()

# #Add GO analysis for all regions ----------------------------------------
dir_out_cluster='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_01/'
folder_out=paste0(dir_out_cluster,'cluster_assigned/')
dir_out_GO='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_001/'
UC_merge=readRDS('../downstream/input/mouse_analysis/UC_only_all_regions.rds')#Define all analyzed regions, were using UC_merge_max_loc_cluster01.rds,4626
cutoff_fn='001'
#Runnning
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
#prepare enhancer background gene list
uc_gr=lapply(UC_merge,function(x) rownames(x))
uc_gr=Reduce(intersect,uc_gr)
uc_gr=convert_GR(uc_gr)
enhancer=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)
#Prepare promoter background gene
tss=get_mm10_tss()

bg_promoter=names(subsetByOverlaps(tss,uc_gr,maxgap = 2000))
enc_type="enhancer"
region_type='all'
GO_out_all=list()
for(ts in tissue_all){
  if(enc_type=="enhancer"){
    bg=bg_enhancer
  }else 
    if(enc_type=="promoter"){
      bg=bg_promoter
      
    }
  GO_out_all[[region_type]][[ts]]=GO_run_tissue(ts,folder_out,enc_type=enc_type,region_type_sel=region_type,bg=bg,DNase=F)
  GO_out_all[[region_type]][[ts]]=lapply(GO_out_all[[region_type]][[ts]],function(x){
    return(list(GO_out_cluster_all=x$GO_out_cluster_all,
                csv_in_ts_clu=cbind(x$csv_in_ts_clu,as.data.table(UC_merge[[ts]][x$csv_in_ts_clu$region,!grepl('max',colnames(UC_merge[[ts]]))]))))
    
    
  })
  
}
saveRDS(GO_out_all,paste0(dir_out_GO,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_',enc_type,'.rds'))
#Non of GO terms is significant
plot_GO_heatmap_all(tissue_all,GO_out_all[[region_type]],region_type=region_type,enc_type="enhancer",ptcount=0,FDR_cutoff=0.2,
                    dir_plot=paste0('../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/UC_',cutoff_fn,'/'))

