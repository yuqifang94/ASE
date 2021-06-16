rm(list=ls())
source('mainFunctions_sub.R')
# clustering K-means --------------------------------------------------------------
total_run=10
#Convert into df with major
cutoffs=0.1
cluster_out=list()
for(fn in c(paste0('uc_',cutoffs,'_',1:10,'.rds'))){
  cluster_in=readRDS(paste0(dir_cluster_in_01,fn))
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




pdf(paste0(figure_path,'proportion_run_kmeans_10_all_regions_mm10.pdf'),width=3,height=3)
for(ts in names(cluster_out)){
  hist(cluster_out[[ts]]$percent_cluster_in,xlab="Proportion of runs in major cluster",main=ts)
  
  
}
dev.off()


# Find regions belong to major cluster ------------------------------------
UC_merge=readRDS(UC_merge_max_loc_file)
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
  write.csv(region_out,paste0(dir_out_cluster01,ts,'.csv'))
}

saveRDS(cluster_region_out,cluster_region_out_fn)
# Put regions with other info ----------------------------------------------------

lapply(names(cluster_out),function(x){
  cluster_out_ts=cluster_out[[x]]
  UC_merge_max_loc_sub_ts=UC_merge_max_loc_sub[[x]]
  
})
# Percent left for: EFP 0.9209754 
# Percent left for: NT 0.9990584 
# Percent left for: forebrain 0.9823054 
# Percent left for: heart 0.9347036 
# Percent left for: hindbrain 0.9234645 
# Percent left for: limb 0.9862634 
# Percent left for: liver 0.8056983 
# Percent left for: midbrain 0.9104696 

cluster_result=readRDS(cluster_region_out_fn)

# Plot heatmap ------------------------------------------------------------

library(RColorBrewer)
library(pheatmap)
library(gplots)
source('mainFunctions_sub.R')
UC_merge=readRDS(UC_merge_file)
UC_merge=lapply(UC_merge,function(x) x[,!grepl('max',colnames(x))])
d <- lapply(UC_merge,function(x) x[,grepl('UC-',colnames(x))])
names(d)=names(UC_merge)
dmml <-readRDS(dmml_cor_file)
dnme <-readRDS(dnme_cor_file)
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

mat_out=matrix(ncol=39)
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
sum(rowSums(is.na(mat_out))!=0)/nrow(mat_out)
#Refine plotting parameters
colann <- data.frame(time=sub('.*:','',colnames(mat_out)),tissue=sub(':.*','',colnames(mat_out)),stringsAsFactors = F)
rownames(colann) <- colnames(mat_out)
c1 <- mouse_color()
c2 <- brewer.pal(10,'Set3')
names(c2) <- 1:10
c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
names(c4) <- sort(unique(colann[,1]))
tiff('../downstream/output/mouse_analysis/clustering/heatmap_acrosstissue/all_sc_N17_ft_kmeans_10run_filtered_all.tiff',width=5000,height=5000,res=300)
pheatmap(scalematrix(mat_out),cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
         annotation_col = colann,show_colnames = F,show_rownames = F,
         gaps_row = row_gap,gaps_col = cumsum(rle(colann[,2])$lengths),
         annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4
                                  #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)
                                  ))
dev.off()

# created merged object for all UC, dMML and dNME ----------------------------------------
#in cpelasm
mml <- readRDS(MML_matrix_file)
mml=convert_GR(mml,direction="matrix")
nme <- readRDS(NME_matrix_file)
nme=convert_GR(nme,direction="matrix")
uc=readRDS(UC_in_matrix_ls_file)
uc=lapply(uc,convert_GR,direction="matrix")
UC_merge=lapply(names(uc),function(x){
  uc_in=uc[[x]]
  mml_in=mml[,grepl(x,colnames(mml))]
  nme_in=nme[,grepl(x,colnames(nme))]
  regions=intersect(intersect(rownames(uc_in),rownames(mml_in)),rownames(nme_in))
  uc_in=uc_in[regions,]
  
  colnames(nme_in)=gsub('.*-','',gsub("-all","",colnames(nme_in)))
  colnames(mml_in)=gsub('.*-','',gsub("-all","",colnames(mml_in)))
  mml_in=mml_in[regions,]
  nme_in=nme_in[regions,]
  time_series=gsub(paste0(x,'-'),'',gsub("-all","",colnames(uc_in)))
  dnme=do.call(cbind,lapply(time_series,function(x){
    return(abs(nme_in[,gsub('-.*','',x)]-nme_in[,gsub('.*-','',x)]))
    
  }))
  dmml=do.call(cbind,lapply(time_series,function(x){
    return(abs(mml_in[,gsub('-.*','',x)]-mml_in[,gsub('.*-','',x)]))
    
  }))
  colnames(uc_in)=paste0("UC-",colnames(uc_in))
  colnames(dmml)=paste0("dMML-",time_series)
  colnames(dnme)=paste0("dNME-",time_series)
  return(cbind(uc_in,dmml,dnme))
})
names(UC_merge)=names(uc)
saveRDS(UC_merge,UC_merge_file)
UC_merge_max_loc=lapply(UC_merge,function(x){
  cat("Percent all data:",sum(rowSums(is.na(x))==0)/nrow(x),'\n')
  x=as.data.frame(x[rowSums(is.na(x))==0,])
  uc_dt=  x[,grepl("UC-",colnames(x))]
  dNME_dt=  x[,grepl("dNME-",colnames(x))]
  dMML_dt=  x[,grepl("dMML-",colnames(x))]
  
  x$dMML_max_pair=apply(dMML_dt,1,max)
  x$dNME_max_pair=apply(dNME_dt,1,max)
  x$UC_max_pair=apply(uc_dt,1,max)
  x$dMML_max_time=gsub('dMML-','',colnames(dMML_dt)[apply(dMML_dt,1,which.max)])
  x$dNME_max_time=gsub('dNME-','',colnames(dNME_dt)[apply(dNME_dt,1,which.max)])
  
  uc_max=apply(uc_dt,1,which.max)
  x$UC_max_time=gsub('UC-','',colnames(uc_dt)[uc_max])
  x$dNME_max_UC_pair=dNME_dt[cbind(seq_along(uc_max), uc_max)]
  #x$UC_max_UC_pair=uc_dt[cbind(seq_along(uc_max), uc_max)]
  x$dMML_max_UC_pair=dMML_dt[cbind(seq_along(uc_max), uc_max)]
  adj_time=paste0(paste0("E",10:15,'.5'),'-',paste0("E",11:16,'.5'))
  uc_max_adj=unlist(apply(x[,(grep(paste0('UC-.*',adj_time,collapse="|",sep=''),colnames(x)))],1,which.max))
  
  x$UC_max_time_adj=gsub('UC-','',colnames(x))[(grepl(paste0('UC-.*',adj_time,collapse="|",sep=''),colnames(x)))][uc_max_adj]
  x$dNME_max_UC_pair_adj=x[,(grepl(paste0('dNME-.*',adj_time,collapse="|",sep=''),colnames(x)))][cbind(seq_along(uc_max_adj), uc_max_adj)]
  x$UC_max_UC_pair_adj=x[,(grepl(paste0('UC-.*',adj_time,collapse="|",sep=''),colnames(x)))][cbind(seq_along(uc_max_adj), uc_max_adj)]
  x$dMML_max_UC_pair_adj=x[,(grepl(paste0('dMML-.*',adj_time,collapse="|",sep=''),colnames(x)))][cbind(seq_along(uc_max_adj), uc_max_adj)]
  return(x)
  
})
names(UC_merge_max_loc)=names(UC_merge)
saveRDS(UC_merge_max_loc,UC_merge_max_loc_file)

cluster=readRDS(paste0(dir_cluster_in_01,'uc_0.1_1.rds'))
UC_merge_max_loc_sub=lapply(names(UC_merge_max_loc),function(x) {
  print(x)
  return(UC_merge_max_loc[[x]][names(cluster[[x]]),])
  
})
names(UC_merge_max_loc_sub)=names(UC_merge_max_loc)
saveRDS(UC_merge_max_loc_sub,UC_merge_max_loc_01_file)