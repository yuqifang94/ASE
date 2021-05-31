source('mainFunctions_sub.R')
# Merge UC,dMML,dNME ------------------------------------------------------


mml <- readRDS('../downstream/output/mouse_analysis/CPEL_outputs/mml_matrix_DNase.rds')
nme <- readRDS('../downstream/output/mouse_analysis/CPEL_outputs/nme_matrix_DNase.rds')

UC_merge=lapply(tissue_all,function(x){
  uc_in=uc[[x]]
  mml_in=mml[,grepl(x,colnames(mml))]
  nme_in=nme[,grepl(x,colnames(nme))]
  regions=intersect(intersect(rownames(uc_in),rownames(mml_in)),rownames(nme_in))
  uc_in=uc_in[regions,]
  
  colnames(nme_in)=gsub('.*-','',gsub("-all","",colnames(nme_in)))
  colnames(mml_in)=gsub('.*-','',gsub("-all","",colnames(mml_in)))
  mml_in=mml_in[regions,]
  nme_in=nme_in[regions,]
  time_series=colnames(uc_in)
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
names(UC_merge)=tissue_all
saveRDS(UC_merge,'../downstream/output/UC_dMML_dNME.rds')
UC_merge_max_loc=lapply(UC_merge,function(x){
  x=as.data.frame(x,keep.rownames = T)
  x$dMML_max_pair=apply(x[,grepl("dMML-",colnames(x))],1,max)
  x$dNME_max_pair=apply(x[,grepl("dNME-",colnames(x))],1,max)
  x$UC_max_pair=apply(x[,grepl("UC-",colnames(x))],1,max)
  x$dMML_max_time=gsub('dMML-','',colnames(x)[grepl("dMML-",colnames(x))][apply(x[,grepl("dMML-",colnames(x))],1,which.max)])
  x$dNME_max_time=gsub('dNME-','',colnames(x)[grepl("dNME-",colnames(x))][apply(x[,grepl("dNME-",colnames(x))],1,which.max)])
  x$UC_max_time=gsub('UC-','',colnames(x)[grepl("UC-",colnames(x))][apply(x[,grepl("UC-",colnames(x))],1,which.max)])
  adj_time=paste0('UC-',paste0("E",10:15,'.5'),'-',paste0("E",11:16,'.5'))
  x$UC_max_time_adj=gsub('UC-','',colnames(x)[(colnames(x)%in%adj_time)][apply(x[,(colnames(x)%in%adj_time)],1,which.max)])
  return(x)
  
})
saveRDS(UC_merge_max_loc,'../downstream/output/UC_merge_max_loc.rds')


motif_Ken_dir='../downstream/input/Ken_motif_binding_site/'
for(fn in dir(motif_Ken_dir,pattern='_TF_motif_site.rds')){
  motif_in=readRDS(paste0(motif_Ken_dir,fn))
  motif_in_dt=do.call(rbind,lapply(names(motif_in), function(x){
    motif=convert_GR(motif_in[[x]],direction="DT")
    
    motif$motif=gsub('.*_','',x)
    motif$motif_full_name=x
    return(motif)
  }))
  motif_in_dt=motif_in_dt[,list(gene=gene,dist=dist,motif=paste(unique(motif),collapse=";"),
                                motif_full_name=paste(unique(motif_full_name),collapse=";")),by=list(region)]
  motif_in_gr=convert_GR(motif_in_dt$region)
  mcols(motif_in_gr)=motif_in_dt[,list(gene,dist,motif,motif_full_name)]
  saveRDS(motif_in_gr,paste0(motif_Ken_dir,gsub("_TF_motif_site.rds","_TF_motif_site_merged.rds",fn)))
}

motif_Ken_dir='../downstream/input/Ken_motif_result/'
for(fn in dir(motif_Ken_dir,pattern='_motif_site_dNME.rds')){
  motif_in=readRDS(paste0(motif_Ken_dir,fn))
  if(length(motif_in)>0){
  motif_in_dt=do.call(rbind,lapply(names(motif_in), function(x){
    motif=convert_GR(motif_in[[x]],direction="DT")
    
    motif$motif=gsub('.*_','',x)
    motif$motif_full_name=x
    return(motif)
  }))
  motif_in_dt=motif_in_dt[,list(cluster=cluster,region_type=region_type,motif=paste(unique(motif),collapse=";"),
                                motif_full_name=paste(unique(motif_full_name),collapse=";")),by=list(region)]
  motif_in_gr=convert_GR(motif_in_dt$region)
  mcols(motif_in_gr)=motif_in_dt[,list(region_type,motif,motif_full_name)]
  saveRDS(motif_in_gr,paste0(motif_Ken_dir,gsub("_motif_site_dNME.rds","_motif_site_dNME_merged.rds",fn)))
  }
}
