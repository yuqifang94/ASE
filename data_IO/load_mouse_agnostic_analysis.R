rm(list=ls())
source("mainFunctions_sub.R")

read.agnostic.mouse<-function(fn,in_dir,replicate="all"){
  fn_sub=gsub('mm10_|.bedGraph|_all_allele_agnostic','',fn)

  #tissue,stage,stat_type,replicate
  tissue=sub('_.*','',fn_sub)
  stat_type=sub('.*_','',fn_sub)
  stage=gsub(paste0(tissue,"_|_",stat_type),'',fn_sub)
  
  file_in=paste0(in_dir,fn)

  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    informME_in$tissue=tissue
    stage=gsub('_5','.5',stage)
    stage=gsub('day','E',stage)
    stage=gsub('E0','P0',stage)
    informME_in$stage=stage
    informME_in$bioreplicate=replicate
    informME_in$Sample=paste(tissue,stage,replicate,sep='-')
    informME_in$stat_type=stat_type
    return(informME_in)
  }
}
#UC
read.agnostic.mouse.uc<-function(file_in,matrix=FALSE,fileter_N=1,gff_in=NA){

  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    #process file name
    file_in=strsplit(file_in,'\\/')[[1]]
    file_in=file_in[length(file_in)]
    comp= strsplit(strsplit(file_in,'_uc.bedGraph')[[1]],'-vs-')[[1]]
    strain=unlist(lapply(strsplit(comp,'_'),function(x) x[1]))
    #if  contain BL6DBA, use ref is BL6DBA
    strain=ifelse('BL6DBA'%in%strain,'BL6DBA','mm10')
    comp=unlist(lapply(strsplit(comp,'_'),function(x) paste(x[-1],collapse = '_')))
    comp=comp[comp!='']
    comp_stage=unlist(lapply(comp,function(x) {x_split=strsplit(x,'_')[[1]]
    x_split=x_split[-length(x_split)][-1]
    x_split=paste(x_split,collapse = '_')
    return(x_split)}))
    tissue1=strsplit(comp[1],'_')[[1]][1]
    tissue2=strsplit(comp[2],'_')[[1]][1]
    #if BL6DBA, the 1st comp_stage is empty

    comp_stage=gsub('_5','.5',comp_stage)
    comp_stage=gsub('day','E',comp_stage)
    comp_stage=gsub('E0','P0',comp_stage)
    replicate=strsplit(comp[1],'_')[[1]][length(strsplit(comp[1],'_')[[1]])]
    replicate=gsub('merged','',replicate)
    informME_in$Sample=paste0(tissue1,'-',comp_stage[1],'-',tissue2,'-',comp_stage[2],'-',replicate)
    informME_in=informME_in[informME_in$N>=fileter_N]
    informME_in$tissue=tissue1
    informME_in$stage=paste0(comp_stage[1],'-',comp_stage[2])
    informME_in$replicate=replicate
    cat('Minimum N:',min(informME_in$N),'\n')
    #informME_in$Ref=strain
    if(matrix){
      informME_in_dt=as.data.table(mcols(informME_in))[,c("score","Sample")]
      colnames(informME_in_dt)=c("UC","Sample")
      informME_in_dt$UC=as.numeric(informME_in_dt$UC)
      informME_in_dt$region=paste0(seqnames(informME_in),":",start(informME_in),"-",end(informME_in))
      informME_in_dt=informME_in_dt[match(gff_in,region),"UC"]
      colnames(informME_in_dt)=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
      return(informME_in_dt)
    }
    else{return(informME_in)}
  }
}
agnostic_matrix_conversion<-function(gr_in,stat='NME'){
  gr_out=granges(unique(gr_in))
  olap=findOverlaps(gr_in,gr_out,type='equal')
  stat_in_df=elementMetadata(gr_in[queryHits(olap)])[c(stat,'Sample')]
  stat_in_df$idx=NA
  stat_in_df$idx[queryHits(olap)]=subjectHits(olap)
  stat_in_df=as.data.table(stat_in_df)
  
  stat_in_df_stat=dcast.data.table(data=stat_in_df,formula=idx~Sample,value.var = stat,fun.aggregate=mean)#remove agg.fun for new run
  gr_out=gr_out[stat_in_df_stat$idx]
  mcols(gr_out)=stat_in_df_stat[,-1]
  return(gr_out)
  
}
#filtering UC region based on available data
UC_filtering<-function(d){
  d <- sapply(d,function(am) {
    am <- am[,!grepl('P0',colnames(am))]
    am <- am[complete.cases(am),]
  })
  k <- table(unlist(sapply(d,rownames)))
  id <- names(k)[k==length(d)]
  d <- sapply(d,function(i) i[id,],simplify = F)
  return(d)
}
# reading in mouse MML and NME --------------------------------------------
#Complimentary regions
dir_comp='../downstream/data/compliment_MML_NME_model_mouse/'
MML_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*mml"),
                               read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
MML_in$MML=MML_in$score
MML_in$score=NULL
NME_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*nme"),read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
NME_in$NME=NME_in$score
NME_in$score=NULL
#Analyzed PRC, DNase and control
dir_analyzed='../downstream/data/DNase_control_PRC_MML_NME_model_mouse/'
MML_in_analyzed=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*mml"),
                               read.agnostic.mouse,in_dir=dir_analyzed,mc.cores=20))
MML_in_analyzed$MML=MML_in_analyzed$score
MML_in_analyzed$score=NULL
NME_in_analyzed=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*nme"),read.agnostic.mouse,in_dir=dir_analyzed,mc.cores=20))
NME_in_analyzed$NME=NME_in_analyzed$score
NME_in_analyzed$score=NULL
#Combine two datasets
NME_in=c(NME_in,NME_in_analyzed)
MML_in=c(MML_in,MML_in_analyzed)
#Convert to matrix: note here I didn't filter out the N<=17 since all NME and MML are intersect with UC regions in later analysis, the filtering in only done in UC
NME_in_matrix=agnostic_matrix_conversion(NME_in[NME_in$N>=2])
MML_in_matrix=agnostic_matrix_conversion(MML_in[MML_in$N>=2],'MML')

saveRDS(NME_in_matrix,NME_matrix_file)
saveRDS(MML_in_matrix,MML_matrix_file)
rm(MML_in)
rm(NME_in)
rm(MML_in_matrix)
rm(NME_in_matrix)
gc()
# reading in mouse UC --------------------------------------------
#Complimentary regions
UC_in_dir='../downstream/data/compliment_UC_non_MDS_mouse/'
UC_in=fastDoCall('c',mclapply(dir(UC_in_dir,pattern = '.*uc.bedGraph'),function(x){UC_in=read.agnostic.mouse.uc(paste(UC_in_dir,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=20))
UC_in$tissue=sub('-.*','',UC_in$Sample)
UC_in$Sample=sub('.5-.*-E1','.5-E1',UC_in$Sample)
#DNase,control,PRC
UC_in_dir_analyzed='../downstream/data/DNase_control_PRC_non_MDS_mouse/'
UC_in_analyzed=fastDoCall('c',mclapply(dir(UC_in_dir,pattern = paste0(paste(unique(UC_in$tissue),collapse='|'),'.*uc.bedGraph')),
                                       function(x){
  UC_in=read.agnostic.mouse.uc(paste(UC_in_dir_analyzed,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=20))
UC_in_analyzed$Sample=sub('.5-.*-E1','.5-E1',UC_in_analyzed$Sample)
UC_in_analyzed=UC_in_analyzed[UC_in_analyzed$Sample %in% unique(UC_in$Sample)]
#Merging data
UC_in=c(UC_in_analyzed[UC_in_analyzed$N>=2&UC_in_analyzed$N<=17],UC_in[UC_in$N>=2&UC_in$N<=17])
#Convert to matrix
UC_in_matrix_ls=mclapply(unique(UC_in$tissue),function(x) agnostic_matrix_conversion(UC_in[UC_in$tissue==x],'UC'),mc.cores=20)
names(UC_in_matrix_ls)=unique(UC_in$tissue)
saveRDS(UC_in_matrix_ls,UC_in_matrix_ls_file)


# UC for mouse MDS comparison ---------------------------------------------
gff_in_compliment=import.gff3(mouse_compliment_gff_file)
gff_in_compliment=paste0(seqnames(gff_in_compliment),':',start(gff_in_compliment),'-',end(gff_in_compliment))
UC_in_MDS_comp=data.table(region=gff_in_compliment)
compliment_MDS_dir='../downstream/data/compliment_UC_MDS_mouse/'
UC_in_MDS_comp_UC=fastDoCall('cbind',
                             mclapply(dir(compliment_MDS_dir,pattern = '.*uc.bedGraph'),function(x){
                               read.agnostic.mouse.uc(paste(compliment_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_compliment)},mc.cores=10))

UC_in_MDS_comp=cbind(UC_in_MDS_comp,UC_in_MDS_comp_UC)
#Filter based on N first to save space
DNase_conrol_MDS_dir='../downstream/data/DNase_control_PRC_MDS_mouse/'
gff_in_DNase=import.gff3(mouse_DNase_control_gff_file)
gff_in_DNase=paste0(seqnames(gff_in_DNase),':',start(gff_in_DNase),'-',end(gff_in_DNase))
UC_in_analyzed_MDS=data.table(region=gff_in_DNase)
UC_in_analyzed_MDS_UC=fastDoCall('cbind',
                                 mclapply(dir(compliment_MDS_dir,pattern = '.*uc.bedGraph'),function(x){
                                   read.agnostic.mouse.uc(paste(DNase_conrol_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_DNase)},mc.cores=10))
UC_in_analyzed_MDS=cbind(UC_in_analyzed_MDS,UC_in_analyzed_MDS_UC)
UC_in_MDS_all=rbind(UC_in_MDS_comp,UC_in_analyzed_MDS)
saveRDS(UC_in_MDS_all,UC_in_MDS_all_file)

# UC for mouse MDS comparison with P0---------------------------------------------
#Note old and new run have different names before and after -vs-
#Fixing this issue
for(fn in dir(compliment_MDS_dir_P0)){
    sample_name=unlist(strsplit(gsub('_uc.bedGraph','',fn),'-vs-'))
    sample_name_rev=paste0(DNase_conrol_MDS_dir,sample_name[2],'-vs-',sample_name[1],'_uc.bedGraph')
    if(file.exists(sample_name_rev)){
         cat('Reversing file name:',fn,' to ',sample_name_rev)
         file.rename(sample_name_rev,paste0(DNase_conrol_MDS_dir,fn))

    }
}
gff_in_compliment=import.gff3(mouse_compliment_gff_file)
gff_in_compliment=paste0(seqnames(gff_in_compliment),':',start(gff_in_compliment),'-',end(gff_in_compliment))
UC_in_MDS_comp_P0=data.table(region=gff_in_compliment)

UC_in_MDS_comp_P0_UC=fastDoCall('cbind',
                             mclapply(dir(compliment_MDS_dir_P0,pattern = '.*uc.bedGraph'),function(x){
                               read.agnostic.mouse.uc(paste(compliment_MDS_dir_P0,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_compliment)},mc.cores=20))
UC_in_MDS_comp_P0_UC=cbind(UC_in_MDS_comp_P0,UC_in_MDS_comp_P0_UC)


gff_in_DNase=import.gff3(mouse_DNase_control_gff_file)
gff_in_DNase=paste0(seqnames(gff_in_DNase),':',start(gff_in_DNase),'-',end(gff_in_DNase))
UC_in_analyzed_MDS_P0=data.table(region=gff_in_DNase)
UC_in_analyzed_MDS_P0_UC=fastDoCall('cbind',
                                 mclapply(dir(compliment_MDS_dir_P0,pattern = '.*uc.bedGraph'),function(x){
                                   read.agnostic.mouse.uc(paste(DNase_conrol_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_DNase)},mc.cores=20))
UC_in_analyzed_MDS_P0=cbind(UC_in_analyzed_MDS_P0,UC_in_analyzed_MDS_P0_UC)
UC_in_MDS_all_P0=rbind(UC_in_MDS_comp_P0_UC,UC_in_analyzed_MDS_P0)

saveRDS(UC_in_MDS_all_P0,UC_in_MDS_all_P0_file)
UC_in_all=readRDS(UC_in_MDS_all_file)
UC_in_MDS_all_P0_all=cbind(UC_in_MDS_all_P0,UC_in_all[,-1])
saveRDS(UC_in_MDS_all_P0_all, UC_in_MDS_all_P0_all_file)
# created merged object for all UC, dMML and dNME ----------------------------------------

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
#Filtering complete data and regions having all data: should be 4876367
UC_merge=readRDS(UC_merge_file)
UC_merge_max_loc=readRDS(UC_merge_max_loc_file)
UC_merge=UC_filtering(UC_merge)
UC_merge_max_loc=UC_merge=UC_filtering(UC_merge_max_loc)
saveRDS(UC_merge_max_loc,UC_merge_max_loc_file)
saveRDS(UC_merge,UC_merge_file)
#Subset regions by UC>0.1
cluster=readRDS(paste0(dir_cluster_in_01,'uc_0.1_1.rds'))
UC_merge_max_loc_sub=lapply(names(UC_merge_max_loc),function(x) {
  print(x)
  return(UC_merge_max_loc[[x]][names(cluster[[x]]),])
  
})
names(UC_merge_max_loc_sub)=names(UC_merge_max_loc)
saveRDS(UC_merge_max_loc_sub,UC_merge_max_loc_01_file)

#Read in mouse NME and scRNA
NME_in=readRDS(NME_matrix_file)
#From JASON

mcols(NME_in)=mcols(NME_in)[,grepl('limb',colnames(mcols(NME_in)))]

gtf <- fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
genes <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
genes$gene_name <- gn
NME_in=dist_calc(NME_in,genes)
#Percent gene covered?
length(unique(NME_in[abs(NME_in$dist)<=3000]$gene))/length(genes[seqnames(genes)!="chrM"])#96%
NME_in_dt=convert_GR(NME_in,dir='DT')
NME_in_dt=melt.data.table(NME_in_dt,id.var=c('dist','gene','region'),value.name='NME',variable.name='stage')

NME_in_dt$hyper_var=-100
NME_in_dt$var=-100
NME_in_dt$mean=-100
for(st in unique(NME_in_dt$stage)){
  tt1=proc.time()[[3]]
  if(file.exists(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))){
    scRNA_in=readRDS(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))
    scRNA_in=scRNA_in[rownames(scRNA_in)%in% unique(c(NME_in_dt[(stage==st)]$gene)),]
    if(nrow(scRNA_in)>0){
      #Add hypervar to TSS 
      NME_in_dt[(stage==st)]$hyper_var=scRNA_in[NME_in_dt[(stage==st)]$gene,"hypervar_logvar"]
      NME_in_dt[(stage==st)]$var=scRNA_in[NME_in_dt[(stage==st)]$gene,"var"]
      NME_in_dt[(stage==st)]$mean=scRNA_in[NME_in_dt[(stage==st)]$gene,"mean"]
      
    }
  }else{cat("File not exist for ",st,'\n')}
  cat('Finish processing ',sub('E','',st),'in: ',proc.time()[[3]]-tt1,'\n')
  
}
saveRDS(NME_in_dt,NME_mouse_MAV_fn)
