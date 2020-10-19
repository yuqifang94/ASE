rm(list=ls())
library(rtracklayer)
library(Gmisc)
library(data.table)
library(ggfortify)
library(pheatmap)
library(matrixStats)
library(RColorBrewer)
library(parallel)


tissue=c(rep('kidney',4),rep('Lung',4),rep('forebrain',8),rep('liver',7),
         rep('heart',8),rep('hindbrain',8),rep('midbrain',8),
         rep('limb',6),rep('EFP',6),rep('NT',5),rep('intestine',4),rep('stomach',4))
stage=c('day14_5','day15_5','day16_5','day0',
        'day14_5','day15_5','day16_5','day0',
        'day10_5','day11_5','day12_5','day13_5','day14_5','day15_5','day16_5','day0',
        'day11_5','day12_5','day13_5','day14_5','day15_5','day16_5','day0',
        'day10_5','day11_5','day12_5','day13_5','day14_5','day15_5','day16_5','day0',
        'day10_5','day11_5','day12_5','day13_5','day14_5','day15_5','day16_5','day0',
        'day10_5','day11_5','day12_5','day13_5','day14_5','day15_5','day16_5','day0',
        'day10_5','day11_5','day12_5','day13_5','day14_5','day15_5',
        'day10_5','day11_5','day12_5','day13_5','day14_5','day15_5',
        'day11_5','day12_5','day13_5','day14_5','day15_5',
        'day14_5','day15_5','day16_5','day0',
       'day14_5','day15_5','day16_5', 'day0')

sample=paste(tissue,stage,sep='-')
meta_df=data.frame(tissue=tissue,stage=stage,sample=sample,stringsAsFactors = F)
read.agnostic.mouse<-function(in_dir,tissue,stage,stat_type,replicate){
  file_in=paste(in_dir,'mm10_',tissue,'_',stage,'_',replicate,'_allele_agnostic_',stat_type,'.bedGraph',sep='')
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
    file_in=gsub('_uc.bedGraph','',file_in)
    file_in=gsub('_jsd.bedGraph','',file_in)
    comp= strsplit(file_in,'-vs-')[[1]]
    comp= strsplit(file_in,'-vs-')[[1]]
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
    informME_in$Sample=paste0(tissue1,'-',comp_stage[1],'-',comp_stage[2],'-',replicate)
    informME_in=informME_in[informME_in$N>=fileter_N]
    cat('Minimum N:',min(informME_in$N),'\n')
    #informME_in$Ref=strain
    if(matrix){
      informME_in$Samplepaste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
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
in_dir=''
read_bed_out=mclapply(1:nrow(meta_df),function(i){
  NME_in=read.agnostic.mouse(in_dir,meta_df$tissue[i],meta_df$stage[i],'nme',replicate='all')
  
  MML_in=read.agnostic.mouse(in_dir,meta_df$tissue[i],meta_df$stage[i],'mml',replicate='all')

  NME_in$NME=NME_in$score
  MML_in$MML=MML_in$score
return(list(MML_in,NME_in))},mc.cores=20)

MML_in=fastDoCall('c',lapply(read_bed_out,function(x) x[[1]]))
NME_in=fastDoCall('c',lapply(read_bed_out,function(x) x[[2]]))

saveRDS(NME_in,'NME_agnostic_mouse_all_merged.rds')
saveRDS(MML_in,'MML_agnostic_mouse_all_merged.rds')

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

NME_in_matrix=agnostic_matrix_conversion(NME_in[NME_in$N>=2])# 0.735 region have all data
MML_in_matrix=agnostic_matrix_conversion(MML_in[MML_in$N>=2],'MML')#0.735 region have all data
saveRDS(NME_in_matrix,'NME_matrix_mouse_all_dedup_N2.rds')
saveRDS(MML_in_matrix,'MML_matrix_mouse_all_dedup_N2.rds')
rm(MML_in)
rm(NME_in)
rm(MML_in_matrix)
rm(NME_in_matrix)
gc()
#UC_run_before_MDS
in_dir='./'
UC_in=GRanges()
UC_in_ls=mclapply(dir(in_dir,pattern = 'mm10.*uc.bedGraph'),function(x){UC_in=read.agnostic.mouse.uc(paste(in_dir,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=24)
UC_in=fastDoCall('c',UC_in_ls)
UC_in=UC_in[UC_in$N>=2]
#saveRDS(UC_in,'UC_agnostic_mouse_dedup_N2_all_time_fix_UC.rds')#74% regiOn have all data
#UC_in_matrix=agnostic_matrix_conversion(UC_in,'UC')#duplicated regions due to error,waiting for one more to finish
UC_in$tissue=sub('-.*','',UC_in$Sample)
UC_in_matrix_ls=mclapply(unique(UC_in$tissue),function(x) agnostic_matrix_conversion(UC_in[UC_in$tissue==x],'UC'),mc.cores=12)
names(UC_in_matrix_ls)=unique(UC_in$tissue)
saveRDS(UC_in_matrix_ls,'UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')#74% regiOn have all data

#read in UC for MDS
in_dir='./'
MDS_file=dir(in_dir,pattern = 'mm10.*uc.bedGraph')

gff_in=import.gff3('../mm10_allele_agnostic_analysis.gff')
gff_in=paste0(seqnames(gff_in),':',start(gff_in),'-',end(gff_in))
UC_out=data.table(region=gff_in)

UC_in=fastDoCall('cbind',
                mclapply(MDS_file,function(x){
                  read.agnostic.mouse.uc(paste(in_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in)},mc.cores=24))
#saveRDS(UC_in,paste0('UC_agnostic_mouse_dedup_MDS_',i,'.rds'))#74% regiOn have all data
UC_out=cbind(UC_out,UC_in)


#read in JSD for MDS
in_dir='./'
JSD_in=GRanges()
JSD_in_ls=mclapply(dir(in_dir,pattern = '.*jsd.bedGraph'),function(x){JSD_in=read.agnostic.mouse.uc(paste(in_dir,x,sep=''))
JSD_in$JSD=JSD_in$score
return(JSD_in)},mc.cores=24)
JSD_in=fastDoCall('c',JSD_in_ls)
JSD_in$tissue=sub('_.*','',JSD_in$Sample)

JSD_in_matrix_ls=mclapply(unique(JSD_in$tissue),function(x) agnostic_matrix_conversion(JSD_in[JSD_in$tissue==x],'JSD'),mc.cores=12)
saveRDS(JSD_in,'JSD_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')#74% regiOn have all data
