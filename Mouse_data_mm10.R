rm(list=ls())
source('mainFunctions_sub.R')
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
in_dir=''
#Read in MML and NME
read_bed_out=mclapply(1:nrow(meta_df),function(i){
  NME_in=read.agnostic.mouse(in_dir,meta_df$tissue[i],meta_df$stage[i],'nme',replicate='all')
  
  MML_in=read.agnostic.mouse(in_dir,meta_df$tissue[i],meta_df$stage[i],'mml',replicate='all')

  NME_in$NME=NME_in$score
  MML_in$MML=MML_in$score
return(list(MML_in,NME_in))},mc.cores=10)

MML_in=fastDoCall('c',lapply(read_bed_out,function(x) x[[1]]))
NME_in=fastDoCall('c',lapply(read_bed_out,function(x) x[[2]]))

saveRDS(NME_in,'NME_agnostic_mouse_all_merged.rds')
saveRDS(MML_in,'MML_agnostic_mouse_all_merged.rds')
#Convert to matrix
NME_in_matrix=agnostic_matrix_conversion(NME_in[NME_in$N>=2])# 0.735 region have all data
MML_in_matrix=agnostic_matrix_conversion(MML_in[MML_in$N>=2],'MML')#0.735 region have all data
saveRDS(NME_in_matrix,'NME_matrix_mouse_all_dedup_N2.rds')
saveRDS(MML_in_matrix,'MML_matrix_mouse_all_dedup_N2.rds')
rm(MML_in)
rm(NME_in)
rm(MML_in_matrix)
rm(NME_in_matrix)
gc()
#UC_run_before_MDS folder
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
saveRDS(UC_out,'UC_agnostic_mouse_N2_MDS_fix.rds')#74% regiOn have all data

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

#Filter out the region with N >=18, UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds
UC_in=readRDS('UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')

gff_in=readGFFAsGRanges('../mm10_allele_agnostic_analysis.gff')
gff_in_sub=gff_in[as.numeric(gff_in$N)<=17]#98.3%
UC_in_sub=mclapply(UC_in,function(x) {return(subsetByOverlaps(x,gff_in_sub,type='equal'))})
# EFP : 0.9786391
# Lung : 0.9792112
# NT : 0.9789671
# forebrain : 0.9785839
# heart : 0.9785712
# hindbrain : 0.978633
# intestine : 0.9788505
# kidney : 0.9789867
# limb : 0.9786628
# liver : 0.978591
# midbrain : 0.9786526
# stomach : 0.9788745
saveRDS(UC_in_sub,'UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix_less_equal_17CG.rds')


jsd_in=readRDS('JSD_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')
names(jsd_in)=unlist(lapply(jsd_in,function(x) sub('_.*','',colnames(mcols(x))[1])))
for(ts in names(jsd_in)){
  colnames(mcols(jsd_in[[ts]]))=gsub('_','-',sub(paste0('-',ts),'',gsub("_all","",colnames(mcols(jsd_in[[ts]])))))
  mcols(jsd_in[[ts]])=mcols(jsd_in[[ts]])[,which(!grepl("P0",colnames(mcols(jsd_in[[ts]]))))]
  jsd_in[[ts]]=jsd_in[[ts]][which(rowSums(is.na(mcols(jsd_in[[ts]])))==0)]

}
saveRDS(jsd_in,'JSD_agnostic_mouse_matrix_dedup_N2_all_merged_ls_ft.rds')
gff_in=readGFFAsGRanges('../mm10_allele_agnostic_analysis.gff')
gff_in_sub=gff_in[as.numeric(gff_in$N)<=17]#98.3%
jsd_in_sub=mclapply(jsd_in,function(x) {return(subsetByOverlaps(x,gff_in_sub,type='equal'))})
for(ts in names(jsd_in_sub)){
  
  
  cat("top 10% quantile for",ts,'is',quantile(as.vector(as.matrix(mcols(jsd_in_sub[[ts]]))),probs=0.9),'\n')
  cat("Proportion regions letft for",ts,"is",length(jsd_in_sub[[ts]])/length(jsd_in[[ts]]),'\n')
}
# top 10% quantile for EFP is 0.2912969
# Proportion regions letft for EFP is 0.9874118
# top 10% quantile for Lung is 0.2613509
# Proportion regions letft for Lung is 0.9849098
# top 10% quantile for NT is 0.2768564
# Proportion regions letft for NT is 0.9861501
# top 10% quantile for forebrain is 0.2898555
# Proportion regions letft for forebrain is 0.9855473
# top 10% quantile for heart is 0.2799574
# Proportion regions letft for heart is 0.9858041
# top 10% quantile for hindbrain is 0.2813578
# Proportion regions letft for hindbrain is 0.9860615
# top 10% quantile for intestine is 0.2584561
# Proportion regions letft for intestine is 0.9845821
# top 10% quantile for kidney is 0.2568008
# Proportion regions letft for kidney is 0.985014
# top 10% quantile for limb is 0.2953809
# Proportion regions letft for limb is 0.9870893
# top 10% quantile for liver is 0.3264372
# Proportion regions letft for liver is 0.9849181
# top 10% quantile for midbrain is 0.2777929
# Proportion regions letft for midbrain is 0.985423
# top 10% quantile for stomach is 0.2534778
# Proportion regions letft for stomach is 0.9843143

saveRDS(jsd_in_sub,'JSD_agnostic_mouse_matrix_dedup_N2_all_merged_ls_less_equal_17CG.rds')