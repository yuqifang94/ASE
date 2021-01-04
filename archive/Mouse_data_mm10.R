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
return(list(MML_in,NME_in))},mc.cores=20)

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
