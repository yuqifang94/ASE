library(data.table)
#generating sbatching files
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

sample=paste(tissue,stage,sep='_')
meta_df=data.table(tissue=tissue,stage=stage,sample=sample,stringsAsFactors = F)
time_comp=c()
sink("../downstream/output/UC_submission_all_merged")
cat('#!/bin/bash\n')
for(ts in unique(meta_df$tissue)){
  meta_df_ts=meta_df[tissue==ts]
  stage_num=nrow(meta_df_ts)
  for(i in 1:stage_num){
    #cat i+1 through stage num
    if(i<stage_num){
      for(j in (i+1):stage_num){
        time_out=paste0('sbatch cpelasm_allele_agnostic_uc.slrm ',"mm10_",meta_df_ts$sample[i],'_all ',"mm10_",meta_df_ts$sample[j],'_all\n')
        time_comp=c(time_comp,time_out)
        cat(time_out)
      }
    }
    
  }
  
}
sink()
#MDS plot require JSD between tissue in addition to previous ones
# comparison_output=data.table(sample1=c(),sample2=c())
# for(ts in unique(meta_df$tissue)){
#   meta_df_ts=meta_df[meta_df$tissue==ts]
#   meta_df_ts_others=meta_df[meta_df$tissue!=ts]
#   for(stage in meta_df_ts$stage){
#     sample1=paste0(ts,'-',stage,'_all')
#     for(sample2 in meta_df_ts_others$sample){
#       sample2=paste0(sample2,"_all")
#     if((!paste(sample1,sample2,sep=';') %in% paste(comparison_output$sample1,comparison_output$sample2,sep=';'))&
#        (!paste(sample2,sample1,sep=';') %in% paste(comparison_output$sample1,comparison_output$sample2,sep=';'))){
#       comparison_output=rbind(comparison_output,data.table(sample1=sample1,sample2=sample2))
#     }
#     }
#   }
# }
samples_all=combn(meta_df$sample,2)
samples_all=data.table(sample1=samples_all[1,],sample2=samples_all[2,])
samples_write=paste0('sbatch cpelasm_allele_agnostic_uc.slrm ',"mm10_",samples_all$sample1,'_all ',"mm10_",samples_all$sample2,'_all\n')

sink("../downstream/output/UC_submission_all_MDS")
cat('#!/bin/bash\n')
cat(samples_write[!unlist(lapply(samples_write,function(x) x%in%time_comp))])
sink()
#the 2 combined should be same as 
