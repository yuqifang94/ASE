#histone marker
read.histone<-function(in_dir,tissue,stage,replicates,marker){
  file_name=paste(in_dir,stage,'_',tissue,'_',replicates,marker,'.bed.gz',sep='')
  cat('Processing:',file_name,'\n')
  if(file.exists(file_name)){
    extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                              qValue = "numeric", peak = "integer")
    bed_in=import(file_name,format='BED',extraCols= extraCols_narrowPeak)
    bed_in$tissue=tissue
    bed_in$stage=stage
    bed_in$Sample=paste(tissue,stage,replicates,sep='-')
    bed_in$marker=marker
    bed_in$replicates=replicates
    return(bed_in)
  }else{cat(file_name,'Not exsit\n')}
  
}

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
meta_df=data.frame(tissue=tissue,stage=stage,sample=sample)
#in_dir='../downstream/data/histone_H3K4me3_peaks/'
in_dir=''
histone_marker=GRanges()
histone_marker=do.call('c',lapply(c('H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K27ac','H3K27me3','H3K36me3','H3K9me3'),function(mk){
  do.call('c',lapply(1:nrow(meta_df),function(i) {
    histone_marker=c(histone_marker,read.histone(in_dir,meta_df$tissue[i],meta_df$stage[i],replicates=NULL,mk))
    
  }))
}))
histone_marker=histone_marker[seqnames(histone_marker)%in%seqlevels(histone_marker)[1:21]]
saveRDS(histone_marker,'histone_marker_mm10_merged.rds')

DNase_mm10=readRDS('../downstream/output/mm10_DNase.rds')
histone_marker=readRDS('../downstream/output/histone_marker_mm10_merged.rds')
histone_marker$singnalValue[histone_marker$qValue>=-log10(0.05)]=0
histone_marker$Sample_marker=paste(histone_marker$Sample,histone_marker$marker,sep='')
DNase_mm10_output=lapply(unique(histone_marker$marker), function(mk,DNase_mm10,sig){
  cat('Processing:',mk,'\n')
  chromatin_sub=histone_marker[histone_marker$marker==mk]
  sp_olap=findOverlaps(DNase_mm10,chromatin_sub)
  sp_df=data.table(qt=queryHits(sp_olap),sample=chromatin_sub$Sample_marker[subjectHits(sp_olap)],
                   signal =mcols(chromatin_sub)[[sig]][subjectHits(sp_olap)])
  sp_df_dc=dcast.data.table(sp_df,qt~sample,value.var = "signal",fun.aggregate = mean,fill = NA)
  for(sp in unique(colnames(sp_df_dc))[-1]){
    elementMetadata(DNase_mm10)[[sp]]=NA
    elementMetadata(DNase_mm10)[[sp]][sp_df_dc$qt]=sp_df_dc[[sp]]
    
  }
  return(DNase_mm10)
},DNase_mm10=DNase_mm10,sig='singnalValue')
saveRDS(DNase_mm10_output,'../downstream/output/DNase_mm10_histone_signal.rds')

DNase_mm10_output=lapply(unique(histone_marker$marker), function(mk,DNase_mm10,sig){
  cat('Processing:',mk,'\n')
  chromatin_sub=histone_marker[histone_marker$marker==mk]
  sp_olap=findOverlaps(DNase_mm10,chromatin_sub)
  sp_df=data.table(qt=queryHits(sp_olap),sample=chromatin_sub$Sample_marker[subjectHits(sp_olap)],
                   signal =mcols(chromatin_sub)[[sig]][subjectHits(sp_olap)])
  sp_df_dc=dcast.data.table(sp_df,qt~sample,value.var = "signal",fun.aggregate = mean,fill = NA)
  for(sp in unique(colnames(sp_df_dc))[-1]){
    elementMetadata(DNase_mm10)[[sp]]=NA
    elementMetadata(DNase_mm10)[[sp]][sp_df_dc$qt]=sp_df_dc[[sp]]
    
  }
  return(DNase_mm10)
},DNase_mm10=DNase_mm10,sig='qValue')
#ATAC peaks
read.ATAC<-function(in_dir,tissue,stage,replicates){
  FC_name=paste(in_dir,stage,'_',tissue,'_',replicates,'_FC','.bigWig',sep='')
  pval_name=paste(in_dir,stage,'_',tissue,'_',replicates,'_FC','.bigWig',sep='')
  cat('Processing:',FC_name,'\n')
  if(file.exists(FC_name)&file.exists(pval_name)){
    FC_in=import.bw(FC_name)
    pval_in=import.bw(FC_name)
    #ATAC-seq are not having same genomic regions
    return(bed_in)
  }else{cat(file_name,'Not exsit\n')}
  
}
in_dir=''

ATAC_marker=GRanges()
for(i in 1:nrow(meta_df)){
  ATAC_marker=c(ATAC_marker,read.ATAC(in_dir,meta_df$tissue[i],meta_df$stage[i],'1',mk))
  
  
}

