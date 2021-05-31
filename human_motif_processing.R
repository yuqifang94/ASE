# Ken motif processing ----------------------------------------------------
source('mainFunctions_sub.R')

#hg19
GR_merge=readRDS(GR_merge_file)
DNase=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp_control.rds')
JASPAR_motif=readRDS('../downstream/output/motif_JASPAR_hg19.rds')

#DNase
GR_merge_DNase=subsetByOverlaps(GR_merge,DNase)
GR_merge_control=subsetByOverlaps(GR_merge,control)


# Use allele-specific analysis -------------------------------
#/scratch/users/yfang27@jhu.edu/yfang/allele/Running_Code
#mean NME
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_DNase,stat_in='NME',mc.cores=12)
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/allelic_motif/JASPAR_motif_hg19_DNase_allelic_NME_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_control,stat_in='NME',mc.cores=12)
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/allelic_motif/JASPAR_motif_hg19_control_allelic_NME_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#mean MML
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=lapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_DNase,stat_in='MML')
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/allelic_motif/JASPAR_motif_hg19_DNase_allelic_MML_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=lapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_control,stat_in='MML')
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/allelic_motif/JASPAR_motif_hg19_control_allelic_MML_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}

# DNase vs control --------------------------------------------------------
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/DNase_hg19_250bp_control.rds')
NME_in=readRDS('../downstream/output/human_analysis/CPEL_outputs/allele_agnostic_hg19_DNase_NME_homogeneous_excluding_dMML.rds')
JASPAR_motif=readRDS('../downstream/output/motif_JASPAR_hg19.rds')

#DNase NME
NME_in_DNase=subsetByOverlaps(NME_in,DNase,type='equal')
NME_in_control=subsetByOverlaps(NME_in,control,type='equal')
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_DNase,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/homogeneous/JASPAR_motif_hg19_NME_',i,'_agnostic_DNase.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#Control NME
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_control,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/homogeneous/JASPAR_motif_hg19_NME_',i,'_agnostic_Control.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}


#DNase MML
MML_in=readRDS('../downstream/output/allele_agnostic_hg19_DNase_MML_homogeneous.rds')
MML_in_DNase=subsetByOverlaps(MML_in,DNase,type='equal')
MML_in_control=subsetByOverlaps(MML_in,control,type='equal')

split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=lapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_DNase,stat_in='MML')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/homogeneous/JASPAR_motif_hg19_MML_',i,'_agnostic_DNase.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#control MML
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_control,stat_in='MML')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/homogeneous/JASPAR_motif_hg19_MML_',i,'_agnostic_Control.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
# DNase region using allelic region -----------------------------------------
GR_merge=readRDS(GR_merge_file)
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
GR_merge_DNase=NME_dNME_ken(DNase,GR_merge,"NME")
GR_merge_DNase_mt=as.matrix(mcols(GR_merge_DNase))
rownames(GR_merge_DNase_mt)=paste0(seqnames(GR_merge_DNase),':',start(GR_merge_DNase),'-',end(GR_merge_DNase))
saveRDS(GR_merge_DNase_mt,'../downstream/output/DNase_mt_SNP_allelic.rds')
#Dnase region using allele-agnostic model at SNP

NME_in=readRDS('../downstream/output/NME_agnostic_ASM.rds')
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
GR_merge_DNase=NME_dNME_ken(DNase,NME_in,"NME")
GR_merge_DNase_mt=as.matrix(mcols(GR_merge_DNase))
rownames(GR_merge_DNase_mt)=paste0(seqnames(GR_merge_DNase),':',start(GR_merge_DNase),'-',end(GR_merge_DNase))
saveRDS(GR_merge_DNase_mt,'../downstream/output/DNase_mt_SNP_agnostic.rds')
#DNase region at DNase regions all
NME_in=readRDS('../downstream/output/allele_agnostic_hg19_DNase_NME.rds')
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
NME_in=subsetByOverlaps(NME_in,DNase,type='equal')
GR_merge_DNase=NME_dNME_ken(DNase,NME_in,"NME")
GR_merge_DNase_mt=as.matrix(mcols(GR_merge_DNase))
rownames(GR_merge_DNase_mt)=paste0(seqnames(GR_merge_DNase),':',start(GR_merge_DNase),'-',end(GR_merge_DNase))
saveRDS(GR_merge_DNase_mt,'../downstream/output/DNase_mt_all_agnostic.rds')

#Add NA column to non_NA thing
samples_in=readRDS('../downstream/output/huamn_samples.rds')
folder_in='../downstream/output/human_motif/allelic_motif/'
folder_in='../downstream/output/human_motif/homo'
for(fn in dir(folder_in,pattern="MML")){
  MML_in=readRDS(paste0(folder_in,fn))
  for(idx in which(unlist(lapply(MML_in,function(x) ncol(mcols(x))!=49)))){
    
    for (sp in samples_in[!samples_in %in% colnames(mcols(MML_in[[idx]]))]){
      
      mcols(MML_in[[idx]])[[sp]]=as.numeric(NA)
      
    }
    
    MML_in[[idx]]=MML_in[[idx]][,samples_in]
    
    
  }
  saveRDS(MML_in,paste0(folder_in,gsub('.rds','.complete.rds',fn)))
}

#Prepare for GWAS analysis
NME_in=readRDS('../downstream/output/human_analysis/CPEL_outputs/allele_agnostic_hg19_DNase_NME_homogeneous_excluding_dMML.rds')
mcols(NME_in)=mcols(NME_in)[,c('Sample','NME')]
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp_control.rds')
NME_in$DNAase="NA"
olap_DNase=findOverlaps(NME_in,DNase,type='equal')
olap_control=findOverlaps(NME_in,control,type='equal')
NME_in$DNAase[queryHits(olap_DNase)]="DNAase"
NME_in$DNAase[queryHits(olap_control)]="control"
saveRDS(NME_in,'../downstream/output/human_analysis/CPEL_outputs/NME_DNAase_control_hg19.rds')
