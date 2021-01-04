# Ken motif processing ----------------------------------------------------
source('mainFunctions_sub.R')

#hg19
GR_merge=readRDS(GR_merge_file)
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/DNase_hg19_250bp_control.rds')
JASPAR_motif=readRDS('../downstream/output/motif_JASPAR_hg19.rds')

#DNase
GR_merge_DNase=subsetByOverlaps(GR_merge,DNase)
GR_merge_control=subsetByOverlaps(GR_merge,control)


# Use mean NME for allele-specific analysis -------------------------------
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

# Use allele-agnostic NME -------------------------------
genomic_features=readRDS(genomic_features_file)
TSS=genomic_features$TSS
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
for(i in 1:12){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/JASPAR_motif_hg19_NME_',i,'_agnostic.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
for (n in 1:3){
  JASPAR_out=list()
  gc()
  for(i in 1:4+(n-1)*4){
    cat(i)
    JASPAR_out=append(JASPAR_out,mclapply(readRDS(paste0('../downstream/output/human_motif/JASPAR_motif_hg19_NME_',i,'_agnostic.rds',sep='')),
                                        function(x) {
                                          olap=findOverlaps(x,TSS,maxgap = 10000)
                                          NA_id=1:length(x)
                                          NA_id=NA_id[!NA_id%in%queryHits(olap)]
                                          mcols(x)[NA_id,]=NA
                                          return(x)},mc.cores=24))
    
  }
  saveRDS(JASPAR_out,paste('../downstream/output/JASPAR_motif_hg19_NME_',n,'_agnostic_merged_10k.rds',sep=''))
}


# DNase vs control --------------------------------------------------------
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/DNase_hg19_250bp_control.rds')
NME_in=readRDS('../downstream/output/human_motif/allele_agnostic_hg19_DNase_NME_homogeneous_excluding_dMML.rds')
JASPAR_motif=readRDS('../downstream/output/motif_JASPAR_hg19.rds')

#DNase
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

for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_control,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_motif/homogeneous/JASPAR_motif_hg19_NME_',i,'_agnostic_Control.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}


# DNase region using allelic ones -----------------------------------------
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


