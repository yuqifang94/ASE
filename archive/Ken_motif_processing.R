source('mainFunctions_sub.R')
#hg19

GR_merge=readRDS(GR_merge_file)
NME_in=readRDS(NME_agnostic_file)
mcols(NME_in)=mcols(NME_in)[,-1]
JASPAR_motif=readRDS('../downstream/output/motif_JASPAR_hg19.rds')
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
for(i in 1:12){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in,stat_in='NME',mc.cores=8)
  saveRDS(JASPAR_motif_sp,paste('../downstream/output/JASPAR_motif_hg19_NME_',i,'_agnostic.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}

for (n in 1:3){
  JASPAR_out=list()
  gc()
  for(i in 1:4+(n-1)*4){
    cat(i)
    JASPAR_out=append(JASPAR_out,lapply(readRDS(paste('../downstream/output/JASPAR_motif_hg19_NME_',i,'_agnostic.rds',sep='')),
                                        function(x) {
                                          olap=findOverlaps(x,genomic_features,maxgap = 10000)
                                          NA_id=1:length(x)
                                          NA_id=NA_id[!NA_id%in%queryHits(olap)]
                                          mcols(x)[NA_id,]=NA
                                          return(x)}))
    
  }
  saveRDS(JASPAR_out,paste('../downstream/output/JASPAR_motif_hg19_NME_',n,'_agnostic_merged_10k.rds',sep=''))
}
#mm10
NME_in=readRDS('../downstream/output/NME_agnostic_mouse_3kTSS_FANTOM.rds')
mcols(NME_in)=mcols(NME_in)[,-c(1,2,3,4,5)]
JASPAR_motif=readRDS('../downstream/input/motif_JASPAR_mm10.rds')
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
for(i in 1:12){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in,stat_in='NME',mc.cores=1)
  saveRDS(JASPAR_motif_sp,paste('../downstream/output/JASPAR_motif_mm10_NME_',i,'_agnostic_bio1.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}

MML_in=readRDS('../downstream/output/MML_agnostic_mouse_3kTSS_FANTOM.rds')
MML_in=MML_in[MML_in$bioreplicate==1]
mcols(MML_in)=mcols(MML_in)[,-c(1,2,3,4,5,6)]
JASPAR_motif=readRDS('../downstream/input/motif_JASPAR_mm10.rds')
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
for(i in 5:12){
  i=12
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in,stat_in='MML',mc.cores=1)
  saveRDS(JASPAR_motif_sp,paste('../downstream/output/JASPAR_motif_mm10_MML_',i,'_agnostic_bio1.rds',sep=''))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}


for (n in 1:3){
  JASPAR_out=list()
  gc()
  for(i in 1:4+(n-1)*4){
    cat(i)
    JASPAR_out=append(JASPAR_out,lapply(readRDS(paste('../downstream/output/JASPAR_motif_hg19_NME_',i,'_agnostic.rds')),
                                        function(x) {
                                          olap=findOverlaps(x,genomic_features,maxgap = 10000)
                                          NA_id=1:length(x)
                                          NA_id=NA_id[!NA_id%in%queryHits(olap)]
                                          mcols(x)[NA_id,]=NA
                                          return(x)}))
    
  }
  saveRDS(JASPAR_out,paste('../downstream/output/JASPAR_motif_hg19_NME_',n,'_agnostic_merged_10k.rds',sep=''))
}
#subset to 10 kb
genomic_features=readRDS(genomic_features_file)
genomic_features=genomic_features$TSS
for (i in 1:3){
  JASPAR_motif=readRDS(paste('../downstream/output/JASPAR_motif_hg19_NME_',i,'_agnostic_merged.rds',sep=''))
  
  
  saveRDS(JASPAR_motif,paste('../downstream/output/JASPAR_motif_hg19_NME_',i,'_agnostic_merged_10k.rds',sep=''))
  JASPAR_motif=NA
  gc()
}
JASPAR_motif_sp=lapply(JASPAR_motif[1:210],NME_dNME_ken,GR_in=NME_in,stat_in='NME')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_NME_1_agnostic.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[211:420],NME_dNME_ken,GR_in=NME_in,stat_in='NME')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_NME_2_agnostic.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[421:630],NME_dNME_ken,GR_in=NME_in,stat_in='NME')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_NME_3_agnostic.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[1:210],NME_dNME_ken,GR_in=GR_merge,stat_in='dNME')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_dNME_1.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[211:420],NME_dNME_ken,GR_in=GR_merge,stat_in='dNME')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_dNME_2.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[421:630],NME_dNME_ken,GR_in=GR_merge,stat_in='dNME')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_dNME_3.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[1:210],NME_dNME_ken,GR_in=GR_merge,stat_in='dNME_pval')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_dNME_pval_1.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[211:420],NME_dNME_ken,GR_in=GR_merge,stat_in='dNME_pval')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_dNME_pval_2.rds')
JASPAR_motif_sp=NA
gc()
JASPAR_motif_sp=lapply(JASPAR_motif[421:630],NME_dNME_ken,GR_in=GR_merge,stat_in='dNME_pval')
saveRDS(JASPAR_motif_sp,'../downstream/output/JASPAR_motif_hg19_dNME_pval_3.rds')

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
