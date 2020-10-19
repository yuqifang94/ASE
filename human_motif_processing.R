# Ken motif processing ----------------------------------------------------
source('mainFunctions_sub.R')
#hg19

GR_merge=readRDS(GR_merge_file)
NME_in=readRDS(NME_agnostic_file)
mcols(NME_in)=mcols(NME_in)[,-1]

JASPAR_motif=readRDS('../downstream/output/motif_JASPAR_hg19.rds')

# Use mean NME for allele-specific analysis -------------------------------
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge,stat_in='NME',mc.cores=8)
  saveRDS(JASPAR_motif_sp,paste('../downstream/output/human_motif/JASPAR_motif_hg19_NME_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}


# Use allele-agnostic NME -------------------------------
