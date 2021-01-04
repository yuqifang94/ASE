source('mainFunctions_sub.R')

# DNase vs control --------------------------------------------------------
NME_in=readRDS('../downstream/output/NME_agnostic_mouse_all_merged.rds')
NME_in=NME_in[NME_in$N>=2]
mcols(NME_in)=mcols(NME_in)[,c("Sample","NME")]
DNase=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp.rds')
control=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp_control.rds')
NME_in_DNase=subsetByOverlaps(NME_in,DNase,type='equal')
NME_in_control=subsetByOverlaps(NME_in,control,type='equal')
JASPAR_motif=readRDS('../downstream/input/motif_JASPAR_mm10.rds')
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
rm(NME_in)
for(i in 1:12){
  tt1=proc.time()[[3]]
  JASPAR_motif_DNase_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_DNase,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_DNase_sp,paste('../downstream/output/mouse_motif/JASPAR_motif_mm10_NME_',i,'_agnostic_merged_DNase.rds'))
  JASPAR_motif_DNase_sp=NA
  JASPAR_motif_control_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_control,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_control_sp,paste('../downstream/output/mouse_motif/JASPAR_motif_mm10_NME_',i,'_agnostic_merged_control.rds'))
  JASPAR_motif_control_sp=NA
  
  proc.time()[[3]]-tt1
}

# assign region -----------------------------------------------------------
nme_cor <- readRDS('../downstream/input/dnmecor.rds')
mml_cor <- readRDS('../downstream/input/dmmlcor.rds')
cor_out=data.table()
for(ts in names(nme_cor)){
  cor_in=data.table(region_all=names(nme_cor[[ts]]),
                    nme_cor=nme_cor[[ts]])
  
  cor_in$mml_cor=mml_cor[[ts]][cor_in$region_all]
  cor_in$tissue=ts
  cor_out=rbind(cor_out,cor_in)
}
cor_out$region_type=NA
cor_out[nme_cor>=0.7&mml_cor<0.7]$region_type="dNME_only"
cor_out[nme_cor>=0.7&mml_cor>=0.7]$region_type="dNME_dMML"
cor_out[nme_cor<0.7&mml_cor>=0.7]$region_type="dMML_only"
cor_out[nme_cor<0.7&mml_cor<0.7]$region_type="non_dNME_non_dMML"
