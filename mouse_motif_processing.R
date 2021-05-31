source('mainFunctions_sub.R')
# DNase vs control tissue-specific--------------------------------------------------------
NME_in=readRDS('../downstream/output/mouse_analysis/CPEL_output/NME_matrix_mouse_all_dedup_N2_all_regions.rds')
DNAase_in=readRDS('../downstream/input/mouse_analysis/motif_site_tissue_specific/DNase_mm10_peak_sample_250bp.rds')
names(DNAase_in)=gsub('C57BL/6 | tissue embryo| tissue male embryo |\\(| days\\)','',names(DNAase_in))
names(DNAase_in)=gsub('embryonic facial prominence','EFP', names(DNAase_in))
names(DNAase_in)=gsub('neural tube','NT', names(DNAase_in))
names(DNAase_in)=gsub(' ','', names(DNAase_in))
names(DNAase_in)=sub('1','-E1', names(DNAase_in))
#NME
colnames(mcols(NME_in))=sub('-all','',colnames(mcols(NME_in)))
mcols(NME_in)=mcols(NME_in)[,which(colnames(mcols(NME_in))%in% names(DNAase_in))]
DNAase_in=DNAase_in[which(names(DNAase_in) %in% colnames(mcols(NME_in)))]
JASPAR_motif=readRDS('../downstream/input/mouse_analysis/motif_site_tissue_specific/motif_JASPAR_mm10.rds')
for(ts in names(DNAase_in)){
  tt1=proc.time()[[3]]
  NME_in_ts=NME_in
  mcols(NME_in_ts)=mcols(NME_in_ts)[,which(colnames(mcols(NME_in))==ts)]
  NME_in_ts=subsetByOverlaps(NME_in_ts,DNAase_in[[ts]])
  colnames(mcols(NME_in_ts))='NME'
  NME_in_ts$Sample=ts
  NME_in_ts=NME_in_ts[!grepl('chrX|chrY',seqnames(NME_in_ts))]
  JASPAR_motif_DNase_sp=mclapply(JASPAR_motif,NME_dNME_ken,GR_in=NME_in_ts,stat_in='NME',mc.cores=24)
  
  saveRDS(JASPAR_motif_DNase_sp,paste('../downstream/output/mouse_analysis/mouse_motif_tissue_specific/JASPAR_motif_mm10_NME_',ts,'_agnostic_merged_DNase.rds'))
  JASPAR_motif_DNase_sp=NA
  # JASPAR_motif_control_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_control,stat_in='NME',mc.cores=24)
  # saveRDS(JASPAR_motif_control_sp,paste('../downstream/output/mouse_motif/mouse_motif_tissue_specific/JASPAR_motif_mm10_NME_',i,'_agnostic_merged_control.rds'))
  # JASPAR_motif_control_sp=NA
  
  proc.time()[[3]]-tt1
}
#MML
MML_in=readRDS('../downstream/output/mouse_analysis/CPEL_output/MML_matrix_mouse_all_dedup_N2_all_regions.rds')
colnames(mcols(MML_in))=sub('-all','',colnames(mcols(MML_in)))
mcols(MML_in)=mcols(MML_in)[,which(colnames(mcols(MML_in))%in% names(DNAase_in))]
DNAase_in=DNAase_in[which(names(DNAase_in) %in% colnames(mcols(MML_in)))]
JASPAR_motif=readRDS('../downstream/input/mouse_analysis/motif_site_tissue_specific/motif_JASPAR_mm10.rds')
for(ts in names(DNAase_in)){
  tt1=proc.time()[[3]]
  MML_in_ts=MML_in
  mcols(MML_in_ts)=mcols(MML_in_ts)[,which(colnames(mcols(MML_in))==ts)]
  MML_in_ts=subsetByOverlaps(MML_in_ts,DNAase_in[[ts]])
  colnames(mcols(MML_in_ts))='MML'
  MML_in_ts$Sample=ts
  MML_in_ts=MML_in_ts[!grepl('chrX|chrY',seqnames(MML_in_ts))]
  JASPAR_motif_DNase_sp=mclapply(JASPAR_motif,NME_dNME_ken,GR_in=MML_in_ts,stat_in='MML',mc.cores=24)
  
  saveRDS(JASPAR_motif_DNase_sp,paste('../downstream/output/mouse_analysis/mouse_motif_tissue_specific/JASPAR_motif_mm10_MML_',ts,'_agnostic_merged_DNase.rds'))
  JASPAR_motif_DNase_sp=NA
  # JASPAR_motif_control_sp=mclapply(JASPAR_motif[split_data==i],MML_dMML_ken,GR_in=MML_in_control,stat_in='NME',mc.cores=24)
  # saveRDS(JASPAR_motif_control_sp,paste('../downstream/output/mouse_motif/mouse_motif_tissue_specific/JASPAR_motif_mm10_NME_',i,'_agnostic_merged_control.rds'))
  # JASPAR_motif_control_sp=NA
  
  proc.time()[[3]]-tt1
}
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
MML_in=readRDS('../downstream/input/MML_agnostic_mouse_all_merged.rds')
MML_in=MML_in[MML_in$N>=2]
mcols(MML_in)=mcols(MML_in)[,c("Sample","MML")]
DNase=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp.rds')
control=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp_control.rds')
MML_in_DNase=subsetByOverlaps(MML_in,DNase,type='equal')
MML_in_control=subsetByOverlaps(MML_in,control,type='equal')
JASPAR_motif=readRDS('../downstream/input/motif_JASPAR_mm10.rds')
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
rm(MML_in)
for(i in 1:12){
  tt1=proc.time()[[3]]
  JASPAR_motif_DNase_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_DNase,stat_in='MML',mc.cores=10)
  
  saveRDS(JASPAR_motif_DNase_sp,paste('../downstream/output/mouse_motif/JASPAR_motif_mm10_MML_',i,'_agnostic_merged_DNase.rds'))
  JASPAR_motif_DNase_sp=NA
  JASPAR_motif_control_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_control,stat_in='MML',mc.cores=10)
  saveRDS(JASPAR_motif_control_sp,paste('../downstream/output/mouse_motif/JASPAR_motif_mm10_MML_',i,'_agnostic_merged_control.rds'))
  JASPAR_motif_control_sp=NA
  
  proc.time()[[3]]-tt1
}


