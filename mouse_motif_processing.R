source('mainFunctions_sub.R')

# DNase vs control --------------------------------------------------------
output_dir='../downstream/output/mouse_motif_Ken/'
NME_in=readRDS(NME_matrix_file)
NME_in=convert_GR(NME_in,direction='DT')
NME_in=melt.data.table(NME_in,id.vars = 'region',variable.name = 'Sample',value.name = 'NME')
NME_in=NME_in[!is.na(NME_in$NME)]
NME_in$Sample=gsub('\\.','-',NME_in$Sample)
NME_in_gr=convert_GR(NME_in$region,direction='GR')
mcols(NME_in_gr)=NME_in[,list(Sample,NME)]
NME_in=NME_in_gr
rm(NME_in_gr)
DNase=readRDS(DNase_mm10_file)
control=readRDS(control_mm10_file)
NME_in_DNase=subsetByOverlaps(NME_in,DNase,type='equal')
NME_in_control=subsetByOverlaps(NME_in,control,type='equal')
#This is from ken need large memory node?
JASPAR_motif=readRDS(JASPAR_motif_mm10_file)
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
rm(NME_in)
#Need lrgmem
for(i in 1:12){
  tt1=proc.time()[[3]]
  JASPAR_motif_DNase_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_DNase,stat_in='NME',mc.cores=24)
  
  saveRDS(JASPAR_motif_DNase_sp,paste0(output_dir,'JASPAR_motif_mm10_NME_',i,'_agnostic_merged_DNase.rds'))
  JASPAR_motif_DNase_sp=NA
  JASPAR_motif_control_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_control,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_control_sp,paste0(output_dir,'JASPAR_motif_mm10_NME_',i,'_agnostic_merged_control.rds'))
  JASPAR_motif_control_sp=NA
  
  proc.time()[[3]]-tt1
}
#MML
MML_in=readRDS(MML_matrix_file)
MML_in=convert_GR(MML_in,direction='DT')
MML_in=melt.data.table(MML_in,id.vars = 'region',variable.name = 'Sample',value.name = 'MML')
MML_in=MML_in[!is.na(MML_in$MML)]
MML_in$Sample=gsub('\\.','-',MML_in$Sample)
MML_in_gr=convert_GR(MML_in$region,direction='GR')
mcols(MML_in_gr)=MML_in[,list(Sample,MML)]
MML_in=MML_in_gr
rm(MML_in_gr)
DNase=readRDS(DNase_mm10_file)
DNase=readRDS(DNase_mm10_file)
control=readRDS(control_mm10_file)
MML_in_DNase=subsetByOverlaps(MML_in,DNase,type='equal')
MML_in_control=subsetByOverlaps(MML_in,control,type='equal')
JASPAR_motif=readRDS('../downstream/input/motif_JASPAR_mm10.rds')
split_data=cut(1:length(JASPAR_motif),breaks=12,label=FALSE)
rm(MML_in)
for(i in 1:12){
  tt1=proc.time()[[3]]
  JASPAR_motif_DNase_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_DNase,stat_in='MML',mc.cores=10)
  
  saveRDS(JASPAR_motif_DNase_sp,paste0(output_dir,'JASPAR_motif_mm10_MML_',i,'_agnostic_merged_DNase.rds'))
  JASPAR_motif_DNase_sp=NA
  JASPAR_motif_control_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_control,stat_in='MML',mc.cores=10)
  saveRDS(JASPAR_motif_control_sp,paste0(output_dir,'JASPAR_motif_mm10_MML_',i,'_agnostic_merged_control.rds'))
  JASPAR_motif_control_sp=NA
  
  proc.time()[[3]]-tt1
}


