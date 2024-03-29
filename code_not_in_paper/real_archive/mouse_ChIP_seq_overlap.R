source('mainFunctions_sub.R')
# Make Bedfile for Ken's motif binding site -------------------------------
Ken_motif_locus=readRDS(tissue_region_motif_all_regions_fn)
for(stat_in in c("dMML","dNME")){
    for(tissue in names(Ken_motif_locus)){
      motif_sig=fread(paste0(Ken_motif_folder,tissue,'_OR_residual_',stat_in,'.csv'))
      bed_file_path=paste0(motif_locus_bed_dir,stat_in,'/',tissue,'/')
      motif_locus_in=Ken_motif_locus[[tissue]][motif_sig$motif]
      ifelse(!dir.exists(file.path(bed_file_path)), dir.create(file.path(bed_file_path)), FALSE)
      lapply(names(motif_locus_in),function(x) {export.bed(motif_locus_in[[x]],paste0(bed_file_path,
                                                                                      gsub("::","_",gsub('.*_','',x)),'.bed'))})
    
    }
}
# Seleccted from Chip-atlas ---------------------------------------------
#embyro_heart: mm10->Embyro->Embyronic heart
#adult_heart: mm10->Cardiovascular->all
#embyro_all: mm10->Embyro->all
#Adult Neuro:  mm10->Neural->all

# Check mosue ChiP-seq result ---------------------------------------------
enhancer_regions_motif_dNME_all=readRDS(enhancer_motif_all_dNME_fn)
enhancer_regions_motif_dMML_all=readRDS(enhancer_motif_all_dMML_fn)

#embyro Heart
factor_in_heart_embyro=fread(paste0(chip_atlas_dir,'mm10_embyro_heart_allTF.bed'),skip = 1,sep='\t')
colnames(factor_in_heart_embyro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_heart_embyro=factor_in_heart_embyro[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_heart_embyro=factor_olap('heart',factor_in_heart_embyro,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="embyro_heart")
ChIP_olap_heart_embyro=ChIP_olap(factor_olap_heart_embyro$factor_in_dNME,factor_olap_heart_embyro$factor_in_dMML,
                                 enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(ChIP_olap_heart_embyro,paste0(ChiP_motif_dir,'ChIP_olap_heart_embyro_heart_all.rds'))


#Adult heart
factor_in_heart=fread(paste0(chip_atlas_dir,'mm10_heart_all_TF.bed'),skip = 1,sep='\t')
colnames(factor_in_heart)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_heart=factor_in_heart[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_heart_adult=factor_olap('heart',factor_in_heart,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="adult")
ChIP_olap_heart_adult=ChIP_olap(factor_olap_heart_adult$factor_in_dNME,factor_olap_heart_adult$factor_in_dMML,
                                 enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(ChIP_olap_heart_adult,paste0(ChiP_motif_dir,'ChIP_olap_heart_adult_heart_all.rds'))

#embyro all:
factor_in_embyro=fread(paste0(chip_atlas_dir,'mm10_embyro_all_TF.bed'),skip = 1,sep='\t')
colnames(factor_in_embyro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_embyro=factor_in_embyro[,list(seqnames,start,end,metadata,log10qval)]

#Heart embyro all
factor_olap_heart_embyro_all=factor_olap('heart',factor_in_embyro,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="embyro_all")
ChIP_olap_heart_embyro_all=ChIP_olap(factor_olap_heart_embyro_all$factor_in_dNME,factor_olap_heart_embyro_all$factor_in_dMML,
                                enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(ChIP_olap_heart_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_heart_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_all_consistent_dMML.rds'))
#Example filtering
ChIP_olap_heart_embyro_all$dNME_region[!is.na(predict_in_ChIP)&predict_in_ChIP!=""&region_type=='NME only']
#Forebrain embyro all
factor_olap_forebrain_embyro_all=factor_olap('forebrain',factor_in_embyro,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="embyro_all")
ChIP_olap_forebrain_embyro_all=ChIP_olap(factor_olap_forebrain_embyro_all$factor_in_dNME,factor_olap_forebrain_embyro_all$factor_in_dMML,
                                     enhancer_regions_motif_dNME_all$forebrain,enhancer_regions_motif_dMML_all$forebrain)
saveRDS(ChIP_olap_forebrain_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_forebrain_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_embyro_all_consistent_dMML.rds'))
#EFP embyro all
factor_olap_EFP_embyro_all=factor_olap('EFP',factor_in_embyro,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="embyro_all")
ChIP_olap_EFP_embyro_all=ChIP_olap(factor_olap_EFP_embyro_all$factor_in_dNME,factor_olap_EFP_embyro_all$factor_in_dMML,
                                         enhancer_regions_motif_dNME_all$EFP,enhancer_regions_motif_dMML_all$EFP)
saveRDS(ChIP_olap_EFP_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_EFP_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_EFP_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_EFP_embyro_all_consistent_dMML.rds'))
#limb embyro all
factor_olap_limb_embyro_all=factor_olap('limb',factor_in_embyro,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="embyro_all")
ChIP_olap_limb_embyro_all=ChIP_olap(factor_olap_limb_embyro_all$factor_in_dNME,factor_olap_limb_embyro_all$factor_in_dMML,
                                   enhancer_regions_motif_dNME_all$limb,enhancer_regions_motif_dMML_all$limb)
saveRDS(ChIP_olap_limb_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_limb_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_limb_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_limb_embyro_all_consistent_dMML.rds'))
#Adult neuro
factor_in_neuro=fread(paste0(chip_atlas_dir,'mm10_neuro_all_TF.bed'),skip = 1,sep='\t')
colnames(factor_in_neuro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_neuro=factor_in_neuro[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_forebrain_adult=factor_olap('forebrain',factor_in_embyro,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="adult")
ChIP_olap_forebrain_adult=ChIP_olap(factor_olap_forebrain_adult$factor_in_dNME,factor_olap_forebrain_adult$factor_in_dMML,
                                         enhancer_regions_motif_dNME_all$forebrain,enhancer_regions_motif_dMML_all$forebrain)
saveRDS(ChIP_olap_forebrain_adult$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_adult_consistent_dNME.rds'))
saveRDS(ChIP_olap_forebrain_adult$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_adult_consistent_dMML.rds'))


