rm(list=ls())
source('mainFunctions_sub.R')

ChiP_motif_dir='../downstream/output/mouse_analysis/motif_analysis/motif_Chip_rds/'
region_in=readRDS(tissue_out_filtered_fn)
#Getting Bin enhancer
enhancer=readRDS(bin_enhancer_rds)
UC_merge=readRDS(UC_merge_max_loc_file)
#Getting target genes
region_in_enhancer=lapply(region_in,function(x) {
  olap=findOverlaps(convert_GR(x$region),enhancer)
  x=x[queryHits(olap)]
  x$gene=enhancer[subjectHits(olap)]$`Target Gene`
  return(x)
})

region_in_enhancer=lapply(region_in_enhancer,function(x){
  x=x[,list(region,tissue,cluster,gene,region_type,dMML_cor,dNME_cor)]
  UC_merge_ts=UC_merge[[unique(x$tissue)]]
  print(head(UC_merge_ts[x$region,'dNME_max_UC_pair_adj']))
  x$dNME_max_pair=UC_merge_ts[x$region,'dNME_max_UC_pair_adj']
  x$dMML_max_pair =UC_merge_ts[x$region,'dMML_max_UC_pair_adj']
  x$UC_max_pair=UC_merge_ts[x$region,'UC_max_UC_pair_adj']
  x$UC_max_time =UC_merge_ts[x$region,'UC_max_time_adj']
  x$UC_max_time=gsub(paste0(unique(x$tissue),'-|-all'),'',x$UC_max_time)
  return(x)
  
})

# #Mouse gene example
GO_out_all=readRDS(GO_01_enhancer_fn)

tissue_sel=names(region_in_enhancer)
for(ts in tissue_sel){
  select_top_GO_out=select_top_GO(GO_out_all$all,ts,ptcount=0,FDR_cutoff=0.1,FC_cutoff=1.5)
  ts_sel_genes=unique(unlist(strsplit(select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5]$genes,';')))
  enhancer_regions_ts=region_in_enhancer[[ts]][gene %in% ts_sel_genes][order(dNME_max_pair,decreasing=T)]
  enhancer_regions_ts$GO_terms=unlist(lapply(enhancer_regions_ts$gene,function(x) paste(select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5][grepl(x,genes)]$Term,collapse = ';')))
  enhancer_regions_ts$GO_ID=unlist(lapply(enhancer_regions_ts$gene,function(x) paste(gsub('GO:','',select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5][grepl(x,genes)]$`GO.ID`),collapse = ';')))
  region_in_enhancer[[ts]]=enhancer_regions_ts
  write.csv(enhancer_regions_ts,paste0(gene_example_dir,ts,'.csv'),row.names = F)
}
#EFP:before 2941, after 2140
#This file contains all regions we analzyed and which motif is their binding site
motif_locus_ken=readRDS(tissue_region_motif_all_regions_fn)

motif_enhancer_dNME=list()
motif_enhancer_dMML=list()
for(ts in names(region_in_enhancer)){
  enhancer_regions_ts=region_in_enhancer[[ts]]
  dMML_motif=motif_sig_Ken(ts,"dMML",motif_locus_ken,enhancer_regions_ts,dir_in_Ken=Ken_motif_folder)
  dNME_motif=motif_sig_Ken(ts,"dNME",motif_locus_ken,enhancer_regions_ts,dir_in_Ken=Ken_motif_folder)
 #motif_region relationship
  motif_enhancer_dNME[[ts]]=dNME_motif$motif_locus
  motif_enhancer_dMML[[ts]]=dMML_motif$motif_locus
  #Region-motif relationship 
  enhancer_regions_ts$dMML_motif=dMML_motif$motif_locus_dt[match(enhancer_regions_ts$region,region)]$motif
  enhancer_regions_ts$dNME_motif=dNME_motif$motif_locus_dt[match(enhancer_regions_ts$region,region)]$motif

  region_in_enhancer[[ts]]=enhancer_regions_ts
  write.csv(enhancer_regions_ts,paste0(region_motif_dir,ts,'.csv'))
}
saveRDS(region_in_enhancer,enhancer_region_fn)


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
# Selecct motif in factor book --------------------------------------------

# Check mosue ChiP-seq result ---------------------------------------------
enhancer_regions_motif_dNME_all=readRDS(enhancer_motif_all_dNME_fn)
enhancer_regions_motif_dMML_all=readRDS(enhancer_motif_all_dMML_fn)

#embyro Heart
factor_in_heart_embyro=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_embyro_heart_allTF.bed',skip = 1,sep='\t')
colnames(factor_in_heart_embyro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_heart_embyro=factor_in_heart_embyro[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_heart_embyro=factor_olap('heart',factor_in_heart_embyro,Ken_motif_folder,stage="embyro_heart")
ChIP_olap_heart_embyro=ChIP_olap(factor_olap_heart_embyro$factor_in_dNME,factor_olap_heart_embyro$factor_in_dMML,
                                 enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(enhancer_regions_heart_embyro$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_consistent_dNME.rds'))


#Adult heart
factor_in_heart=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_heart_all_TF.bed',skip = 1,sep='\t')
colnames(factor_in_heart)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_heart=factor_in_heart[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_heart_adult=factor_olap('heart',factor_in_heart,Ken_motif_folder,stage="adult")
ChIP_olap_heart_adult=ChIP_olap(factor_olap_heart_adult$factor_in_dNME,factor_olap_heart_adult$factor_in_dMML,
                                 enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(ChIP_olap_heart_adult$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_adult_consistent_dNME.rds'))

#embyro all:
factor_in_embyro=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_embyro_all_TF.bed',skip = 1,sep='\t')
colnames(factor_in_embyro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_embyro=factor_in_embyro[,list(seqnames,start,end,metadata,log10qval)]

#Heart embyro all
factor_olap_heart_embyro_all=factor_olap('heart',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_heart_embyro_all=ChIP_olap(factor_olap_heart_embyro_all$factor_in_dNME,factor_olap_heart_embyro_all$factor_in_dMML,
                                enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(ChIP_olap_heart_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_heart_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_all_consistent_dMML.rds'))
#Forebrain embyro all
factor_olap_forebrain_embyro_all=factor_olap('forebrain',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_forebrain_embyro_all=ChIP_olap(factor_olap_forebrain_embyro_all$factor_in_dNME,factor_olap_forebrain_embyro_all$factor_in_dMML,
                                     enhancer_regions_motif_dNME_all$forebrain,enhancer_regions_motif_dMML_all$forebrain)
saveRDS(ChIP_olap_forebrain_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_forebrain_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_embyro_all_consistent_dMML.rds'))
#EFP embyro all
factor_olap_EFP_embyro_all=factor_olap('EFP',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_EFP_embyro_all=ChIP_olap(factor_olap_EFP_embyro_all$factor_in_dNME,factor_olap_EFP_embyro_all$factor_in_dMML,
                                         enhancer_regions_motif_dNME_all$EFP,enhancer_regions_motif_dMML_all$EFP)
saveRDS(ChIP_olap_EFP_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_EFP_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_EFP_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_EFP_embyro_all_consistent_dMML.rds'))
#limb embyro all
factor_olap_limb_embyro_all=factor_olap('limb',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_limb_embyro_all=ChIP_olap(factor_olap_limb_embyro_all$factor_in_dNME,factor_olap_limb_embyro_all$factor_in_dMML,
                                   enhancer_regions_motif_dNME_all$limb,enhancer_regions_motif_dMML_all$limb)
saveRDS(ChIP_olap_limb_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_limb_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_limb_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_limb_embyro_all_consistent_dMML.rds'))
#Adult neuro
factor_in_neuro=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_neuro_all_TF.bed',skip = 1,sep='\t')
colnames(factor_in_neuro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_neuro=factor_in_neuro[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_forebrain_adult=factor_olap('forebrain',factor_in_embyro,Ken_motif_folder,stage="adult")
ChIP_olap_forebrain_adult=ChIP_olap(factor_olap_forebrain_adult$factor_in_dNME,factor_olap_forebrain_adult$factor_in_dMML,
                                         enhancer_regions_motif_dNME_all$forebrain,enhancer_regions_motif_dMML_all$forebrain)
saveRDS(ChIP_olap_forebrain_adult$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_adult_consistent_dNME.rds'))
saveRDS(ChIP_olap_forebrain_adult$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_adult_consistent_dMML.rds'))

# #Liver
# #Embyro
# factor_in_liver=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_embyro_liver.bed',skip = 1,sep='\t')
# colnames(factor_in_liver)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
# factor_in_liver=factor_in_liver[,list(seqnames,start,end,metadata,log10qval)]
# factor_in_liver$metadata=gsub('%20','_',factor_in_liver$metadata)
# factor_in_liver$metadata=gsub('%3','_',factor_in_liver$metadata)
# forebrain_ken_dNME=fread(paste0(Ken_motif_folder,'liver_OR_residual_dNME.csv'))
# forebrain_ken_dMML = fread(paste0(Ken_motif_folder,'liver_OR_residual_dMML.csv'))
# factor_in_liver_dNME=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dNME$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# factor_in_liver_dMML=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dMML$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# 
# for(motif in gsub('.*_','',liver_ken_dNME$motif)){
#   motif_ChIP=factor_in_liver_dNME[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_liver_embyro_dNME_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }
# for(motif in gsub('.*_','',liver_ken_dMML$motif)){
#   motif_ChIP=factor_in_liver_dMML[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dMML/mm10_liver_embyro_dMML_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }
# #Adult
# factor_in_liver=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_liver.bed',skip = 1,sep='\t')
# colnames(factor_in_liver)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
# factor_in_liver=factor_in_liver[,list(seqnames,start,end,metadata,log10qval)]
# factor_in_liver$metadata=gsub('%20','_',factor_in_liver$metadata)
# factor_in_liver$metadata=gsub('%3','_',factor_in_liver$metadata)
# liver_ken_dNME=fread('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_result_20210406/liver_OR_residual_dNME.csv')
# liver_ken_dMML = fread('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_result_20210406/liver_OR_residual_dMML.csv')
# factor_in_liver_dNME=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dNME$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# factor_in_liver_dMML=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dMML$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# 
# for(motif in gsub('.*_','',liver_ken_dNME$motif)){
#   motif_ChIP=factor_in_liver_dNME[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_liver_adult_dNME_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }
# for(motif in gsub('.*_','',liver_ken_dMML$motif)){
#   motif_ChIP=factor_in_liver_dMML[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dMML/mm10_forebrain_adult_dMML_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }


# Find examples by looking at overlap between motif and ChIP data -----------
enhancer_regions_motif_dNME_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dNME_all.rds')
enhancer_regions_motif_dMML_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dMML_all.rds')
#Heart 
factor_in_embyro_heart_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_embyro_heart_dNME.rds')))
factor_in_heart_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_heart_dNME.rds')))
heart_dNME_pubmed=enhancer_regions_motif_dNME_all$heart
heart_dNME_pubmed=heart_dNME_pubmed[,list(region,tissue,cluster,gene,UC_max_time,region_type,PMID,dNME_max_pair,dMML_max_pair,UC_max_pair,dMML_motif,dNME_motif)]
olap=findOverlaps(convert_GR(heart_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_embyro_heart_dNME))
heart_dNME_pubmed_motif=heart_dNME_pubmed[queryHits(olap)][order(dNME_max_pair,decreasing=T)]

olap=findOverlaps(convert_GR(heart_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_heart_dNME))
heart_dNME_pubmed_motif=heart_dNME_pubmed[queryHits(olap)][order(dNME_max_pair,decreasing=T)]
#Forebrain
factor_in_embyro_forebrain_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_embyro_forebrain_dNME.rds')))
factor_in_forebrain_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_adult_forebrain_dNME.rds')))
forebrain_dNME_pubmed=enhancer_regions_motif_dNME_all$forebrain
forebrain_dNME_pubmed=forebrain_dNME_pubmed[,list(region,tissue,cluster,gene,UC_max_time,region_type,PMID,dNME_max_pair,dMML_max_pair,UC_max_pair,dMML_motif,dNME_motif)]
olap=findOverlaps(convert_GR(forebrain_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_embyro_forebrain_dNME))
forebrain_dNME_pubmed_motif=forebrain_dNME_pubmed[queryHits(olap)]
forebrain_dNME_pubmed_motif$motif_ChIP=factor_in_embyro_forebrain_dNME[subjectHits(olap)]$metadata
forebrain_dNME_pubmed_motif$motif_in_ChIP=forebrain_dNME_pubmed_motif[,list(motif_in_ChIP=grepl(gsub(';','|',dNME_motif),motif_ChIP,ignore.case = T)),by=list(region,motif_ChIP)]$motif_in_ChIP

olap=findOverlaps(convert_GR(forebrain_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_forebrain_dNME))
forebrain_dNME_pubmed_motif=forebrain_dNME_pubmed[queryHits(olap)]
forebrain_dNME_pubmed_motif$motif_ChIP=factor_in_forebrain_dNME[subjectHits(olap)]$metadata
forebrain_dNME_pubmed_motif$motif_in_ChIP=forebrain_dNME_pubmed_motif[,list(motif_in_ChIP=grepl(gsub(';','|',dNME_motif),motif_ChIP,ignore.case = T)),by=list(region,motif_ChIP)]$motif_in_ChIP
# archived ----------------------------------------------------------------



liver_regions=fread('../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dNME.csv')
olap=findOverlaps(convert_GR(liver_regions$region),GRanges(seqname=factor_in_liver_dNME$seqnames,IRanges(start=factor_in_liver_dNME$start,end=factor_in_liver_dNME$end)))

write.csv(liver_regions[queryHits(olap)],'../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dNME_ChiP.csv')

liver_regions=fread('../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dMML.csv')
olap=findOverlaps(convert_GR(liver_regions$region),GRanges(seqname=factor_in_liver_dMML$seqnames,IRanges(start=factor_in_liver_dMML$start,end=factor_in_liver_dMML$end)))

write.csv(liver_regions[queryHits(olap)],'../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dMML_ChiP.csv')


#Subset the motif
#Ken motif
heart_Ken_binding_site=readRDS('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_binding_site/tissue_region_motif_all.rds')
heart_enhancer_GO_motif=

olap_motif=findOverlaps(convert_GR(region_in_enhancer_GO$heart$region),factor_in_GATA4_gr)

unique(region_in_enhancer$heart[queryHits(olap_motif)])[order(dNME_max_pair,decreasing = T)]
unique(region_in_enhancer_GO$heart[queryHits(olap_motif)][grepl("heart|cardi|aort|ventricu|vascula|muscle|actin",GO_result),
                                   list(region,cluster,dNME_max_pair,dNME_max_time,UC_max_pair,dMML_max_pair,gene),with=T])[order(dNME_max_pair,decreasing=T)]




#Examples for correlation
both_example = data.frame(value=as.numeric(UC_merge$heart["chr4:97739493-97739742",!grepl("max",colnames(UC_merge$heart))]),
                          stage_stat=colnames(UC_merge$heart["chr4:97739493-97739742",!grepl("max",colnames(UC_merge$heart))]))

both_example$stat=gsub('-.*','',both_example$stage_stat)
both_example$stage=sub('dMML-|dNME-|UC-','',both_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(both_example[both_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()

dNME_only_example = data.frame(value=as.numeric(UC_merge$heart["chr10:80117344-80117593",!grepl("max",colnames(UC_merge$heart))]),
                               stage_stat=colnames(UC_merge$heart["chr10:80117344-80117593",!grepl("max",colnames(UC_merge$heart))]))

dNME_only_example$stat=gsub('-.*','',dNME_only_example$stage_stat)
dNME_only_example$stage=sub('dMML-|dNME-|UC-','',dNME_only_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(dNME_only_example[dNME_only_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()

dMML_only_example = data.frame(value=as.numeric(UC_merge$heart["chr10:62313679-62313928",!grepl("max",colnames(UC_merge$heart))]),
                               stage_stat=colnames(UC_merge$heart["chr10:62313679-62313928",!grepl("max",colnames(UC_merge$heart))]))

dMML_only_example$stat=gsub('-.*','',dMML_only_example$stage_stat)
dMML_only_example$stage=sub('dMML-|dNME-|UC-','',dMML_only_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(dMML_only_example[dMML_only_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()


neither_example = data.frame(value=as.numeric(UC_merge$heart["chr1:151764174-151764423",!grepl("max",colnames(UC_merge$heart))]),
                             stage_stat=colnames(UC_merge$heart["chr1:151764174-151764423",!grepl("max",colnames(UC_merge$heart))]))

neither_example$stat=gsub('-.*','',neither_example$stage_stat)
neither_example$stage=sub('dMML-|dNME-|UC-','',neither_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(neither_example[neither_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()
#IGV examples

neuro_motif_ChIP=fread('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_forebrain_adult_selected_dNME_ChiPatlas.bed')
forebrain_motif=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/forebrain_dNME_Pubmed_annotated_region_only.csv')
olap=findOverlaps(convert_GR(forebrain_motif$region),makeGRangesFromDataFrame(neuro_motif_ChIP,seqnames.field = "V1",start.field = "V2",end.field = "V3"))
forebrain_motif=forebrain_motif[unique(queryHits(olap))][order(dNME_max_UC_pair,decreasing = T)]
write.csv(forebrain_motif,'../downstream/output/mouse_analysis/motif_analysis/region_motif/forebrain_dNME_Pubmed_annotated_region_only_ChIp_olap.csv')

heart_motif_ChIP=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_heart_GATA4.bed',skip = 1,sep='\t')
heart_motif_ChIP=fread('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_heart_embyro_selected_dNME_ChiPatlas.bed',sep=' ')
heart_motif=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/heart_motif_region_only.csv')
olap=findOverlaps(convert_GR(heart_motif$region),makeGRangesFromDataFrame(heart_motif_ChIP,seqnames.field = "V1",start.field = "V2",end.field = "V3"))
heart_motif=heart_motif[unique(queryHits(olap))][order(dNME_max_UC_pair,decreasing = T)]
write.csv(heart_motif,'../downstream/output/mouse_analysis/motif_analysis/region_motif/heart_dNME_Pubmed_annotated_region_only_ChIp_olap2.csv')


