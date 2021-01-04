rm(list=ls())
source("mainFunctions_sub.R")

#Define ggplot theme

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
#Density analysis
GR_merge=readRDS(GR_merge_file)
#####Subsetting by DNase region will not give us enough power to do it

###Reading in data
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
# motif_gene=fastDoCall('c',lapply(dir('../downstream/input/JASPAR_out/'),function(x){
#   cat(x,'\n')
#   readRDS(paste0('../downstream/input/JASPAR_out/',x))
#   
# }))
# saveRDS(motif_gene,motif_gene_file)
motif_gene <- readRDS(motif_gene_file)#See motif_break_array.R, default setting
#All regions
#NME
motif_dir=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="NME")
colnames(motif_dir)[c(1,6,7)]=c('TF','Pvalue','Proportion')
motif_dir$FDR=p.adjust(motif_dir$Pvalue,method="BH")
write.csv(motif_dir[FDR<=0.1&Proportion>0.5,list(TF,Proportion,Pvalue,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
write.csv(motif_dir[FDR<=0.1&Proportion<0.5,list(TF,Proportion=1-Proportion,`Pvalue`,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/table2_motif_prefer_low_NME.csv')

#MML
motif_dir=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="MML")
colnames(motif_dir)[c(1,6,7)]=c('TF','Pvalue','Proportion')
motif_dir$FDR=p.adjust(motif_dir$Pvalue,method="BH")
write.csv(motif_dir[FDR<=0.1&Proportion<0.5,list(TF,Proportion=1-Proportion,Pvalue,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/tablS1_motif_prefer_low_MML.csv')
write.csv(motif_dir[FDR<=0.1&Proportion>0.5,list(TF,Proportion,Pvalue,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/tableS2_motif_prefer_high_MML.csv')

# DNase analysis currently not in use -------------------------------------
# #DNase
# #NME
# DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
# motif_dir=direction_calc_enriched_subj(motif_gene,subsetByOverlaps(variant_HetCpG_meta,DNase),
#                                        unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="NME")
# colnames(motif_dir)[c(1,6,7)]=c('TF','Pvalue','Proportion')
# motif_dir$FDR=p.adjust(motif_dir$Pvalue,method="BH")
# write.csv(motif_dir[FDR<=0.1&Proportion<0.5,list(TF,Proportion=1-Proportion,Pvalue,FDR)], row.names =F,
#           '../downstream/output/graphs/motif_preference_table/DNase/table1_motif_prefer_low_NME_DNase.csv')
# write.csv(motif_dir[FDR<=0.1&Proportion>0.5,list(TF,Proportion,Pvalue,FDR)], row.names =F,
#           '../downstream/output/graphs/motif_preference_table/DNase/table2_motif_prefer_high_NME_DNase.csv')
# #MML
# motif_dir=direction_calc_enriched_subj(motif_gene,subsetByOverlaps(variant_HetCpG_meta,DNase),
#                                        unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="MML")
# colnames(motif_dir)[c(1,6,7)]=c('TF','Pvalue','Proportion')
# motif_dir$FDR=p.adjust(motif_dir$Pvalue,method="BH")
# write.csv(motif_dir[FDR<=0.1&Proportion<0.5,list(TF,Proportion=1-Proportion,Pvalue,FDR)], row.names =F,
#           '../downstream/output/graphs/motif_preference_table/DNase/table3_motif_prefer_low_MML_DNase.csv')
# write.csv(motif_dir[FDR<=0.1&Proportion>0.5,list(TF,Proportion,Pvalue,FDR)], row.names =F,
#           '../downstream/output/graphs/motif_preference_table/DNase/table4_motif_prefer_high_MML_DNase.csv')



# get OMIM annotation of those high NME motif only or high MML only -------
low_MML_motif=fread('../downstream/output/graphs/motif_preference_table/All_regions/tablS1_motif_prefer_low_MML.csv')
high_NME_motif=fread('../downstream/output/graphs/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
OMIM=fread('../downstream/input/genemap2.txt',skip=3)
OMIM$`Gene Symbols`=strsplit(as.character(OMIM$`Gene Symbols`),', ')
OMIM=OMIM[Phenotypes!=""]
motif_ent_OMIM=OMIM_annotation(high_NME_motif,OMIM)
motif_ent_only_OMIM=OMIM_annotation(high_NME_motif[!(TF%in%low_MML_motif$TF)],OMIM)
motif_low_MML_only_OMIM=OMIM_annotation(low_MML_motif[!(TF%in%high_NME_motif$TF)],OMIM)
motif_shared_OMIM=OMIM_annotation(low_MML_motif[TF%in%high_NME_motif$TF],OMIM)
write.csv(motif_low_MML_only_OMIM, row.names =F,
          "../downstream/output/graphs/motif_preference_table/All_regions/tableX_motif_prefer_low_MML_only.csv")
write.csv(motif_ent_only_OMIM, row.names =F,
          "../downstream/output/graphs/motif_preference_table/All_regions/table3_motif_prefer_high_NME_only.csv")
write.csv(motif_shared_OMIM, row.names =F,
          "../downstream/output/graphs/motif_preference_table/tableXmotif_shared_only.csv")
# write.csv(unlist(strsplit(gsub('\\(var.3\\)','',gsub('\\(var.2\\)','',motif_low_MML_only_OMIM$TF)),"::")),
#           "../downstream/output/motif_prefer_low_MML_only.csv")
# write.csv(unlist(strsplit(gsub('\\(var.3\\)','',gsub('\\(var.2\\)','',motif_ent_only_OMIM$TF)),"::")),
#           "../downstream/output/motif_prefer_high_NME_only.csv")
# write.csv(unlist(strsplit(gsub('\\(var.3\\)','',gsub('\\(var.2\\)','',motif_shared_OMIM$TF)),"::")),
#           "../downstream/output/motif_shared.csv")

#Ken's analysis preprocessing is in human_motif_processing.R


# TBD  --------------------------------------------------------------------

#find gnomic features enriched 
#OR: features vs non-feature, target vs non-targe, in dNME-ASM
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene <- readRDS(motif_gene_file)#See motif_break_array.R, default setting
genomic_features=readRDS(genomic_features_file)
selected_features=c("CpG island","CpG shore","CpG shelf","CpG open sea","gene body","exon","intron","intergenic","promoter","TSS")
motif_prefer_high_NME=fread('../downstream/output/graphs/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
motif_gene_high_NME=motif_gene[motif_gene$geneSymbol%in% motif_prefer_high_NME$TF]
variant_HetCpG_meta_dNME_ASM=variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff]
variant_HetCpG_meta_dNME_ASM$region=paste0(seqnames(variant_HetCpG_meta_dNME_ASM),":",start(variant_HetCpG_meta_dNME_ASM))
motif_gene_high_NME$region=paste0(seqnames(motif_gene_high_NME),':',start(motif_gene_high_NME))
# variant_HetCpG_meta_dNME_ASM_motif=variant_HetCpG_meta_dNME_ASM[variant_HetCpG_meta_dNME_ASM$region%in%motif_gene_high_NME_region]
# dNME_ASM_non_motif=variant_HetCpG_meta_dNME_ASM[!variant_HetCpG_meta_dNME_ASM$region%in%motif_gene_high_NME_region]
#Extract motif binding regions
ah = AnnotationHub()
chromHMM_state=ENCODE_to_sample(unique(variant_HetCpG_meta$Sample))
chromHMM_list=list()
suppressMessages({for (sp in chromHMM_state$sample[!is.na(chromHMM_state$ENCODE)]){
  ah_num=names(query(ah, c("chromhmmSegmentations", chromHMM_state$ENCODE[chromHMM_state$sample==sp])))
  chromHMM_list[[sp]]=ah[[ah_num]]
  # }
}
})
cont_out_df=data.table()
chromHMM_motif_all_TF=data.table()
genomic_features_OR_out=data.table()
OR_cont_table=list()
for(TF in unique(motif_gene_high_NME$geneSymbol)){
  OR_cont_table[[TF]]=data.table()
  #ALT-REF
  motif_gene_high_NME_TF=motif_gene_high_NME[motif_gene_high_NME$geneSymbol==TF]
  variant_dNME_ASM_TF=variant_HetCpG_meta_dNME_ASM
  variant_dNME_ASM_TF$allele_diff=0
  olap=findOverlaps(variant_dNME_ASM_TF,motif_gene_high_NME_TF)
  variant_dNME_ASM_TF[queryHits(olap)]$allele_diff=motif_gene_high_NME_TF[subjectHits(olap)]$alleleDiff  

  variant_dNME_ASM_TF$allele_diff_NME=variant_dNME_ASM_TF$altNME-variant_dNME_ASM_TF$refNME
  #Fit the function
  variant_dNME_ASM_TF$ASM="No"
  variant_dNME_ASM_TF$ASM[sign(variant_dNME_ASM_TF$allele_diff_NME)==sign(variant_dNME_ASM_TF$allele_diff)]="Yes"
  genomic_features_OR=data.table()
  for(ft in selected_features){
    OR_feature=NULL
    OR_feature=testEnrichmentFeature_stat(variant_dNME_ASM_TF,genomic_features[[ft]],output_ft=1)
    if(!is.null(OR_feature)){
      OR_feature[[1]]$feature=ft
      OR_cont_table[[TF]]=rbind(OR_cont_table[[TF]],OR_feature[[1]])
      OR_feature=OR_feature[[2]]
    genomic_features_OR=rbind(genomic_features_OR,
                              data.table(feature=ft,OR_feature=OR_feature$estimate,p_value=OR_feature$p.value,
                                                             lower_CI=OR_feature$conf.int[1],upper_CI=OR_feature$conf.int[2])
                              )
    }
  }
  if(nrow(genomic_features_OR)>0){
    genomic_features_OR$TF=TF
    genomic_features_OR_out=rbind(genomic_features_OR_out,genomic_features_OR)
  }
  #Get samples for chromHMM analysis
  sample_chromHMM=names(table(variant_dNME_ASM_TF[variant_dNME_ASM_TF$ASM=="Yes"]$Sample))[table(variant_dNME_ASM_TF[variant_dNME_ASM_TF$ASM=="Yes"]$Sample)>=5]
  chromHMM_ls=list()
  for(sp in sample_chromHMM){
    chromHMM_in=chromHMM_list[[sp]]
    count_table=list()
    out_df=data.table()
    for(st in unique(chromHMM_in$name)){
      OR_chromHMM=NULL
      OR_chromHMM=testEnrichmentFeature_stat(variant_dNME_ASM_TF,chromHMM_in[chromHMM_in$name==st],output_ft=1)
      if(!is.null(OR_chromHMM)){
      count_table[[st]]=OR_chromHMM[[1]]
      OR_chromHMM=OR_chromHMM[[2]]
      out_df=rbind(out_df,
                   data.table(state=st,OR=OR_chromHMM$estimate,p_value=OR_chromHMM$p.value,
                              lower_CI=OR_chromHMM$conf.int[1],upper_CI=OR_chromHMM$conf.int[2]))
      }
    
      
    }
    out_df$Sample=sp
    
    chromHMM_ls[[sp]]=list(out_df,count_table)
  }
  if(length(chromHMM_ls)>0){chromHMM_motif_all=chromHMM_combine(chromHMM_ls)}
  if(nrow(chromHMM_motif_all)>0){
    chromHMM_motif_all$TF=TF
    out_df$TF=TF
    cont_out_df=rbind(cont_out_df,out_df)
    chromHMM_motif_all_TF=rbind(chromHMM_motif_all_TF,chromHMM_motif_all)
  }
  cat("Finish TF:",TF,'\n')
}

chromHMM_motif_all_TF$FDR=p.adjust(chromHMM_motif_all_TF$p_value,method="BH")
genomic_features_OR$FDR=p.adjust(genomic_features_OR$pvalue,method='BH')
saveRDS(chromHMM_motif_all_TF,'../downstream/output/chromHMM_motif_all_TF.rds')
saveRDS(genomic_features_OR_out,'../downstream/output/genomic_features_OR_motif.rds')
saveRDS(OR_cont_table,'../downstream/output/OR_cont_table_motif.rds')
chromHMM_motif_all_TF_sig=chromHMM_motif_all_TF[FDR<=0.1&OR>=2]
chromHMM_motif_all_TF_sig_heatmap=dcast.data.table(chromHMM_motif_all_TF[TF%in%chromHMM_motif_all_TF_sig$TF],states~TF,value.var = "OR")
chromHMM_motif_all_TF_sig_FDR=as.matrix(dcast.data.table(chromHMM_motif_all_TF[TF%in%chromHMM_motif_all_TF_sig$TF],states~TF,value.var = "FDR")[,-1])
chromHMM_motif_all_TF_sig_FDR[chromHMM_motif_all_TF_sig_FDR <= 0.1]="*"
chromHMM_motif_all_TF_sig_FDR[chromHMM_motif_all_TF_sig_FDR > 0.1]=""
heatmap_mt=as.matrix(chromHMM_motif_all_TF_sig_heatmap[,-1])
rownames(heatmap_mt)=chromHMM_motif_all_TF_sig_heatmap$states
breaksList=seq(0,3,0.05)
heatmap_col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
pdf('../downstream/output/chromHMM_SNP_feature.pdf')
pheatmap(t(heatmap_mt),cluster_cols = F,cluster_rows = F,display_numbers=t(chromHMM_motif_all_TF_sig_FDR),
         color=heatmap_col,breaks = breaksList)
dev.off()
#Find CTCF example

TF="CTCF"
motif_gene_high_NME_TF=motif_gene_high_NME[motif_gene_high_NME$geneSymbol==TF]
variant_dNME_ASM_TF=variant_HetCpG_meta_dNME_ASM
variant_dNME_ASM_TF$allele_diff=0
olap=findOverlaps(variant_dNME_ASM_TF,motif_gene_high_NME_TF)
variant_dNME_ASM_TF[queryHits(olap)]$allele_diff=motif_gene_high_NME_TF[subjectHits(olap)]$alleleDiff  
variant_dNME_ASM_TF$pctRef=NA
variant_dNME_ASM_TF[queryHits(olap)]$pctRef=motif_gene_high_NME_TF[subjectHits(olap)]$pctRef  
variant_dNME_ASM_TF$pctAlt=NA
variant_dNME_ASM_TF[queryHits(olap)]$pctAlt=motif_gene_high_NME_TF[subjectHits(olap)]$pctAlt    
variant_dNME_ASM_TF$allele_diff_NME=variant_dNME_ASM_TF$altNME-variant_dNME_ASM_TF$refNME
#Fit the function
variant_dNME_ASM_TF$ASM="No"
variant_dNME_ASM_TF$ASM[sign(variant_dNME_ASM_TF$allele_diff_NME)==sign(variant_dNME_ASM_TF$allele_diff)]="Yes"
cTCF_loc=variant_dNME_ASM_TF[variant_dNME_ASM_TF$ASM=="Yes"]
cTCF_loc=cTCF_loc[cTCF_loc$N>=4]
cTCF_loc=cTCF_loc[order(cTCF_loc$dNME,decreasing=T)]

source('plotMB.R')
pdf('../downstream/output/mouse_examples/motif_CTCF.pdf')
plotMB(subsetByOverlaps(motif_gene_high_NME_TF,cTCF_loc[2]),'STL003-874104')
dev.off()
# #Motif vs mouse enhancer enrichment
# motif_ent_only_OMIM=fread("../downstream/output/graphs/table4_motif_prefer_high_NME.csv")
# tissue=c('limb','forebrain','heart')
# motif_sig_out=data.table()
# motif_sig_TF_mouse=c()
# for(ts in tissue){
#   folder_in=paste0('../downstream/input/mouse_motif_enrichment_enhancer/',ts,'/')
#   for(fn in dir(folder_in,pattern='csv')){
#     mouse_motif_in=fread(paste0(folder_in,fn))
#     motif_sig=sub('.*_','',mouse_motif_in[FDR<=0.1]$motif)
#     motif_sig_TF_mouse=c(motif_sig_TF_mouse,motif_sig)
#      motif_ent_only_OMIM_sig=motif_ent_only_OMIM[TF %in% motif_sig]
#       if(length(motif_sig)>0){
#       motif_ent_only_OMIM_sig$mouse=gsub('.csv','',fn)
#       motif_sig_out=rbind(motif_sig_out,motif_ent_only_OMIM_sig)}
#   }
#   
# }
# 
# motif_sig_out=motif_sig_out[,list(mouse_enrich=paste(sort(gsub('motif_','',mouse)),collapse = ",")),
#                             by=list(TF,same_dir,opposite_dir,Proportion,Pvalue,FDR,OMIM)]
# write.csv(motif_sig_out,"../downstream/output/graphs/tableS2_motif_prefer_high_NME_only_mouse.csv")
# motif_ent_only_OMIM$mouse_enrich=NA
# motif_sig_out$V1=NA
# motif_ent_only_OMIM=rbind(motif_sig_out,motif_ent_only_OMIM[!TF%in%motif_sig_out$TF])
# write.csv(motif_ent_only_OMIM,"../downstream/output/graphs/tableS2_motif_prefer_high_NME_only_mouse_info.csv")
# #All JASPAR
# JASPAR_2020=readRDS('../downstream/input/JASPAR_2020_human_PFM.rds')
# JASPAR_name=fread('../downstream/input/JASPAR_human_redundant_2020.csv')
# JASPAR_name=JASPAR_name[ID %in% names(JASPAR_2020)]
# rm(JASPAR_2020)
# #OR or enrichment
# motif_sig_TF_mouse=unique(motif_sig_TF_mouse)
# ent_only_mouse=sum(!is.na(motif_ent_only_OMIM$mouse_enrich))
# ent_only_non_mouse=sum(is.na(motif_ent_only_OMIM$mouse_enrich))
# non_ent_only_mouse=length(motif_sig_TF_mouse)-ent_only_mouse
# non_ent_only_non_mouse=sum(!(JASPAR_name$Name%in%motif_ent_only_OMIM$TF)&!(JASPAR_name$Name%in%motif_sig_TF_mouse))
# 
# fisher.test(matrix(c(ent_only_mouse,ent_only_non_mouse,non_ent_only_mouse,non_ent_only_non_mouse),nrow=2))

# human_mono_motif_TSV=as.data.frame(read.table('../downstream/input/JASPAR_human_redundant_2020.csv',sep=',',header=T,stringsAsFactors = F))
# more_ent_enrich=motif_family_enrich(unique(motif_dir$TF[motif_dir$FDR<=0.1 & motif_dir$prob>0.5]),
#                                     unique(motif_gene$geneSymbol),human_mono_motif_TSV)
# more_less_ent_enrich=motif_family_enrich(unique(motif_dir$TF[motif_dir$FDR<=0.1 & motif_dir$prob<0.5]),
#                                     unique(motif_gene$geneSymbol),human_mono_motif_TSV)
# write.csv(unique(motif_dir$TF[motif_dir$FDR<=0.1 & motif_dir$prob>0.5],
#                  '../downstream/output/graphs/tableS1_motif_prefer_ent.csv'))
# write.csv(more_ent_enrich,'../downstream/output/graphs/table1_motif_family_prefer_ent.csv')
# write.csv(more_less_ent_enrich,'../downstream/output/graphs/table2_motif_family_not_prefer_ent.csv')
# OMIM=read.csv('../downstream/input/genemap2.txt',fill = T,sep='\t',header=T,stringsAsFactors =F,na.strings = "NA",skip=3)
# OMIM$Gene.Symbols=strsplit(as.character(OMIM$Gene.Symbols),', ')
# 
# OMIM_annotation<-function(motif_in,OMIM){
# motif_in=motif_in[order(motif_in$prob,decreasing=T),]
# motif_in$OMIM=NA
# for(tf in unique(motif_in$TF)){
#     tf_in=gsub('\\(var.2\\)','',tf)
#     tf_in=gsub('\\(var.3\\)','',tf_in)
#     
#     tf_in=unlist(strsplit(tf_in,'::'))
#     #print(tf_in)
#     OMIM_disease=OMIM$Phenotypes[which(unlist(lapply(OMIM$Gene.Symbols,function(x) any(x%in% tf_in))))]
#     if(length(OMIM_disease)>0){motif_in$OMIM[motif_in$TF==tf]=as.list(OMIM_disease)}
#   }
#   motif_in$OMIM=unlist(lapply(motif_in$OMIM,function(x) paste(x,collpase=',')))
#   return(motif_in)
# }
# motif_ent_OMIM=OMIM_annotation(motif_dir[motif_dir$FDR<=0.1 & motif_dir$prob>0.5,],OMIM)
# write.csv(motif_ent_OMIM[c('TF','total_data','same_dir','opposite_dir','prob','FDR','OMIM')],
#           '../downstream/output/graphs/tableS1_motif_prefer_ent_OMIM.csv')
# motif_non_ent_OMIM=OMIM_annotation(motif_dir[motif_dir$FDR<=0.1 & motif_dir$prob<0.5,],OMIM)
# write.csv(motif_non_ent_OMIM[c('TF','total_data','same_dir','opposite_dir','prob','FDR','OMIM')],
#           '../downstream/output/graphs/table3_motif_not_prefer_ent_OMIM.csv')
#LiftOver tools

