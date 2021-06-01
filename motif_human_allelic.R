rm(list=ls())
source("mainFunctions_sub.R")

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()

# get all variant ---------------------------------------------------------

variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)

motif_gene <- readRDS(motif_gene_file)#See motif_break_array.R, default setting
#All regions
#NME
motif_dir_dNME=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="NME")
colnames(motif_dir_dNME)[c(1,6,7)]=c('TF','Pvalue','Proportion')
motif_dir_dNME$FDR=p.adjust(motif_dir_dNME$Pvalue,method="BH")
saveRDS(motif_dir_dNME,'../downstream/output/human_analysis/motif_analysis/dNME_all.rds')
motif_dir_dNME=readRDS('../downstream/output/human_analysis/motif_analysis/dNME_all.rds')
write.csv(motif_dir_dNME[FDR<=0.1&Proportion>0.5,list(TF,Proportion,Pvalue,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
write.csv(motif_dir_dNME[FDR<=0.1&Proportion<0.5,list(TF,Proportion=1-Proportion,`Pvalue`,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/table2_motif_prefer_low_NME.csv')

#MML
motif_dir_dMML=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="MML")
colnames(motif_dir_dMML)[c(1,6,7)]=c('TF','Pvalue','Proportion')
motif_dir_dMML$FDR=p.adjust(motif_dir_dMML$Pvalue,method="BH")
saveRDS(motif_dir_dMML,'../downstream/output/human_analysis/motif_analysis/dMML_all.rds')
motif_dir_dMML=readRDS('../downstream/output/human_analysis/motif_analysis/dMML_all.rds')
write.csv(motif_dir_dMML[FDR<=0.1&Proportion<0.5,list(TF,Proportion=1-Proportion,Pvalue,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/tablS1_motif_prefer_low_MML.csv')
write.csv(motif_dir_dMML[FDR<=0.1&Proportion>0.5,list(TF,Proportion,Pvalue,FDR)], row.names =F,
          '../downstream/output/graphs/motif_preference_table/All_regions/tableS2_motif_prefer_high_MML.csv')


# Check against Sharp's data ----------------------------------------------
#Supp data 5 from Sharp paper
motif_sharp=as.data.table(readxl::read_xlsx('../downstream/input/Sharp_motif_list.xlsx',skip=1))
motif_sharp_sig=motif_sharp$TFBS[which(as.numeric(motif_sharp$`Bonferroni-corrected P-value`)<=0.05)]
length(motif_sharp_sig)#46
sum(motif_sharp_sig%in%motif_dir$TF)#21
sum(motif_sharp_sig%in%motif_dir[FDR<=0.1]$TF)#4
sum(motif_sharp$TFBS[which(as.numeric(motif_sharp$`Bonferroni-corrected P-value`)>0.05)]%in%motif_dir[FDR<=0.1]$TF)#9
sum(motif_sharp$TFBS[which(as.numeric(motif_sharp$`Bonferroni-corrected P-value`)>0.05)]%in%motif_dir[FDR>0.1]$TF)#14
fisher.test(matrix(c(4,17,9,14),nrow=2))#About random chance, 44 motif they analyze we also analyze

#Find low MML only and high NME only

low_MML=motif_dir_dMML[FDR<=0.1&Proportion<0.5]$TF
high_NME=motif_dir_dNME[FDR<=0.1&Proportion>0.5]$TF

motif_dir_dMML$Proportion_high_MML=motif_dir_dMML$Proportion
motif_dir_dMML$Proportion_low_MML=1-motif_dir_dMML$Proportion_high_MML
motif_dir_dMML$dMML_pvalue=motif_dir_dMML$Pvalue
motif_dir_dMML$dMML_FDR=motif_dir_dMML$FDR

motif_dir_dNME$proportion_high_NME=motif_dir_dNME$Proportion
motif_dir_dNME$proportion_low_NME=1-motif_dir_dNME$Proportion
motif_dir_dNME$dNME_pvalue=motif_dir_dNME$Pvalue
motif_dir_dNME$dNME_FDR=motif_dir_dNME$FDR

low_MML_only=cbind(data.table(TF=low_MML[!low_MML%in%high_NME]),
                   motif_dir_dMML[TF%in%low_MML[!low_MML%in%high_NME]],
                   motif_dir_dNME[TF%in%low_MML[!low_MML%in%high_NME]])


high_NME_only=cbind(data.table(high_NME[!high_NME%in%low_MML]),
                    motif_dir_dNME[TF%in%high_NME[!high_NME%in%low_MML]],
                    motif_dir_dMML[TF%in%high_NME[!high_NME%in%low_MML]])

write.csv(low_MML_only[order(Proportion_low_MML),list(TF,Proportion_low_MML,dMML_pvalue,dMML_FDR,proportion_high_NME,dNME_pvalue,dNME_FDR)],
          '../downstream/output/human_analysis/motif_analysis/motif_low_MML_only.csv')
write.csv(high_NME_only[order(Proportion_low_MML),list(TF,proportion_high_NME,dNME_pvalue,dNME_FDR,Proportion_low_MML,dMML_pvalue,dMML_FDR)],
          '../downstream/output/human_analysis/motif_analysis/motif_high_NME_only.csv')
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
low_MML_motif=fread('../downstream/output/human_analysis/motif_analysis/motif_low_MML_only.csv')
high_NME_motif=fread('../downstream/output/human_analysis/motif_analysis/motif_high_NME_only.csv')
OMIM=fread('../downstream/input/human_analysis/SNP_biology/genemap2.txt',skip=3)
OMIM$`Gene Symbols`=strsplit(as.character(OMIM$`Gene Symbols`),', ')
OMIM=OMIM[Phenotypes!=""]
motif_ent_OMIM=OMIM_annotation(high_NME_motif,OMIM)
motif_ent_only_OMIM=OMIM_annotation(high_NME_motif[!(TF%in%low_MML_motif$TF)],OMIM)
motif_low_MML_only_OMIM=OMIM_annotation(low_MML_motif[!(TF%in%high_NME_motif$TF)],OMIM)
motif_shared_OMIM=OMIM_annotation(low_MML_motif[TF%in%high_NME_motif$TF],OMIM)
write.csv(motif_low_MML_only_OMIM, row.names =F,
          "../downstream/output/human_analysis/motif_analysis/tableX_motif_prefer_low_MML_only.csv")
write.csv(motif_ent_only_OMIM, row.names =F,
          "../downstream/output/human_analysis/motif_analysis/table3_motif_prefer_high_NME_only.csv")
write.csv(motif_shared_OMIM, row.names =F,
          "../downstream/output/human_analysis/motif_analysis/tableXmotif_shared_only.csv")
# write.csv(unlist(strsplit(gsub('\\(var.3\\)','',gsub('\\(var.2\\)','',motif_low_MML_only_OMIM$TF)),"::")),
#           "../downstream/output/motif_prefer_low_MML_only.csv")
# write.csv(unlist(strsplit(gsub('\\(var.3\\)','',gsub('\\(var.2\\)','',motif_ent_only_OMIM$TF)),"::")),
#           "../downstream/output/motif_prefer_high_NME_only.csv")
# write.csv(unlist(strsplit(gsub('\\(var.3\\)','',gsub('\\(var.2\\)','',motif_shared_OMIM$TF)),"::")),
#           "../downstream/output/motif_shared.csv")

#Ken's analysis preprocessing is in human_motif_processing.R


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

