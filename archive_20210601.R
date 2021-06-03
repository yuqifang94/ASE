# counting ASM ------------------------------------------------------------

#Examples:
#dNME example1:
subsetByOverlaps(GR_merge[GR_merge$Sample=="ectoderm_paired - HUES64"],GRanges(seqnames="chr14",IRanges(start=104552150,end=104552495)))
#dNM example2:
subsetByOverlaps(GR_merge[GR_merge$Sample=="Psoas_Muscle_single - STL003"],GRanges(seqnames="chr17",IRanges(start=33750066,end=33770266)))
#dMML_example1
subsetByOverlaps(GR_merge[GR_merge$Sample=="stem_27_undifferentiated_paired - HUES64"],GRanges(seqnames="chr11",IRanges(start=2720817,end=2721033)))
#dMML_example2
subsetByOverlaps(GR_merge[GR_merge$Sample=="endoerm_27_paired - HUES64"],GRanges(seqnames="chr20",IRanges(start=32308087,end=32308287)))
#Find examples for dNME
GR_merge_dNME=GR_merge[GR_merge$N>=4&GR_merge$dNME_pval<=pval_cutoff]
GR_merge_dNME=GR_merge_dNME[order(GR_merge_dNME$dNME,decreasing=T)]
GR_merge_dNME_dt=convert_GR(GR_merge_dNME,direction="DT")
#Example 1: dNME
#chr8:3,770,843-3,771,336
#HUES64 stem
#chr14:104,566,670-104,566,870
#chr4:7,308,694-7,309,119
#chr12:125,583,449-125,583,649
#chr11:103,242,820-103,243,199

#Examples mpre dNME
#chr6:20,481,578-20,481,925, edoderm
#     chr2 100438919-100439371,ectoderm  
#chr3   76733036-76733603 ectoderm
#   chr1   18993859-18994059 , endoderm
# chr9     4023395-4023595   esc

Imprinted_Genes <- as.data.frame(read_excel("../downstream/input/human_analysis/imprinting_ASE/Imprinted Genes.xlsx"))
GR_merge_dMML=GR_merge[GR_merge$genes_promoter%in%Imprinted_Genes$Gene&GR_merge$dMML_pval<=pval_cutoff]
#chr6 144329267-144329363 149
#chr6:144,329,572-144,329,772 150
#chr6 144329572-144329772

#Mouse analysis
NME_in=readRDS('../downstream/input/mouse_analysis/NME_agnostic_mouse_all_merged.rds')
MML_in=readRDS('../downstream/input/mouse_analysis/MML_agnostic_mouse_all_merged.rds')
gtf=fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table=F)
promoter_in=gtf <- gtf[gtf[,3]=='gene',]
type <- sub('\".*','',sub('.*gene_type \"','',gtf[,9]))
gtf <- gtf[type=='protein_coding',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gene_body <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
gene_body$gene=gn
olap=findOverlaps(MML_in,gene_body)
MML_in$gene=NA
MML_in$gene[queryHits(olap)]=gene_body$gene
olap=findOverlaps(NME_in,gene_body)
NME_in$gene=NA
NME_in$gene[queryHits(olap)]=gene_body$gene

MML_in_gata4=subsetByOverlaps(MML_in,gene_body[gene_body$gene=="Gata4"])
MML_in_gata4=sort(MML_in_gata4)

NME_in_gata4=subsetByOverlaps(NME_in,gene_body[gene_body$gene=="Gata4"])
NME_in_gata4=sort(NME_in_gata4)
MML_in_gata4[MML_in_gata4$N>=4&MML_in_gata4$MML>=0.96&MML_in_gata4$tissue%in%c("heart","forebrain","limb")]#Forebrain E16.5: chr14 63212488-63212737

NME_in_gata4[NME_in_gata4$N>=4&NME_in_gata4$NME>=0.97&NME_in_gata4$tissue%in%c("heart","forebrain","limb")]#forebrain E11.5


# GR_merge=readRDS('GR_merge_final12_ls.rds')
# in_dir='../allele_agnostic_hg19_cov5_3kb_FANTOM/'
# NME_in=GRanges()
# MML_in=GRanges()
# for(fn in  dir(in_dir)){
#   cat('Reading in',fn,'\n')
#   stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
#   sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
#   subject_in=sub('_.*','',sample_in)
#   tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
#   sample_in=paste0(tissue_in,' - ',subject_in)
#   
#   if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
#   if(stat_in=="NME"){
#     NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$Sample==sample_in]))}
#   else if(stat_in=="MML"){
#     MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$Sample==sample_in]))}else
#     {cat("Error stat_in:", stat_in,'\n')}
# }
# 
# NME_in$NME=NME_in$score
# MML_in$MML=MML_in$score
# NME_in=NME_in[!(NME_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
# MML_in=MML_in[!(MML_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
# NME_in=NME_in[NME_in$N>=2]
# MML_in=MML_in[MML_in$N>=2]
# saveRDS(NME_in,"allele_agnostic_hg19_cov10_3kb_FANTOM_NME.rds")
# saveRDS(MML_in,"allele_agnostic_hg19_cov10_3kb_FANTOM_MML.rds")
# NME_in$hyper_var_fn[NME_in$Sample %in% c('42_embryonic_stem_cell_single - H9','stem_27_undifferentiated_paired - HUES64',"ESC - H1")]=
#   paste('../downstream/input/scRNA/HESC_1.rds',sep = ';')
# archive from motif human allelic ----------------------------------------

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


# gff_file_generation -----------------------------------------------------


##Enrichment of enhancer in DNase region
DNAase=readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp.rds')
control=c(readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp_control.rds'),
          readRDS('../downstream/input/mouse_analysis/mm10_PRC.rds'),
          readRDS('../downstream/output/mouse_analysis/mm10_allele_agnostic_analysis_compliment.rds'))
enhancer=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')

enhancer_DNase=length(subsetByOverlaps(DNAase,enhancer))
non_enhancer_DNase=length(DNAase)-enhancer_DNase
enhancer_control=length(subsetByOverlaps(control,enhancer))
non_enhancer_control=length(control)-enhancer_control
OR=fisher.test(matrix(c(enhancer_DNase,non_enhancer_DNase,enhancer_control,non_enhancer_control),nrow=2))

TSS=get_mm10_tss()

enhancer_DNase=length(subsetByOverlaps(DNAase,TSS))
non_enhancer_DNase=length(DNAase)-enhancer_DNase
enhancer_control=length(subsetByOverlaps(control,TSS))
non_enhancer_control=length(control)-enhancer_control
OR=fisher.test(matrix(c(enhancer_DNase,non_enhancer_DNase,enhancer_control,non_enhancer_control),nrow=2))


# CPEL imprinting analysis ------------------------------------------------

NME_in=readRDS('../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_ASM.rds')
MML_in=readRDS('../downstream/output/human_analysis/CPEL_outputs/MML_agnostic_ASM.rds')
NME_in_tb=as.data.table(mcols(NME_in))
NME_in_tb$region=paste0(seqnames(NME_in),':',start(NME_in),'-',end(NME_in))
MML_in_tb=as.data.table(mcols(MML_in))
MML_in_tb$region=paste0(seqnames(MML_in),':',start(MML_in),'-',end(MML_in))
GR_merge_tb=rbind(NME_in_tb,MML_in_tb)
rm(NME_in_tb)
rm(MML_in_tb)
rm(NME_in)
rm(MML_in)
GR_merge_tb$K=NULL
GR_merge_tb$N=NULL
GR_merge=readRDS(GR_merge_file)
GR_merge_tb_asm=rbind(
  data.table(score=GR_merge$dNME,Sample=GR_merge$Sample,statistics='dNME',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dNME_pval,Sample=GR_merge$Sample,statistics='dNME_pval',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dMML,Sample=GR_merge$Sample,statistics='dMML',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dMML_pval,Sample=GR_merge$Sample,statistics='dMML_pval',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=(GR_merge$MML1+GR_merge$MML2)/2,Sample=GR_merge$Sample,statistics='MML_ASM',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=(GR_merge$NME1+GR_merge$NME2)/2,Sample=GR_merge$Sample,statistics='NME_ASM',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge)))
)
rm(GR_merge)

GR_merge_tb=rbind(GR_merge_tb,GR_merge_tb_asm)
GR_merge_tb=dcast.data.table(GR_merge_tb,Sample+region~statistics,value.var = "score")
saveRDS(GR_merge_tb,'../downstream/output/human_analysis/imprinting/GR_merge_ASM_comp.rds')
