
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


# Ken_motif_prep ----------------------------------------------------------



#Prepare for GWAS analysis
NME_in=readRDS('../downstream/output/human_analysis/CPEL_outputs/allele_agnostic_hg19_DNase_NME_homogeneous_excluding_dMML.rds')
mcols(NME_in)=mcols(NME_in)[,c('Sample','NME')]
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp_control.rds')
NME_in$DNAase="NA"
olap_DNase=findOverlaps(NME_in,DNase,type='equal')
olap_control=findOverlaps(NME_in,control,type='equal')
NME_in$DNAase[queryHits(olap_DNase)]="DNAase"
NME_in$DNAase[queryHits(olap_control)]="control"
saveRDS(NME_in,'../downstream/output/human_analysis/CPEL_outputs/NME_DNAase_control_hg19.rds')

# DNase region using allelic region -----------------------------------------
GR_merge=readRDS(GR_merge_file)
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
GR_merge_DNase=NME_dNME_ken(DNase,GR_merge,"NME")
GR_merge_DNase_mt=as.matrix(mcols(GR_merge_DNase))
rownames(GR_merge_DNase_mt)=paste0(seqnames(GR_merge_DNase),':',start(GR_merge_DNase),'-',end(GR_merge_DNase))
saveRDS(GR_merge_DNase_mt,'../downstream/output/DNase_mt_SNP_allelic.rds')
#Dnase region using allele-agnostic model at SNP

NME_in=readRDS('../downstream/output/NME_agnostic_ASM.rds')
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
GR_merge_DNase=NME_dNME_ken(DNase,NME_in,"NME")
GR_merge_DNase_mt=as.matrix(mcols(GR_merge_DNase))
rownames(GR_merge_DNase_mt)=paste0(seqnames(GR_merge_DNase),':',start(GR_merge_DNase),'-',end(GR_merge_DNase))
saveRDS(GR_merge_DNase_mt,'../downstream/output/DNase_mt_SNP_agnostic.rds')
#DNase region at DNase regions all
NME_in=readRDS('../downstream/output/allele_agnostic_hg19_DNase_NME.rds')
DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
NME_in=subsetByOverlaps(NME_in,DNase,type='equal')
GR_merge_DNase=NME_dNME_ken(DNase,NME_in,"NME")
GR_merge_DNase_mt=as.matrix(mcols(GR_merge_DNase))
rownames(GR_merge_DNase_mt)=paste0(seqnames(GR_merge_DNase),':',start(GR_merge_DNase),'-',end(GR_merge_DNase))
saveRDS(GR_merge_DNase_mt,'../downstream/output/DNase_mt_all_agnostic.rds')

#Add NA column to non_NA thing
samples_in=readRDS('../downstream/output/huamn_samples.rds')
folder_in='../downstream/output/human_analysis/Ken_motif/allelic_motif_hg19/'
folder_in='../downstream/output/human_analysis/Ken_motif/homo'
for(fn in dir(folder_in,pattern="MML")){
  MML_in=readRDS(paste0(folder_in,fn))
  for(idx in which(unlist(lapply(MML_in,function(x) ncol(mcols(x))!=49)))){
    
    for (sp in samples_in[!samples_in %in% colnames(mcols(MML_in[[idx]]))]){
      
      mcols(MML_in[[idx]])[[sp]]=as.numeric(NA)
      
    }
    
    MML_in[[idx]]=MML_in[[idx]][,samples_in]
    
    
  }
  saveRDS(MML_in,paste0(folder_in,gsub('.rds','.complete.rds',fn)))
}

# CpG_density_NME_MML.R ---------------------------------------------------
genomic_features=readRDS(genomic_features_file)
variant_HetCpG_meta$CpG_island=FALSE
variant_HetCpG_meta$CpG_island[queryHits(olap)]=TRUE

#mutation_tri_unique[gainCG_idx]=paste0(sub('.*->','', mutation_tri_unique[gainCG_idx]),'->',sub('->.*','', mutation_tri_unique[gainCG_idx]))

# # #Lumping all
# variant_HetCpG_meta_dt$CpG_change="No change"
# variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri ))]$CpG_change='Lose CG'
# variant_HetCpG_meta_dt[(!grepl('CG',REF_tri )) & (grepl('CG',ALT_tri ))]$CpG_change='Gain CG'
# variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (grepl('CG',ALT_tri))]$CpG_change='No change'

# #Calculate relative dNME
# variant_HetCpG_meta_dt$dNME_relative=variant_HetCpG_meta_dt$refNME-variant_HetCpG_meta_dt$altNME
# #convert all to less CG - more CG, in this case, a minus sign is added in the ones "Lose CG" since in that case altNME has fewer CG than refNME
# variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dNME_relative=-variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dNME_relative

#substr(variant_HetCpG_meta_dt$tri_SNP,2,2)="X"
#For C->G SNP, change background
# if(sn != 'C->G'){
#   #Write in method
#   variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
# variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
# }else{
#   #Separate G->C in SNP and C->G SNP
#   #G->C SNP
#   tri_l=gsub('->.*','',tri)
#   tri_r=gsub('.*->','',tri)
#   if(paste0(substring(tri_l, 2, 2),'->',substring(tri_r, 2, 2))=='G->C'){
#     
#     variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[paste0(substring(gsub('->.*','',tri_SNP_unique), 2, 2),'->',substring(gsub('.*->','',tri_SNP_unique), 2, 2))=='G->C' &
#                                                         dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
#     
#   }else{
#     variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[paste0(substring(gsub('->.*','',tri_SNP_unique), 2, 2),'->',substring(gsub('.*->','',tri_SNP_unique), 2, 2))=='C->G' &
#                                                         dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
#     
#     
#   }
#}
# variant_SNP_tri$SNP=paste0(gsub('.*->','',variant_SNP_tri$SNP),'->',gsub('->.*','',variant_SNP_tri$SNP))
# sn=  paste0(gsub('.*->','',sn),'->',gsub('->.*','',sn))
# SNP_all[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR))+geom_bar(fill=color_theme[sn],stat="identity")+ylab('Enrichment')+
#   geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.1,position=position_dodge(.9))+ggtitle(sn)+xlab("")+
# theme_glob+theme(legend.position = "none")+
#   #theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+geom_text(aes(label=significant,y=upperCI*1+sig_h),vjust = sig_v,hjust=sig_h)+coord_flip()
# SNP_box[[sn]]=ggplot(variant_HetCpG_meta_dt[SNP==sn & dNME_pval<=pval_cutoff],aes(x=reorder(tri_SNP_unique, -NME_relative, FUN = median) ,y=NME_relative,fill=CpG_change))+
#   geom_boxplot()+
#   ylab('refNME-altNME')+theme_glob+theme(legend.position = "bottom")+ ylim(c(-1,1))+ 
#   scale_fill_manual(values=c("No CG change"="grey","Gain CG"="light blue"))+
#   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   xlab('')+coord_flip()

# pdf('../downstream/output/graphs/Figure3/variant_OR_tri2.pdf',width=14,height=14)
# SNP_all=SNP_all[c("C>G", names(SNP_all)[names(SNP_all)!="C>G"])]
# ggarrange(plotlist=SNP_all, nrow=4,ncol=3)
# dev.off()
#After import_font, use loadfont(device='pdf')
# font_import('C:/Users/vince/Downloads/dejavu-fonts-ttf-2.37')
# loadfonts(device='pdf')
#cairo_pdf('../downstream/output/graphs/Figure2/Figure-4B-variant_OR_tri3_two_cat_greater_CG_bg.pdf',width=10,height=7,family = "DejaVu Sans")
# pdf('../downstream/output/graphs/Figure3/variant_box_tri.pdf',width=14,height=14)
# SNP_box=SNP_box[c("C>G", names(SNP_box)[names(SNP_box)!="C>G"])]
# ggarrange(plotlist=SNP_box, nrow=4,ncol=3,common.legend = T,legend="bottom")
# dev.off()



OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_2cat.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

#For dMML

#Calculate relative dMML
variant_HetCpG_meta_dt=readRDS('../downstream/output/human_analysis/CPEL_outputs/variant_HetCpG_meta_dt.rds')
single_SNP_unique=unique_mutation(unique(variant_HetCpG_meta_dt$SNP))
variant_HetCpG_meta_dt$SNP=single_SNP_unique[variant_HetCpG_meta_dt$SNP]
mutation_tri_unique=unique_mutation(unique(variant_HetCpG_meta_dt$tri_SNP))
#Finding the ones losing CG
gainCG_idx=which(!grepl("CG",sub('.*-','',mutation_tri_unique))&grepl("CG",sub('-.*','',mutation_tri_unique)))
mutation_tri_unique[gainCG_idx]=paste0(sub('.*-','', mutation_tri_unique[gainCG_idx]),'-',sub('-.*','', mutation_tri_unique[gainCG_idx]))
variant_HetCpG_meta_dt$tri_SNP_unique=mutation_tri_unique[variant_HetCpG_meta_dt$tri_SNP]
#Lumping all
variant_HetCpG_meta_dt$CpG_change="No change"
variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri ))]$CpG_change='Lose CG'
variant_HetCpG_meta_dt[(!grepl('CG',REF_tri )) & (grepl('CG',ALT_tri ))]$CpG_change='Gain CG'
variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (grepl('CG',ALT_tri))]$CpG_change='No change'
variant_HetCpG_meta_dt$dMML_relative=variant_HetCpG_meta_dt$refMML-variant_HetCpG_meta_dt$altMML
#convert all to less CG - more CG, in this case, a minus sign is added in the ones "Lose CG" since in that case altNME has fewer CG than refNME
variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dMML_relative=-variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dMML_relative
#Convert everything to gain CG
variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$CpG_change='Gain CG'
SNP_all=list()
SNP_het=list()
SNP_box=list()
sig_v=0.75
sig_h=-0.05
color_theme=c(rainbow(length(unique(variant_HetCpG_meta_dt$SNP))))
variant_SNP_tri=data.table()
names(color_theme)=unique(variant_HetCpG_meta_dt$SNP)
for (sn in unique(variant_HetCpG_meta_dt$SNP)){
  for(tri in unique(variant_HetCpG_meta_dt[SNP==sn]$tri_SNP_unique)){
    variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dMML_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff,stat_in="MML")
    variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
    variant_SNP_tri=rbind(variant_SNP_tri,variant_SNP_tri_OR)
  }
  #get HetCpG
  variant_SNP_tri$CpG_change=factor(variant_SNP_tri$CpG_change,levels=c('Gain CG','No change','Lose CG'))
  variant_SNP_tri=variant_SNP_tri[order(OR,decreasing=F)]
  variant_SNP_tri$SNP=factor(variant_SNP_tri$SNP,levels = variant_SNP_tri$SNP)
  
  variant_SNP_tri$FDR=p.adjust(variant_SNP_tri$pvalue,method='BH')
  variant_SNP_tri$significant=add.significance.stars(variant_SNP_tri$FDR, cutoffs = c(0.05, 0.01, 0.001))
  variant_SNP_tri=variant_SNP_tri[!is.infinite(OR)]
  SNP_all[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR))+geom_bar(fill=color_theme[sn],stat="identity")+ylab('Enrichment')+
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.1,position=position_dodge(.9))+ggtitle(sn)+xlab("")+
    theme_glob+theme(legend.position = "none")+
    #theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+geom_text(aes(label=significant,y=upperCI*1+sig_h),vjust = sig_v,hjust=sig_h)+coord_flip()
  
  SNP_het[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR,fill=CpG_change))+geom_bar(stat="identity")+ylab('Enrichment')+xlab("")+
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.4,position=position_dodge(.9),size=0.25)+ggtitle(sn)+ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+
    theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+geom_text(aes(label=significant,y=upperCI*1),vjust =sig_v,hjust=sig_h)+coord_flip()
  
  
  
  
  variant_SNP_tri=data.table()
}
# pdf('../downstream/output/graphs/Figure3/variant_OR_tri2.pdf',width=14,height=14)
# SNP_all=SNP_all[c("C>G", names(SNP_all)[names(SNP_all)!="C>G"])]
# ggarrange(plotlist=SNP_all, nrow=4,ncol=3)
# dev.off()

pdf('../downstream/output/graphs/Figure2/Figure-S3-variant_OR_tri3_dMML_cat.pdf',width=10,height=7)
#SNP_het=SNP_het[c("C>G", names(SNP_het)[names(SNP_het)!="C>G"])]
ggarrange(plotlist=SNP_het, nrow=2,ncol=3,common.legend = T,legend="bottom")
dev.off()
# pdf('../downstream/output/graphs/Figure3/variant_box_tri.pdf',width=14,height=14)
# SNP_box=SNP_box[c("C>G", names(SNP_box)[names(SNP_box)!="C>G"])]
# ggarrange(plotlist=SNP_box, nrow=4,ncol=3,common.legend = T,legend="bottom")
# dev.off()
#All
OR_all_SNP_change=data.table()

for(sn in unique(variant_HetCpG_meta_dt$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc(variant_HetCpG_meta_dt[dMML_pval<=pval_cutoff],sn,"CpG_change",stat_in="MML"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_dMML.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

#CpG island: -0.189
OR_all_SNP_change=data.table()

for(sn in unique( variant_HetCpG_meta_dt[CpG_island==TRUE]$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc( variant_HetCpG_meta_dt[CpG_island==TRUE][dMML_pval<=pval_cutoff],sn,"CpG_change",stat_in="MML"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_dMML_island.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

#CpG island: 0.333
OR_all_SNP_change=data.table()

for(sn in unique(variant_HetCpG_meta_dt[CpG_island==FALSE]$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc(variant_HetCpG_meta_dt[CpG_island==FALSE][dMML_pval<=pval_cutoff],sn,"CpG_change",stat_in="MML"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_dMML_non_island.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

#Region  analysis
#plot_CpG_number(GR_merge[GR_merge$dNME_pval<=pval_cutoff])

olap=findOverlaps(GR_merge,genomic_features$`CpG island`)

#GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1/CG_allele_extend_g2)]
#dNME
#ratio
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+

# ggplot(GR_merge_dt_sig_density_diff,aes(x=as.factor(round(density_diff,digits = 2)),y=dNME_relative))+geom_violin(fill='light blue')+
#   xlab("CpG density ratio")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
#   theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
#stat_summary(fun=median, geom="point")+
# ggplot(GR_merge_dt_sig_density_diff,aes(x=as.factor(round(density_diff,digits = 2)),y=dNME_relative))+geom_violin(fill='light blue')+
#   xlab("CpG density ratio")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
#   theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# ylim(c(-0.5,0.3))+


pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dNME_quantile.pdf',width=7,height=7)
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
GR_merge_dt_S4=GR_merge_dt[dNME_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_S4$density_diff=findInterval(GR_merge_dt_S4$density_diff,quantile(GR_merge_dt_S4$density_diff,prob=seq(0,0.9,0.1)))
GR_merge_dt_S4$density_diff=factor(paste0(GR_merge_dt_S4$density_diff*10,'%'),levels=paste0( seq(0,100,10),"%"))
ggplot(GR_merge_dt_S4,aes(x=density_diff,y=dNME_relative))+geom_violin(fill='light blue')+
  xlab("CpG density difference quantile")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#dMML

#0.214, significant
cor.test(GR_merge_dt[dMML_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$dMML_relative, 
         GR_merge_dt[dMML_pval<=pval_cutoff&dMML_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$density_diff)

GR_merge_dt$CpG_stat="No difference"
GR_merge_dt[CpGdiff!=0]$CpG_stat="With CpG difference"
GR_merge_dt$CpG_stat=factor(GR_merge_dt$CpG_stat,levels = c("With CpG difference","No difference"))
#Make this one better
#Figure 3A
GR_merge_dt$dMML_relative_more_less=GR_merge_dt$dMML_relative
GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dMML_relative_more_less=GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dMML_relative*sign(GR_merge_dt[GR_merge_dt$CpGdiff!=0]$CpGdiff)
t.test(GR_merge_dt[CpGdiff!=0&dMML_pval<=pval_cutoff]$dMML_relative_more_less,alternative="less")
pdf('../downstream/output/graphs/Figure3/Figure3A_CpG_number_MML.pdf',width=7,height=7)
ggplot(GR_merge_dt[dMML_pval<=pval_cutoff],aes(y=dMML_relative_more_less,x=CpG_stat,fill=CpG_stat))+
  geom_violin()+xlab("")+
  theme_glob+ylab('relative dMML')+theme(legend.position = "none")

dev.off()
pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dMML_ratio.pdf',width=7,height=7)
ggplot(GR_merge_dt[dMML_pval<=pval_cutoff&density_diff!=1],aes(x=as.factor(round(density_diff,digits = 2)),y=dMML_relative))+geom_violin(fill='light blue')+
  xlab("CpG density ratio")+ylab("relative dMML")+stat_summary(fun=median, geom="point")+theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dMML_quantile.pdf',width=7,height=7)
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
GR_merge_dt_S4=GR_merge_dt[dMML_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_S4$density_diff=findInterval(GR_merge_dt_S4$density_diff,quantile(GR_merge_dt_S4$density_diff,prob=seq(0,0.9,0.1)))
GR_merge_dt_S4$density_diff=factor(paste0(GR_merge_dt_S4$density_diff*10,'%'),levels=paste0( seq(0,100,10),"%"))
ggplot(GR_merge_dt_S4,aes(x=density_diff,y=dMML_relative))+geom_violin(fill='light blue')+
  xlab("CpG density difference quantile")+ylab("relative dMML")+stat_summary(fun=median, geom="point")+theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# ggplot(GR_merge_dNME[dNME_pval<=pval_cutoff&CpGdiff!=0],aes(x=dNME_relative_more_less))+geom_density()+xlab("relative dNME")+theme_glob
# ggplot(GR_merge_dNME[dNME_pval<=pval_cutoff&CpGdiff!=0],aes(y=dNME_relative_more_less,x=as.factor(density_diff)))+
#   geom_boxplot()+xlab("CpG difference")+
#   theme_glob+ylab('relative dNME')

# # check NME vs NME from allele specific way -------------------------------
GR_merge_tb=readRDS('../downstream/output/GR_merge_ASM_comp.rds')
GR_merge_tb=GR_merge_tb[!is.na(MML_ASM)&!is.na(MML)&!is.na(NME)]#3284912/52263042,total ASM regions=3332744
pdf('../downstream/output/graphs/Figure3/Fig-S5A-dNME_NME_all_pt.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb,aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb$NME_ASM,GR_merge_tb$NME)
pdf('../downstream/output/graphs/Figure3/Fig-S5B-dNME_NME_all_dMML_ASM_pt.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb[dMML_pval<=pval_cutoff],aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb[dMML_pval<=pval_cutoff]$NME_ASM,GR_merge_tb[dMML_pval<=pval_cutoff]$NME)
pdf('../downstream/output/graphs/Figure3/Fig-S5C-dNME_NME_non_dMML.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb[dMML_pval>pval_cutoff],aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+ xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+  
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb[dMML_pval>pval_cutoff]$NME_ASM,GR_merge_tb[dMML_pval>pval_cutoff]$NME)
#quantile(NME_in$density,prob=seq(0,1,0.2),na.rm=T)
#NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
png('../downstream/output/graphs/Figure3/FigureS_CpG_density_NME_line.png',width=1080,height=1080)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in[sample(1:length(NME_in),round(length(NME_in)/5))])),aes(x=density, y=NME))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot_CG_exp_non_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in[-queryHits(olap_islands)])),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()   

pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot_CG_exp_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in[queryHits(olap_islands)])),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off() 
png('../downstream/output/graphs/FigureS6/NME_density_feature_box.png')
# my_comparisons <- list( c("islands", "shores"), c("islands", "shelf"), c("islands", "open sea"),
#                         c("shores", "shelf"),c("shores", "open sea"),c("shelf", "open sea"))
ggplot(CpG_density_NME,aes(x=feature,y=NME))+
  geom_boxplot()+
  #geom_violin()+ 
  theme_glob+xlab("genomic features")
#stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()

# #same analysis on dMML-ASM
MML_in=readRDS(MML_agnostic_file)
MML_in$CG_hg19=countOverlaps(MML_in,CpG_hg19)
MML_in_gr=unique(granges(MML_in))
gr_seq=getSeq(Hsapiens,MML_in_gr,as.character=T)
MML_in_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
#Do it for MML
MML_in_olap=findOverlaps(MML_in,MML_in_gr,type='equal')
MML_in$CGcont_exp[queryHits(MML_in_olap)]=MML_in_gr$CGcont_exp[subjectHits(MML_in_olap)]
MML_in$density=MML_in$CG_hg19/MML_in$CGcont_exp
MML_in=readRDS('../downstream/input/human_analysis/MML_allele_agnostic_merge_20k_homogeneous_excluding_dMML2_CG.rds')
# MML_in$density=MML_in$CG_hg19/width(MML_in)
# MML_in=MML_in[seqnames(MML_in)%in%paste0("chr",1:22)]
# MML_in$density_quant=findInterval(MML_in$density,seq(0,0.1,0.01))
# MML_in$density_quant=factor(quant_conv[MML_in$density_quant],levels=quant_conv)
MML_in$density_quant=findInterval(MML_in$density,seq(0,1,0.1))
#NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
MML_in$density_quant=factor(quant_conv[MML_in$density_quant],levels=quant_conv)
olap=findOverlaps(MML_in,genomic_features$`CpG island`)
cor.test(MML_in$density,MML_in$MML,method='pearson')#-0.346
cor.test(MML_in[queryHits(olap)]$density,MML_in[queryHits(olap)]$MML,method='pearson')#-0.395
cor.test(MML_in[-queryHits(olap)]$density,MML_in[-queryHits(olap)]$MML,method='pearson')#0.0059
#Inverse correlation between MML an d CpG content: 19325872 [-0.46]
pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in)),aes(x=density_quant, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in[queryHits(olap)])),aes(x=density_quant, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot_non_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in[-queryHits(olap)])),aes(x=density_quant, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_line.png',width=1080,height=1080)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in[MML_in$density<=1][sample(1:length(MML_in[MML_in$density<=1]),round(length(MML_in)/3))])),aes(x=density, y=MML))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Figure S5
olap_islands=findOverlaps(MML_in,genomic_features$`CpG island`)
olap_shores=findOverlaps(MML_in,genomic_features$`CpG shore`)
olap_shelf=findOverlaps(MML_in,genomic_features$`CpG shelf`)
olap_open_sea=findOverlaps(MML_in,genomic_features$`CpG open sea`)

CpG_density_MML=rbind(data.table(MML=MML_in$MML[queryHits(olap_islands)],feature='islands'),
                      data.table(MML=MML_in$MML[queryHits(olap_shores)],feature='shores'),
                      data.table(MML=MML_in$MML[queryHits(olap_shelf)],feature='shelf'),
                      data.table(MML=MML_in$MML[queryHits(olap_open_sea)],feature='open sea'))



pdf('../downstream/output/graphs/FigureS6/MML_density_feature_box.pdf')
# my_comparisons <- list( c("islands", "shores"), c("islands", "shelf"), c("islands", "open sea"),
#                         c("shores", "shelf"),c("shores", "open sea"),c("shelf", "open sea"))
ggplot(CpG_density_MML,aes(x=feature,y=MML))+
  geom_boxplot(outlier.shape = NA)+
  #geom_violin()+ 
  theme_glob+xlab("genomic features")
#stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()

enhancer=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')
tss=get_mm10_tss()
pdf('../downstream/output/graphs/Figure3/mouse_MML_density_enhancer_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(mml,enhancer),stat_name="MML")
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_MML_density_promoter_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(mml,tss,maxgap = 2000),stat_name="MML")
dev.off()

pdf('../downstream/output/graphs/Figure3/mouse_NME_density_enhancer_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(nme,enhancer),stat_name="NME")
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_NME_density_promoter_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(nme,tss,maxgap = 2000),stat_name="NME")
dev.off()

#Looking only looking at ones in cluster result
cluster_out=readRDS('../downstream/input/mouse_analysis/clustering/tissue_specific/kmeans_10run/uc_11.rds')

nme_dt=convert_GR(nme,direction="DT")
nme_dt=melt.data.table(nme_dt,id.vars = c("CG_mm10","CGcont_exp","density","region","density_quant"),variable.name = "Sample",value.name="stat_in")
nme_dt=nme_dt[!is.na(stat_in)]
mml_dt=convert_GR(mml,direction="DT")
mml_dt=melt.data.table(mml_dt,id.vars = c("CG_mm10","CGcont_exp","density","region","density_quant"),variable.name = "Sample",value.name="stat_in")
mml_dt=mml_dt[!is.na(stat_in)]
nme_dt_clu=data.table()
mml_dt_clu=data.table()
for(ts in names(cluster_out)){
  nme_dt_clu=rbind(nme_dt_clu,nme_dt[grepl(ts,Sample) & (region %in%names(cluster_out[[ts]]))])
  mml_dt_clu=rbind(mml_dt_clu,mml_dt[grepl(ts,Sample) & (region %in%names(cluster_out[[ts]]))])
  
}


pdf('../downstream/output/graphs/Figure3/mouse_MML_density_enhancer_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
enhancer_olap=findOverlaps(convert_GR(mml_dt_clu$region),enhancer)
print(ggplot(mml_dt_clu[queryHits(enhancer_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_MML_density_promoter_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
promoter_olap=findOverlaps(convert_GR(mml_dt_clu$region),tss,maxgap=2000)
print(ggplot(mml_dt_clu[queryHits(promoter_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

pdf('../downstream/output/graphs/Figure3/mouse_NME_density_enhancer_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
enhancer_olap=findOverlaps(convert_GR(nme_dt_clu$region),enhancer)
print(ggplot(nme_dt_clu[queryHits(enhancer_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_NME_density_promoter_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
promoter_olap=findOverlaps(convert_GR(nme_dt_clu$region),tss,maxgap=2000)
print(ggplot(nme_dt_clu[queryHits(promoter_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

# pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
# ggplot(as.data.frame(mcols(MML_in[queryHits(olap)])),aes(x=density_quant, y=MML))+
#   ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
#   ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()  
# pdf('../downstream/output/graphs/Figure3/FigureS6_CpG_density_MML_boxplot_non_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
# 
# ggplot(as.data.frame(mcols(MML_in[-queryHits(olap)])),aes(x=density_quant, y=MML))+
#   ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
#   ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()  


# Plot density vs NME with points -----------------------------------------
# NME_in_dt=data.table(NME=NME_in$NME,density=NME_in$density)
# digits_round=4
# NME_in_dt_agg=NME_in_dt[, list(NME=median(NME),Bottom25=quantile(NME,probs=0.25),
#                                    top25=quantile(NME,probs=0.75)), 
#                             by = list(density = round(density,digits=digits_round))]
# NME_in_dt_agg$Bottom25= predict(loess(Bottom25~density,NME_in_dt_agg),newdata=NME_in_dt_agg$density)
# NME_in_dt_agg$top25= predict(loess(top25~density,NME_in_dt_agg),newdata=NME_in_dt_agg$density)
# png('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME.png',width=1080,height=1080)#Totally having 69530406 points
# ggplot(NME_in_dt_agg,aes(x=density, y=NME))+
#   ylim(c(0,1))+geom_smooth(method="loess",se=FALSE)+theme_glob+xlab("CpG density")+
#   ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)+
#   theme(axis.title.x=element_text(hjust=0.5,size=48,face="bold"),
#         axis.title.y=element_text(hjust=0.5,size=48,face="bold"),
#         axis.text.x=element_text(size=46),
#         axis.text.y=element_text(size=46))+
#   geom_point(data=NME_in_dt,size=0.1,aes(x=density,y=NME),alpha=0.1)
# dev.off()              


# #LiftOver
# GR_merge=readRDS(GR_merge_file)
# ch = import.chain('../downstream/data/hg19ToMm10.over.chain')
# cur=granges(unique(GR_merge[GR_merge$dMML_pval<=pval_cutoff]))
# seqlevelsStyle(cur) = "UCSC"  # necessary
# cur19 = unlist(liftOver(cur, ch))
# overlap=list()
# overlap_dat=data.table()
# for(fn in dir('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/',pattern='all.csv')){
#   csv_in=fread(paste0('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/',fn))
#   region_in=convert_GR(unique(csv_in$region))
#   ts=gsub('_all.csv','',fn)
#   region_olap=subsetByOverlaps(region_in,cur19)
#   region_olap=paste0(seqnames(region_olap),':',start(region_olap),'-',end(region_olap))
#   overlap[[ts]]=csv_in[region%in%region_olap]
#   overlap_dat=rbind(overlap_dat,data.table(ts=ts,total_regions=length(region_in),overlap=length(region_olap),
#                                            olap_GO=sum(region_olap_GO%in%csv_in[GO_result!=""]$region)))
# }

#Find which one lose CG increase Entropy and motif prefer higher ent
motif_gene <- readRDS(motif_gene_file)

olap=findOverlaps(variant_HetCpG_meta,motif_gene)
variant_HetCpG_meta_motif=variant_HetCpG_meta[queryHits(olap)]
variant_HetCpG_meta_motif$ref_score=motif_gene$scoreRef[subjectHits(olap)]
variant_HetCpG_meta_motif$alt_score=motif_gene$scoreRef[queryHits(olap)]
variant_HetCpG_meta_lose_CG=variant_HetCpG_meta_motif[grepl('CG',variant_HetCpG_meta_motif$REF_tri)&!(grepl('CG',variant_HetCpG_meta_motif$ALT_tri))]

binom.test(sum(variant_HetCpG_meta_lose_CG$ref_score<variant_HetCpG_meta_lose_CG$alt_score),length(variant_HetCpG_meta_lose_CG),
           p=sum(variant_HetCpG_meta_motif$ref_score<variant_HetCpG_meta_motif$alt_score)/length(variant_HetCpG_meta_motif))
#MML and NME curve for 
MML_in=readRDS('../downstream/input/MML_allele_agnostic_merge_20k_homogeneous_excluding_dMML2_CG.rds')
NME_in=readRDS('../downstream/input/NME_allele_agnostic_merge_20k_homogeneous_excluding_dMML2_CG.rds')
NME_in$stat_type="NME"
MML_in$stat_type="MML"
plot_dt=as.data.table(rbind(mcols(NME_in)[,c("score","density","stat_type")],mcols(MML_in)[,c("score","density","stat_type")]))
saveRDS(plot_dt,'../downstream/output/plot_dt_NME_MML.rds')

rm(MML_in)
rm(NME_in)
plot_dt=readRDS('../downstream/output/plot_dt_NME_MML.rds')
pdf('../downstream/output/FigureS_CpG_density_MML_NME_line_density_small1.pdf',width=7,height=7)#Totally having 69530406 points
ggplot(plot_dt[density<=1],aes(x=density, y=score,color=stat_type))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="bottom")
dev.off()
pdf('../downstream/output/FigureS_CpG_density_MML_NME_line_density_all.pdf',width=7,height=7)#Totally having 69530406 points
ggplot(plot_dt,aes(x=density, y=score,color=stat_type))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="bottom")
dev.off()


# preprocess CPEL ---------------------------------------------------------


#cut into pieces to do fastDoCall
UC_in_MDS_comp_analyzed_dt=fastDoCall('rbind',lapply(UC_in_MDS_comp_analyzed,convert_GR,direction="DT"))
UC_in_MDS_analyzed$tissue=sub('-.*','',UC_in_MDS_analyzed$Sample)
UC_in_MDS_analyzed$Sample=sub('.5-.*-E1','.5-E1',UC_in_MDS_analyzed$Sample)
UC_in_MDS_analyzed=UC_in_MDS_analyzed[UC_in_MDS_analyzed$Sample %in% unique(UC_in_MDS_comp$Sample)]
UC_in_MDS_analyzed$tissue=sub('-.*','',UC_in_MDS_analyzed$Sample)
UC_in_MDS_analyzed$Sample=sub('.5-.*-E1','.5-E1',UC_in_MDS_analyzed$Sample)
saveRDS(UC_in_MDS_analyzed,'UC_in_MDS_analyzed.rds')


#UC_run_before_MDS folder
in_dir='../downstream/data/mouse_analysis/UC_run_before_MDS/'
UC_in=GRanges()
UC_in_ls=mclapply(dir(in_dir,pattern = 'mm10.*uc.bedGraph'),function(x){UC_in=read.agnostic.mouse.uc(paste(in_dir,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=24)
UC_in=fastDoCall('c',UC_in_ls)
UC_in=UC_in[UC_in$N>=2]
#saveRDS(UC_in,'UC_agnostic_mouse_dedup_N2_all_time_fix_UC.rds')#74% regiOn have all data
#UC_in_matrix=agnostic_matrix_conversion(UC_in,'UC')#duplicated regions due to error,waiting for one more to finish
UC_in$tissue=sub('-.*','',UC_in$Sample)
UC_in_matrix_ls=mclapply(unique(UC_in$tissue),function(x) agnostic_matrix_conversion(UC_in[UC_in$tissue==x],'UC'),mc.cores=12)
names(UC_in_matrix_ls)=unique(UC_in$tissue)
saveRDS(UC_in_matrix_ls,'../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')#74% regiOn have all data

#read in UC for MDS
in_dir='../downstream/data/mouse_analysis/'
MDS_file=dir(in_dir,pattern = 'mm10.*uc.bedGraph')

gff_in=import.gff3('../mm10_allele_agnostic_analysis.gff')
gff_in=paste0(seqnames(gff_in),':',start(gff_in),'-',end(gff_in))
UC_out=data.table(region=gff_in)

UC_in=fastDoCall('cbind',
                 mclapply(MDS_file,function(x){
                   read.agnostic.mouse.uc(paste(in_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in)},mc.cores=24))
#saveRDS(UC_in,paste0('UC_agnostic_mouse_dedup_MDS_',i,'.rds'))#74% regiOn have all data
UC_out=cbind(UC_out,UC_in)
saveRDS(UC_out,'../downstream/input/UC_agnostic_mouse_N2_MDS_fix.rds')#74% regiOn have all data


# from clustering.R -------------------------------------------------------

#Plot non-repeats

re_web=readRDS('../downstream/output/mouse_analysis/repeats/re_web.rds')
dir_in='../downstream/input/mouse_analysis/clustering/tissue_specific/currently_in_use/ts_cluster_0_1_non_repeats/'
lapply(dir(dir_in,pattern='csv'),function(x){
  
  csv_in=fread(paste0(dir_in,x))
  olap=findOverlaps(convert_GR(csv_in$regions,direction='GR'),re_web,minoverlap =0.5*mean(width(convert_GR(csv_in$regions,direction='GR'))))
  cat('percent removed for ',unique(csv_in$tisue),'is ',length(unique(queryHits(olap)))/nrow(csv_in),'\n')
  write.csv(csv_in[-unique(queryHits(olap))],paste0(dir_in,x),row.names = F)
  
  
  
})

# summary_stat=data.table()
# for(fn in dir(dir_in,pattern='csv')){
#   csv_in=fread(paste0(dir_in,fn))
#   tissue=gsub('.csv','',fn)
#   remain_region=csv_in$region%in%cluster_result[[tissue]]$regions
#   total_region=nrow(csv_in)
#   cluster_result_ts=cluster_result[[tissue]]$cluster
#   names(cluster_result_ts)=cluster_result[[tissue]]$regions
#   csv_in=csv_in[remain_region]
#   #Looking into the detail
#   csv_in$old_cluster=csv_in$cluster
#   csv_in$cluster=cluster_result_ts[csv_in$region]
#   for(clu in 1:10){
#     #Find the major cluster goes into
#     csv_in_clu=csv_in[old_cluster==clu]  
#     major_cluster=names(which.max(table(csv_in_clu$cluster)))
#     
#     summary_stat=rbind(summary_stat,
#                        data.table(remained_region=sum(remain_region),
#                                   total_region=total_region,
#                                   consistent_region=sum(csv_in_clu$cluster==major_cluster),
#                                   tissue=tissue,
#                                   cluster=clu,
#                                   total_region_cluster=nrow(csv_in_clu)))
#   }
#   csv_in$old_cluster=NULL
#   write.csv(csv_in,paste0(dir_out,fn),row.names = F)
# }
# summary_stat$percent_consistent=summary_stat$consistent_region/summary_stat$total_region_cluster
# saveRDS(summary_stat,'../downstream/output/CpG_N17_filter_clustering_consistencyv2.rds')
# write.csv(summary_stat,'../downstream/output/percent_consistent_N17_kmeans.csv')
# hist(summary_stat$percent_consistent,main="kmeans",xlab="proportion of region having consistent assignment")
# # Prepare folder of desired regions: ../downstream/input/ts_cluster_0_1/ ----------------------------------------------
# uc=readRDS('../downstream/output/uc_matrix_DNase.rds')
# analyzed_regions=readRDS('../downstream/output/mm10_allele_agnostic_analysis.rds')
# analyzed_regions_N17_N2=analyzed_regions[analyzed_regions$N>=2&analyzed_regions$N<=17]
# analyzed_regions_N17_N2=paste0(seqnames(analyzed_regions_N17_N2),":",start(analyzed_regions_N17_N2),'-',end(analyzed_regions_N17_N2))
# uc_ft=lapply(uc,function(x) x[rownames(x) %in%analyzed_regions_N17_N2, ])
# for(x in names(uc_ft)){
#   
#   cat(x,':',nrow(uc_ft[[x]])/nrow(uc[[x]]),'\n')
# }
# 
# # EFP : 0.9766308 
# # NT : 0.9742735 
# # forebrain : 0.9732378 
# # heart : 0.9737051 
# # hindbrain : 0.9741982 
# # limb : 0.97602 
# # liver : 0.972034 
# # midbrain : 0.9730241 
# saveRDS(uc_ft,'../downstream/output/uc_matrix_DNase_ft_N17.rds')
# UC_merge=readRDS('../downstream/output/UC_merge_max_loc.rds')
# UC_merge=lapply(UC_merge,function(x) x[rownames(x) %in%analyzed_regions_N17_N2, ])
# saveRDS(UC_merge,'../downstream/output/UC_merge_max_loc_ft_N17.rds')

# tiff(paste0('../downstream/output/heatmap_acrosstissue/',n,'.tiff'),width=3000,height=3000,res=300)
# #png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
# pheatmap(mat,cluster_rows = F,annotation_row = rowann,cluster_cols = F,
#          annotation_col = colann,show_colnames = F,show_rownames = F,
#          gaps_row = row_gap,gaps_col = cumsum(rle(colann[,1])$lengths),
#          annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4,dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)))
# dev.off()
#png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)

# Put regions with other info ----------------------------------------------------

lapply(names(cluster_out),function(x){
  cluster_out_ts=cluster_out[[x]]
  UC_merge_max_loc_sub_ts=UC_merge_max_loc_sub[[x]]
  
})
# Percent left for: EFP 0.9209754 
# Percent left for: NT 0.9990584 
# Percent left for: forebrain 0.9823054 
# Percent left for: heart 0.9347036 
# Percent left for: hindbrain 0.9234645 
# Percent left for: limb 0.9862634 
# Percent left for: liver 0.8056983 
# Percent left for: midbrain 0.9104696 

cluster_result=readRDS(cluster_region_out_fn)

#Pre-process CPEL
# 
# rm(list=ls())
# library(data.table)
# library(rtracklayer)
# library(Gmisc)
# #Read in FeDMR 
# dir_in='../downstream/input/FeDMR/'
# FeDMR_fn=dir(dir_in,pattern='tsv')
# FeDMR_in=fastDoCall('c',lapply(dir(dir_in,pattern='tsv'),function(fn){
#   sp=gsub('.tsv','',fn)
#   sp=gsub("feDMR_","",sp)
#   tissue=gsub(".*_","",sp)
#   stage=gsub(paste0("_",tissue),"",sp)
#   stage=gsub("_5",".5",stage)
#   tissue=gsub("craniofacial","EFP",tissue)
#   tissue=gsub("tube","NT",tissue)
#   fn_in=fread(paste0(dir_in,fn))
#   fn_in=makeGRangesFromDataFrame(fn_in,seqnames.field = "chrom",keep.extra.columns = T)
#   fn_in$tissue=tissue
#   fn_in$stage=stage
#   fn_in$sample=paste(tissue,stage,sep='-')
#   return(fn_in)
#   
# })
# )
# FeDMR_in_mcols=as.data.table(mcols(FeDMR_in))
# FeDMR_in_mcols$tissue_exist=TRUE
# FeDMR_in_gr=unique(FeDMR_in)
# FeDMR_in_gr=FeDMR_in_gr[order(FeDMR_in_gr$dmr_id,decreasing=T)]
# mcols(FeDMR_in_gr)=mcols(FeDMR_in_gr)[,c("dmr_id","score","tissue_specificity","stage_specificity")]
# FeDMR_in_gr$dmr_id_original=FeDMR_in_gr$dmr_id
# FeDMR_in_gr$dmr_id=NULL
# FeDMR_in_mcols_dc=dcast.data.table(FeDMR_in_mcols,dmr_id~sample,value.var = "tissue_exist",fill=FALSE)
# FeDMR_in_mcols_dc=FeDMR_in_mcols_dc[order(FeDMR_in_mcols_dc$dmr_id,decreasing = T)]
# mcols(FeDMR_in_gr)=cbind(mcols(FeDMR_in_gr),FeDMR_in_mcols_dc)
# print(which(FeDMR_in_gr$dmr_id_original!=FeDMR_in_gr$dmr_id))
# FeDMR_in_gr$dmr_id_original=NULL
# saveRDS(FeDMR_in_gr,'../downstream/output/FeDMR.rds')
# # percent_cov=NULL
# # UC_in=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')
# # for(tissue in unique(FeDMR_in_mcols$tissue)){
# # percent_cov=rbind(percent_cov,data.frame(tissue=tissue,
# #   percent=length(subsetByOverlaps(FeDMR_in_gr[apply(mcols(FeDMR_in_gr)[,gsub('-.*','',colnames(mcols(FeDMR_in_gr)))==tissue],1,any)],UC_in[[tissue]]))/
# #     length(FeDMR_in_gr[apply(mcols(FeDMR_in_gr)[,gsub('-.*','',colnames(mcols(FeDMR_in_gr)))==tissue],1,any)])
# #   ))
# # }
# 
# #Read in chromHMM
# 
# #ChromHMM
# 
# chromHMM_pooled= read_chromHMM_bed('../downstream/input/chromHMM_mm10/pooled/','pooled')
# saveRDS(chromHMM_pooled,'../downstream/output/chromHMM.rds')
# chromHMM_enhancer=chromHMM_pooled[sub('-.*','',as.character(chromHMM_pooled$chrom_state))=="En"]
# saveRDS(chromHMM_enhancer,'../downstream/output/chromHMM_enhancer.rds')



# GO analysis -------------------------------------------------------------

# #filter by repeats
# 
# dir_in='../downstream/input/mouse_analysis/clustering/tissue_specific/currently_in_use/ts_cluster_0_1_non_repeats/region_assigned/'
# re_web=readRDS('../downstream/output/mouse_analysis/repeats/re_web.rds')
# cor_dt_filtered=readRDS('../downstream/output/mouse_analysis/correlation/correlation_dt_N17_kmeans_10run_filtered_all_regions_non_repeats.rds')
# names(cor_dt_filtered)=NULL
# olap_enhancer=findOverlaps(enhancer,re_web,minoverlap = mean(width(convert_GR(do.call('rbind',cor_dt_filtered)$region))))
# enhancer_bg=subsetByOverlaps(enhancer[unique(queryHits(olap_enhancer))],uc_gr)
# bg_enhancer=unique(enhancer_bg$`Target Gene`)
# olap_tss=findOverlaps(enhancer,re_web,maxgap = 2000)
# bg_promoter=names(subsetByOverlaps(tss[unique(queryHits(olap_tss))],uc_gr,maxgap = 2000))

# # Write motif analysis result ---------------------------------------------
# motif_dir='../downstream/input/mouse_analysis/motif_analysis/mouse_motif_enrichment_0510/'
# motif_all=data.table()
# for(fn in dir(motif_dir,pattern='csv')){
#   stat=gsub('.csv','',gsub('.*_','',fn))
#   csv_in=fread(paste0(motif_dir,fn))
#   csv_in$motif=gsub('.*_','',csv_in$motif)
#   csv_out=csv_in[,grepl("motif|FDR",colnames(csv_in)),with=F]
#   colnames(csv_out)=c('motif','FDR')
#   csv_out$stat=stat
#   csv_out$tissue=gsub('_.*','',fn)
#   motif_all=rbind(motif_all,csv_out)
#   }
# write.csv(motif_all[FDR<=0.1,list(tissue,motif,FDR,stat)],'../downstream/output/mouse_analysis/motif_analysis/motif_all_Ken.csv')


# Mouse motif preprocessing.R ---------------------------------------------

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


# promoter_vs_enhancer_mouse ----------------------------------------------

#max dMML/dNME in non-adjacent vs adjacent 

csv_out$dMML_max_time=gsub('E','',gsub('\\.5','',csv_out$dMML_max_time))
csv_out$dMML_max_time_diff=abs(as.numeric(gsub('-.*','',csv_out$dMML_max_time))-as.numeric(gsub('.*-','',csv_out$dMML_max_time)))
csv_out$dNME_max_time=gsub('E','',gsub('\\.5','',csv_out$dNME_max_time))
csv_out$dNME_max_time_diff=abs(as.numeric(gsub('-.*','',csv_out$dNME_max_time))-as.numeric(gsub('.*-','',csv_out$dNME_max_time)))
csv_out$UC_max_time=gsub('E','',gsub('\\.5','',csv_out$UC_max_time))
csv_out$UC_max_time_diff=abs(as.numeric(gsub('-.*','',csv_out$UC_max_time))-as.numeric(gsub('.*-','',csv_out$UC_max_time)))
#boxplot time of maximum value
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dMML_max_time_diff,fill=states))+geom_boxplot()+ylab("dMML time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dNME_max_time_diff,fill=states))+geom_boxplot()+ylab("dNME time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=UC_max_time_diff,fill=states))+geom_boxplot()+ylab("UC time change")

ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=dMML_max_time_diff,color=states))+geom_density()+ylab("dMML time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=dNME_max_time_diff,color=states))+geom_density()+ylab("dNME time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=UC_max_time_diff,color=states))+geom_density()+ylab("dNME time change")
#From mouse example.R
# mouse_motif_overlap -----------------------------------------------------
# tissue_sel=c("forebrain","heart",'limb','midbrain','hindbrain','liver','EFP')
# GO_out_all=readRDS('../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run/GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_enhancer.rds')
# GO_out_all=lapply(GO_out_all$all,function(x) do.call(rbind,lapply(x,function(y) 
#   cbind(y$csv_in_ts_clu[,gsub('UC-','',colnames(y$csv_in_ts_clu)) %in% paste0('E',10:20,'.5-E',11:21,'.5'),with=FALSE],
#         y$csv_in_ts_clu[,gsub('dMML-','',colnames(y$csv_in_ts_clu)) %in% paste0('E',10:20,'.5-E',11:21,'.5'),with=FALSE],
#         y$csv_in_ts_clu[,gsub('dNME-','',colnames(y$csv_in_ts_clu)) %in% paste0('E',10:20,'.5-E',11:21,'.5'),with=FALSE],
#    #y$csv_in_ts_clu[,grepl("UC-|dNME-|dMML-",colnames(y$csv_in_ts_clu)),with=FALSE],
#     y$csv_in_ts_clu[,list(cluster,region,region_type)]))))
# cor_in=readRDS('../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_ref11.rds')
# GO_out_all=lapply(GO_out_all,function(x) {
#   uc_dt=  x[,grepl("UC-",colnames(x)),with=FALSE]
#   dNME_dt=  x[,grepl("dNME-",colnames(x)),with=FALSE]
#   dMML_dt=  x[,grepl("dMML-",colnames(x)),with=FALSE]
#   uc_max=apply(uc_dt,1,which.max)
#   x$UC_max_time=gsub('UC-','',colnames(uc_dt)[uc_max])
#   x$dNME_max_UC_pair=as.data.frame(dNME_dt)[cbind(seq_along(uc_max), uc_max)]
#   x$UC_max_UC_pair=as.data.frame(uc_dt)[cbind(seq_along(uc_max), uc_max)]
#   x$dMML_max_UC_pair=as.data.frame(dMML_dt)[cbind(seq_along(uc_max), uc_max)]
#   return(x)
# })
# enhancer_regions=readRDS('../downstream/output/mouse_analysis/enhancers/GO_regions_enchancer.rds')

# CpG_mm10=CpG_mm10=getCpgSitesmm10()
# end(CpG_mm10)=start(CpG_mm10)+1
# motif_locus_ken_CG=lapply(motif_locus_ken,function(motifs) {
#   for(mt in names(motifs)){
#     if(length(motifs[[mt]])>0){
#     motifs[[mt]]$mouse_region=motifs[[mt]]$region
#     motifs[[mt]]$region=NULL
#     motifs[[mt]]$motif=mt
#     motifs[[mt]]$NCG=countOverlaps(motifs[[mt]],CpG_mm10,minoverlap = 2)
#     }
#     
#   }
#   return(do.call(c,motifs))
#   
# })
# motif_locus_ken_CG=lapply(motif_locus_ken_CG,function(x){
#   
#   names(x)=NULL
#   return(do.call('c',x))
#   
# }
#   
#   )
# motif_locus_ken_CG_olap=lapply(motif_locus_ken_CG,function(x) x[x$NCG>0])
# motif_locus_ken_CG_olap=lapply(motif_locus_ken_CG_olap,convert_GR,dir="DT")
# 
# 
# saveRDS(motif_locus_ken_CG,'../downstream/output/mouse_analysis/motif_analysis/tissue_region_motif_all_locus_CG.rds')
  # enhancer_regions_ts$dNME_max_UC_pair =GO_out_all[[ts]][match(enhancer_regions_ts$region,region)]$dNME_max_UC_pair
  # enhancer_regions_ts$dMML_max_UC_pair =GO_out_all[[ts]][match(enhancer_regions_ts$region,region)]$dMML_max_UC_pair
  # enhancer_regions_ts$dNME_UC_cor=cor_in[[ts]][match(enhancer_regions_ts$region,region)]$dNME_cor
  # enhancer_regions_ts$dMML_UC_cor=cor_in[[ts]][match(enhancer_regions_ts$region,region)]$dMML_cor

      
    #if(length(PMID_unique)>300){
      #Counting number of citations
    # breaks=ceiling(length(PMID_unique)/300)
    # PMID_cut=cut(1:length(PMID_unique),breaks,label=FALSE)
    # #Sys.sleep(300)
    # pmc_ref_count=c()
    # for(idx in 1:breaks){
    # 
    #   pmc_ref_count=c(pmc_ref_count,do.call(c,lapply(entrez_summary(db="pubmed", id=PMID_unique[PMID_cut==idx]),function(x) x$pmcrefcount)))
    #   Sys.sleep(120)
    # }
    # }else( pmc_ref_count=do.call(c,lapply(entrez_summary(db="pubmed", id=PMID_unique),function(x) x$pmcrefcount)))
    # 
    # gene_pub_med$pmc_count=pmc_ref_count[gene_pub_med$PMID]
    #gene_pub_med_collapse=gene_pub_med[,list(PMID=paste(PMID,collapse=";"),pmc_count=paste(pmc_count,collapse = ';'),title=paste(title,collapse=' ; ')),by=gene]
    # Writing output ----------------------------------------------------------
# 
# 
# for(ts in names(enhancer_regions_motif_dNME_all)){
#   write.csv(enhancer_regions_motif_dNME_all[[ts]][,list(region,tissue,dNME_motif,cluster,gene,region_type,
#                                               dNME_max_UC_pair ,dMML_max_UC_pair,UC_max_time,
#                                               UC_max_time,UC_max_pair,dNME_UC_cor,
#                                               dMML_UC_cor,dNME_UC_cor,PMID,title,GO_ID)],
#             paste0(region_motif_dir,ts,'_dNME_Pubmed_annotated.csv'),
#             row.names = F,quote=F)
#   
#   
# }
# write.csv(do.call(rbind,lapply(enhancer_regions_motif_dNME_all,function(x) 
#   x[region_type=="dNME_only",list(tissue,dNME_motif,gene,UC_max_time,
#                                   dNME_max_UC_pair,dNME_UC_cor,
#                                   dMML_max_UC_pair,dMML_UC_cor,
#                                   UC_max_pair,
#                                   PMID,GO_ID)])),
#   '../downstream/output/mouse_analysis/motif_analysis/region_motif/all_regions_dNME_only.csv',row.names = F,quote=F)
# write.csv(do.call(rbind,lapply(enhancer_regions_motif_dMML_all,function(x) 
#   x[region_type=="dMML_only",list(tissue,dMML_motif,gene,UC_max_time,
#                                   dNME_max_UC_pair,dNME_UC_cor,
#                                   dMML_max_UC_pair,dMML_UC_cor,
#                                   UC_max_pair,
#                                   PMID,GO_ID)])),
#   '../downstream/output/mouse_analysis/motif_analysis/region_motif/all_regions_dMML_only.csv',row.names = F,quote=F)
# #Read in and rank the dNME motif and dMML motif separately
# dir_motif_gene='../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/'
# #dNME only motif, need to have annotated GO terms
# motif_dNME_out=data.table()
# for(fn in dir(dir_motif_gene,pattern='dNME')){
#   csv_in=fread(paste0(dir_motif_gene,fn))
#   # need to have annotated GO terms
#   GO_in=fread(paste0('../downstream/output/mouse_analysis/GO_analysis/gene_examples/',unique(csv_in$tissue),'.csv'))
#   csv_in=csv_in[region %in% GO_in$region&region_type=="dNME_only"]
#   csv_in$GO_terms=GO_in[match(csv_in$region,region)]$GO_terms
#   #Change to dNME accordingly
#   csv_in$rank1=csv_in[,list(rank=rank(dNME_UC_cor)),by=motif]$rank
#   csv_in$rank2=csv_in[,list(rank=rank(dNME_max_UC_pair)),by=motif]$rank
#   motif_dNME_out=rbind(motif_dNME_out,csv_in[,.SD[order(rank1+rank2,decreasing=T)],by=motif])
# 
# }
# write.csv(motif_dNME_out[,list(motif,region,gene,tissue,
#                                 UC_max_time ,
#                                 dNME_max_UC_pair,UC_max_pair,dMML_max_UC_pair,
#                                 dNME_UC_cor,dMML_UC_cor,
#                                 cluster,GO_terms)],'../downstream/output/mouse_analysis/motif_analysis/mouse_motif_example_GO_dNME.csv',
#           row.names = F)
# 
# 
# #Check the distribution in forebrain
# 
# forebrain_dMML=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/Forebrain_dMML_temp.csv')
# forebrain_dNME=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/Forebrain_dNME_temp.csv')
# 
# # archived ----------------------------------------------------------------
# 
# 
# #dNME region
# region_sample=list()
# set.seed(123)
# for(ts in names(GO_out_all)){
# 
#   motif_dNME=Ken_motif_merge(readRDS(paste0('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_locus/',ts,'_motif_site_dNME.rds')))
#   if(!is.null(motif_dNME)){
#     GO_out_all_ts=GO_out_all[[ts]]
#     GO_out_all_ts_dNME=GO_out_all_ts[region_type=="dNME_only"]
# 
#     GO_out_all_ts_dNME$motif_dNME=convert_GR(motif_dNME,direction="DT")[match(GO_out_all_ts_dNME$region,region)]$motif
#     GO_out_all_ts_dNME$tissue=ts
#     GO_out_all_ts_dNME_motif=GO_out_all_ts_dNME[!is.na(motif_dNME)]
#     region_sample[[ts]]=GO_out_all_ts_dNME_motif[sample(1:nrow(GO_out_all_ts_dNME_motif),20)]
#   }
# }
# #dMML region
# region_sample=list()
# set.seed(123)
# for(ts in names(GO_out_all)){
# 
#   motif_dMML=Ken_motif_merge(readRDS(paste0('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_locus/',ts,'_motif_site_dMML.rds')))
#   if(!is.null(motif_dMML)){
#     GO_out_all_ts=GO_out_all[[ts]]
#     GO_out_all_ts_dMML=GO_out_all_ts[region_type=="dMML_only"]
# 
#     GO_out_all_ts_dMML$motif_dMML=convert_GR(motif_dMML,direction="DT")[match(GO_out_all_ts_dMML$region,region)]$motif
#     GO_out_all_ts_dMML$tissue=ts
#     GO_out_all_ts_dMML_motif=GO_out_all_ts_dMML[!is.na(motif_dMML)]
#     if(nrow(GO_out_all_ts_dMML_motif)>20){
#     region_sample[[ts]]=GO_out_all_ts_dMML_motif[sample(1:nrow(GO_out_all_ts_dMML_motif),20)]
#     }else(region_sample[[ts]]=GO_out_all_ts_dMML_motif)
#   }
# }
# write.csv(region_sample$heart[,list(region,UC_max_time,UC_max_UC_pair,dMML_max_UC_pair,dNME_max_UC_pair,motif_dMML,tissue)],
#           '../downstream/output/mouse_analysis/examples/dMML_heart_GO_motif.csv')


# Writing output ----------------------------------------------------------
# 
# 
# for(ts in names(enhancer_regions_motif_dNME_all)){
#   write.csv(enhancer_regions_motif_dNME_all[[ts]][,list(region,tissue,dNME_motif,cluster,gene,region_type,
#                                               dNME_max_UC_pair ,dMML_max_UC_pair,UC_max_time,
#                                               UC_max_time,UC_max_pair,dNME_UC_cor,
#                                               dMML_UC_cor,dNME_UC_cor,PMID,title,GO_ID)],
#             paste0(region_motif_dir,ts,'_dNME_Pubmed_annotated.csv'),
#             row.names = F,quote=F)
#   
#   
# }
# write.csv(do.call(rbind,lapply(enhancer_regions_motif_dNME_all,function(x) 
#   x[region_type=="dNME_only",list(tissue,dNME_motif,gene,UC_max_time,
#                                   dNME_max_UC_pair,dNME_UC_cor,
#                                   dMML_max_UC_pair,dMML_UC_cor,
#                                   UC_max_pair,
#                                   PMID,GO_ID)])),
#   '../downstream/output/mouse_analysis/motif_analysis/region_motif/all_regions_dNME_only.csv',row.names = F,quote=F)
# write.csv(do.call(rbind,lapply(enhancer_regions_motif_dMML_all,function(x) 
#   x[region_type=="dMML_only",list(tissue,dMML_motif,gene,UC_max_time,
#                                   dNME_max_UC_pair,dNME_UC_cor,
#                                   dMML_max_UC_pair,dMML_UC_cor,
#                                   UC_max_pair,
#                                   PMID,GO_ID)])),
#   '../downstream/output/mouse_analysis/motif_analysis/region_motif/all_regions_dMML_only.csv',row.names = F,quote=F)
# #Read in and rank the dNME motif and dMML motif separately
# dir_motif_gene='../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/'
# #dNME only motif, need to have annotated GO terms
# motif_dNME_out=data.table()
# for(fn in dir(dir_motif_gene,pattern='dNME')){
#   csv_in=fread(paste0(dir_motif_gene,fn))
#   # need to have annotated GO terms
#   GO_in=fread(paste0('../downstream/output/mouse_analysis/GO_analysis/gene_examples/',unique(csv_in$tissue),'.csv'))
#   csv_in=csv_in[region %in% GO_in$region&region_type=="dNME_only"]
#   csv_in$GO_terms=GO_in[match(csv_in$region,region)]$GO_terms
#   #Change to dNME accordingly
#   csv_in$rank1=csv_in[,list(rank=rank(dNME_UC_cor)),by=motif]$rank
#   csv_in$rank2=csv_in[,list(rank=rank(dNME_max_UC_pair)),by=motif]$rank
#   motif_dNME_out=rbind(motif_dNME_out,csv_in[,.SD[order(rank1+rank2,decreasing=T)],by=motif])
# 
# }
# write.csv(motif_dNME_out[,list(motif,region,gene,tissue,
#                                 UC_max_time ,
#                                 dNME_max_UC_pair,UC_max_pair,dMML_max_UC_pair,
#                                 dNME_UC_cor,dMML_UC_cor,
#                                 cluster,GO_terms)],'../downstream/output/mouse_analysis/motif_analysis/mouse_motif_example_GO_dNME.csv',
#           row.names = F)
# 
# 
# #Check the distribution in forebrain
# 
# forebrain_dMML=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/Forebrain_dMML_temp.csv')
# forebrain_dNME=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/Forebrain_dNME_temp.csv')
# 

# # archived ----------------------------------------------------------------
# 
# 
# #dNME region
# region_sample=list()
# set.seed(123)
# for(ts in names(GO_out_all)){
# 
#   motif_dNME=Ken_motif_merge(readRDS(paste0('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_locus/',ts,'_motif_site_dNME.rds')))
#   if(!is.null(motif_dNME)){
#     GO_out_all_ts=GO_out_all[[ts]]
#     GO_out_all_ts_dNME=GO_out_all_ts[region_type=="dNME_only"]
# 
#     GO_out_all_ts_dNME$motif_dNME=convert_GR(motif_dNME,direction="DT")[match(GO_out_all_ts_dNME$region,region)]$motif
#     GO_out_all_ts_dNME$tissue=ts
#     GO_out_all_ts_dNME_motif=GO_out_all_ts_dNME[!is.na(motif_dNME)]
#     region_sample[[ts]]=GO_out_all_ts_dNME_motif[sample(1:nrow(GO_out_all_ts_dNME_motif),20)]
#   }
# }
# #dMML region
# region_sample=list()
# set.seed(123)
# for(ts in names(GO_out_all)){
# 
#   motif_dMML=Ken_motif_merge(readRDS(paste0('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_locus/',ts,'_motif_site_dMML.rds')))
#   if(!is.null(motif_dMML)){
#     GO_out_all_ts=GO_out_all[[ts]]
#     GO_out_all_ts_dMML=GO_out_all_ts[region_type=="dMML_only"]
# 
#     GO_out_all_ts_dMML$motif_dMML=convert_GR(motif_dMML,direction="DT")[match(GO_out_all_ts_dMML$region,region)]$motif
#     GO_out_all_ts_dMML$tissue=ts
#     GO_out_all_ts_dMML_motif=GO_out_all_ts_dMML[!is.na(motif_dMML)]
#     if(nrow(GO_out_all_ts_dMML_motif)>20){
#     region_sample[[ts]]=GO_out_all_ts_dMML_motif[sample(1:nrow(GO_out_all_ts_dMML_motif),20)]
#     }else(region_sample[[ts]]=GO_out_all_ts_dMML_motif)
#   }
# }
# write.csv(region_sample$heart[,list(region,UC_max_time,UC_max_UC_pair,dMML_max_UC_pair,dNME_max_UC_pair,motif_dMML,tissue)],
#           '../downstream/output/mouse_analysis/examples/dMML_heart_GO_motif.csv')
#Get all motif binding site for motifs in Ken's list and factor book and find examples
#Need to have few brown SNPs 
#Need to have large examples
#Looking for genes related to those tissue
#Start with heart and forebrain
#Factorbook
# #CHip_atlas: GATA4 mouse embyro
# factor_in=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_heart_GATA4.bed',skip = 1,sep='\t')
# colnames(factor_in)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
# factor_in=factor_in[,list(seqnames,start,end,metadata,log10qval)]
# factor_in$metadata=gsub('%20','_',factor_in$metadata)
# factor_in$metadata=gsub('%3','_',factor_in$metadata)
# #gata4 in heart
# factor_in_GATA4=factor_in[grepl("GATA",factor_in$metadata)]
# factor_in_GATA4_gr=makeGRangesFromDataFrame(factor_in_GATA4)
#pre_process_CPEL
# MDS for DNase only regions ----------------------------------------------

DNase_conrol_MDS_dir='../downstream/data/DNase_control_PRC_MDS_mouse/'
gff_in_DNase=import.gff3(mouse_DNase_control_gff_file)
gff_in_DNase=paste0(seqnames(gff_in_DNase),':',start(gff_in_DNase),'-',end(gff_in_DNase))
UC_in_analyzed_MDS=data.table(region=gff_in_DNase)
UC_in_analyzed_MDS_UC=fastDoCall('cbind',
                                 mclapply(dir(DNase_conrol_MDS_dir,pattern = '.*uc.bedGraph'),function(x){
                                   read.agnostic.mouse.uc(paste(DNase_conrol_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_DNase)},mc.cores=10))
UC_in_analyzed_MDS=cbind(UC_in_analyzed_MDS,UC_in_analyzed_MDS_UC)
saveRDS(UC_in_analyzed_MDS,'../downstream/output/mouse_analysis/CPEL_outputs/UC_MDS_N2_DNase_control_PRC.rds')
DNase_mm10=readRDS(DNase_mm10_file)
DNase_mm10=convert_GR(DNase_mm10,dir='DT')
UC_in_analyzed_MDS_DNase=UC_in_analyzed_MDS[region%in% DNase_mm10$region]
saveRDS(UC_in_analyzed_MDS_DNase,'../downstream/output/mouse_analysis/CPEL_outputs/UC_MDS_N2_DNase_only.rds')
############From mouse_ChIP_seq_overlap.R
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

# Find examples by looking at overlap between motif and ChIP data -----------
enhancer_regions_motif_dNME_all=readRDS(enhancer_motif_all_dNME_fn)
enhancer_regions_motif_dMML_all=readRDS(enhancer_motif_all_dMML_fn)
#Heart 
factor_in_embyro_heart_dNME=readRDS(paste0(ChiP_motif_dir,'ChIP_olap_heart_embyro_heat_all.rds'))
factor_in_embyro_heart_dNME=factor_in_embyro_heart_dNME$dNME_region
factor_in_heart_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'ChIP_olap_heart_adult_heart_all.rds')))
heart_dNME_pubmed=enhancer_regions_motif_dNME_all$heart
heart_dNME_pubmed=heart_dNME_pubmed[,list(region,tissue,cluster,gene,UC_max_time,region_type,PMID,dNME_max_pair,dMML_max_pair,UC_max_pair,dMML_motif,dNME_motif)]
olap=findOverlaps(convert_GR(heart_dNME_pubmed$region),convert_GR(factor_in_embyro_heart_dNME$region))
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
forebrain_dNME_pubmed_motif$motif_in_ChIP=forebrain_dNME_pubmed_motif[,list(motif_in_ChIP=grepl(gsub(';','|',dNME_motif),motif_ChIP,ignore.case = T)),
by=list(region,motif_ChIP)]$motif_in_ChIP

olap=findOverlaps(convert_GR(forebrain_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_forebrain_dNME))
forebrain_dNME_pubmed_motif=forebrain_dNME_pubmed[queryHits(olap)]
forebrain_dNME_pubmed_motif$motif_ChIP=factor_in_forebrain_dNME[subjectHits(olap)]$metadata
forebrain_dNME_pubmed_motif$motif_in_ChIP=forebrain_dNME_pubmed_motif[,list(motif_in_ChIP=grepl(gsub(';','|',dNME_motif),motif_ChIP,ignore.case = T)),
by=list(region,motif_ChIP)]$motif_in_ChIP

############From mouse_ChIP_seq_overlap.R###########end############
#NME_MAV#########################Start######################

# # NME vs VMR currently not in use--------------------------------------------------------------
# NME_in=readRDS(NME_agnostic_file)
# #Brain
# load("../downstream/input/vmrs_hg19_brain.rda")
# vmr_HC2=vmrs_hg19$HC2
# vmr_HC1=vmrs_hg19$HC1
# names(vmr_HC2)=NULL
# names(vmr_HC1)=NULL
# #Do HC2
# vmr=do.call(c,vmr_HC2)
# saveRDS(vmr,'../downstream/output/vmr_HC2.rds')
# vmr=readRDS('../downstream/output/vmr_HC2.rds')
# NME_in_brain=NME_in[NME_in$Sample%in%c('Brain_Hippocampus_middle_paired - 149','Brain_Hippocampus_middle_paired - 150')]
# OR_quant=data.frame()
# for(percent in unique(NME_in_brain$quant_score)){
#   OR=OR_VMR(NME_in_brain,vmr,percent,NME_quant='quant_score')
#   OR_quant=rbind(OR_quant,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
#   
# }
# write.csv(OR_quant,'../downstream/output/brain_VMR.csv')
# 
# #chromHMM state get enhancer
# ah = AnnotationHub()
# ENCODE_name=ENCODE_to_sample(unique(NME_in_brain$Sample))
# ah_num=names(query(ah, c("chromhmmSegmentations", unique(ENCODE_name$ENCODE))))
# chromHMM=ah[[ah_num]]
# chromHMM_enc=chromHMM[chromHMM$abbr%in%c("7_Enh","6_EnhG")]
# NME_in_brain_enc=subsetByOverlaps(NME_in_brain,chromHMM_enc)
# vmr_enc=subsetByOverlaps(vmr,chromHMM_enc)
# #Find nearest genes
# genomic_features=readRDS(genomic_features_file)
# genomic_features=genomic_features$TSS
# NME_in_brain_enc=dist_calc(NME_in_brain_enc,genomic_features)
# vmr_enc=dist_calc(vmr_enc,genomic_features)
# #Find overlap between nearest genes
# NME_in_brain_enc_highNME=NME_in_brain_enc[NME_in_brain_enc$NME>=quantile(NME_in_brain$NME,prob=0.99)]#0.87 cutoff
# NME_in_brain_enc_highNME=as.data.table(mcols(NME_in_brain_enc_highNME))
# NME_in_brain_enc_highNME=NME_in_brain_enc_highNME[,list(NME=mean(NME),dist=mean(abs(dist))),by=list(gene)]
# NME_in_brain_enc_highNME=NME_in_brain_enc_highNME[order(NME,decreasing = T)]$gene
# bg=unique(c(vmr_enc$gene,NME_in_brain_enc$gene))
# NME_VMR=sum((bg %in%NME_in_brain_enc_highNME)&(bg %in%vmr_enc$gene))
# NME_nonVMR=sum((bg %in%NME_in_brain_enc_highNME)&!(bg %in%vmr_enc$gene))
# VMR_nonNME=sum(!(bg %in%NME_in_brain_enc_highNME)&(bg %in%vmr_enc$gene))
# nonVMR_nonNME=sum(!(bg %in%NME_in_brain_enc_highNME)&!(bg %in%vmr_enc$gene))
# 
# fisher.test(matrix(c(NME_VMR,NME_nonVMR,VMR_nonNME,nonVMR_nonNME),nrow=2))
# 
# write(unique(NME_in_brain_enc_highNME[NME_in_brain_enc_highNME%in%vmr_enc$gene]),'../downstream/output/VMR_NME.txt')
#Hypervariance for expression entropy pre-processing
GR_merge=readRDS(GR_merge_file)
#Only use merged data for H1
GR_merge=GR_merge[!(GR_merge$Sample%in%c("rep1 - H1","rep2 - H1"))]
genomic_features=readRDS(genomic_features_file)
