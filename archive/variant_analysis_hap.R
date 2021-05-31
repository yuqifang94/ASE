# Genomics
# Source main functions
#setwd("~/code/HASM-MetaAnalysis/")
source("mainFunctions_sub.R")

# Find number of overlapped regions ---------------------------------------
GR_merge=readRDS(GR_merge_file)
#GR_merge=readRDS('../downstream/output/GR_merge_CpG_CGcont.rds')
dMML=sum(GR_merge$dMML_pval<=pval_cutoff)
dNME=sum(GR_merge$dNME_pval<=pval_cutoff)
dMML_dNME=sum(GR_merge$dNME_pval<=pval_cutoff&GR_merge$dMML_pval<=pval_cutoff)
#Run it tonight
olap_merge=c()
for(i in 1:50000){olap_merge=c(olap_merge,length(subsetByOverlaps(GR_merge[sample(1:length(GR_merge),dMML,replace=F)],
                                          GR_merge[sample(1:length(GR_merge),dNME,replace=F)],type='equal')))}
sum(olap_merge<=dMML_dNME)/length(olap_merge)#=0
saveRDS(olap_merge,'../downstream/output/olap_merge.rds')

# Plotting dNME vs dMML and NME vs MML ------------------------------------
#NME vs MML
GR_merge=readRDS(GR_merge_file)
qt_005<-function(x){return(quantile(x,probs=0.25))}
qt_095<-function(x){return(quantile(x,probs=0.75))}
#GR dNME vs dMML at ASM
GR_merge_df=rbind(data.frame(MML=GR_merge$MML1,NME=GR_merge$NME1,dNME_pval=GR_merge$dNME_pval,dMML_pval=GR_merge$dMML_pval),
                  data.frame(MML=GR_merge$MML2,NME=GR_merge$NME2,dNME_pval=GR_merge$dNME_pval,dMML_pval=GR_merge$dMML_pval))
#Aggregate NME, using quantiles, 0.05 and 0.95
GR_NME_df_agg=aggregate(GR_merge_df$NME,by=list(round(GR_merge_df$MML,digits=2)),FUN=median)
GR_NME_df_agg$low_CI=aggregate(round(GR_merge_df$NME,digits=2),by=list(round(GR_merge_df$MML,digits=2)),FUN=qt_005)$x
GR_NME_df_agg$high_CI=aggregate(round(GR_merge_df$NME,digits=2),by=list(round(GR_merge_df$MML,digits=2)),FUN=qt_095)$x
names(GR_NME_df_agg)=c("MML","median","Bottom25","top25")
GR_NME_df_agg$Bottom25= predict(loess(Bottom25~MML,GR_NME_df_agg),newdata=GR_NME_df_agg$MML)
GR_NME_df_agg$top25= predict(loess(top25~MML,GR_NME_df_agg),newdata=GR_NME_df_agg$MML)
#GR_NME_df_agg=melt(GR_NME_df_agg,id="MML")
pdf('../downstream/output/graphs/Figure1/NME_vs_MML_with_quantile.pdf')
#Plotting
ggplot(GR_NME_df_agg,aes(x=MML, y=median))+
  xlim(c(0,1))+ylim(c(0,1))+ggtitle("MML and NME relationship")+geom_smooth(method="loess",se=FALSE)+
  ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)
dev.off()
#dNME vs dMML
GR_merge_df=data.frame(dMML=GR_merge$dMML,dNME=GR_merge$dNME,dNME_pval=GR_merge$dNME_pval,dMML_pval=GR_merge$dMML_pval,N=GR_merge$N)
GR_merge_df$MML_Hap=GR_merge_df$dMML_pval<=pval_cutoff
ggplot(GR_merge_df,aes(x=dNME, color=MML_Hap))+
  xlim(c(0,1))+ggtitle("dNME distribution at dMML region")+geom_density( kernel = "gaussian")
GR_df_agg=aggregate(GR_merge_df$dNME,by=list(round(GR_merge_df$dMML,digits=4)),FUN=median)
GR_df_agg$low_CI=aggregate(round(GR_merge_df$dNME,digits = 4),by=list(round(GR_merge_df$dMML,digits=4)),FUN=qt_005)$x
GR_df_agg$high_CI=aggregate(round(GR_merge_df$dNME,digits=4),by=list(round(GR_merge_df$dMML,digits=4)),FUN=qt_095)$x
names(GR_df_agg)=c("dMML","median","Bottom25","top25")
GR_df_agg$Bottom25= predict(loess(Bottom25~dMML,GR_df_agg),newdata=GR_df_agg$dMML)
GR_df_agg$top25= predict(loess(top25~dMML,GR_df_agg),newdata=GR_df_agg$dMML)
#GR_df_agg=melt(GR_df_agg,id="dMML")
###plotting
pdf('../downstream/output/graphs/Figure1/dNME_vs_dMML_quantile_differential_median.pdf')
ggplot(GR_df_agg,aes(x=dMML, y=median))+
  xlim(c(0,1))+ylim(c(0,0.7))+ggtitle("dMML and dNME relationship")+geom_smooth(method="loess",se=FALSE)+
  ylab("dNME")+theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)+
  scale_linetype_manual(values=c("solid","twodash", "twodash"))+scale_color_manual(values=c("Blue","Blue","Blue"))
dev.off()

# allele-agnositic NME vs dMML and dNME --------------------------------------------
###Test for dNME vs anagonisitc NME
GR_merge=readRDS(GR_merge_file)
GR_calc=data.frame()
in_dir='../downstream/data/allele_agnostic/'
GR_merge$fn=paste(GR_merge$Subject,GR_merge$tissue,sep='_')
for (fn in unique(GR_merge$fn)){
  
  fn_mml=paste(in_dir,fn,'_phased_allele_agnostic_mml.bedGraph',sep='')
  fn_nme=paste(in_dir,fn,'_phased_allele_agnostic_nme.bedGraph',sep='')
  
  if(file.exists(fn_mml)&file.exists(fn_nme)){
    cat('Processing',fn,'\n')
    GR_calc=rbind(GR_calc,allele_agnostic_merge(GR_merge[GR_merge$fn==fn],fn_nme,fn_mml))
    
  }
}
pdf('../downstream/output/dMML_effect_NME_all.pdf')
for (sp in unique(GR_calc$sample)){
  print(ggplot(GR_calc[GR_calc$sample==sp,],aes(x=NME,color=dMML_ASM))+geom_density(size=1)+labs(color='MML-Hap')+
          theme(legend.position = 'bottom')+ggtitle(paste(sp,sum(GR_calc$dMML_ASM[GR_calc$sample==sp]))))
}
dev.off()
pdf('../downstream/output/dMML_effect_NME_all_regions.pdf')
  print(ggplot(GR_calc,aes(x=NME,color=dMML_ASM))+geom_density(size=1)+labs(color='MML-Hap')+
          theme(legend.position = 'bottom')+ggtitle('NME near dMML regions'))
dev.off()


# Variant frequency -------------------------------------------------------

variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
variant_HetCpG_meta$variants=variants_collapase(variant_HetCpG_meta)
variant_freq=table(variant_HetCpG_meta$variants)
variant_freq=as.data.frame(variant_freq/sum(variant_freq))
colnames(variant_freq)=c('SNP','frequency')
#Plot SNP frequency
title='SNP frequency'
pdf('../downstream/output/SNP_frequency.pdf')
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
ggplot(variant_freq,aes(x=SNP,y=frequency,fill=SNP)) +ylim(c(0,0.5))+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+theme_bar
dev.off()

# Finding the overlap between monoallelic expressed gene ------------------
#Need to check *
GR_merge=readRDS(GR_merge_file)
MAE_BAE_data_Gimelbrant <- as.data.frame(read_excel("../downstream/input/MAE_BAE_data_Gimelbrant.xlsx"),stringsAsFactors=F)
#Enrichment in MAE dMML , not quiet enriched
MAE=MAE_BAE_data_Gimelbrant$Gene[ MAE_BAE_data_Gimelbrant$`MAE=1_BAE=0`==1]
#Find GR overlap the MAE:NME gene promoter
MAE_GR_dNME=GR_merge[GR_merge$dNME_pval<=pval_cutoff][GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff] %in% MAE]
#write(unique(MAE_GR_dNME$genes_promoter),"../downstream/output/dNME_sig_gene_body_MAE.txt",pval_cutoff=0.1)

#dNME enrichment promoter
MAE_enrich(GR_merge[GR_merge$genes_promoter%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_promoter',stat='dNME_pval',MAE=MAE)
#dMML enrichment promoter
MAE_enrich(GR_merge[GR_merge$genes_promoter%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_promoter',stat='dMML_pval',MAE=MAE)

#dNME enrichment body
MAE_enrich(GR_merge[GR_merge$genes_body%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_body',stat='dNME_pval',MAE=MAE)
#dMML enrichment body
MAE_enrich(GR_merge[GR_merge$genes_body%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_body',stat='dMML_pval',MAE=MAE)


# Processing bulk RNA data ------------------------------------------------
#HUES64:SRP000996 

###Reading in necessary variables
GR_merge=readRDS(GR_merge_file)
###Reading in RNA-seq data
#"../downstream/input/RNA_seq/"
fds_H1=c("H1_rep1_genome1","H1_rep2_genome1","H1_rep3_genome1","H1_rep1_genome2","H1_rep2_genome2","H1_rep3_genome2")
res_RNA_H1=RNA_seq_process("../downstream/input/RNA_seq/",fds_H1)
H1_df=RNA_df(GR_merge[GR_merge$Subject=='H1'],res_RNA_H1[!is.na(res_RNA_H1$log2FoldChange),])

fds_HUES64_ectoderm=c("HUES64_ectoderm_genome1","HUES64_ectoderm_genome2")
res_RNA_HUES64_ectoderm=RNA_seq_process("../downstream/input/RNA_seq/",fds_HUES64_ectoderm,condition_rep = 1)
HUES64_ectoderm_df=RNA_df(GR_merge[GR_merge$Sample=='ectoderm_paired - HUES64'],res_RNA_HUES64_ectoderm[!is.na(res_RNA_HUES64_ectoderm$log2FoldChange),])

fds_HUES64_endoerm=c("HUES64_endoderm_genome1","HUES64_endoderm_genome2")
res_RNA_HUES64_endoerm=RNA_seq_process("../downstream/input/RNA_seq/",fds_HUES64_endoerm,condition_rep = 1)
HUES64_endoerm_df=RNA_df(GR_merge[GR_merge$Sample=='endoerm_27_paired - HUES64'],res_RNA_HUES64_endoerm[!is.na(res_RNA_HUES64_endoerm$log2FoldChange),])

fds_HUES64_stem=c("HUES64_stem_genome1","HUES64_stem_genome2")
res_RNA_HUES64_stem=RNA_seq_process("../downstream/input/RNA_seq/",fds_HUES64_stem,condition_rep = 1)
HUES64_stem=RNA_df(GR_merge[GR_merge$Sample=='stem_27_undifferentiated_paired - HUES64'],res_RNA_HUES64_stem[!is.na(res_RNA_HUES64_stem$log2FoldChange),])

fds_HUES64_mesoderm=c("HUES64_mesoderm_genome1","HUES64_mesoderm_genome2")
res_RNA_HUES64_mesoderm=RNA_seq_process("../downstream/input/RNA_seq/",fds_HUES64_mesoderm,condition_rep = 1)
HUES64_mesoderm=RNA_df(GR_merge[GR_merge$Sample=='mesoderm_23_paired - HUES64'],res_RNA_HUES64_mesoderm[!is.na(res_RNA_HUES64_mesoderm$log2FoldChange),])
###At dMML hap, large correlation with RNA-seq, data too few not significant
RNA_seq_df=rbind(H1_df,HUES64_ectoderm_df,HUES64_endoerm_df,HUES64_stem,HUES64_mesoderm)
cor.test(RNA_seq_df$dMML_relative[RNA_seq_df$dMML_pval<=pval_cutoff],RNA_seq_df$ASE_log2FC[RNA_seq_df$dMML_pval<=pval_cutoff])
ggplot(RNA_seq_df[RNA_seq_df$dMML_pval<=pval_cutoff,],aes(x=dMML_relative,y=ASE_log2FC))+geom_smooth(method='lm')+geom_point()+
  xlab('relative dMML (genome2-genome1)')+ylab('log2Fc (genome2/genome1)')

cor.test(RNA_seq_df$dNME[RNA_seq_df$dNME_pval<=pval_cutoff],abs(RNA_seq_df$ASE_log2FC[RNA_seq_df$dNME_pval<=pval_cutoff]))
ggplot(RNA_seq_df[RNA_seq_df$dNME_pval<=pval_cutoff,],aes(x=dNME_relative,y=ASE_log2FC))+geom_smooth(method='lm')+geom_point()+
  xlab('relative dNME (genome2-genome1)')+ylab('log2Fc (genome2/genome1)')
# Imprinting enrichment ---------------------------------------------------
GR_merge=readRDS(GR_merge_file)

Imprinted_Genes <- as.data.frame(read_excel("../downstream/input/Imprinted Genes.xlsx"))
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',Imprinted_Genes$Gene)
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',
           Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Paternal'])
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',
           Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Maternal'])
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',Imprinted_Genes$Gene)
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',
           Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Paternal'])
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',
           Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Maternal'])
#ecdf
GR_merge_genes_promoter=GR_merge[!is.na(GR_merge$genes_promoter)]
GR_merge_genes_promoter$imprinted=FALSE
GR_merge_genes_promoter$imprinted[GR_merge$genes_promoter%in% Imprinted_Genes$Gene]=TRUE
GR_merge_genes_df=GR_merge_genes_df=data.frame(dMML=GR_merge_genes_promoter$dMML,imprinted=GR_merge_genes_promoter$imprinted)
ecdf_df_out=data.frame()
for(type in unique(GR_merge_genes_df$imprinted)){
  GR_merge_genes_promoter_estimate=ecdf(GR_merge_genes_promoter$dMML[GR_merge_genes_df$imprinted==type])
  imprinted_uq=seq(0,1,0.001)
  ecdf_df_out=rbind(ecdf_df_out,data.frame(dMML=imprinted_uq,quant=GR_merge_genes_promoter_estimate(imprinted_uq),imprinted=type,stringsAsFactors = F))
}
pdf('../downstream/output/graphs/imprinted_ecdf.pdf')
ggplot(ecdf_df_out,aes(x=dMML,y=quant,group=imprinted,color=imprinted))+geom_line(size=0.5)+xlab('dMML')+theme(legend.position = 'bottom')+
  ylab('cumulative probability')
dev.off()
# ChromHMM annotations ----------------------------------------------------


###Maybe having error in shiny package, install dev version

ah = AnnotationHub()
#d <- display(ah)
#18 state also put rep and zinc fingure together but have more detailed specification in TSS
#H1: AH46858, search 15 state chromatin, E003, 
#H9: AH46863
#HUES64: AH46871, E016
GR_merge=readRDS(GR_merge_file)
#Maybe looking at 
ENCODE_name=ENCODE_to_sample(unique(GR_merge$Sample))
# motif_gene=readRDS(motif_gene_file)
#GR_merge=GR_merge[GR_merge$N>1]
#checking if the motif is enriched in active region
# olap_motif=findOverlaps(GR_merge,motif_gene)
# GR_merge$motif=1
# GR_merge$motif[unique(queryHits(olap_motif))]=0
# NME_in=readRDS(NME_agnostic_file)
# MML_in=readRDS(MML_agnostic_file)
# MML_in$high_MML=1
# MML_in$high_MML[MML_in$quant_score=='75%-100%']=0
# MML_in$low_MML=1
# MML_in$low_MML[MML_in$quant_score=='0-25%']=0
# 
# NME_in$high_NME=1
# NME_in$high_NME[NME_in$quant_score=='75%-100%']=0
# NME_in$low_NME=1
# NME_in$low_NME[NME_in$quant_score=='0-25%']=0

#0.2
chromHMM_dNME_all_ls=list()
chromHMM_dMML_all_ls=list()
# chromHMM_motif_all_ls=list()
chromHMM_region_all=list()
# chromHMM_MML_high_all_ls=list()
# chromHMM_MML_low_all_ls=list()
# chromHMM_NME_high_all_ls=list()
# chromHMM_NME_low_all_ls=list()
ah_gr=GRanges()
#GR_merge_1=GR_merge[GR_merge$N==1]
#Do it for all available data
for (sp in ENCODE_name$sample[!is.na(ENCODE_name$ENCODE)]){
  ah_num=names(query(ah, c("chromhmmSegmentations", ENCODE_name$ENCODE[ENCODE_name$sample==sp])))
  chromHMM=ah[[ah_num]]
  #chromHMM=chromHMM_region_all[[sp]]
  chromHMM_dNME_all_ls[[sp]]=chromHMM_OR(GR_merge, chromHMM,sp)
  chromHMM_dMML_all_ls[[sp]]=chromHMM_OR(GR_merge, chromHMM,sp,stat="dMML_pval")
  #chromHMM_motif_all_ls[[sp]]=chromHMM_OR(GR_merge, chromHMM,sp,stat="motif")
  # if(!sp %in%c("rep1 - H1","rep2 - H1")){
  # chromHMM_MML_high_all_ls[[sp]]=chromHMM_OR(MML_in, chromHMM,sp,stat="high_MML")
  # chromHMM_MML_low_all_ls[[sp]]=chromHMM_OR(MML_in, chromHMM,sp,stat="low_MML")
  # chromHMM_NME_high_all_ls[[sp]]=chromHMM_OR(NME_in, chromHMM,sp,stat="high_NME")
  # chromHMM_NME_low_all_ls[[sp]]=chromHMM_OR(NME_in, chromHMM,sp,stat="low_NME")
 chromHMM_region_all[[sp]]=chromHMM
  # }
}


#####Summing up cont table
chromHMM_dNME_all=chromHMM_combine(chromHMM_dNME_all_ls)
chromHMM_dMML_all=chromHMM_combine(chromHMM_dMML_all_ls)
#chromHMM_motif_all=chromHMM_combine(chromHMM_motif_all_ls)
# chromHMM_MML_high_all=chromHMM_combine(chromHMM_MML_high_all_ls)
# chromHMM_MML_low_all=chromHMM_combine(chromHMM_MML_low_all_ls)
# chromHMM_NME_high_all=chromHMM_combine(chromHMM_NME_high_all_ls)
# chromHMM_NME_low_all=chromHMM_combine(chromHMM_NME_low_all_ls)

#####Subset regions according to the chromHMM states
# GR_merge_states=list()
# for(chrom_states in unique(chromHMM$name)){
# GR_merge_sp_states=GRanges()
# 
# chrom_states= chromHMM_dNME_all$states[chromHMM_dNME_all$OR<=0.5]
# for(sp in names(chromHMM_region_all)){
#   
#   GR_merge_sp_states=c(GR_merge_sp_states,subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$Sample==sp],
#                                                            chromHMM_region_all[[sp]][chromHMM_region_all[[sp]]$name%in%chrom_states]))
# }
# #GR_merge_states[[chrom_states]]=GR_merge_sp_states
# 
# #}

#Plot OR with error bar
pdf('../downstream/output/dNME_chromHMM.pdf')  
ggplot(chromHMM_dNME_all,aes(y=OR,x=states,fill=states))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+ coord_flip()+
  theme(axis.text.x = element_text(hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("dNME-ASM states annotation")+geom_text(aes(label=round(OR,digits = 2)), hjust=-3, color="black", size=3.5)
dev.off()
pdf('../downstream/output/dMML_chromHMM.pdf')
chromHMM_dMML_all=chromHMM_dMML_all[order(chromHMM_dMML_all$OR,decreasing=F),]
chromHMM_dMML_all$states=factor(chromHMM_dMML_all$states,levels=chromHMM_dMML_all$states)
ggplot(chromHMM_dMML_all,aes(x=states,y=OR,fill=states))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+ylim(0,11)+ coord_flip()+
  theme(axis.text.x = element_text(hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("dMML-ASM states annotation")+geom_text(aes(label=round(OR,digits = 2)), hjust=-3, color="black", size=3.5)
dev.off()

ggplot(chromHMM_motif_all,aes(x=states,y=OR,fill=states))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("motif states annotation")+geom_text(aes(label=round(OR,digits = 2)), vjust=3, color="black", size=3.5)

pdf('../downstream/output/high_MML_chromHMM.pdf')
ggplot(chromHMM_MML_high_all,aes(x=states,y=OR,fill=states))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("motif states annotation")+geom_text(aes(label=round(OR,digits = 2)), vjust=3, color="black", size=3.5)
dev.off()
pdf('../downstream/output/low_MML_chromHMM.pdf')
ggplot(chromHMM_MML_low_all,aes(x=states,y=OR,fill=states))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("motif states annotation")+geom_text(aes(label=round(OR,digits = 2)), vjust=3, color="black", size=3.5)
dev.off()


# Processing bulk ATAC-seq result --------------------------------------------

######reading in data
GR_merge=readRDS(GR_merge_file)
GR_merge_H1=GR_merge[GR_merge$Subject=='H1']
ATAC_H1=as.data.frame(read.table("../downstream/input/H1_ATAC_allele_count.txt",header=T,stringsAsFactors = F))
colnames(ATAC_H1)=c("seqname","start","end","rep1_g1","rep1_g2","rep2_g1","rep2_g2")
sampleTable <- data.frame(condition = factor(c("genome1","genome2","genome1","genome2")))
rownames(sampleTable) <- colnames(ATAC_H1)[4:7]
rownames(ATAC_H1)=paste(ATAC_H1$seqname,ATAC_H1$start,ATAC_H1$end,sep="_")
###Perfrom RNA-like analysis on ATAC-seq data
dds_ATAC = DESeqDataSetFromMatrix(countData=as.matrix(ATAC_H1[,4:7]), colData=sampleTable,design=~condition)
#Here, it uses piror distribution to aviod zero if betaPrior is set TRUE and use shrunkenLFC, does not seems change result
#MAP: maximum a posteriori, MLE: maximum likelihood##Learn more about this later
dds_ATAC_lfc<-DESeq(dds_ATAC,betaPrior=FALSE)
res_ATAC_lfc<-results(dds_ATAC_lfc,name=c("condition_genome2_vs_genome1"))
###lfc works better
##############read in ATAC-seq CPM data they are same in regions
ATAC_H1_CPM=as.data.frame(read.table("../downstream/input/H1_ATAC_allele_CPM.txt",header=T,stringsAsFactors = F))
colnames(ATAC_H1_CPM)=c("seqname","start","end","rep1_g1_CPM","rep1_g2_CPM","rep2_g1_CPM","rep2_g2_CPM")
ATAC_H1_CPM=makeGRangesFromDataFrame(as.data.frame(ATAC_H1_CPM),keep.extra.columns = T)
elementMetadata(ATAC_H1_CPM)=cbind(res_ATAC_lfc,elementMetadata(ATAC_H1_CPM))
#############calcualte log2FC using CPM and pseudocount
ATAC_H1_CPM$log2FC_CPM=log2((ATAC_H1_CPM$rep1_g2_CPM+ATAC_H1_CPM$rep2_g2_CPM+1)/(ATAC_H1_CPM$rep1_g1_CPM+ATAC_H1_CPM$rep2_g1_CPM+1))

ATAC_H1_CPM$genome2_CPM=(ATAC_H1_CPM$rep1_g2_CPM+ATAC_H1_CPM$rep2_g2_CPM)/2
ATAC_H1_CPM$genome1_CPM=(ATAC_H1_CPM$rep1_g1_CPM+ATAC_H1_CPM$rep2_g1_CPM)/2
############Putting data into GR merge object
olap_ATAC=findOverlaps(GR_merge_H1,ATAC_H1_CPM,select="first",maxgap =500)
#####Check overlapped regions
sum(!is.na(olap_ATAC))#24783
GR_merge_H1_ATAC=GR_merge_H1#[queryHits(olap_ATAC)]
GR_merge_H1_ATAC$ATAC_FC=NA
GR_merge_H1_ATAC$ATAC_FC=ATAC_H1_CPM$log2FoldChange[olap_ATAC]
GR_merge_H1_ATAC$ATAC_FC_pval=ATAC_H1_CPM$pvalue[olap_ATAC]
GR_merge_H1_ATAC$ATAC_FC_pval_adj=p.adjust(GR_merge_H1_ATAC$ATAC_FC_pval,method="BH")
GR_merge_H1_ATAC$ATAC_FC_CPM=ATAC_H1_CPM$log2FC_CPM[olap_ATAC]
GR_merge_H1_ATAC$g1_CPM=ATAC_H1_CPM$genome1_CPM[olap_ATAC]
GR_merge_H1_ATAC$g2_CPM=ATAC_H1_CPM$genome2_CPM[olap_ATAC]#126027 
GR_merge_H1_ATAC$log2FC_CPM=ATAC_H1_CPM$log2FC_CPM[olap_ATAC]
GR_merge_H1_ATAC=GR_merge_H1_ATAC[!is.na(GR_merge_H1_ATAC$ATAC_FC)]#17692
GR_merge_H1_ATAC$dMML_relative=GR_merge_H1_ATAC$MML2-GR_merge_H1_ATAC$MML1
GR_merge_H1_ATAC$dNME_relative=GR_merge_H1_ATAC$NME2-GR_merge_H1_ATAC$NME1
####test correlation
cor.test(GR_merge_H1_ATAC$dMML_relative[GR_merge_H1_ATAC$dMML_pval<=pval_cutoff],
         GR_merge_H1_ATAC$log2FC_CPM[GR_merge_H1_ATAC$dMML_pval<=pval_cutoff])#-0.6302
cor(GR_merge_H1_ATAC$dMML_relative,GR_merge_H1_ATAC$log2FC_CPM)#-0.6302
dMML_plot=data.frame(dMML=GR_merge_H1_ATAC$dMML_relative,ATAC_FC_CPM=GR_merge_H1_ATAC$log2FC_CPM,dMML_pval=GR_merge_H1_ATAC$dMML_pval,
                     dNME=GR_merge_H1_ATAC$dNME_relative,dNME_pval=GR_merge_H1_ATAC$dNME_pval)
pdf('../downstream/output/dMML/ATAC-seq_change.pdf')
ggplot(dMML_plot[dMML_plot$dMML_pval<=pval_cutoff,],aes(x=dMML,y=ATAC_FC_CPM))+geom_smooth()+geom_point(alpha=0.5)+
  ylab('Allelic accessibility change')
dev.off()
cor.test(GR_merge_H1_ATAC$dNME_relative[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff],
         GR_merge_H1_ATAC$log2FC_CPM[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff])
ggplot(dMML_plot[dMML_plot$dNME_pval<=pval_cutoff,],aes(x=dNME,y=ATAC_FC_CPM))+geom_smooth()+geom_point(alpha=0.5)+
  ylab('Allelic accessibility change')


# Hypervariance vs allele-agnostic NME ------------------------------------
GR_merge=readRDS(GR_merge_file)
genomic_features=readRDS(genomic_features_file)
#GR_merge_H1=GR_merge[GR_merge$Sample=='merged - H1']
#GR_merge_H1$score=GR_merge_H1$dNME
#H1_hyper_var=readRDS('../downstream/output/H1_hypervar.rds')

###Ploting hypervar vs distance to TSS
in_dir='../downstream/data/allele_agnostic_2kb/'
NME_hypervar_calc=GRanges()
MML_hypervar_calc=GRanges()
NME_meanvar_calc=GRanges()
MML_meanvar_calc=GRanges()
GR_calc=data.frame()
GR_merge$fn=paste(GR_merge$Subject,GR_merge$tissue,sep='_')
scRNA_result=data.frame()
#For the mean expression: log2 or not?

for (fn in unique(GR_merge$fn)){
  
  fn_mml=paste(in_dir,fn,'_phased_allele_agnostic_mml.bedGraph',sep='')
  fn_nme=paste(in_dir,fn,'_phased_allele_agnostic_nme.bedGraph',sep='')
  hyper_var_file=unique(unlist(GR_merge$hyper_var_fn[GR_merge$fn==fn]))
  cat('Processing',fn,'\n')
  if(file.exists(fn_mml)&file.exists(fn_nme)&length(hyper_var_file)>0){
    sp_hyper_var=read_hypervar(hyper_var_file)
    #scRNA_result=rbind(scRNA_result,sp_hyper_var)
    sp_hyper_var$log2mean=log2(sp_hyper_var$mean)
    #Add hypervaribility inforamtion
    NME_hypervar_calc=c(NME_calc,dist_plot_calc(fn_nme,sp_hyper_var,genomic_features,GR_merge[GR_merge$fn==fn]))
    MML_hypervar_calc=c(NME_calc,dist_plot_calc(fn_nme,sp_hyper_var,genomic_features,GR_merge[GR_merge$fn==fn],stat='MML'))
    MML_meanvar_calc=c(MML_calc,dist_plot_calc(fn_mml,sp_hyper_var,genomic_features,GR_merge[GR_merge$fn==fn],var_stat=quote(log2mean),stat='MML'))
    NME_meanvar_calc=c(MML_calc,dist_plot_calc(fn_mml,sp_hyper_var,genomic_features,GR_merge[GR_merge$fn==fn],var_stat=quote(log2mean),stat='NME'))
  }
}

cor_out=data.frame()
NME_calc_agg_out=data.frame()
#heatmap for hypervaribility vs NME,500bp
for(sp in unique(NME_calc$Sample)){
  NME_calc_sp=NME_calc[NME_calc$Sample==sp&abs(NME_calc$dist)<=500]
  
  NME_calc_agg=aggregate(NME_calc_sp$score,by=list(round(NME_calc_sp$hypervarquant0001,digits =2)),FUN=median)
  colnames(NME_calc_agg)=c('hyper_var','NME')
  NME_calc_agg$Sample=sp
  NME_calc_agg$cor=cor.test(NME_calc_sp$score,NME_calc_sp$hyper_var)$estimate
  NME_calc_agg_out=rbind(NME_calc_agg_out,NME_calc_agg)
}
NME_calc_agg_out=NME_calc_agg_out[order(NME_calc_agg_out$cor,decreasing = F),]
NME_calc_agg_out$Sample=factor(NME_calc_agg_out$Sample,levels = unique(NME_calc_agg_out$Sample))
pdf('../downstream/output/NME_heatmap_mean_hypervar.pdf')
ggplot(NME_calc_agg_out,aes(hyper_var,Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
  xlab('Hypervaribility quantile')+ylab('Sample')+theme(legend.position = 'bottom')
dev.off()
cor.test(NME_calc$score[abs(NME_calc$dist)<=500],NME_calc$hyper_var[abs(NME_calc$dist)<=500])
pdf('../downstream/output/correlation_quant_raw_sample_all_NME_hyper_var_0.pdf')
for(sp in unique(NME_calc$Sample)){
  NME_calc_sp=NME_calc[NME_calc$dist<=300&NME_calc$dist>=0&NME_calc$Sample==sp&NME_calc$hyper_var>0]
  NME_plot_sp_df=data.frame(hyper_varibility=NME_calc_sp$hyper_var,NME=NME_calc_sp$score)
  # print(ggplot(NME_plot_sp_df,aes(x=factor(round(hyper_varibility,digits=1)),y=NME))+ ggtitle(sp)+
  #  geom_violin(trim=FALSE)+stat_summary(fun=median, geom="point", size=2, color="red")+xlab('rounded hypervaribility'))
  print(ggplot(NME_plot_sp_df,aes(x=factor(round(hyper_varibility,digits=0)),y=NME))+ ggtitle(sp)+
          geom_boxplot(outlier.shape = NA)+stat_summary(fun=mean, geom="point", size=2, color="red")+xlab('rounded hypervaribility'))
  print(ggplot(NME_plot_sp_df,aes(x=hyper_varibility,y=NME))+ ggtitle(sp)+
          geom_point(alpha=0.1)+xlab('rounded hypervaribility')+ylab('NME'))
  cor_out=rbind(cor_out,data.frame(sample=sp,cor=cor.test(NME_plot_sp_df$NME,NME_plot_sp_df$hyper_varibility)$estimate[[1]]))
}
dev.off()

hypo_var_genes=NME_calc$gene[abs(NME_calc$dist)<=500&NME_calc$quant_score=='0-25%'&NME_calc$quant=='0-25%']
hyper_var_genes=NME_calc$gene[abs(NME_calc$dist)<=500&NME_calc$hypervarquant0001>=0.75&NME_calc$scorequant0001>=0.75]
background_genes=NME_calc$gene[abs(NME_calc$dist)<=500]
GO_hyper_var=GO_anno(hyper_var_genes,background_genes,topNodes = 11618)
GO_hyper_var_sig=GO_hyper_var[[1]][GO_hyper_var[[1]]$FC>=2&GO_hyper_var[[1]]$Significant>=10,]
GO_hyper_var_sig$qval=p.adjust(GO_hyper_var_sig$classicFisher)
write.csv(GO_hyper_var_sig[,c(1,2,7,8)],'../downstream/output/high_NME_GO.csv')
#cat to Overleaf format
GO_hyper_var_sig$FC=round(GO_hyper_var_sig$FC,digits = 2)

GO_hyper_var_sig=GO_hyper_var_sig[order(GO_hyper_var_sig$qval),]
GO_hyper_var_sig$qval=format(GO_hyper_var_sig$qval,digits = 2)
for(i in 1:20){cat(paste(GO_hyper_var_sig[i,c(1,2,7,8)],collapse = ' & '),'\\\\','\n')}
# background_genes=unique(GR_merge$genes_promoter)
# hyper_var_genes=unique(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff])
# GO_hyper_var=GO_anno(hyper_var_genes,background_genes)
# GO_hyper_var_sig=GO_hyper_var[[1]][GO_hyper_var[[1]]$FC>=2&GO_hyper_var[[1]]$Significant>=10,]
# GO_hyper_var_sig$qval=p.adjust(GO_hyper_var_sig$classicFisher)
# write.csv(GO_hyper_var_sig[,c(1,2,7,8)],'../downstream/output/high_dNME_GO.csv')


# #Violin plot
# 
# ggplot(data.frame(hyper_varibility=factor(round(NME_calc$hyper_var[abs(NME_calc$dist)<=500])),NME=NME_calc$score[abs(NME_calc$dist)<=500]),aes(x=hyper_varibility,y=NME))+ 
#   geom_violin(trim=FALSE)+stat_summary(fun=median, geom="point", size=2, color="red")+xlab('rounded hypervaribility')
# 
# pdf('../downstream/output/graphs/Figure2/NME_correlation_agnostic2.pdf')
# dist_plot_run(NME_calc)
# dev.off()
# pdf('../downstream/output/graphs/Figure2/MML_correlation_agnostic2.pdf')
# dist_plot_run(MML_calc)
# dev.off()



# Gene ontology -----------------------------------------------------------
#high NME genes in dNME
GR_merge=readRDS(GR_merge_file)
GR_merge$highNME=GR_merge$NME1>=NME_quant|GR_merge$NME2>=NME_quant
dNME_GO=GO_anno(unique(GR_merge$genes_promoter[GR_merge$highNME]),unique(GR_merge$genes_promoter))
dNME_GO=GO_anno(unique(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff]),unique(GR_merge$genes_promoter))
dNME_GO_out=dNME_GO[[1]]

dNME_GO[[1]][dNME_GO[[1]]$Significant>=10,c(2,7)]
dMML_GO=GO_anno(unique(GR_merge$genes_promoter[GR_merge$dMML_pval<=pval_cutoff]),unique(GR_merge$genes_promoter))
dMML_GO[[1]][,c(2,7)]
#high dNME genes
write(unique(unlist(GR_merge$genes_promoter[GR_merge$dNME<=pval_cutoff])),'../downstream/output/dNME_promoter_sig.txt')
write(unique(unlist(GR_merge$genes_promoter)),'../downstream/output/dNME_promoter_background.txt')
#high NME and hypervaribile genes
NME_in=readRDS(NME_agnostic_file)
genomic_features=readRDS(genomic_features_file)
NME_in=dist_calc(NME_in,genomic_features$TSS)
GO=GO_anno(NME_in$gene[NME_in$quant_score=="75%-100%"&abs(NME_in$dist)<=500],NME_in$gene[abs(NME_in$dist)<=500])
GO_out=GO[[1]]
GREAT_output=list()
for (sp in unique(NME_in$Sample)){
  NME_sp=NME_in[NME_in$quant_score=="75%-100%"&abs(NME_in$dist)<=500&NME_in$Sample==sp]
  NME_sp=NME_sp[order(NME_sp$NME,decreasing=T)][1:5000]
  GREAT_output[[sp]]= submitGreatJob(NME_sp, species = 'hg19')
  # export.bed(NME_sp,paste('../downstream/output/high_NME_sp/high_NME_TSS_500_sig_',sp,'.bed',sep=''))
  # NME_in_sp=NME_in[NME_in$Sample==sp&abs(NME_in$dist)<=500]
  # NME_in_sp_df=mcols(NME_in_sp)
  # NME_in_sp_df_agg=aggregate(NME_in_sp_df$NME,by=list(NME_in_sp_df$gene),FUN=mean)
  # colnames(NME_in_sp_df_agg)=c('gene','mean NME')
  # write.table(matrix(unlist(NME_in_sp_df_agg$gene[order(NME_in_sp_df_agg$`mean NME`,decreasing = T)]),nrow=1),quote=FALSE,
  #             paste('../downstream/output/high_NME_sp/high_NME_TSS_500_sig_',sp,'.txt',sep=''),sep="\n",row.names = FALSE,col.names = F)
}
GREAT_output_tb=lapply(GREAT_output,getEnrichmentTables)
GREAT_output_tb_BP=lapply(GREAT_output_tb,function(x) x$`GO Biological Process`[x$`GO Biological Process`$Binom_Adjp_BH<=0.1&
                                                                                  x$`GO Biological Process`$Binom_Observed_Region_Hits>=10,])
BP_sp_out_df=data.frame()
for(sp in names(GREAT_output_tb)){
  BP_sp=GREAT_output_tb_BP[[sp]]
  BP_sp$Sample=sp
  BP_sp_out_df=rbind(BP_sp_out_df,BP_sp)
}
BP_sp_out_df$neg_log10_BH=-log10(BP_sp_out_df$Binom_Adjp_BH)
BP_sp_out_df_cast=dcast(BP_sp_out_df,name ~Sample,value.var="neg_log10_BH")
BP_sp_out_df_cast=BP_sp_out_df_cast[order(rowMeans(BP_sp_out_df_cast[,-1],na.rm = T),decreasing=T),]
BP_sp_out_df_cast_unique=BP_sp_out_df_cast[which(apply(BP_sp_out_df_cast,1,function(x) sum(is.na(x)))>=30),]
BP_sp_out_df_sub=BP_sp_out_df[BP_sp_out_df$name%in%BP_sp_out_df_cast$name[1:30],]
BP_sp_out_df_sub$name=factor(BP_sp_out_df_sub$name,levels=rev(BP_sp_out_df_cast$name[1:30]))
ggplot(BP_sp_out_df_sub,aes(y=name,x=Sample,fill=neg_log10_BH))+geom_tile()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
write(unique(NME_in$gene[NME_in$quant_score=="75%-100%"&abs(NME_in$dist)<=500]),'../downstream/output/NME_TSS_500_sig.txt')
write(unique(NME_in$gene[abs(NME_in$dist)<=500]),'../downstream/output/NME_TSS_500_background.txt')


# NME distribution at different density -----------------------------------

#Using GR_merge to do it
GR_merge=readRDS(GR_merge_file)
GR_merge$CpGdiff=GR_merge$g1CG-GR_merge$g2CG
#Looking for Het CpG effect over ASM different regions 
NME_allele_ASM_calc=allele_calc_plot(GR_merge[GR_merge$dNME_pval<=pval_cutoff],'NME','../downstream/output/dNME/',
                                     'NME distribution at NME-Hap')
MML_allele_ASM_calc=allele_calc_plot(GR_merge[GR_merge$dMML_pval<=pval_cutoff],'MML','../downstream/output/dMML/',
                                     'MML distribution at MML-Hap')
NME_allele_ASM_calc=allele_diff_merge(GR_merge[GR_merge$dNME_pval<=pval_cutoff])
NME_allele_all_calc=allele_diff_merge(GR_merge)
saveRDS(NME_allele_ASM_calc,'../downstream/output/dNME/NME_ASM_het_calc_run3.rds')
saveRDS(NME_allele_all_calc,'../downstream/output/dNME/NME_all_het_calc_run3.rds')

# Local density vs dNME ---------------------------------------------------------
GR_merge=readRDS(GR_merge_file)
GR_merge$density=GR_merge$CG_hg19_extend/GR_merge$gff_size_extend
#GR_merge$density=GR_merge$CG_hg19_extend/GR_merge$CGcont_exp
GR_merge$CG_diff=GR_merge$g1CG-GR_merge$g2CG
GR_merge$dNME=GR_merge$NME1-GR_merge$NME2
GR_merge$dMML=GR_merge$MML1-GR_merge$MML2
density_df=data.frame(density=GR_merge$density[GR_merge$CG_diff!=0],dNME=GR_merge$dNME[GR_merge$CG_diff!=0],
                      density_diff=GR_merge$CG_diff[GR_merge$CG_diff!=0]/GR_merge$gff_size_extend[GR_merge$CG_diff!=0],
                      dMML=GR_merge$dMML[GR_merge$CG_diff!=0],
                      dNME_pval=GR_merge$dNME_pval[GR_merge$CG_diff!=0],dMML_pval=GR_merge$dMML_pval[GR_merge$CG_diff!=0])
#This is currently in use
# cor.test(density_df$dNME[abs(density_df$density_diff)<=0.002&density_df$dNME_pval<=pval_cutoff],
#          density_df$density_diff[abs(density_df$density_diff)<=0.002&density_df$dNME_pval<=pval_cutoff])#-0.426
# cor.test(density_df$dMML[abs(density_df$density_diff)<=0.002&density_df$dMML_pval<=pval_cutoff],
#          density_df$density_diff[abs(density_df$density_diff)<=0.002&density_df$dMML_pval<=pval_cutoff])#-0.426

cor.test(density_df$dNME[density_df$dNME_pval<=pval_cutoff],
         density_df$density_diff[density_df$dNME_pval<=pval_cutoff])
cor.test(density_df$dMML[density_df$dMML_pval<=pval_cutoff],
         density_df$density_diff[density_df$dMML_pval<=pval_cutoff])
# 
# density_df_agg=aggregate(abs(density_df$dNME),by=list(round(density_df$density,digits=4)),FUN=median)
# colnames(density_df_agg)=c('density','dNME')
# cor.test(density_df_agg$dNME,density_df_agg$density)
# cor.test(density_df$dNME,log2(density_df$density),method="spearman")#-0.002
# 
# #Violin plot
# #bin based on log2 density
# density_quant=quantile(log2(density_df$density_diff),prob=(seq(0,1,0.1)))
# density_df$density_quant=findInterval(log2(density_df$density),density_quant)
# quant=paste(seq(0,100,10),'%',sep='')
# quant=quant[1:10]
# density_df$density_quant=as.factor(quant[density_df$density_quant])
# ggplot(density_df,aes(x=density_quant,y=abs(dNME)))+ geom_violin(trim=FALSE)+stat_summary(fun.y=median, geom="point", size=2, color="red")
# density_df$density_round=as.factor(quant[density_df$density_quant])
# #dMML have more significant negative correlation than dNME
# pdf('../downstream/output/dNME/local_density_dNME.pdf')
# ggplot(density_df_agg,aes(x=density, y=abs(dNME)))+
#   ylim(c(0,0.15))+ggtitle("All regions")+geom_point(alpha=0.1)
#   ylab("dNME")+xlab("local CpG density")
# dev.off()
# jpeg('../downstream/output/dNME/local_density_dNME.jpg')
# ggplot(density_df[density_df$dNME_pval<=pval_cutoff,],aes(x=log2(density), y=abs(dNME)))+
#   ggtitle("All regions")+geom_smooth()+ylim(0,1)+geom_point(alpha=0.01)+
#   ylab("dNME")+xlab("log2(local density)")
# dev.off()
# pdf('../downstream/output/dNME/local_density_diff_dNME.pdf')
# ggplot(density_df[abs(density_df$density_diff)<=0.002&density_df$dNME_pval<=pval_cutoff,],aes(x=density_diff, y=dNME))+
#   ggtitle("All regions")+geom_smooth()+ylim(-1,1)+geom_point(alpha=0.1)+
#   ylab("dNME")+xlab("local CpG density difference")
# dev.off()


# allele-agnostic analysis on density ------------------------------------------------
NME_in=readRDS(NME_agnostic_file)
NME_in=NME_in[!is.na(NME_in$score)]
NME_in=NME_in[!seqnames(NME_in) %in%c('chrX','chrY','chrMT')]

#genomic_features=readRDS(genomic_features_file)
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
#NME vs different density
#NME_in_uniq=unique(NME_in)
#gr_seq=getSeq(Hsapiens,NME_in_uniq,as.character=T)
#NME_in_uniq$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
#saveRDS(NME_in_uniq,'../downstream/input/CG_Exp_NME_in.rds')
#NME_in_uniq=readRDS('../downstream/input/CG_Exp_NME_in.rds')
NME_in$CG_hg19=countOverlaps(NME_in,CpG_hg19)
#olap=findOverlaps(NME_in,NME_in_uniq,type='equal')
#NME_in$CG_exp[queryHits(olap)]=NME_in_uniq$CGcont_exp[subjectHits(olap)]
NME_in$density=NME_in$CG_hg19/width(NME_in)
cor.test(NME_in$density,NME_in$NME,method='spearman')#or pearson,0.15
#NME_in$density=NME_in$CG_hg19/NME_in$CG_exp
MML_in=readRDS(MML_agnostic_file)
MML_in$CG_hg19=countOverlaps(MML_in,CpG_hg19)
MML_in$density=MML_in$CG_hg19/width(MML_in)
cor.test(MML_in$density,MML_in$MML,method='spearman')#or pearson,-0.08551124 

density_quant=quantile(NME_in$density,prob=(seq(0,1,0.1)))
NME_in$density_quant=findInterval(NME_in$density,density_quant)
NME_in$density_quant[NME_in$density_quant==11]=10
quant=paste(seq(0,100,10),'%',sep='')
quant=paste(quant[1:10],paste(seq(10,100,10),'%',sep=''),sep='-')
NME_in$density_quant=as.factor(quant[NME_in$density_quant])
# ggplot(data.frame(density_quant=NME_in$density_quant,NME=NME_in$NME),aes(x=density_quant,y=NME))+ 
#   geom_violin(trim=FALSE)+stat_summary(fun.y=median, geom="point", size=2, color="red")+xlab('quantile of (density)')
qt_025<-function(x){return(quantile(x,probs=0.25))}
qt_075<-function(x){return(quantile(x,probs=0.75))}
qt_050<-function(x){return(quantile(x,probs=0.5))}
NME_in_df=data.frame(NME=NME_in$NME,density=NME_in$density)
NME_in_df_agg=aggregate(NME_in_df$NME,by=list(round(NME_in_df$density,digits=2)),FUN=qt_050)
NME_in_df_agg$low_CI=aggregate(NME_in_df$NME,by=list(round(NME_in_df$density,digits=2)),FUN=qt_025)$x
NME_in_df_agg$high_CI=aggregate(NME_in_df$NME,by=list(round(NME_in_df$density,digits=2)),FUN=qt_075)$x
names(NME_in_df_agg)=c("density","median","Bottom25","top25")
NME_in_df_agg$Bottom25= predict(loess(Bottom25~density,NME_in_df_agg),newdata=NME_in_df_agg$density)
NME_in_df_agg$top25= predict(loess(top25~density,NME_in_df_agg),newdata=NME_in_df_agg$density)
#GR_NME_df_agg=melt(GR_NME_df_agg,id="MML")
pdf('../downstream/output/graphs/Figure1/NME_vs_density_with_quantile.pdf')
###Plotting
ggplot(NME_in_df_agg,aes(x=density, y=median))+
   ylim(c(0,1))+ggtitle("NME density relationship")+geom_smooth(method="loess",se=FALSE)+
  ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)
dev.off()

MML_in_df=data.frame(NME=MML_in$MML,density=log10(MML_in$density))
MML_in_df_agg=aggregate(MML_in_df$NME,by=list(round(MML_in_df$density,digits=2)),FUN=qt_050)
MML_in_df_agg$low_CI=aggregate(MML_in_df$NME,by=list(round(MML_in_df$density,digits=2)),FUN=qt_025)$x
MML_in_df_agg$high_CI=aggregate(MML_in_df$NME,by=list(round(MML_in_df$density,digits=2)),FUN=qt_075)$x
names(MML_in_df_agg)=c("density","median","Bottom25","top25")
MML_in_df_agg$Bottom25= predict(loess(Bottom25~density,MML_in_df_agg),newdata=MML_in_df_agg$density)
MML_in_df_agg$top25= predict(loess(top25~density,MML_in_df_agg),newdata=MML_in_df_agg$density)
#GR_NME_df_agg=melt(GR_NME_df_agg,id="MML")
pdf('../downstream/output/graphs/MML_vs_density_with_quantile.pdf')
###Plotting
ggplot(MML_in_df_agg,aes(x=density, y=median))+
  ylim(c(0,1))+ggtitle("MML density relationship")+geom_smooth(method="loess",se=FALSE)+
  ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)
dev.off()

pdf('../downstream/output/NME_vs_density_raw.pdf')
ggplot(data.frame(density_round=as.factor(NME_in$density_round),NME=NME_in$NME),aes(x=density_round,y=NME))+ 
  geom_violin(trim=FALSE)+stat_summary(fun=median, geom="point", size=2, color="red")+xlab('density')
dev.off()
#CG ratio is not working
jpeg('../downstream/output/graphs/Figure2/NME_density_agnostic_density_pt.jpg')
ggplot(data.frame(NME=NME_in$NME,density=NME_in$density),aes(x=density,y=NME))+geom_smooth(size=1)+
  xlab('density')+ylab('NME')+geom_point()
dev.off()
pdf('../downstream/output/graphs/Figure2/NME_correlation_agnostic2.pdf')
dist_plot_run(NME_calc)
dev.off()
# #attemp bootstrap test
# for (i in 1:10000){
#   idx=sample(1:length(NME_in),1000000,replace = F)
#   cor_NME=c(cor_NME,cor(NME_in$score[idx],log2(NME_in$density)[idx]))
# }
ggplot(data.frame(NME=NME_in$score,density=NME_in$density),aes(y=NME,x=density))+geom_smooth()
###NME at different genomic features: promoters
olap_islands=findOverlaps(NME_in,genomic_features$`CpG island`)
olap_shores=findOverlaps(NME_in,genomic_features$`CpG shore`)
olap_shelf=findOverlaps(NME_in,genomic_features$`CpG shelf`)
olap_open_sea=findOverlaps(NME_in,genomic_features$`CpG open sea`)

CpG_density_NME=rbind(data.frame(NME=NME_in$NME[queryHits(olap_islands)],feature='islands'),
                      data.frame(NME=NME_in$NME[queryHits(olap_shores)],feature='shores'),
                      data.frame(NME=NME_in$NME[queryHits(olap_shelf)],feature='shelf'),
                      data.frame(NME=NME_in$NME[queryHits(olap_open_sea)],feature='open sea'))
aggregate(CpG_density_NME$NME,by=list(CpG_density_NME$feature),FUN=median)
pdf('../downstream/output/NME_density_feature.pdf')
my_comparisons <- list( c("islands", "shores"), c("islands", "shelf"), c("islands", "open sea"),
                        c("shores", "shelf"),c("shores", "open sea"),c("shelf", "open sea"))
ggplot(CpG_density_NME,aes(x=feature,y=NME))+geom_boxplot(outlier.shape = NA)+ stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()
###Different 
###NME at different genomic features: promoters
olap_promoter=findOverlaps(NME_in,genomic_features$`promoter`)
olap_exon=findOverlaps(NME_in,genomic_features$`exon`)
olap_intron=findOverlaps(NME_in,genomic_features$`intron`)
olap_intergenic=findOverlaps(NME_in,genomic_features$`intergenic`)

gene_NME=rbind(data.frame(NME=NME_in$NME[queryHits(olap_promoter)],feature='promoter'),
               data.frame(NME=NME_in$NME[queryHits(olap_exon)],feature='exon'),
               data.frame(NME=NME_in$NME[queryHits(olap_intron)],feature='intron'),
               data.frame(NME=NME_in$NME[queryHits(olap_intergenic)],feature='intergenic'))
aggregate(gene_NME$NME,by=list(gene_NME$feature),FUN=median)
pdf('../downstream/output/NME_density_gene.pdf')
my_comparisons <- list( c("promoter", "intergenic"), c("exon", "intergenic"), c("intron", "intergenic"),
                        c("promoter", "exon"),c("promoter", "intron"),c("exon", "intron"))
ggplot(gene_NME,aes(x=feature,y=NME))+geom_boxplot(outlier.shape = NA)+xlab('gene feature')+ stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()

# Enrichment of high NME region in high hypervaribility region ------------
NME_enrich_cont_table=list()
MML_enrich_cont_table=list()
for(tissue in unique(GR_merge$tissue)){
  print(tissue)
  NME_enrich_cont_table[[tissue]]=stat_hyper_enrich(GR_merge[GR_merge$tissue==tissue])
  MML_enrich_cont_table[[tissue]]=stat_hyper_enrich(GR_merge[GR_merge$tissue==tissue],stat='MML')
}
CMH_test(do.call(rbind,NME_enrich_cont_table))
CMH_test(do.call(rbind,MML_enrich_cont_table))

# Motif_enrichment analysis -----------------------------------------------

# ###Loading DNAse information
#https://egg2.wustl.edu/roadmap/data/byDataType/dnase/
# load("../downstream/input/state_calls_prom.RData")
# prom_DNase_gr=do.call('c',lapply(rownames(max_states),make_gr_dnase))
# elementMetadata(prom_DNase_gr)=max_states
# saveRDS(prom_DNase_gr,"../downstream/input/prom_DNase.rds")
# load("../downstream/input/state_calls_enh.RData")
# enchancer_DNase_gr=do.call('c',lapply(rownames(max_states),make_gr_dnase))
# elementMetadata(enchancer_DNase_gr)=max_states
# saveRDS(enchancer_DNase_gr,"../downstream/input/enchancer_DNase.rds")

enchancer_DNase_gr=readRDS("../downstream/input/enchancer_DNase.rds")
prom_DNase_gr=readRDS("../downstream/input/prom_DNase.rds")
DNase_all=c(prom_DNase_gr,enchancer_DNase_gr)
#export.bed(DNase_all,'../downstream/output/DNAseI_all.bed')
###Reading in data
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene <- readRDS(motif_gene_file)
motif_gene=subsetByOverlaps(motif_gene,variant_HetCpG_meta,type='equal')
olap=findOverlaps(variant_HetCpG_meta,motif_gene)
# variant_HetCpG_meta=variant_HetCpG_meta[queryHits(olap)]
# variant_HetCpG_meta$alleleDiff=motif_gene$alleleDiff[subjectHits(olap)]
# variant_HetCpG_meta$NME_diff=variant_HetCpG_meta$altNME-variant_HetCpG_meta$refNME
# variant_HetCpG_meta$MML_diff=variant_HetCpG_meta$altMML-variant_HetCpG_meta$refMML
# # #If need to only analyze DNase data
# variant_HetCpG_meta=subsetByOverlaps(variant_HetCpG_meta,DNase_all)
# motif_gene=subsetByOverlaps(motif_gene,DNase_all)
###Running enrichment analysis, currently only use direction_calc_enriched_subj
motif_sig_df=motif_enrich(motif_gene,variant_HetCpG_meta,pval_cutoff =pval_cutoff,dist=500)
# subset allelediff >1.5 for log
motif_dir=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1)
motif_gene_sig=motif_gene[motif_gene$geneSymbol%in% motif_dir_sig_ent$TF]
variant_HetCpG_meta$motif=FALSE
olap=findOverlaps(variant_HetCpG_meta,motif_gene_sig)
variant_HetCpG_meta$motif[queryHits(olap)]=TRUE
GO_motif_dNME=GO_anno(unique(unlist(variant_HetCpG_meta$genes_promoter[variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$motif])),
                      unique(unlist(variant_HetCpG_meta$genes_promoter)))
GO_motif_dNME_out=GO_motif_dNME[[1]]

##getting Ken's motif


motif_dir$qval_binom=p.adjust(motif_dir$binom.pval,method='BH')
saveRDS(motif_dir,'../downstream/output/motif_dirction_all_JASPAR_default.rds')
###Running directionality analysis
#motif_dir_sig=motif_dir[motif_dir$qval_OR<=0.1,]
#motif_dir_sig=readRDS('../downstream/output/motif_dirction_all.rds')
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
#motif_dir=readRDS('../downstream/output/motif_dir.rds')
motif_dir$qvalue=p.adjust(motif_dir$binom.pval,method='BH')
motif_dir_sig_ent=motif_dir[motif_dir$qval_binom<=0.1&motif_dir$prob>0.5,]
motif_dir_sig_ent_select=motif_dir[motif_dir$qval_binom<=0.1&motif_dir$prob>=0.6,]
motif_dir_sig_non_ent=motif_dir[motif_dir$qval_binom<=0.1&motif_dir$prob<0.5,]
OMIM=read.csv('../downstream/input/genemap2.txt',fill = T,sep='\t',header=T,stringsAsFactors =F,na.strings = "NA",skip=3)
OMIM$Gene.Symbols=strsplit(as.character(OMIM$Gene.Symbols),', ')
motif_dir_sig_ent=motif_dir_sig_ent[order(motif_dir_sig_ent$prob,decreasing=T),]
motif_dir_sig_ent$OMIM=NA
for(tf in unique(motif_dir_sig_ent$TF)){
  tf_in=gsub('\\(var.2\\)','',tf)
  tf_in=gsub('\\(var.3\\)','',tf_in)
  
  tf_in=unlist(strsplit(tf_in,'::'))
  #print(tf_in)
  OMIM_disease=OMIM$Phenotypes[which(unlist(lapply(OMIM$Gene.Symbols,function(x) any(x%in% tf_in))))]
  if(length(OMIM_disease)>0){motif_dir_sig_ent$OMIM[motif_dir_sig_ent$TF==tf]=as.list(OMIM_disease)}
}
motif_dir_sig_ent$OMIM=unlist(lapply(motif_dir_sig_ent$OMIM,function(x) paste(x,collapase=',')))
write.csv(motif_dir_sig_ent[c('TF','total_data','same_dir','opposite_dir','prob','qvalue','OMIM')],
          '../downstream/output/motif_prefer_ent_OMIM.csv')
#Motif family
human_mono_motif_TSV=as.data.frame(read.table('../downstream/input/JASPAR_human_redundant_2020.csv',sep=',',header=T,stringsAsFactors = F))
motif_family_enrich(unique(motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1]),unique(motif_gene$geneSymbol),human_mono_motif_TSV)
more_ent_enrich=motif_family_enrich(unique(motif_dir$TF[motif_dir$qval_binom<=0.1 & motif_dir$prob>0.5]),
                                    unique(motif_gene$geneSymbol),human_mono_motif_TSV)
motif_family_enrich(motif_dir_sig_ent_select$TF, unique(motif_gene$geneSymbol),human_mono_motif_TSV)
##testing for adding pseducount
# test_pse=data.frame()
# for (ps in seq(1,10000,100)){
#   more_ent_enrich=motif_family_enrich(unique(motif_dir$TF[motif_dir$qval_binom<=0.1 & motif_dir$prob>0.5]),
#                                       unique(motif_gene$geneSymbol),human_mono_motif_TSV,pse_count=ps)
#   test_pse=rbind(test_pse,data.frame(OR=more_ent_enrich$OR,pse=ps,family=more_ent_enrich$family,var=more_ent_enrich$variance))
# }
# ggplot(test_pse[test_pse$pse>=10,],aes(x=pse,y=var,color=family))+geom_point()+xlab('Pseduocount')+ylab('var(OR)')+
#   theme(legend.position = "none")
# ggplot(test_pse[test_pse$pse>=10,],aes(x=pse,y=var,color=family))+geom_point()+xlab('Pseduocount')+ylab('OR')+theme(legend.position = "none")

more_ent_enrich=more_ent_enrich[order(more_ent_enrich$OR,decreasing = F),]
more_ent_enrich$family=factor(more_ent_enrich$family,levels=more_ent_enrich$family)
pdf('../downstream/output/motif_family_100.pdf',width=20)
ggplot(more_ent_enrich,aes(x=family,y=OR))+
geom_bar(stat="identity",color="black",position=position_dodge(0.9),fill='lightblue')+
  #geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.2,position=position_dodge(0.9))+ 
  coord_flip()+ theme(axis.text.x = element_text(hjust = 1),legend.position = "none")+xlab("Motif family")+
  ggtitle("Motif family enrichment")+geom_text(aes(label=round(OR,digits = 2)), hjust=-0.5, color="black", size=3.5)
dev.off()
more_ent_enrich=more_ent_enrich[more_ent_enrich$OR>1,]
more_ent_enrich=more_ent_enrich[order(more_ent_enrich$OR,decreasing = T),]
write.csv(more_ent_enrich[,c(1,2,5,6,11)],'../downstream/output/motif_family_prefer_high_entropy.csv')
human_mono_motif_TSV_ent=human_mono_motif_TSV[human_mono_motif_TSV$Name%in%motif_dir$TF[motif_dir$qval_binom<=0.1 & motif_dir$prob>0.5],c(2,4,5)]
human_mono_motif_TSV_ent=unique(human_mono_motif_TSV_ent)
write.csv(human_mono_motif_TSV_ent[order(human_mono_motif_TSV_ent$Family),],
          '../downstream/output/motif_family_prefer_high_entropy_motif_metadata.csv')

less_ent_enrich=motif_family_enrich(unique(motif_dir$TF[motif_dir$qval_binom<=0.1& motif_dir$prob<0.5]),
                                    unique(motif_gene$geneSymbol),human_mono_motif_TSV)
less_ent_enrich=less_ent_enrich[less_ent_enrich$OR>1,]
less_ent_enrich=less_ent_enrich[order(less_ent_enrich$qvalue,decreasing = F),]
write.csv(less_ent_enrich[,c(1,2,5,6,11)],'../downstream/output/motif_family_prefer_low_entropy.csv')
less_ent_enrich=less_ent_enrich[order(less_ent_enrich$OR,decreasing = F),]
less_ent_enrich$family=factor(less_ent_enrich$family,levels=less_ent_enrich$family)
pdf('../downstream/output/motif_less_family_100.pdf')
ggplot(less_ent_enrich[less_ent_enrich$qvalue<=0.1,],aes(x=family,y=OR))+
  geom_bar(stat="identity",color="black",position=position_dodge(0.9),fill='lightblue')+
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.2,position=position_dodge(0.9))+ 
  coord_flip()+ theme(axis.text.x = element_text(hjust = 1),legend.position = "none")+xlab("GWAS traits")+
  ggtitle("GWAS trait enrichment")+geom_text(aes(label=round(OR,digits = 2)), hjust=-0.5, color="black", size=3.5)
dev.off()

# Motif vs literature -----------------------------------------------------

variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene <- readRDS(motif_gene_file)
motif_gene=subsetByOverlaps(motif_gene,variant_HetCpG_meta,type='equal')
motif_high_ent= read.csv('../downstream/output/motif_prefer_ent_OMIM.csv')
variant_motif_ent=subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],
                                   motif_gene[motif_gene$geneSymbol%in%motif_high_ent$TF])
JASPAR_bindingsite=readRDS('../downstream/output/motif_JASPAR_hg19.rds')
names(JASPAR_bindingsite)=unlist(lapply(strsplit(names(JASPAR_bindingsite),'_'),function(x) x[2]))
JASPAR_bindingsite_ent=JASPAR_bindingsite[motif_high_ent$TF]
#Nat comm
snv_literature1=readRDS('../downstream/output/OC_snv.rds')
snv_literature1=readRDS('../downstream/output/CHIP_SNP.rds')
snv_olap=GRanges()
for (TF in names(JASPAR_bindingsite_ent)){
  variant_motif_ent=subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],
                                     motif_gene[motif_gene$geneSymbol==TF])
  snv_olap=c(snv_olap,subsetByOverlaps(subsetByOverlaps(JASPAR_bindingsite[[TF]],variant_motif_ent),
                                       snv_literature1))
  if(length(snv_olap)>0){print(snv_olap)}
  
}
variant=subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],snv_olap)
GR_merge=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],variant)
motif=subsetByOverlaps(motif_gene,snv_olap)
source('plotMB.R')
plotMB(motif,'HUES64-1265282')
subsetByOverlaps(subsetByOverlaps(JASPAR_bindingsite,variant_motif_ent),snv_literature1)
subsetByOverlaps(variant_motif_ent,snv_literature1,maxgap = 10)#136217 vs 63929, none
subsetByOverlaps(snv_literature1,variant_motif_ent,maxgap = 10)
subsetByOverlaps(variant_motif_ent,CHIP_SNP_out)#136217 vs 63929, none
subsetByOverlaps(variant_HetCpG_meta,snv_literature1)#136217 vs 63929, none
#Nat genetics
CTCF_cluster=as.data.frame(read_xlsx('../downstream/input/Nature_gentic_supp4.xlsx',sheet = 'Mutation clusters'))
CTCF_cluster_cord=strsplit(CTCF_cluster$`Site coordinates (GRCh37)`,':')
CTCF_cluster_df=data.frame(chr=paste('chr',unlist(lapply(CTCF_cluster_cord,function(x) x[1])),sep=''),
                           region=unlist(lapply(CTCF_cluster_cord,function(x) x[2])),stringsAsFactors = F)
CTCF_cluster_df$start=unlist(lapply(strsplit(CTCF_cluster_df$region,'-'),function(x) x[1]))
CTCF_cluster_df$end=unlist(lapply(strsplit(CTCF_cluster_df$region,'-'),function(x) x[2]))
CTCF_cluster_gr=makeGRangesFromDataFrame(CTCF_cluster_df)
JASPAR_CTCF=JASPAR_bindingsite_ent[['CTCF']]
dNME_CTCF=subsetByOverlaps(JASPAR_CTCF,variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff])
dNME_CTCF_cluster=subsetByOverlaps(dNME_CTCF,CTCF_cluster_gr)

# NME_VMR -----------------------------------------------------------------
NME_in=readRDS(NME_agnostic_file)
#Brain
load("../downstream/input/vmrs_hg19_brain.rda")
vmr_HC2=vmrs_hg19$HC2
vmr_HC1=vmrs_hg19$HC1
names(vmr_HC2)=NULL
names(vmr_HC1)=NULL
#Do HC2
vmr=do.call(c,vmr_HC2)
saveRDS(vmr,'../downstream/output/vmr_HC2.rds')
vmr=readRDS('../downstream/output/vmr_HC2.rds')
NME_in_brain=NME_in[NME_in$Sample%in%c('Brain_Hippocampus_middle_paired - 149','Brain_Hippocampus_middle_paired - 150')]
olap=findOverlaps(NME_in_brain,vmr)
NME_in_brain$VMR=NA
NME_in_brain$VMR[queryHits(olap)]=vmr$meanSDS[subjectHits(olap)]
#cor.test(NME_in_brain$VMR,NME_in_brain$NME,method="spearman")#-0.01201045, 0.00784
NME_in_brain=NME_in_brain[order(NME_in_brain$NME,decreasing = T)]
NME_in_brain$quant_score_1p=findInterval(NME_in_brain$NME,quantile(NME_in_brain$NME,prob=seq(0,1,0.01)))
NME_in_brain$quant_score_1p=NME_in_brain$quant_score_1p-1
# OR_quant_1p=data.frame()
# 
# for(percent in unique(NME_in_brain$quant_score_1p)){
#   OR=OR_VMR(NME_in_brain,vmr,percent,NME_quant='quant_score_1p')
#   OR_quant_1p=rbind(OR_quant_1p,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
#   
# }
# ggplot(OR_quant_1p,aes(x=quant,y=OR))+geom_point()+geom_smooth()+xlab('NME quantile')+ylab('OR of VMR')
OR_quant=data.frame()
for(percent in unique(NME_in_brain$quant_score)){
  OR=OR_VMR(NME_in_brain,vmr,percent,NME_quant='quant_score')
  OR_quant=rbind(OR_quant,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
  
}

#Lung
load("../downstream/input/List_of_VMRs_lung.rda")
VMR_Lung=list_of_VMRs$Lung_all
NME_in_lung=NME_in[NME_in$Sample%in%c("Lung_single - STL001","Lung_single - STL002")]
OR_quant=data.frame()
for(percent in unique(NME_in_lung$quant_score)){
  OR=OR_VMR(NME_in_lung,VMR_Lung,percent,NME_quant='quant_score')
  OR_quant=rbind(OR_quant,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
  
}
OR_quant[order(OR_quant$OR,decreasing=T),]


##corSIV analysis
Waterland_CorSIV <- as.data.frame(read_excel("../downstream/input/Waterland_CorSIV.xls", sheet = "S3"),stringsAsFactors=F)
CorSIV_loc=strsplit(Waterland_CorSIV$USCS_Coordinates_CoRSIV,':')
CorSIV_loc_start=unlist(lapply(CorSIV_loc,function(x) x[1]))
CorSIV_loc_ranges=lapply(CorSIV_loc,function(x) IRanges(as.numeric(strsplit(x[2],'-')[[1]][1]),width=601))
CorSIV_gr=unique(GRanges(seqnames=CorSIV_loc_start,ranges=do.call('c',CorSIV_loc_ranges)))

OR_quant_corSIV_sp=list()
for(sp in unique(NME_in$Sample)){
  OR_quant_corSIV=data.frame()
  tt1=proc.time()[[3]]
  for(percent in unique(NME_in$quant_score)){
    OR=OR_VMR(NME_in[NME_in$Sample==sp],CorSIV_gr,percent,NME_quant='quant_score')
    OR_quant_corSIV=rbind(OR_quant_corSIV,
                          data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2],Sample=sp))
    OR_quant_corSIV_sp[[sp]]=OR_quant_corSIV
  }
  cat('Finish processing',sp,'in',proc.time()[[3]]-tt1,'\n')
}
saveRDS(OR_quant_corSIV_sp,'../downstream/output/OR_quant_corSIV_sp.rds')



# chekcing variant age ----------------------------------------------------
age_out=readRDS('../downstream/input/age_atlas_sub.rds')
qual_score=0.5
#10055361, 3575815 have SNP,867652 SNP, check source of variant id
age_out=age_out[age_out$QualScore_Jnt>=qual_score&age_out$QualScore_Mut>=qual_score&age_out$QualScore_Rec>=qual_score]
age_out=age_out[as.character(age_out$AlleleRef)==as.character(age_out$AlleleAnc)]
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
genome_freq=readRDS('../downstream/input/genome_1k_variant.rds')
olap=findOverlaps(variant_HetCpG_meta,age_out)
variant_HetCpG_meta$variant_age=NA
variant_HetCpG_meta$variant_age[queryHits(olap)]=age_out$AgeMedian_Jnt[subjectHits(olap)]
#variant_HetCpG_meta$Ref_age[queryHits(olap)]=gsub(' ','',as.character(age_out$AlleleRef[subjectHits(olap)]))
olap=findOverlaps(variant_HetCpG_meta,genome_freq)
variant_HetCpG_meta$DAF=NA
variant_HetCpG_meta$DAF[queryHits(olap)]=genome_freq$MAF[subjectHits(olap)]
variant_HetCpG_meta=variant_HetCpG_meta[!is.na(variant_HetCpG_meta$variant_age)]
variant_HetCpG_meta$dNME=variant_HetCpG_meta$altNME-variant_HetCpG_meta$refNME
variant_HetCpG_meta=subsetByOverlaps(variant_HetCpG_meta,genomic_features$promoter)
variant_age_out=data.frame()
all_df_out=data.frame()
for (sp in unique(variant_HetCpG_meta$Sample)){
  cat('processing sample:',sp,'\n')
  #dNME_sp_smaller=which(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$Sample==sp&variant_HetCpG_meta$dNME<0)
  #dNME_sp_larger=which(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$Sample==sp&variant_HetCpG_meta$dNME>0)
  dNME_sp=which(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$Sample==sp)
  dMML_sp=which(variant_HetCpG_meta$dMML_pval<=pval_cutoff&variant_HetCpG_meta$Sample==sp)
  non_dNME_dMML_sp=variant_HetCpG_meta$dMML_pval>pval_cutoff& variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$Sample==sp
  all_df_sp=rbind(data.frame(age=variant_HetCpG_meta$variant_age[dNME_sp],
                             DAF=variant_HetCpG_meta$DAF[dNME_sp],
                             type='dNME-Hap',stringsAsFactors = F),
                  data.frame(age=ifelse(length(dMML_sp)==0,NA,variant_HetCpG_meta$variant_age[dMML_sp]),
                             DAF=ifelse(length(dMML_sp)==0,NA,variant_HetCpG_meta$DAF[dMML_sp]),
                             type=ifelse(length(dMML_sp)==0,NA,'dMML-Hap'),stringsAsFactors = F),
                  data.frame(age=variant_HetCpG_meta$variant_age[non_dNME_dMML_sp], 
                             DAF=variant_HetCpG_meta$DAF[non_dNME_dMML_sp],
                             type='non-dMML-dNME-Hap',stringsAsFactors = F)
                  # data.frame(age=variant_HetCpG_meta$variant_age[non_dNME_dMML_sp], 
                  #            DAF=variant_HetCpG_meta$DAF[non_dNME_dMML_sp],
                  #            type='non-dMML-dNME-Hap',stringsAsFactors = F),
                  # data.frame(age=variant_HetCpG_meta$variant_age[dNME_sp_larger], 
                  #            DAF=variant_HetCpG_meta$DAF[dNME_sp_larger],
                  #            type='dNME-Hap-alt-larger',stringsAsFactors = F)
                  
  )
  
  all_df_sp=all_df_sp[!is.na(all_df_sp$DAF)&!is.na(all_df_sp$age),]
  all_df_sp$sample=sp
  all_df_sp$tissue=unique(variant_HetCpG_meta$tissue[dNME_sp])
  all_df_sp$subject=strsplit(sp,' - ')[[1]][2]
  all_df_out=rbind(all_df_out,all_df_sp)
  cat('plotting\n')
  # plot for each sample
  #  print(ggplot(all_df_sp,aes(x=DAF,y=age,color=type))+geom_smooth(size=1)+xlab('DAF')+ylab('variant age')+theme(legend.position = 'bottom')+
  #  ggtitle(sp))
  #   # ggplot(all_df)+geom_histogram(aes(x=round(DAF,digits = 1),y = ..density..,fill=type,color=type),binwidth=0.1,position="dodge")+xlab('DAF')+
  #   #   theme(legend.position = 'bottom')
  #   print(ggplot(all_df_sp,aes(x=age,color=type))+geom_density(size=1)+xlab('variant age')+theme(legend.position = 'bottom')+ggtitle(sp))
  #   print(ggplot(all_df_sp,aes(x=age,color=type))+stat_ecdf(size=1)+xlab('variant age')+theme(legend.position = 'bottom')+ggtitle(sp))
}
all_df_out$old=all_df_out$age<=12500
all_df_out_sp_old=aggregate(all_df_out$old,by=list(all_df_out$sample,all_df_out$type),FUN=mean)
names(all_df_out_sp_old)=c('Sample','type','Young_sample')
wilcox.test(all_df_out_sp_old$`Young_sample`[all_df_out_sp_old$type=='dNME-Hap'],all_df_out_sp_old$`Young_sample`[all_df_out_sp_old$type=='dMML-Hap'])
wilcox.test(all_df_out_sp_old$`Young_sample`[all_df_out_sp_old$type=='dNME-Hap'],all_df_out_sp_old$`Young_sample`[all_df_out_sp_old$type=='non-dMML-dNME-Hap'])
my_comparisons <- list( c("dNME-Hap-alt-larger", "non-dMML-dNME-Hap"), c("dNME-Hap-alt-larger", "dNME-Hap-alt-smaller"), c("dNME-Hap-alt-larger", "dMML-Hap"))
my_comparisons <- list( c("dNME-Hap", "non-dMML-dNME-Hap"), c("dNME-Hap", "dMML-Hap"),c("dMML-Hap", "non-dMML-dNME-Hap"))
my_comparisons <- list( c("dNME-Hap", "non-dMML-dNME-Hap"))
pdf('../downstream/output/age_percent12500.pdf')
ggplot(all_df_out_sp_old[all_df_out_sp_old$type%in%c('dNME-Hap','non-dMML-dNME-Hap'),],aes(x=type,y=Young_sample))+
  geom_boxplot(outlier.shape = NA)+ylab('Percent of regions <=12500 generations')+ 
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()
ggplot(all_df_out,aes(x=age,group=type))+geom_density(size=1,aes(color=type))+xlab('Age')+theme(legend.position = 'bottom')
pdf('../downstream/output/age_ecdf.pdf')
ecdf_df_out=data.frame()
for(type in unique(all_df_out$type)){
  all_df_out_estimate=ecdf(all_df_out$age[all_df_out$type==type])
  age_uq=seq(1,max(unique(all_df_out$age)),20)
  ecdf_df_out=rbind(ecdf_df_out,data.frame(age=age_uq,quant=all_df_out_estimate(age_uq),type=type,stringsAsFactors = F))
}
ggplot(ecdf_df_out,aes(x=age,y=quant,group=type,color=type))+geom_line(size=0.5)+xlab('Age')+theme(legend.position = 'bottom')
dev.off()
#Gene ontology analysis for those young SNPs
genes_sig=unique(variant_HetCpG_meta$genes_promoter[variant_HetCpG_meta$variant_age<=12500&variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$N>=1])
genes_all=unique(variant_HetCpG_meta$genes_promoter[variant_HetCpG_meta$variant_age<=12500&variant_HetCpG_meta$N>=1])
GO_young=GO_anno(genes_sig,genes_all)
GO_young_genes=GO_young[[1]]
GO_young_genes_sig=GO_young_genes[GO_young_genes$Significant>=10&GO_young_genes$FC>=2,]
GO_young_genes_sig$qval=p.adjust(GO_young_genes_sig$classicFisher)
GO_young_genes_sig=GO_young_genes_sig[GO_young_genes_sig$qval<=0.1,]
GO_young_obj=GO_young[[2]]
gene_out=genesInTerm(GO_young_obj, GO_young_genes_sig$GO.ID)
names(gene_out)=GO_young_genes_sig$Term
gene_out=lapply(gene_out,function(x) x[x%in%genes_sig])


# GWAS analysis -----------------------------------------------------------
#Check GWAS traits have hg19 coordinates
variant_trait=readRDS('../downstream/input/human_analysis/variant_traits.rds')
# variant_trait_gr=do.call(c,variant_trait)
# variant_trait_gr=unique(variant_trait_gr)
# genome_freq=readRDS('../downstream/input/genome_1k_variant.rds')
# genome_freq$variant_id=gsub('dbSNP_153\\:','',genome_freq$Dbxref)
# olap=findOverlaps(variant_trait,genome_freq,type='equal')
#Check proportion of dNME SNP from GWAS
#get traits location

variant_trait_gr=do.call(c,variant_trait)
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
#variant_HetCpG_meta=readRDS(GR_merge_file)
elementMetadata(variant_HetCpG_meta)=elementMetadata(variant_HetCpG_meta)[,c('dNME','dMML','dNME_pval','dMML_pval')]
length(subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],variant_trait_gr,maxgap = 0))
olap=findOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],variant_trait_gr)
variant_HetCpG_meta$GWAS=NA
variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff][queryHits(olap)]$GWAS=variant_trait_gr[subjectHits(olap)]$`DISEASE/TRAIT`


#  23827/380526= 0.06261596
#dMML:1236/14006=0.088
length(subsetByOverlaps(variant_HetCpG_meta,variant_trait_gr,maxgap = 1000))# 601586/7967588=0.076
length(subsetByOverlaps(variant_HetCpG_meta[!is.na(variant_HetCpG_meta$genes_promoter)],variant_trait,maxgap = 1000))# 7264/230407=0.0315, #ntrait>1-> 
#735/230407=0.0032
length(subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff&!is.na(variant_HetCpG_meta$genes_promoter)],
                        variant_trait))#180/6351=0.028,ntrait>1->11/6351=0.0017
#check gwas near motif
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
#motif_dir=readRDS('../downstream/output/motif_dir.rds')
motif_dir$qvalue=p.adjust(motif_dir$binom.pval,method='BH')
motif_dir_sig_ent=motif_dir[motif_dir$qval_binom<=0.1&motif_dir$prob>0.5,]
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene=readRDS(motif_gene_file)
motif_gene=motif_gene[motif_gene$geneSymbol%in%motif_dir_sig_ent$TF]
olap=findOverlaps(variant_HetCpG_meta,motif_gene)
variant_HetCpG_meta$Hign_ent_motif=FALSE
variant_HetCpG_meta$Hign_ent_motif[queryHits(olap)]=TRUE
variant_HetCpG_meta$dNME_pval[!variant_HetCpG_meta$Hign_ent_motif]=1
#GWAS enrichment at different ranges, mc.core settings
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=500)
saveRDS(dNME_traits,'../downstream/output/dNME_traits_500_gr.rds')
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=1000)
saveRDS(dNME_traits,'../downstream/output/dNME_traits_1000_motif_SNP.rds')
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=5000)
saveRDS(dNME_traits,'../downstream/output/dNME_traits_5k_gr.rds')
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=10000)
saveRDS(dNME_traits,'../downstream/output/dNME_traits_10k_gr.rds')
dNME_traits=readRDS('../downstream/output/dNME_traits_1000_gr2.rds')
#dNME_traits=readRDS('../downstream/output/dNME_traits_1000_motif_SNP.rds')
#downstream analysis
dNME_traits_df=do.call(rbind,lapply(dNME_traits,function(x) x[[1]]))
dNME_traits_gr=do.call(c,lapply(dNME_traits,function(x) x[[2]]))
dNME_traits_df=dNME_traits_df[dNME_traits_df$dNME_trait>=10&dNME_traits_df$OR>1,]
dNME_traits_df$qvalue=p.adjust(dNME_traits_df$p_value,method='BH')
dNME_traits_sig=dNME_traits_df[dNME_traits_df$qvalue<=0.1,]
dNME_traits_sig=dNME_traits_sig[order(dNME_traits_sig$OR,decreasing=F),]
write.csv(dNME_traits_sig,'../downstream/output/GWAS_dNME.csv')
#gene ontology on genes showing GWAS traits
dNME_traits_sig$trait=factor(dNME_traits_sig$trait,levels=dNME_traits_sig$trait)
dNME_traits_gr_sig=dNME_traits_gr[dNME_traits_gr$traits%in%dNME_traits_sig$trait]
dNME_traits_gr_sig_gene=unique(unlist(strsplit(unlist(strsplit(dNME_traits_gr_sig$genes,', ')),'- ')))
variant_trait_gr_background=subsetByOverlaps(variant_trait_gr,GR_merge,maxgap = 1000)
dNME_traits_gr_gene=unique(unlist(strsplit(unlist(strsplit(
  variant_trait_gr[variant_trait_gr$`DISEASE/TRAIT`%in% dNME_traits_gr$traits]$MAPPED_GENE,', ')),'- ')))
GO_GWAS=GO_anno(dNME_traits_gr_sig_gene,dNME_traits_gr_gene,topNodes=  10613)
GR_merge=readRDS(GR_merge_file)
GR_merge_traits=subsetByOverlaps(GR_merge,dNME_traits_gr_sig,maxgap = 1000)
GO_GWAS_body=GO_anno(unique(unlist(GR_merge_traits$genes_body)),unique(unlist(GR_merge$genes_body)))
pdf('../downstream/output/GWAS_traits_500_gr.pdf',width=15,height=15)
ggplot(dNME_traits_sig,aes(y=OR,x=trait))+
  geom_bar(stat="identity",color="black",position=position_dodge(0.9),fill='lightblue')+
  #geom_errorbar(aes(ymin=log(lower_CI),ymax=log(upper_CI)),width=0.2,position=position_dodge(0.9))+ 
  coord_flip()+ theme(axis.text.x = element_text(hjust = 1),legend.position = "none")+xlab("GWAS traits")+
  ggtitle("GWAS trait enrichment")+geom_text(aes(label=round(OR,digits = 2)), hjust=-0.5, color="black", size=3.5)
dev.off()


# vQTL and high NME regions -----------------------------------------------
NME_in=readRDS(NME_agnostic_file)
vQTL_SNP=read.csv('../downstream/input/vQTL.csv')
colnames(vQTL_SNP)=c('SNP','trait')
vQTL_loc=snpsById( SNPlocs.Hsapiens.dbSNP151.GRCh38,as.character(vQTL_SNP$SNP))
ch = import.chain('../downstream/input/hg38ToHg19.over.chain')
seqlevelsStyle(vQTL_loc) = 'UCSC'
vQTL_loc = liftOver(vQTL_loc, ch)
vQTL_loc=unlist(vQTL_loc)

vQTL_loc$trait=as.character(vQTL_SNP$trait[match(vQTL_loc$RefSNP_id,vQTL_SNP$SNP)])
NME_in_vQTL_out=GRanges()
for(traits in unique(vQTL_loc$trait)){
NME_in_vQTL=subsetByOverlaps(NME_in,vQTL_loc[vQTL_loc$trait==traits],maxgap =1000)
NME_in_vQTL$traits=traits
NME_in_vQTL_out=c(NME_in_vQTL_out,NME_in_vQTL)
}
NME_in_vQTL_df=rbind(data.frame(NME=NME_in_vQTL_out$NME,traits=NME_in_vQTL_out$traits),
                     data.frame(NME=NME_in$NME[sample(1:length(NME_in),15665,replace = F)],traits='All'))
my_comparisons <- list(c("HC", "All"), c("BMR", "All"), c("BFP", "All"),c("FFR", "All"),c("WC", "All"),c("BMI", "All"),c("BMD", "All"),c("WHR", "All"))
pdf('../downstream/output/vQTL_traits.pdf')
ggplot(NME_in_vQTL_df,aes(x=reorder(traits, NME, FUN = median),y=NME))+geom_boxplot(outlier.shape = NA)+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()




#Archive

#Archived
# Checking cutoff distribution at each N ----------------------------------

cutoff_N_min=data.frame()
for (i in unique(GR_merge$N)){
  cutoff_out=data.frame(dNME_cutoff=min(GR_merge$dNME[GR_merge$dNME_pval<=pval_cutoff & GR_merge$N==i]),
                        dMML_cutoff=min(GR_merge$dMML[GR_merge$dMML_pval<=pval_cutoff & GR_merge$N==i]),
                        N=i)
  cutoff_N_min=rbind(cutoff_N_min,cutoff_out)
}
cutoff_N_min[order(cutoff_N_min$N),]
ggplot(cutoff_N_min,aes(x=N,y=dNME_cutoff))+ geom_bar(stat="identity", position=position_dodge(),fill='light blue')+
  geom_text(aes(label=round(dNME_cutoff,digits=3)), position=position_dodge(1))+theme(legend.position = 'bottom')+ylim(0,1)

ggplot(cutoff_N_min,aes(x=N,y=dMML_cutoff))+ geom_bar(stat="identity", position=position_dodge(),fill='light blue')+
  geom_text(aes(label=round(dMML_cutoff,digits=3)), position=position_dodge(1))+theme(legend.position = 'bottom')+ylim(0,1)
# SNP type enrichment in ASM ----------------------------------------------

variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
#There's a slight depletion in A-T SNP in dMML-ASM
#variant_HetCpG_meta=variant_HetCpG_meta[variant_HetCpG_meta$N>1]
variant_HetCpG_meta$variants=variants_collapase(variant_HetCpG_meta)
#seleted_sample= unique(variant_HetCpG_meta$Subject)[1:9]
#variant_HetCpG_meta=variant_HetCpG_meta[variant_HetCpG_meta$Subject !="GM12878"]

# OR for each SNP ---------------------------------------------------------

#calculate OR for each type of SNP
#OR_df_out_sp=list()
OR_df_out=data.frame()
for (stat_in in c("dMML",'dNME',"UC")){
  OR_output=list()
  for (SNP_type in unique(variant_HetCpG_meta$variants)){OR_output[[SNP_type]]=variants_OR(variant_HetCpG_meta,SNP_type,stat_in)}
  OR_df=do.call(rbind,lapply(OR_output,function(x) data.frame(lower_CI=x$conf.int[1],upper_CI=x$conf.int[2],OR=x$estimate)))
  OR_df$stat=stat_in
  OR_df$SNP=names(OR_output)
  OR_df_out=rbind(OR_df_out,OR_df)
}
SNP_fill=data.frame(SNP=sort(unique(OR_df_out$SNP)),
                    color=c('red','coral','darkturquoise','orchid3','dodgerblue3','seagreen3'),stringsAsFactors = F)
OR_df_out$SNP_fill=SNP_fill$color[match(OR_df_out$SNP,SNP_fill$SNP)]
#Only for dNME and dMML
title='SNP enrichement'
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom",plot.title = element_text(hjust = 0.5))
ggplot(OR_df_out[OR_df_out$stat%in% c('dNME','dMML'),],aes(x=stat,y=OR,fill=SNP)) +ylim(c(0,2))+xlab("")+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+theme_bar+scale_fill_manual(values=SNP_fill$color)
-#dNME, dMML and UC
title='SNP enrichement'
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom",plot.title = element_text(hjust = 0.5))
ggplot(OR_df_out,aes(x=stat,y=OR,fill=SNP)) +ylim(c(0,2))+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+theme_bar+scale_fill_manual(values=SNP_fill$color)
#dNME SNP
title='dNME SNP enrichement'
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom",plot.title = element_text(hjust = 0.5))
ggplot(OR_df_out[OR_df_out$stat=='dNME',],aes(x=SNP,y=OR,fill=SNP)) +ylim(c(0,2))+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+theme_bar+scale_fill_manual(values=SNP_fill$color)

# OR for trinucleotide ----------------------------------------------------


#one trinucleotide
tri_OR_out_SNP_all=list()
for (stat_type in c("dMML",'dNME',"UC")){
  tri_OR_out_SNP=list()
  for(vari in unique(variant_HetCpG_meta$variants)){
    tri_OR_out=list()
    for(tri_mask in unique(variant_HetCpG_meta$mask_tri)){
      tri_OR_out[[tri_mask]]=tri_nucleo_OR(variant_HetCpG_meta
                                           [variant_HetCpG_meta$variants==vari],tri_mask,stat_type)
    }
    tri_OR_df=do.call(rbind,lapply(tri_OR_out,function(x) data.frame(lower_CI=x$conf.int[1],upper_CI=x$conf.int[2],OR=x$estimate)))
    tri_OR_df$tri=rownames(tri_OR_df)
    tri_OR_df$SNP=vari
    tri_OR_df=tri_OR_df[order(tri_OR_df$OR,decreasing = T),]
    tri_OR_out_SNP[[vari]]=tri_OR_df
  }
  tri_OR_out_SNP_all[[stat_type]]=tri_OR_out_SNP
}
#Plot for dMML or dNME or UC
dNME_tri_OR_all=tri_OR_out_SNP_all$dNME
dNME_tri_plot_ls=list()
for(vari in names(dNME_tri_OR_all)){
  title=vari
  
  theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="",plot.title = element_text(hjust = 0.5))
  dNME_tri_plot_ls[[vari]]=ggplot(dNME_tri_OR_all[[vari]],aes(x=reorder(tri,-OR),y=OR,fill=SNP)) +ylim(c(0,4))+
    geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
    geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                  position=position_dodge(.9))+theme_bar+xlab('')+scale_fill_manual(values=SNP_fill$color[SNP_fill$SNP==vari])
}
dNME_tri_plot_ls=dNME_tri_plot_ls[order(names(dNME_tri_plot_ls))]
grid.arrange(grobs=dNME_tri_plot_ls, nrow = 2)


# HetCpG enrichment -------------------------------------------------------


variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
variant_HetCpG_meta$ASM=NA
variant_HetCpG_meta$ASM[variant_HetCpG_meta$dNME_pval<=pval_cutoff]='Yes'
variant_HetCpG_meta$ASM[variant_HetCpG_meta$dNME_pval>pval_cutoff]='No'
variant_HetCpG_meta=variant_HetCpG_meta[!is.na(variant_HetCpG_meta$ASM)]
variant_HetCpG_meta$HetCpG=variant_HetCpG_meta$g1CG!=variant_HetCpG_meta$g2CG
  ASM_enrich=data.frame(sp=NULL,OR=NULL,lower_CI=NULL,upper_CI=NULL)
  for (sp in unique(variant_HetCpG_meta$Sample)){
    ASM_enrich_stat=ASM_het_enrichment(unique(variant_HetCpG_meta[variant_HetCpG_meta$Sample==sp&variant_HetCpG_meta$N>=2]))
    ASM_enrich=rbind(ASM_enrich,data.frame(sp=sp,subjects=strsplit(sp,' - ')[[1]][2],
                                           OR=ASM_enrich_stat$estimate,
                                           lower_CI=ASM_enrich_stat$conf.int[1],
                                           upper_CI=ASM_enrich_stat$conf.int[2]))
  }
#Box plot of all types of statistics, try a CMH test here
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom", plot.title = element_text(hjust=0.5))

ggplot(ASM_enrich,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0,5)+
  ggtitle('dNME Het CpG enrichment')+xlab('Sample name')+ylab('Odds Ratio')+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+theme_bar




# dNME vs local density and ecdf ------------------------------------------


NME_ASM=readRDS('../downstream/output/dNME/NME_ASM_het_calc_run3.rds')
#Plot cumulative distribution of NME with HetCpG and without HetCpG at ASM region
plot(ecdf(NME_ASM$diff_NME[NME_ASM$CpGdiff!=0]),main='cumulative distribution of dNME at ASM',xlab='dNME',xlim=c(-1,1),lwd=1)
lines(ecdf(NME_ASM$diff_NME[NME_ASM$CpGdiff==0]),col='red',lwd=3)
#Plot density difference vs dNME at ASM
density_df=data.frame(density_diff=round(log10(abs(NME_ASM$density_diff[NME_ASM$density_diff!=0])),digits=3),
                      dNME=NME_ASM$diff_NME[NME_ASM$density_diff!=0])
density_df_agg=aggregate(abs(density_df$dNME),by=list(density_df$density_diff),FUN=median)
ggplot(density_df_agg,aes(x=Group.1, y=x))+
  ylim(c(0,1))+ggtitle("dNME change as density difference change")+geom_smooth(method="lm")+
  ylab("dNME")+xlab("log10(density difference)")+geom_point()




# distance to genomic features, find odds ratio instead -------------------------------------------
#OR calculation
GR_merge=readRDS(GR_merge_file)
genomic_features=readRDS(genomic_features_file)
OR_feature_out=data.frame()
for (stat_type in c('dNME_pval','dMML_pval')){
  for(feature in c('CpG island','CpG shore','CpG shelf','CpG open sea')){
    GR_merge$ASM='No'
    GR_merge$ASM[mcols(GR_merge)[[stat_type]]<=pval_cutoff]='Yes'
    OR_out=testEnrichmentFeature_stat(GR_merge,genomic_features[[feature]])[[2]]
    OR_feature_out=rbind(OR_feature_out,data.frame(feature=feature,statistics=stat_type,OR=OR_out$estimate,pvalue=OR_out$p.value,lowerCI=OR_out$conf.int[1],upperCI=OR_out$conf.int[2]))
  }
}
pdf('../downstream/output/genomic_feature_enrich.pdf')
ggplot(OR_feature_out,aes(x=feature,y=OR,fill=statistics))+geom_bar(stat="identity", position=position_dodge())+theme(legend.position = 'bottom')+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9)) 
dev.off()

#filtering CpG islands bed
CpG_islands=readRDS("../downstream/input/cpg_islands_hg19.rds")

GR_merge=readRDS(GR_merge_file)
promoter_island_overlap=findOverlaps(genomic_features$promoter,genomic_features$`CpG island`)
promoter_with_CpG_islands=genomic_features$promoter[queryHits(promoter_island_overlap)]
promoter_without_CpG_islands=genomic_features$promoter[-queryHits(promoter_island_overlap)]
dist_plot=5000
#Find distnace to CpG island for different regions
###Calcualte distance to the feature
GR_merge=dist_calc(GR_merge,genomic_features$`CpG island`)

plot_df=rbind(data.frame(dist=GR_merge$dist[GR_merge$dMML_pval<=pval_cutoff&abs(GR_merge$dist)<=dist_plot],stat='dMML_ASM'),
              data.frame(dist=GR_merge$dist[GR_merge$dNME_pval<=pval_cutoff&abs(GR_merge$dist)<=dist_plot],stat='dNME_ASM'),
              data.frame(dist=GR_merge$dist[abs(GR_merge$dist)<=dist_plot],stat='All'))
plot_df$dist_round=round(plot_df$dist, digits=-2)
plot_df_density=data.frame()
for(stats in unique(plot_df$stat)){
  tb=table(plot_df$dist_round[plot_df$stat==stats])
  plot_df_density=rbind(plot_df_density,
                        data.frame(stat=stats,dist_round=as.numeric(names(tb)),
                                   density_dist=as.numeric(tb/sum(tb))))
}
pdf('../downstream/output/distance_to_CpG_islands_features.pdf')
ggplot(plot_df_density,aes(x=dist_round,y=density_dist,color=stat))+geom_line(size=1)+
  xlab('Distance to CpG islands (bp)')+theme(legend.position = 'bottom')+ylab('density')
dev.off()
plot(density(plot_df$dist_round[plot_df$stat=='dMML_ASM']))
###Promoters with or without islands
GR_merge=dist_calc(GR_merge,promoter_without_CpG_islands)

plot_df=rbind(data.frame(dist=GR_merge$dist[GR_merge$dMML_pval<=pval_cutoff&abs(GR_merge$dist)<=dist_plot],stat='dMML_ASM'),
              data.frame(dist=GR_merge$dist[GR_merge$dNME_pval<=pval_cutoff&abs(GR_merge$dist)<=dist_plot],stat='dNME_ASM'),
              data.frame(dist=GR_merge$dist[abs(GR_merge$dist)<=dist_plot],stat='All'))
ggplot(plot_df,aes(x=dist,color=stat))+geom_density(size=1)+xlab('Distance to CpG islands (bp)')+theme(legend.position = 'bottom')
###Enhancers
enhancers <- readRDS("../downstream/input/enchancer_DNase.rds") #474,004 enhancer region
GR_merge=dist_calc(GR_merge,enhancers)
dist_plot=5000
plot_df=rbind(data.frame(dist=GR_merge$dist[GR_merge$dMML_pval<=pval_cutoff&abs(GR_merge$dist)<=dist_plot],stat='dMML_ASM'),
              data.frame(dist=GR_merge$dist[GR_merge$dNME_pval<=pval_cutoff&abs(GR_merge$dist)<=dist_plot],stat='dNME_ASM'),
              data.frame(dist=GR_merge$dist[abs(GR_merge$dist)<=dist_plot],stat='All'))
ggplot(plot_df,aes(x=dist,color=stat))+geom_density(size=1)+xlab('enhancers (bp)')+theme(legend.position = 'bottom')




# Number of CpG analyzed --------------------------------------------------


CpG_hg19=readRDS("../downstream/input/CpG_hg19.rds")
olap_GR=findOverlaps(CpG_hg19,GR_merge)
length(unique(CpG_hg19[queryHits(olap_GR)]))#3385511 #13%
#Find invidual sample coverage
individual_CG=data.frame(sample=NULL,nCG=NULL)
for (sp in unique(GR_merge$Sample)){
  GR_individual=GR_merge[GR_merge$Sample==sp]
  individual_hg=data.frame(sample=sp,nCG=length(subsetByOverlaps(CpG_hg19,GR_individual)))
  individual_CG=rbind(individual_CG,individual_hg)
}
###CpG island enrichment
genomic_features=readRDS("../downstream/input/genomic_features2020.rds")
CpG_islands_olap=findOverlaps(CpG_hg19,genomic_features$`CpG island`)
total_CpG=length(CpG_hg19)#26752702
GR_islands=sum(queryHits(olap_GR) %in% queryHits(CpG_islands_olap))
GR_not_islands=sum(queryHits(olap_GR) %in% which(!((1:total_CpG) %in% queryHits(CpG_islands_olap))))
not_GR_not_islands=sum((!(1:total_CpG) %in% queryHits(olap_GR)) & !((1:total_CpG) %in% queryHits(CpG_islands_olap)))
not_GR_islands=sum(which(!((1:total_CpG) %in% queryHits(olap_GR))) %in% queryHits(CpG_islands_olap))
#CpG analyzed in island: 215957
#CpG analyzed not in island: 3385511-215957=3169554, 3385511-487832= 2897679
#CpG not analyzed in island: 2742391,4488865
#CpG not analyzed not in island: 26752702-3385511-2742391=20624800,26752702-3385511-4488865=18878326
fisher.test(matrix(c(GR_islands,GR_not_islands,not_GR_islands,not_GR_not_islands),nrow=2)) #OR=0.7572037, new run: 1.1




# DEVG analysis -----------------------------------------------------------
GR_merge=readRDS(GR_merge_file)
DEVG_out=data.frame()
for(tissue in unique(GR_merge$tissue)){
  DEVG_out=rbind(DEVG_out,DEVG_analysis(GR_merge,tissue))
  
}
DEVG_out=DEVG_out[apply(DEVG_out[,1:4],1,function(x) any(x!=0)),]
DEVG_out_percent=as.data.frame(t(apply(DEVG_out[,1:4],1,function(x) x/sum(x)*100)))
DEVG_out_percent$tissue=DEVG_out$tissue

DEVG_gr=DEVG_analysis(GR_merge,NA)
DEVG_genes=lapply(DEVG_gr,function(x) unique(x$TSS))
hign_NME_combined=GO_anno(c(DEVG_genes$high_NME_two_hyprvarible,DEVG_genes$high_NME_one_hyprvarible),unique(unlist(DEVG_genes)))
hign_NME_combined_sig=hign_NME_combined[[1]][hign_NME_combined[[1]]$Significant>=5&hign_NME_combined[[1]]$FC>=2&hign_NME_combined[[1]]$classicFisher<=0.05,]
#Check genes
all_genes=unique(unlist(DEVG_genes))
gene_term=genesInTerm(hign_NME_combined[[2]],hign_NME_combined_sig$GO.ID)
names(gene_term)=hign_NME_combined_sig$Term
gene_term=lapply(gene_term,function(x) x[x%in%all_genes])
cat(unique(unlist(gene_term)),sep=',')

hign_NME_two_GO=GO_anno(DEVG_genes$high_NME_two_hyprvarible,unique(unlist(DEVG_genes)))
head(hign_NME_two_GO[[1]][hign_NME_two_GO[[1]]$Significant>=10&hign_NME_two_GO[[1]]$FC>=1.5,],n=30)
hign_NME_one_GO=GO_anno(DEVG_genes$high_NME_one_hyprvarible,unique(unlist(DEVG_genes)))
head(hign_NME_one_GO[[1]][hign_NME_one_GO[[1]]$Significant>=10&hign_NME_one_GO[[1]]$FC>=1.5,],n=30)
low_dNME_non_dNME_GO=GO_anno(DEVG_genes$low_dNME_non_dNME,unique(unlist(DEVG_genes)))
head(low_dNME_non_dNME_GO[[1]][low_dNME_non_dNME_GO[[1]]$Significant>=5&low_dNME_non_dNME_GO[[1]]$FC>=1.5,])

#DEVG age analysis
DEVG_age=lapply(DEVG_gr,age_out=age_out,function(x,age_out) {
                if(length(x)>0){
                olap=findOverlaps(x,age_out)
                x$age=NA
                x$age[queryHits(olap)]=age_out$AgeMedian_Jnt[subjectHits(olap)]
                return(x)}})
DEVG_age_df=rbind(data.frame(NME_type='high_NME_two_hyper_Variable',age=DEVG_age$high_NME_two_hyprvarible$age),
                  data.frame(NME_type='high_NME_one_hyper_Variable',age=DEVG_age$high_NME_one_hyprvarible$age),
                  data.frame(NME_type='low_dNME_non_dNME',age=DEVG_age$low_dNME_non_dNME$age)
)

ggplot(DEVG_age_df,aes(x=age,color=NME_type))+stat_ecdf(size=1)+theme(legend.position = 'bottom')

#Promoter
DEVG_genes=DEVG_analysis(GR_merge,NA,stat_var='hyper_var_promoter')
DEVG_genes=lapply(DEVG_genes,function(x) unique(x$genes_promoter))
hign_NME_two_GO=GO_anno(DEVG_genes$high_NME_two_hyprvarible,unique(unlist(DEVG_genes)))
head(hign_NME_two_GO[[1]][hign_NME_two_GO[[1]]$Significant>=10&hign_NME_two_GO[[1]]$FC>=1.5,],n=30)
hign_NME_one_GO=GO_anno(DEVG_genes$high_NME_one_hyprvarible,unique(unlist(DEVG_genes)))
head(hign_NME_one_GO[[1]][hign_NME_one_GO[[1]]$Significant>=5&hign_NME_one_GO[[1]]$FC>=1.5,],n=30)
low_dNME_non_dNME_GO=GO_anno(DEVG_genes$low_dNME_non_dNME,unique(unlist(DEVG_genes)))
head(low_dNME_non_dNME_GO[[1]][low_dNME_non_dNME_GO[[1]]$Significant>=5&low_dNME_non_dNME_GO[[1]]$FC>=1.5,])

# Check if new dNME SNP is due to C-T mutation ----------------------------
variant_HetCpG_meta$variants=paste(variant_HetCpG_meta$REF,variant_HetCpG_meta$ALT,sep='-')

#OR=1.26
young_CT=sum(variant_HetCpG_meta$variant_age<=20000&variant_HetCpG_meta$variants=='C-T')
old_CT=sum(variant_HetCpG_meta$variant_age>20000&variant_HetCpG_meta$variants=='C-T')
young_non_CT=sum(variant_HetCpG_meta$variant_age<=20000&variant_HetCpG_meta$variants!='C-T')
old_non_CT=sum(variant_HetCpG_meta$variant_age>20000&variant_HetCpG_meta$variants!='C-T')
fisher.test(matrix(c(young_CT,old_CT,young_non_CT,old_non_CT),nrow=2))

young_CT=sum(variant_HetCpG_meta$variant_age<=20000&variant_HetCpG_meta$variants%in%c('C-T','G-A'))
old_CT=sum(variant_HetCpG_meta$variant_age>20000&variant_HetCpG_meta$variants%in%c('C-T','G-A'))
young_non_CT=sum(variant_HetCpG_meta$variant_age<=20000&!variant_HetCpG_meta$variants%in%c('C-T','G-A'))
old_non_CT=sum(variant_HetCpG_meta$variant_age>20000&!variant_HetCpG_meta$variants%in%c('C-T','G-A'))
fisher.test(matrix(c(young_CT,old_CT,young_non_CT,old_non_CT),nrow=2))

# young_CT=sum(variant_HetCpG_meta$variant_age<=20000&(variant_HetCpG_meta$g1CG>variant_HetCpG_meta$g2CG))
# old_CT=sum(variant_HetCpG_meta$variant_age>20000&variant_HetCpG_meta$g1CG>variant_HetCpG_meta$g2CG)
# young_non_CT=sum(variant_HetCpG_meta$variant_age<=20000&variant_HetCpG_meta$g1CG<=variant_HetCpG_meta$g2CG)
# old_non_CT=sum(variant_HetCpG_meta$variant_age>20000&variant_HetCpG_meta$g1CG<=variant_HetCpG_meta$g2CG)
# fisher.test(matrix(c(young_CT,old_CT,young_non_CT,old_non_CT),nrow=2))

# young_CT=sum(variant_HetCpG_meta$variant_age<=10000&(variant_HetCpG_meta$dNME_pval<=pval_cutoff))
# old_CT=sum(variant_HetCpG_meta$variant_age>10000&variant_HetCpG_meta$dNME_pval<=pval_cutoff)
# young_non_CT=sum(variant_HetCpG_meta$variant_age<=10000&variant_HetCpG_meta$dNME_pval>pval_cutoff)
# old_non_CT=sum(variant_HetCpG_meta$variant_age>10000&variant_HetCpG_meta$dNME_pval>pval_cutoff)
# fisher.test(matrix(c(young_CT,old_CT,young_non_CT,old_non_CT),nrow=2))


# young_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&(variant_HetCpG_meta$g1CG!=variant_HetCpG_meta$g2CG))
# old_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$g1CG!=variant_HetCpG_meta$g2CG)
# young_non_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$g1CG==variant_HetCpG_meta$g2CG)
# old_non_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$g1CG==variant_HetCpG_meta$g2CG)
# fisher.test(matrix(c(young_CT,old_CT,young_non_CT,old_non_CT),nrow=2))
density_dNME_df=data.frame()
GR_merge$Sample2=paste(GR_merge$Subject,'_',GR_merge$tissue,'_phased',sep='')
for (sp in unique(GR_merge$Sample2)){
  density_dNME=sum(GR_merge$dNME_pval<=pval_cutoff&(GR_merge$density_diff!=0) & GR_merge$Sample2==sp)
  density_non_dNME=sum(GR_merge$dNME_pval>pval_cutoff&GR_merge$density_diff!=0 & GR_merge$Sample2==sp)
  nondensity_dNME=sum(GR_merge$dNME_pval<=pval_cutoff&GR_merge$density_diff==0 & GR_merge$Sample2==sp)
  nondensity_nondNME=sum(GR_merge$dNME_pval>pval_cutoff&GR_merge$density_diff==0 & GR_merge$Sample2==sp)
  density_dNME_df=rbind(density_dNME_df,data.frame(sample=sp,
                                                   OR=fisher.test(matrix(c(density_dNME,density_non_dNME,nondensity_dNME,nondensity_nondNME),nrow=2))$estimate))
  
}
CT_dNME_df=data.frame()
sample2_df=data.frame(Sample1=unique(GR_merge$Sample),sample2=unique(GR_merge$Sample2))
variant_HetCpG_meta$Sample2=sample2_df$sample2[match(variant_HetCpG_meta$Sample,sample2_df$Sample1)]
for (sp in unique(variant_HetCpG_meta$Sample2)){
  dNME_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants=='C-T'&variant_HetCpG_meta$Sample2==sp)
  nondNME_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$variants=='C-T'&variant_HetCpG_meta$Sample2==sp)
  dNME_non_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants!='C-T'&variant_HetCpG_meta$Sample2==sp)
  nondNME_non_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$variants!='C-T'&variant_HetCpG_meta$Sample2==sp)
  CT_dNME_df=rbind(CT_dNME_df,data.frame(sample=sp,
                                                   OR=fisher.test(matrix(c(dNME_CT,nondNME_CT,dNME_non_CT,nondNME_non_CT),nrow=2))$estimate))
  
}



CT_dNME_df=data.frame()
sample2_df=data.frame(Sample1=unique(GR_merge$Sample),sample2=unique(GR_merge$Sample2))
variant_HetCpG_meta$Sample2=sample2_df$sample2[match(variant_HetCpG_meta$Sample,sample2_df$Sample1)]
for (sp in unique(variant_HetCpG_meta$Sample2)){
  dNME_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants%in%c('C-T','G-A')&variant_HetCpG_meta$Sample2==sp)
  nondNME_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$variants%in%c('C-T','G-A')&variant_HetCpG_meta$Sample2==sp)
  dNME_non_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&!variant_HetCpG_meta$variants%in%c('C-T','G-A')&variant_HetCpG_meta$Sample2==sp)
  nondNME_non_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&!variant_HetCpG_meta$variants%in%c('C-T','G-A')&variant_HetCpG_meta$Sample2==sp)
  CT_dNME_df=rbind(CT_dNME_df,data.frame(sample=sp,
                                         OR=fisher.test(matrix(c(dNME_CT,nondNME_CT,dNME_non_CT,nondNME_non_CT),nrow=2))$estimate))
  
}

cov_sample_g2=aggregate(coverage_all$genome2_cov[coverage_all$genome1_cov>=5&coverage_all$genome2_cov>=5],
                        by=list(coverage_all$Sample[coverage_all$genome1_cov>=5&coverage_all$genome2_cov>=5]), mean)
cov_sample_g1=aggregate(coverage_all$genome1_cov[coverage_all$genome1_cov>=5&coverage_all$genome2_cov>=5],
                        by=list(coverage_all$Sample[coverage_all$genome1_cov>=5&coverage_all$genome2_cov>=5]), mean)

cov_sample_g1$Group.1[cov_sample_g1$Group.1=="H1_merged"]='H1_merged_phased'
cov_sample_g1$Group.1[cov_sample_g1$Group.1=="STL003_Adipose_Tissue_single_phased"]='STL003_Adipose_single_phased'
density_dNME_df$g1cov=cov_sample_g1$x[match(density_dNME_df$sample,cov_sample_g1$Group.1)]
density_dNME_df$g2cov=cov_sample_g2$x[match(density_dNME_df$sample,cov_sample_g1$Group.1)]
density_dNME_df=density_dNME_df[!is.na(density_dNME_df$g1cov),]#cor=-0.787
cor(density_dNME_df$OR,density_dNME_df$g1cov)
plot(density_dNME_df$OR,density_dNME_df$g1cov,xlab='OR',ylab='coverage')

CT_dNME_df$g1cov=cov_sample_g1$x[match(CT_dNME_df$sample,cov_sample_g1$Group.1)]
CT_dNME_df$g2cov=cov_sample_g2$x[match(CT_dNME_df$sample,cov_sample_g1$Group.1)]
CT_dNME_df=CT_dNME_df[!is.na(CT_dNME_df$g1cov),]#cor=-0.787
cor(CT_dNME_df$OR,CT_dNME_df$g1cov)
plot(CT_dNME_df$OR,CT_dNME_df$g1cov,xlab='OR',ylab='coverage')
#increase entropy
median(c(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME>0&variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants=='C-T']))
median(c(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME>0&variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants!='C-T']))
median(c(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$variants=='C-T']))
median(c(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$variants!='C-T']))
median(c(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME<0&variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants=='C-T']))
median(c(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME<0&variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants!='C-T']))

#OR=0.918
for (var in unique(variant_HetCpG_meta$variants)){
  dNME_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants==var)
  nondNME_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$variants==var)
  dNME_non_CT=sum(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$variants!=var)
  nondNME_non_CT=sum(variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$variants!=var)
  cat(var,':',fisher.test(matrix(c(dNME_CT,nondNME_CT,dNME_non_CT,nondNME_non_CT),nrow=2))$estimate,'\n')
}


dNME_smaller=which(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$dNME<0&!is.na(variant_HetCpG_meta$TSS))
dNME_larger=which(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$dNME>0&!is.na(variant_HetCpG_meta$TSS))
non_dNME=which(variant_HetCpG_meta$dNME_pval>pval_cutoff&!is.na(variant_HetCpG_meta$TSS))
all_df=rbind(data.frame(age=variant_HetCpG_meta$variant_age[dNME_smaller],
                           DAF=variant_HetCpG_meta$DAF[dNME_smaller],
                           type='dNME-Hap-alt-smaller',stringsAsFactors = F),
                data.frame(age=variant_HetCpG_meta$variant_age[non_dNME], 
                           DAF=variant_HetCpG_meta$DAF[non_dNME],
                           type='non-dNME-Hap',stringsAsFactors = F),
                data.frame(age=variant_HetCpG_meta$variant_age[dNME_larger], 
                           DAF=variant_HetCpG_meta$DAF[dNME_larger],
                           type='dNME-Hap-alt-larger',stringsAsFactors = F)
                
)

#compare with ken's result
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_agnostic_ken=read.csv('../downstream/output/motif_enrich_high_NME.csv')
motif_agnostic_ken$TF=unlist(lapply(strsplit(as.character(motif_agnostic_ken$motif),'_'),function(x) x[2]))
motif_agnostic_ken=motif_agnostic_ken[motif_agnostic_ken$TF %in% motif_dir$TF,]
motif_agnostic_ken$rank=rank(-(motif_agnostic_ken$ave_score))
motif_dir=motif_dir[motif_dir$TF %in% motif_agnostic_ken$TF,]
motif_dir$rank=rank(-motif_dir$prob)
cor_df=data.frame(TF=motif_agnostic_ken$TF,ken_rank=motif_agnostic_ken$rank,stringsAsFactors = F)
cor_df$allele_rank=motif_dir$rank[match(cor_df$TF,motif_dir$TF)]#cor=0.11
cor_df$allele_rank_qvalue=motif_dir$qvalue[match(cor_df$TF,motif_dir$TF)]#cor=0.11
#cor_df=cor_df[cor_df$allele_rank_qvalue<=0.1,]
plot(cor_df$ken_rank,cor_df$allele_rank,xlab='allele-specific motif rank',ylab='allele-agnositc motif rank')
cor.test(cor_df$ken_rank,cor_df$allele_rank)
#N>=1 vs N>=2
motif_dir_N2=readRDS('../downstream/output/motif_dir_N2.rds')
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_dir_N2$rank=rank(-motif_dir_N2$prob)
motif_dir$rank=rank(-motif_dir$prob)
motif_dir=motif_dir[motif_dir$TF %in% motif_dir_N2$TF,]
cor_df=data.frame(TF=motif_dir_N2$TF,N2_rank=motif_dir_N2$rank,stringsAsFactors = F)
cor_df$allele_rank=motif_dir$rank[match(cor_df$TF,motif_dir$TF)]#cor=0.11
plot(cor_df$N2_rank,cor_df$allele_rank,xlab='rank N>=1',ylab='rank N>=2')
cor.test(cor_df$N2_rank,cor_df$allele_rank)#Cor =0.5
#Promoter vs non promoter
motif_dir_promoter=motif_dir
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_dir_promoter$rank=rank(-motif_dir_promoter$prob)
motif_dir$rank=rank(-motif_dir$prob)
motif_dir=motif_dir[motif_dir$TF %in% motif_dir_promoter$TF,]
cor_df=data.frame(TF=motif_dir_promoter$TF,promoter_rank=motif_dir_promoter$rank,stringsAsFactors = F)
cor_df$allele_rank=motif_dir$rank[match(cor_df$TF,motif_dir$TF)]
plot(cor_df$promoter_rank,cor_df$allele_rank,xlab='rank promoter',ylab='rank all')
cor.test(cor_df$promoter_rank,cor_df$allele_rank) #Cor =0.1

#Promoter vs non promoter
motif_dir_promoter=motif_dir
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_dir_promoter$rank=rank(-motif_dir_promoter$prob)
motif_dir$rank=rank(-motif_dir$prob)
motif_dir=motif_dir[motif_dir$TF %in% motif_dir_promoter$TF,]
cor_df=data.frame(TF=motif_dir_promoter$TF,promoter_rank=motif_dir_promoter$rank,stringsAsFactors = F)
cor_df$allele_rank=motif_dir$rank[match(cor_df$TF,motif_dir$TF)]#cor=0.11
plot(cor_df$promoter_rank,cor_df$allele_rank,xlab='rank DNAaseI',ylab='rank all')
cor.test(cor_df$promoter_rank,cor_df$allele_rank)#Cor =0.5

#UC vs dMML and dNME
UC_dt=data.table(rbind(
  data.table(data=GR_merge$dNME[GR_merge$UC>=0.5],stats="dNME"),
  data.table(data=GR_merge$dMML[GR_merge$UC>=0.5],stats="dMML")))

ggplot(UC_dt,aes(x=data,fill=stats))+geom_density(alpha=0.5)
