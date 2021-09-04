# Genomics
# Source main functions
#setwd("~/code/HASM-MetaAnalysis/")
rm(list=ls())
source("mainFunctions_sub.R")
#Loading data
genomic_features=readRDS(genomic_features_file)
GR_merge=readRDS(GR_merge_file)
#Only use merged data for H1
#Define ggplot theme
theme_glob=theme(plot.title = element_text(hjust = 0.5,size=36),
                 axis.title.x=element_text(hjust=0.5,size=36,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=36,face="bold"),
                 axis.text.x=element_text(size=32),
                 axis.text.y=element_text(size=32),
                 legend.text = element_text(size=32),
                 legend.title = element_blank())+theme_classic()
#Figure S1:

selected_features=c("CpG island","CpG shore","CpG shelf","CpG open sea","gene body","exon","intron","intergenic","promoter","TSS")

genomic_features_OR=lapply(selected_features,function(x){
  #NME enrichment
  GR_merge$ASM="No"
  GR_merge$ASM[GR_merge$dNME_pval<=pval_cutoff]="Yes"
  NME_enrich=testEnrichmentFeature_stat(GR_merge,genomic_features[[x]],maxgap=0)[[2]]
  #MML_enrichment
  GR_merge$ASM="No"
  GR_merge$ASM[GR_merge$dMML_pval<=pval_cutoff]="Yes"
  MML_enrich=testEnrichmentFeature_stat(GR_merge,genomic_features[[x]],maxgap=0)[[2]]

  return(rbind(data.table(OR=NME_enrich$estimate,pvalue=NME_enrich$p.value,
                          lowerCI=NME_enrich$conf.int[1],upperCI=NME_enrich$conf.int[2],region_type="dNME-ASM",feature=x),
               data.table(OR=MML_enrich$estimate,pvalue=NME_enrich$p.value,
                          lowerCI=MML_enrich$conf.int[1],upperCI=MML_enrich$conf.int[2],region_type="dMML-ASM",feature=x)))
})

genomic_features_OR=do.call(rbind,genomic_features_OR)
genomic_features_OR$star=add.significance.stars(genomic_features_OR$pvalue, cutoffs = c(0.05, 0.01, 0.001))
genomic_features_OR$star=gsub(' ','',genomic_features_OR$star)
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
pdf('../downstream/output/graphs/FigureS1/feature_enrichment.pdf',width=10,height=5)
print(ggplot(genomic_features_OR,aes(x=feature,y=OR,fill=region_type)) + 
  geom_bar(stat="identity", position=position_dodge())+ylim(0, 6.5)+
  ggtitle("")+xlab('Sample name')+ylab('Odds Ratio')+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,
                position=position_dodge(.9))+  
  geom_text(aes(y=upperCI,label=round(OR,digits = 2)),vjust=-0.5,position = position_dodge(1),size=3)+
  geom_text(aes(y=upperCI,label=star),vjust=-1.5,position = position_dodge(1))+
  theme_bar+theme_glob+theme(legend.position = 'bottom'))
dev.off()

#Figure 2A,B
library(readxl)
Imprinted_Genes <- as.data.frame(read_excel("../downstream/input/human_analysis/imprinting_ASE/Imprinted Genes.xlsx"))
Imprinted_Genes$Expressed_allele= Imprinted_Genes$`Expressed Allele`
#ecdf: Figure 2A
GR_merge_genes_promoter=GR_merge[!is.na(GR_merge$genes_promoter)]
GR_merge_genes_promoter$imprinted="Non-Imprinted"
GR_merge_genes_promoter$imprinted[unlist(lapply(GR_merge_genes_promoter$genes_promoter,function(x) any(x%in% Imprinted_Genes$Gene)))]="Imprinted"
GR_merge_genes_df=data.table(regon=paste0(seqnames(GR_merge_genes_promoter),":",start(GR_merge_genes_promoter),'-',end(GR_merge_genes_promoter)),
                                                   dMML=GR_merge_genes_promoter$dMML,imprinted=GR_merge_genes_promoter$imprinted,dMML_pvalue=GR_merge_genes_promoter$dMML_pval,
                             genes=GR_merge_genes_promoter$genes_promoter)
GR_merge_genes_df$expressed_allele=lapply(GR_merge_genes_df$genes,function(x) Imprinted_Genes$Expressed_allele[match(x,Imprinted_Genes$Gene)])
GR_merge_genes_df$genes=unlist(lapply(GR_merge_genes_df$genes,function(x) paste(x,collapse = ";")))
GR_merge_genes_df$expressed_allele=unlist(lapply(GR_merge_genes_df$expressed_allele,function(x) paste(x,collapse = ";")))
write.csv(GR_merge_genes_df[dMML_pvalue<=pval_cutoff&imprinted=="Imprinted"],'../downstream/output/imprinted_regions.csv')
#Getting ecdf
ecdf_df_out=data.frame()
for(type in unique(GR_merge_genes_df$imprinted)){
  GR_merge_genes_promoter_estimate=ecdf(GR_merge_genes_promoter$dMML[GR_merge_genes_promoter$imprinted==type])
  imprinted_uq=seq(0,1,0.001)
  ecdf_df_out=rbind(ecdf_df_out,data.frame(stat_diff=imprinted_uq,quant=GR_merge_genes_promoter_estimate(imprinted_uq),imprinted=type,stringsAsFactors = F))
}
#Getting density ecdf files

pdf('../downstream/output/graphs/Figure2/imprinted_2A.pdf',width=3.5,height=3.5)
print(ggplot(ecdf_df_out,aes(x=dMML,y=quant,group=imprinted,color=imprinted))+geom_line(size=1.5)+xlab('dMML')+
  ylab('cumulative probability')+theme_glob+theme(legend.position = 'bottom'))
 # print(ggplot(GR_merge_genes_df,aes(x=dMML,fill=imprinted,color=imprinted))+geom_density(size=1.5,alpha=0.4)+
 #         xlab('dMML')+theme(legend.position = 'bottom')+theme_glob)
dev.off()

#Getting ecdf for dNME
ecdf_df_out_dNME=data.frame()
for(type in unique(GR_merge_genes_df$imprinted)){
  GR_merge_genes_promoter_estimate=ecdf(GR_merge_genes_promoter$dNME[GR_merge_genes_promoter$imprinted==type])
  imprinted_uq=seq(0,1,0.001)
  ecdf_df_out_dNME=rbind(ecdf_df_out_dNME,data.frame(stat_diff=imprinted_uq,quant=GR_merge_genes_promoter_estimate(imprinted_uq),imprinted=type,stringsAsFactors = F))
}
ecdf_df_out_dNME$stat_type='dNME'
ecdf_df_out$stat_type='dMML'
ecdf_df_out_dNME$imprinted=paste0('dNME-',ecdf_df_out_dNME$imprinted)
ecdf_df_out$imprinted=paste0('dMML-',ecdf_df_out$imprinted)
ecdf_df_out=rbind(ecdf_df_out,ecdf_df_out_dNME)

GR_merge_non_promoter_estimate=ecdf(GR_merge[is.na(GR_merge$genes_promoter)]$dNME)
ecdf_df_out=rbind(rbind(ecdf_df_out,data.frame(dNME=imprinted_uq,quant=GR_merge_non_promoter_estimate(imprinted_uq),imprinted='Non-promoter',stringsAsFactors = F)))
#Getting density ecdf files

pdf('../downstream/output/graphs/Figure2/imprinted_dNME.pdf',width=3.5,height=3.5)
ggplot(ecdf_df_out,aes(x=stat_diff,y=quant,group=imprinted,color=imprinted,alpha=stat_type))+geom_line(size=1.5)+xlab('differential stats')+
        ylab('cumulative probability')+theme_glob+theme(legend.position = 'bottom')+ scale_alpha_manual(name = "stat_type", values = c(1,.25))+
  scale_color_manual(name='imprinted',values=c('red','blue','red','blue'))+
  guides(color=guide_legend(nrow=2,byrow=TRUE),stat_type=guide_legend(nrow=2,byrow=TRUE))
# print(ggplot(GR_merge_genes_df,aes(x=dMML,fill=imprinted,color=imprinted))+geom_density(size=1.5,alpha=0.4)+
#         xlab('dMML')+theme(legend.position = 'bottom')+theme_glob)
dev.off()
ggplot(data=as.data.frame(mcols(GR_merge)),aes(x=dMML,y=dNME))+geom_smooth()+
  geom_density2d(geom = "raster",  aes(fill = after_stat(density)), contour = FALSE,n=200)+
  xlab('dMML')+ylab('dNME')+geom_abline(slope=1,intercept = 0)+xlim(c(0,1))+ylim(c(0,1))


#Bar plot Figure 2B
dMML_OR=rbind(cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',Imprinted_Genes$Gene),
                    data.frame(allele="All Imprinted\ngenes")),
              cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',
           Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Paternal']), data.frame(allele="Paternally\nexpressed")),
            cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',
           Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Maternal']), data.frame(allele="Maternally\nexpressed")))
dMML_OR$stat_type="dMML-ASM"
dMML_OR$star=add.significance.stars(dMML_OR$pvalue, cutoffs = c(0.05, 0.01, 0.001))
dNME_OR=rbind(cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',Imprinted_Genes$Gene),
                    data.frame(allele="All Imprinted genes")),
              cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',
                               Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Paternal']), data.frame(allele="Paternally expressed")),
              cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',
                               Imprinted_Genes$Gene[Imprinted_Genes$`Expressed Allele`=='Maternal']), data.frame(allele="Maternally expressed")))
dNME_OR$stat_type="dNME-ASM"
dNME_OR$star=add.significance.stars(dNME_OR$pvalue, cutoffs = c(0.05, 0.01, 0.001))

OR_barplot=dMML_OR

pdf('../downstream/output/graphs/Figure2/Imprinted_2B.pdf',width=3.5,height=3.5, useDingbats=FALSE)
print(ggplot(OR_barplot,aes(x=allele,y=OR))+
  geom_bar(stat="identity",position=position_dodge(0.9),fill="light blue",width=0.5)+ylim(c(0,105))+
  geom_errorbar(aes(ymin=lowerCI,ymax=upperCI),width=0.2,position=position_dodge(0.9))+ylab("Enrichment")+
  xlab("")+geom_text(aes(label=round(OR,digits = 2),y=upperCI),position=position_dodge(0.9),vjust=-0.5,color="black", size=3)+
  geom_text(aes(y=upperCI,label=star),vjust=-1,position = position_dodge(0.9))+
  theme_glob)
dev.off()
# Finding the overlap between monoallelic expressed gene ------------------
MAE_BAE_data_Gimelbrant <- as.data.frame(read_excel("../downstream/input/human_analysis/imprinting_ASE/MAE_BAE_data_Gimelbrant.xlsx"),stringsAsFactors=F)
#Enrichment in MAE dMML , not quiet enriched
MAE=MAE_BAE_data_Gimelbrant$Gene[ MAE_BAE_data_Gimelbrant$`MAE=1_BAE=0`==1]
#dMML enrichment promoter
write.csv(unique(unlist(GR_merge[GR_merge$dMML_pval<=pval_cutoff]$genes_promoter[GR_merge[GR_merge$dMML_pval<=pval_cutoff]$genes_promoter%in% MAE_BAE_data_Gimelbrant$Gene])),
          '../downstream/output/human_analysis/dMML_analysis/MAE_genes.csv')
dMML_MAE=MAE_enrich(GR_merge[GR_merge$genes_promoter%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_promoter',stat='dMML_pval',MAE=MAE)
dNME_MAE=MAE_enrich(GR_merge[GR_merge$genes_promoter%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_promoter',stat='dNME_pval',MAE=MAE)
cat("MAE enriched in dMML:", dMML_MAE$OR,'with pvalue=',dMML_MAE$pvalue,'\n')
cat("MAE enriched in dNME:", dNME_MAE$OR,'with pvalue=',dNME_MAE$pvalue,'\n')
#Figure 2C chromHMM enrichment
# ChromHMM annotations ----------------------------------------------------
ah = AnnotationHub()
ENCODE_name=ENCODE_to_sample(unique(GR_merge$Sample))
chromHMM_dMML_all_ls=list()
ah_gr=GRanges()
#Do it for all available data, check
suppressMessages({for (sp in ENCODE_name$sample[!is.na(ENCODE_name$ENCODE)]){
ah_num=names(AnnotationHub::query(ah, c("chromhmmSegmentations", ENCODE_name$ENCODE[ENCODE_name$sample==sp])))
  chromHMM=ah[[ah_num]]
  chromHMM_dMML_all_ls[[sp]]=chromHMM_OR(GR_merge, chromHMM,sp,stat="dMML_pval")
  # }
}
})

#####Summing up cont table
chromHMM_dMML_all=chromHMM_combine(chromHMM_dMML_all_ls)

#Plot OR with error bar
chromHMM_dMML_all=chromHMM_dMML_all[order(chromHMM_dMML_all$OR,decreasing=F),]
chromHMM_dMML_all$states=factor(chromHMM_dMML_all$states,levels=chromHMM_dMML_all$states)
chromHMM_dMML_all$FDR=p.adjust(chromHMM_dMML_all$p_value,method='BH')
chromHMM_dMML_all$star=add.significance.stars(chromHMM_dMML_all$FDR, cutoffs = c(0.05, 0.01, 0.001))
chromHMM_dMML_all$states=factor(chromHMM_dMML_all$states,levels=chromHMM_dMML_all$states[order(chromHMM_dMML_all$OR,decreasing=F)])
pdf('../downstream/output/graphs/Figure2/Figure2C_crhomHMM.pdf',width=5,height=5)
print(ggplot(chromHMM_dMML_all,aes(x=states,y=OR,fill=states))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+ coord_flip()+
  theme_glob+xlab("chromHMM states")+theme(legend.position = "")+ylim(c(0,13))+
  geom_text(aes(label=round(OR,digits = 2),y=upper_CI),hjust=-0.5, color="black", size=3)+
  geom_text(aes(label=star,y=upper_CI+0.5),hjust=-1, color="black", size=3))
dev.off()

#Figure 2D chromatin accessibility

# Processing bulk ATAC-seq result --------------------------------------------

######reading in data
GR_merge_H1=GR_merge[GR_merge$Subject=='H1']
#This file from Ken
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
#This file from Ken
ATAC_H1_CPM=as.data.frame(read.table("../downstream/input/H1_ATAC_allele_CPM.txt",header=T,stringsAsFactors = F))
colnames(ATAC_H1_CPM)=c("seqname","start","end","rep1_g1_CPM","rep1_g2_CPM","rep2_g1_CPM","rep2_g2_CPM")
ATAC_H1_CPM=makeGRangesFromDataFrame(as.data.frame(ATAC_H1_CPM),keep.extra.columns = T)
elementMetadata(ATAC_H1_CPM)=cbind(res_ATAC_lfc,elementMetadata(ATAC_H1_CPM))
#############calcualte log2FC using CPM and pseudocount
ATAC_H1_CPM$log2FC_CPM=log2((ATAC_H1_CPM$rep1_g2_CPM+ATAC_H1_CPM$rep2_g2_CPM+1)/(ATAC_H1_CPM$rep1_g1_CPM+ATAC_H1_CPM$rep2_g1_CPM+1))

ATAC_H1_CPM$genome2_CPM=(ATAC_H1_CPM$rep1_g2_CPM+ATAC_H1_CPM$rep2_g2_CPM)/2
ATAC_H1_CPM$genome1_CPM=(ATAC_H1_CPM$rep1_g1_CPM+ATAC_H1_CPM$rep2_g1_CPM)/2
############Putting data into GR merge object
olap_ATAC=findOverlaps(GR_merge_H1,ATAC_H1_CPM,select="all",maxgap =500)
#####Check overlapped regions
sum(!is.na(olap_ATAC))#12696
qt_dt=cbind(data.table(qt=queryHits(olap_ATAC)),as.data.table(mcols(ATAC_H1_CPM[subjectHits(olap_ATAC)])))
qt_dt=qt_dt[,lapply(.SD,mean),by=qt]#.SD: subset of Data, subset each column as list
GR_merge_H1_ATAC=GR_merge_H1[qt_dt$qt]
mcols(GR_merge_H1_ATAC)=cbind(mcols(GR_merge_H1_ATAC),qt_dt[,-1])
GR_merge_H1_ATAC$dMML_relative=GR_merge_H1_ATAC$MML2-GR_merge_H1_ATAC$MML1
GR_merge_H1_ATAC$dNME_relative=GR_merge_H1_ATAC$NME2-GR_merge_H1_ATAC$NME1
####test correlation
cor_dMML=cor.test(GR_merge_H1_ATAC$dMML_relative[GR_merge_H1_ATAC$dMML_pval<=pval_cutoff],
         GR_merge_H1_ATAC$log2FC_CPM[GR_merge_H1_ATAC$dMML_pval<=pval_cutoff])#-0.6302
cor_dNME=cor.test(GR_merge_H1_ATAC$dNME_relative[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff],
         GR_merge_H1_ATAC$log2FC_CPM[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff])
ATAC_plot=data.frame(dMML=GR_merge_H1_ATAC$dMML_relative,ATAC_FC_CPM=GR_merge_H1_ATAC$log2FC_CPM,dMML_pval=GR_merge_H1_ATAC$dMML_pval,
                     dNME=GR_merge_H1_ATAC$dNME_relative,dNME_pval=GR_merge_H1_ATAC$dNME_pval)
pdf('../downstream/output/graphs/Figure2/Figure2D_ATAC-seq_dMML.pdf',width=3.5,height=3.5)
print(ggplot(ATAC_plot[ATAC_plot$dMML_pval<=pval_cutoff,],aes(x=dMML,y=ATAC_FC_CPM))+geom_smooth()+geom_point(alpha=0.5)+
  ylab('Allelic accessibility change')+theme_glob)+xlab('relative dMML')
dev.off()
#Not in use Figure S3
# pdf('../downstream/output/graphs/FigureS3/Figure2D_ATAC-seq_dNME.pdf',width=3.5,height=3.5)
# print(ggplot(ATAC_plot[ATAC_plot$dNME_pval<=pval_cutoff,],aes(x=dNME,y=ATAC_FC_CPM))+geom_smooth()+geom_point(alpha=0.5)+
#   ylab('Allelic accessibility change')+theme_glob)+xlab('relative dNME')
# dev.off()
