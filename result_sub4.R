rm(list=ls())
source("mainFunctions_sub.R")
#Define ggplot theme

theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))
plot_CpG_number<-function(allele_gr_ASM,stat="NME",theme_in=theme_glob){
  allele_het_df=rbind(data.frame(value=c(elementMetadata(allele_gr_ASM)[,paste0(stat,'1')][allele_gr_ASM$CpGdiff>0],
                                         elementMetadata(allele_gr_ASM)[,paste0(stat,'2')][allele_gr_ASM$CpGdiff<0]),type='More CpG'),
                      data.frame(value=c(elementMetadata(allele_gr_ASM)[,paste0(stat,'2')][allele_gr_ASM$CpGdiff>0],
                                         elementMetadata(allele_gr_ASM)[,paste0(stat,'1')][allele_gr_ASM$CpGdiff<0]),type='less CpG'))
  print(ggplot(allele_het_df,aes(x=value,color=type,fill=type))+
          geom_density(alpha=0.25,size=1)+xlab(stat)+
          theme_in+theme(legend.position="bottom",legend.title = element_blank())+ylim(c(0,4))+
          scale_color_manual(values=c("blue","red"))+scale_fill_manual(values=c("blue","red"))+
          geom_hline(yintercept=0, colour="white", size=1))
}

#Density analysis
GR_merge=readRDS(GR_merge_file)
#Only use merged data for H1

GR_merge$CpGdiff=GR_merge$g1CG-GR_merge$g2CG
#Figure 3C * check
pdf('../downstream/output/graphs/Figure3/Figure3C_CpG_number_NME.pdf',width=3.5,height=3.5)
plot_CpG_number(GR_merge[GR_merge$dNME_pval<=pval_cutoff])
dev.off()
#CpG density vs dNME 
GR_merge$dNME_relative=GR_merge$NME1-GR_merge$NME2
cor.test(GR_merge$dNME_relative[GR_merge$dNME_pval<=pval_cutoff], GR_merge$density_diff[GR_merge$dNME_pval<=pval_cutoff])
#Figure 3D 
NME_in=readRDS(NME_agnostic_file)
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
NME_in$CG_hg19=countOverlaps(NME_in,CpG_hg19)
NME_in$density=NME_in$CG_hg19/width(NME_in)
NME_in=NME_in[seqnames(NME_in)%in%paste0("chr",1:22)]
cor.test(NME_in$density,NME_in$NME,method='pearson')#or pearson,0.15
#Make boxplot of this
#quantile(NME_in$density,prob=seq(0,1,0.2),na.rm=T)
NME_in$density_quant=findInterval(NME_in$density,seq(0,0.1,0.01))
#NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
quant_conv=c(paste0(seq(0,0.09,0.01),'-',seq(0.01,0.1,0.01)),'>0.1')
NME_in$density_quant=factor(quant_conv[NME_in$density_quant],levels=quant_conv)
pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in)),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()              


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
genomic_features=readRDS(genomic_features_file)
#Figure S5
olap_islands=findOverlaps(NME_in,genomic_features$`CpG island`)
olap_shores=findOverlaps(NME_in,genomic_features$`CpG shore`)
olap_shelf=findOverlaps(NME_in,genomic_features$`CpG shelf`)
olap_open_sea=findOverlaps(NME_in,genomic_features$`CpG open sea`)

CpG_density_NME=rbind(data.table(NME=NME_in$NME[queryHits(olap_islands)],feature='islands'),
                      data.table(NME=NME_in$NME[queryHits(olap_shores)],feature='shores'),
                      data.table(NME=NME_in$NME[queryHits(olap_shelf)],feature='shelf'),
                      data.table(NME=NME_in$NME[queryHits(olap_open_sea)],feature='open sea'))



pdf('../downstream/output/graphs/FigureS5/NME_density_feature.pdf',width=3.5,height=3.5)
# my_comparisons <- list( c("islands", "shores"), c("islands", "shelf"), c("islands", "open sea"),
#                         c("shores", "shelf"),c("shores", "open sea"),c("shelf", "open sea"))
ggplot(CpG_density_NME,aes(x=feature,y=NME))+geom_boxplot(outlier.shape = NA)+ 
  theme_glob+xlab("Genomic Features")
  #stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()

