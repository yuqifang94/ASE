rm(list=ls())
source("mainFunctions_sub.R")
#Define ggplot theme

theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))

# Check SNP frequency -----------------------------------------------------
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
variant_HetCpG_meta=as.data.table(mcols(variant_HetCpG_meta))
variant_HetCpG_meta$SNP=apply(variant_HetCpG_meta[,list(REF,ALT)],1,function(x) paste(x,collapse = '>'))
variant_HetCpG_meta$tri_SNP=paste0(variant_HetCpG_meta$REF_tri,'>',variant_HetCpG_meta$ALT_tri)
#Lumping all
variant_HetCpG_meta$CpG_change="No change"
variant_HetCpG_meta[(grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri ))]$CpG_change='Lose CG'
variant_HetCpG_meta[(!grepl('CG',REF_tri )) & (grepl('CG',ALT_tri ))]$CpG_change='Gain CG'
variant_HetCpG_meta[(grepl('CG',REF_tri )) & (grepl('CG',ALT_tri))]$CpG_change='No change'
#substr(variant_HetCpG_meta$tri_SNP,2,2)="X"
SNP_all=list()
SNP_het=list()
SNP_box=list()
color_theme=c(rainbow(length(unique(variant_HetCpG_meta$SNP))))
variant_SNP_tri=data.table()
names(color_theme)=unique(variant_HetCpG_meta$SNP)
for (sn in unique(variant_HetCpG_meta$SNP)){
  for(tri in unique(variant_HetCpG_meta[SNP==sn]$tri_SNP)){
    variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta[SNP==sn &dNME_pval<=pval_cutoff],tri,"tri_SNP",pval_cutoff)
    variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta[SNP==sn & tri_SNP==tri]$CpG_change))
    variant_SNP_tri=rbind(variant_SNP_tri,variant_SNP_tri_OR)
  }
  #get HetCpG
  variant_SNP_tri$CpG_change=factor(variant_SNP_tri$CpG_change,levels=c('Gain CG','No change','Lose CG'))
  variant_SNP_tri=variant_SNP_tri[order(OR,decreasing=T)]
  variant_SNP_tri$SNP=factor(variant_SNP_tri$SNP,levels = variant_SNP_tri$SNP)
  
  variant_SNP_tri$FDR=p.adjust(variant_SNP_tri$pvalue,method='BH')
  variant_SNP_tri$significant=add.significance.stars(variant_SNP_tri$FDR, cutoffs = c(0.05, 0.01, 0.001))
  SNP_all[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR))+geom_bar(fill=color_theme[sn],stat="identity")+ylab('OR (REF > ALT)')+
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ggtitle(sn)+xlab("")+
  theme_glob+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0,max(variant_SNP_tri$upperCI)*1.1))+geom_text(aes(label=significant,y=upperCI*1),vjust = -0.5)
 
   SNP_het[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR,fill=CpG_change))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ggtitle(sn)+ylim(c(0,max(variant_SNP_tri$upperCI)*1.1))+
    theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
     scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+geom_text(aes(label=significant,y=upperCI*1),vjust = -0.5)
   
   SNP_box[[sn]]=ggplot(variant_HetCpG_meta[SNP==sn & dNME_pval<=pval_cutoff],aes(x=reorder(tri_SNP, -NME_relative, FUN = median) ,y=NME_relative,fill=CpG_change))+geom_boxplot()+
     ylab('refNME-altNME')+theme_glob+theme(legend.position = "bottom")+ ylim(c(-1,1))+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab('')
   
  variant_SNP_tri=data.table()
}
# pdf('../downstream/output/graphs/Figure3/variant_OR_tri2.pdf',width=14,height=14)
# SNP_all=SNP_all[c("C>G", names(SNP_all)[names(SNP_all)!="C>G"])]
# ggarrange(plotlist=SNP_all, nrow=4,ncol=3)
# dev.off()

pdf('../downstream/output/graphs/Figure2/Figure-S3-variant_OR_tri3.pdf',width=14,height=14)
SNP_het=SNP_het[c("C>G", names(SNP_het)[names(SNP_het)!="C>G"])]
ggarrange(plotlist=SNP_het, nrow=4,ncol=3,common.legend = T,legend="bottom")
dev.off()
# pdf('../downstream/output/graphs/Figure3/variant_box_tri.pdf',width=14,height=14)
# SNP_box=SNP_box[c("C>G", names(SNP_box)[names(SNP_box)!="C>G"])]
# ggarrange(plotlist=SNP_box, nrow=4,ncol=3,common.legend = T,legend="bottom")
# dev.off()

OR_all_SNP_change=data.table()

for(sn in unique(variant_HetCpG_meta$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc(variant_HetCpG_meta[dNME_pval<=pval_cutoff],sn,"CpG_change"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=T)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")
dev.off()
#Density analysis
GR_merge=readRDS(GR_merge_file)

GR_merge$CpGdiff=GR_merge$g1CG-GR_merge$g2CG

#plot_CpG_number(GR_merge[GR_merge$dNME_pval<=pval_cutoff])

#CpG density vs dNME,here uses an extended density 
GR_merge$dNME_relative=GR_merge$NME1-GR_merge$NME2
GR_merge_dt=as.data.table(mcols(GR_merge))
#abs difference
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
#ratio
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1/CG_allele_extend_g2)]
cor.test(GR_merge_dt[dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$dNME_relative, 
         GR_merge_dt[dNME_pval<=pval_cutoff&dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$density_diff)


GR_merge_dt$CpG_stat="No difference"
GR_merge_dt[CpGdiff!=0]$CpG_stat="With CpG difference"
GR_merge_dt$CpG_stat=factor(GR_merge_dt$CpG_stat,levels = c("With CpG difference","No difference"))
#Make this one better
#Figure 3A
GR_merge_dt$dNME_relative_more_less=GR_merge_dt$dNME_relative
GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative_more_less=GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative*sign(GR_merge_dt[GR_merge_dt$CpGdiff!=0]$CpGdiff)
t.test(GR_merge_dt[CpGdiff!=0&dNME_pval<=pval_cutoff]$dNME_relative_more_less,alternative="less")
pdf('../downstream/output/graphs/Figure3/Figure3A_CpG_number_NME.pdf',width=7,height=7)
ggplot(GR_merge_dt[dNME_pval<=pval_cutoff],aes(y=dNME_relative_more_less,x=CpG_stat,fill=CpG_stat))+
  geom_violin()+xlab("")+
  theme_glob+ylab('relative dNME')+theme(legend.position = "none")

dev.off()
pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dNME_ratio.pdf',width=7,height=7)
ggplot(GR_merge_dt[dNME_pval<=pval_cutoff&density_diff!=0],aes(x=as.factor(round(density_diff,digits = 1)),y=dNME_relative))+geom_violin(fill='light blue')+
  xlab("CpG density ratio")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dNME_quantile.pdf',width=7,height=7)
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
GR_merge_dt_S4=GR_merge_dt[dNME_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_S4$density_diff=findInterval(GR_merge_dt_S4$density_diff,quantile(GR_merge_dt_S4$density_diff,prob=seq(0,0.9,0.1)))
GR_merge_dt_S4$density_diff=factor(paste0(GR_merge_dt_S4$density_diff*10,'%'),levels=paste0( seq(0,100,10),"%"))
ggplot(GR_merge_dt_S4,aes(x=density_diff,y=dNME_relative))+geom_violin(fill='light blue')+
  xlab("CpG density quantile")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
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


# Figure 3B: allele-agnotic density ---------------------------------------

NME_in=readRDS(NME_agnostic_file)
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
NME_in$CG_hg19=countOverlaps(NME_in,CpG_hg19)
NME_in_gr=unique(granges(NME_in))
gr_seq=getSeq(Hsapiens,NME_in_gr,as.character=T)
NME_in_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
NME_in_olap=findOverlaps(NME_in,NME_in_gr,type='equal')
NME_in$CGcont_exp[queryHits(NME_in_olap)]=NME_in_gr$CGcont_exp[subjectHits(NME_in_olap)]
NME_in$density=NME_in$CG_hg19/NME_in$CGcont_exp
NME_in=NME_in[seqnames(NME_in)%in%paste0("chr",1:22)]
cor.test(NME_in$density,NME_in$NME,method='pearson')
#Make boxplot of this
#quantile(NME_in$density,prob=seq(0,1,0.2),na.rm=T)
NME_in$density_quant=findInterval(NME_in$density,seq(0,1,0.1))
#NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
NME_in$density_quant=factor(quant_conv[NME_in$density_quant],levels=quant_conv)
pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot_CG_exp.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in)),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()   
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



pdf('../downstream/output/graphs/FigureS6/NME_density_feature.pdf',width=3.5,height=3.5)
# my_comparisons <- list( c("islands", "shores"), c("islands", "shelf"), c("islands", "open sea"),
#                         c("shores", "shelf"),c("shores", "open sea"),c("shelf", "open sea"))
ggplot(CpG_density_NME,aes(x=feature,y=NME))+geom_boxplot(outlier.shape = NA)+ 
  theme_glob+xlab("Genomic Features")
#stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()

# #same analysis on dMML-ASM
# MML_in=readRDS(MML_agnostic_file)
# MML_in$CG_hg19=countOverlaps(MML_in,CpG_hg19)
# MML_in$density=MML_in$CG_hg19/width(MML_in)
# MML_in=MML_in[seqnames(MML_in)%in%paste0("chr",1:22)]
# MML_in$density_quant=findInterval(MML_in$density,seq(0,0.1,0.01))
# MML_in$density_quant=factor(quant_conv[MML_in$density_quant],levels=quant_conv)
# olap=findOverlaps(MML_in,genomic_features$`CpG island`)
# cor.test(MML_in$density,MML_in$MML,method='pearson')#-0.434
# cor.test(MML_in[queryHits(olap)]$density,MML_in[queryHits(olap)]$MML,method='pearson')#-
# cor.test(MML_in[-queryHits(olap)]$density,MML_in[-queryHits(olap)]$MML,method='pearson')#-
# pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
# ggplot(as.data.frame(mcols(MML_in)),aes(x=density_quant, y=MML))+
#   ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
#   ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()      
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

