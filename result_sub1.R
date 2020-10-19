# Genomics
# Source main functions
#setwd("~/code/HASM-MetaAnalysis/")
rm(list=ls())
source("mainFunctions_sub.R")
#Define ggplot theme

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
# Find number of overlapped regions ---------------------------------------
GR_merge=readRDS(GR_merge_file)
#Only use merged data for H1
GR_merge=GR_merge[!(GR_merge$Sample%in%c("rep1 - H1","rep2 - H1"))]
dMML=sum(GR_merge$dMML_pval<=pval_cutoff)
dNME=sum(GR_merge$dNME_pval<=pval_cutoff)
dMML_dNME=sum(GR_merge$dNME_pval<=pval_cutoff&GR_merge$dMML_pval<=pval_cutoff)
cat('Number of dMML:',dMML,'\n')
cat('Number of dNME:',dNME,'\n')
cat('Number of dNME and dMML:',dMML_dNME,'\n')
# #Run it tonight
# olap_merge=c()
# for(i in 1:50000){olap_merge=c(olap_merge,length(subsetByOverlaps(GR_merge[sample(1:length(GR_merge),dMML,replace=F)],
#                                                                   GR_merge[sample(1:length(GR_merge),dNME,replace=F)],type='equal')))}
# sum(olap_merge<=dMML_dNME)/length(olap_merge)#=0
# saveRDS(olap_merge,'../downstream/output/olap_merge.rds')
#Examples:
#dNME example1:
subsetByOverlaps(GR_merge[GR_merge$Sample=="ectoderm_paired - HUES64"],GRanges(seqnames="chr14",IRanges(start=104552150,end=104552495)))
#dNM example2:
subsetByOverlaps(GR_merge[GR_merge$Sample=="Psoas_Muscle_single - STL003"],GRanges(seqnames="chr17",IRanges(start=33750066,end=33770266)))
#dMML_example1
subsetByOverlaps(GR_merge[GR_merge$Sample=="stem_27_undifferentiated_paired - HUES64"],GRanges(seqnames="chr11",IRanges(start=2720817,end=2721033)))
#dMML_example2
subsetByOverlaps(GR_merge[GR_merge$Sample=="endoerm_27_paired - HUES64"],GRanges(seqnames="chr20",IRanges(start=32308087,end=32308287)))
#UC vs dNME or dMML
plot_dt=rbind(
  data.table(value=GR_merge$dNME[GR_merge$UC>=0.5],statistic="dNME"),
  data.table(value=GR_merge$dMML[GR_merge$UC>=0.5],statistic="dMML"))

ggplot(plot_dt,aes(x=value,fill=statistic))+geom_density(alpha=0.5)
plot_dt=rbind(
  data.table(UC=GR_merge$UC[GR_merge$dNME>=0.5],statistic="high dNME"),
  data.table(UC=GR_merge$UC[GR_merge$dMML>=0.5],statistic="high dMML"))

ggplot(plot_dt,aes(x=value,fill=statistic))+geom_density(alpha=0.5)

#Figure S1
digits_round=2
#Fig 1A:dNME vs dMML
GR_merge_dt=data.table(dMML=GR_merge$dMML,dNME=GR_merge$dNME,dNME_pval=GR_merge$dNME_pval,dMML_pval=GR_merge$dMML_pval)
GR_merge_dt=GR_merge_dt[GR_merge_dt$dNME_pval<=pval_cutoff|GR_merge_dt$dMML_pval<=pval_cutoff]
GR_merge_dt_agg=GR_merge_dt[, list(dNME=round(median(dNME),digits=digits_round),
                                   Bottom25=round(quantile(dNME,probs=0.25),digits=digits_round),
                                   top25=round(quantile(dNME,probs=0.75),digits=digits_round)), 
                            by = list(dMML = round(dMML,digits=digits_round))]

GR_merge_dt_agg$Bottom25= predict(loess(Bottom25~dMML,GR_merge_dt_agg),newdata=GR_merge_dt_agg$dMML)
GR_merge_dt_agg$top25= predict(loess(top25~dMML,GR_merge_dt_agg),newdata=GR_merge_dt_agg$dMML)
###plotting
pdf('../downstream/output/graphs/FigureS1/dNME_vs_dMML_quantile_differential_median.pdf',width=5,height=5)
print(ggplot(GR_merge_dt_agg,aes(x=dMML, y=dNME))+
  xlim(c(0,1))+ylim(c(0,0.7))+ggtitle("dMML and dNME relationship")+geom_smooth(method="loess",se=FALSE)+
  ylab("dNME")+theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)+theme_glob+
  scale_linetype_manual(values=c("solid","twodash", "twodash"))+scale_color_manual(values=c("Blue","Blue","Blue")))
dev.off()

# Plotting dNME vs dMML and NME vs MML ------------------------------------
#Figure S1B: MML and NME
GR_merge_dt=rbind(data.table(MML=GR_merge$MML1,NME=GR_merge$NME1,dNME_pval=GR_merge$dNME_pval,dMML_pval=GR_merge$dMML_pval),
                  data.table(MML=GR_merge$MML2,NME=GR_merge$NME2,dNME_pval=GR_merge$dNME_pval,dMML_pval=GR_merge$dMML_pval))
#Aggregate NME, using quantiles, 0.05 and 0.95
GR_merge_dt_agg=GR_merge_dt[, list(NME=round(median(NME),digits=digits_round),
                                   Bottom25=round(quantile(NME,probs=0.25),digits=digits_round),
                                   top25=round(quantile(NME,probs=0.75),digits=digits_round)), 
                            by = list(MML = round(MML,digits=digits_round))]

GR_merge_dt_agg$Bottom25= predict(loess(Bottom25~MML,GR_merge_dt_agg),newdata=GR_merge_dt_agg$MML)
GR_merge_dt_agg$top25= predict(loess(top25~MML,GR_merge_dt_agg),newdata=GR_merge_dt_agg$MML)
pdf('../downstream/output/graphs/FigureS1/NME_vs_MML_with_quantile.pdf',width=5,height=5)
#Plotting
print(ggplot(GR_merge_dt_agg,aes(x=MML, y=NME))+
  xlim(c(0,1))+ylim(c(0,1))+ggtitle("MML and NME relationship")+geom_smooth(method="loess",se=FALSE)+
  ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)+theme_glob)
dev.off()

