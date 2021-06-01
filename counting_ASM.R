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
dMML=sum(GR_merge$dMML_pval<=pval_cutoff)
dNME=sum(GR_merge$dNME_pval<=pval_cutoff)
dMML_dNME=sum(GR_merge$dNME_pval<=pval_cutoff&GR_merge$dMML_pval<=pval_cutoff)
cat('Number of regions:',length(GR_merge),'\n')
cat('Number of dMML:',dMML,'\n')
cat('Number of dNME:',dNME,'\n')
cat('Number of dNME and dMML:',dMML_dNME,'\n')
#Run it tonight
olap_merge=c()
for(i in 1:50000){olap_merge=c(olap_merge,length(subsetByOverlaps(GR_merge[sample(1:length(GR_merge),dMML,replace=F)],
                                                                  GR_merge[sample(1:length(GR_merge),dNME,replace=F)],type='equal')))}
sum(olap_merge<=dMML_dNME)/length(olap_merge)#=0
saveRDS(olap_merge,'../downstream/output/olap_merge.rds')
