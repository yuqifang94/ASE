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
# Find overlap by random change -------------------------------------------
rm(list=ls())
source("mainFunctions_sub.R")
GR_merge=readRDS(GR_merge_file)
GR_merge_dt=convert_GR(GR_merge,direction='DT')
library(mcreplicate)
set.seed(123)
olap_perm_test<-function(region,dMML_pval,dNME_pval,total_perm=5000,pval_cutoff=0.1){
    dMML_sig=sum(dMML_pval<=pval_cutoff)
    dNME_sig=sum(dNME_pval<=pval_cutoff)
    dNME_dMML_olap=sum(dNME_pval<=pval_cutoff&dMML_pval<=pval_cutoff)
    length_all=length(region)
    return(sum(mc_replicate(total_perm,length(intersect(region[sample(1:length_all,dMML_sig)],region[sample(1:length_all,dNME_sig)])),mc.cores=24)>=dNME_dMML_olap))
}
GR_merge_dt_perm=GR_merge_dt[,list(perm_larger=olap_perm_test(region,dMML_pval,dNME_pval,total_perm=10000)),by=list(Sample)]
GR_merge_dt_perm_all=olap_perm_test(GR_merge_dt$region,GR_merge_dt$dMML_pval,GR_merge_dt$dNME_pval,total_perm=10000)
saveRDS(GR_merge_dt_perm,'../downstream/output/human_analysis/CPEL_outputs/GR_merge_dt_perm.rds')
saveRDS(GR_merge_dt_perm_all,'../downstream/output/human_analysis/CPEL_outputs/GR_merge_dt_perm_all.rds')