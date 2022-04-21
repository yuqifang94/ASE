source('mainFunctions_sub.R')
NME_all=readRDS(NME_agnostic_all_file)
GR_merge=readRDS(GR_merge_file)
CpG_hg19=getCpgSitesH19()#26752702
subsetByOverlaps(CpG_hg19,GR_merge[GR_merge$dMML_pval<=0.1])#25105
subsetByOverlaps(CpG_hg19,NME_all)#21234156
MML_all=readRDS(MML_agnostic_all_file)
NME_all_dt=convert_GR(NME_all,direction="DT")
MML_all_dt=convert_GR(MML_all,direction="DT")
# Plotting NME vs MML ------------------------------------

merged_dt=rbind(MML=MML_all_dt$MML,NME=NME_all_dt$NME)
#Aggregate NME, using quantiles, 0.05 and 0.95
merged_dt_agg=merged_dt[, list(NME=round(median(NME),digits=digits_round),
                                   Bottom25=round(quantile(NME,probs=0.25),digits=digits_round),
                                   top25=round(quantile(NME,probs=0.75),digits=digits_round)), 
                            by = list(MML = round(MML,digits=digits_round))]

merged_dt_agg$Bottom25= predict(loess(Bottom25~MML,merged_dt_agg),newdata=merged_dt_agg$MML)
merged_dt_agg$top25= predict(loess(top25~MML,merged_dt_agg),newdata=merged_dt_agg$MML)
pdf('../downstream/output/human_analysis/QC/NME_vs_MML_with_quantile_agnostic.pdf',width=5,height=5)
#Plotting
print(ggplot(merged_dt_agg,aes(x=MML, y=NME))+
        xlim(c(0,1))+ylim(c(0,1))+ggtitle("MML and NME relationship")+geom_smooth(method="loess",se=FALSE)+
        ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)+theme_glob)
dev.off()