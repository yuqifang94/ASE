source("mainFunctions_sub.R")
GR_merge=readRDS(GR_merge_file)
digits_round=2
# Plotting NME vs MML ------------------------------------

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
# Plotting dNME vs dMML ------------------------------------
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
