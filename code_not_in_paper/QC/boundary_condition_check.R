source('mainFunctions_sub.R')
# HUES64_mesoderm_23_paired_phased_tmml_pvals.bedGraph CPEL_version2 4.2M, Apr 28 (filter 5, no boudnary), 98230
# HUES64_mesoderm_23_paired_phased_tmml_pvals.bedGraph newrun3 9.4M, May 7 (filter 5, boundary), 171985
#Contains boundary condition
GR_contain_boundary=import.subject('../downstream/data/ASM_run3/bedGraph_diff/')#2,221,880 unique, total 28164174
GR_contain_boundary=convert_GR(GR_contain_boundary,direction="DT")
GR_no_boundary=readRDS(GR_file) #3,332,744 unique, total 16601685
GR_no_boundary = convert_GR(GR_no_boundary,direction="DT")
GR_no_boundary = GR_no_boundary[!Sample%in%c("rep1 - H1","rep2 - H1")]
gff_in=readRDS(gff_in_file)
gff_in=convert_GR(gff_in,direction="DT")
GR_contain_boundary$NCpG = gff_in[match(GR_contain_boundary$region,region)]$N
GR_no_boundary$NCpG = gff_in[match(GR_no_boundary$region,region)]$N
sum(GR_contain_boundary[NCpG>=2&Statistic=="dMML" & Sample %in% GR_no_boundary$Sample]$pvalue<=0.1)#9696
sum(GR_contain_boundary[NCpG>=2&Statistic=="dNME"& Sample %in% GR_no_boundary$Sample]$pvalue<=0.1)#10497
sum(GR_no_boundary[NCpG>=2&Statistic=="dMML"]$pvalue<=0.1)#6803
sum(GR_no_boundary[NCpG>=2&Statistic=="dNME"]$pvalue<=0.1)#29826
# #Get min dNME or dMML at pval<0.1 with each N
# GR_no_boundary_summary=GR_no_boundary[pvalue<=0.1 & NCpG>1,list(minScore = min(score)),by=list(Statistic,NCpG,Sample)][order(NCpG)]
# #GR_no_boundary_summary=GR_no_boundary_summary[,list(minScore=mean(minScore)),by=list(Statistic,NCpG)]
# GR_no_boundary_summary$CPEL = "Without reads at boundary condition"
# GR_contain_boundary_summary=GR_contain_boundary[pvalue<=0.1& NCpG>1,list(minScore = min(score)),by=list(Statistic,NCpG,Sample)][order(NCpG)]
# #GR_contain_boundary_summary=GR_contain_boundary_summary[,list(minScore=mean(minScore)),by=list(Statistic,NCpG)]
# GR_contain_boundary_summary$CPEL = "With reads at boundary condition"
# dMML_comparison = rbind(GR_no_boundary_summary[Statistic=="dMML"],GR_contain_boundary_summary[Statistic=="dMML"])
# dNME_comparison = rbind(GR_no_boundary_summary[Statistic=="dNME"],GR_contain_boundary_summary[Statistic=="dNME"])
# pdf("../downstream/output/human_analysis/QC/dMML_cutoff_boundary.pdf",width = 6, height=5)
# print(ggplot(dMML_comparison,aes(x = NCpG,y=minScore,color=CPEL))+geom_smooth()+ylim(c(0,1))+
#     ggtitle("dMML")+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title = element_blank()))
# dev.off()
# pdf("../downstream/output/human_analysis/QC/dNME_cutoff_boundary.pdf",width = 6, height=5)
# print(ggplot(dNME_comparison,aes(x = NCpG,y=minScore,color=CPEL))+geom_smooth()+ylim(c(0,1))+
#     ggtitle("dNME")+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title = element_blank()))
# dev.off()

#Get min dNME or dMML at pval<0.1 with each N
GR_no_boundary_summary=GR_no_boundary[pvalue<=0.1 & NCpG>1,list(minScore = min(score)),by=list(Statistic,NCpG)][order(NCpG)]
#GR_no_boundary_summary=GR_no_boundary_summary[,list(minScore=mean(minScore)),by=list(Statistic,NCpG)]
GR_no_boundary_summary$CPEL = "Without reads at boundary condition"
GR_contain_boundary_summary=GR_contain_boundary[pvalue<=0.1& NCpG>1,list(minScore = min(score)),by=list(Statistic,NCpG)][order(NCpG)]
#GR_contain_boundary_summary=GR_contain_boundary_summary[,list(minScore=mean(minScore)),by=list(Statistic,NCpG)]
GR_contain_boundary_summary$CPEL = "With reads at boundary condition"
dMML_comparison = rbind(GR_no_boundary_summary[Statistic=="dMML"],GR_contain_boundary_summary[Statistic=="dMML"])
dNME_comparison = rbind(GR_no_boundary_summary[Statistic=="dNME"],GR_contain_boundary_summary[Statistic=="dNME"])
pdf("../downstream/output/human_analysis/QC/dMML_cutoff_boundary.pdf",width = 6, height=5)
print(ggplot(dMML_comparison,aes(x = NCpG,y=minScore,fill=CPEL))+geom_bar(stat="identity", position=position_dodge())+ylim(c(0,1))+
    ggtitle("dMML")+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title = element_blank()))
dev.off()
pdf("../downstream/output/human_analysis/QC/dNME_cutoff_boundary.pdf",width = 6, height=5)
print(ggplot(dNME_comparison,aes(x = NCpG,y=minScore,fill=CPEL))+geom_bar(stat="identity", position=position_dodge())+ylim(c(0,1))+
    ggtitle("dNME")+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title = element_blank()))
dev.off()