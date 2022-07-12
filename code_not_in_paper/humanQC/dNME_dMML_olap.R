source('mainFunctions_sub.R')
GR_merge=readRDS(GR_merge_file)
merged_dt=readRDS("../downstream/output/human_analysis/CPEL_outputs/NME_MML_merged_dt.rds")
GR_merge_dt=convert_GR(GR_merge,direction="DT")
#Separate tissue GO
#write.csv(data.table(Sample=unique(GR_merge$Sample)),"../downstream/output/human_analysis/dMML_dNME_GO/human_sample_tissue.csv")
sample_tissue=fread("../downstream/output/human_analysis/dMML_dNME_GO/human_sample_tissue.csv")
GR_merge_dt$tissue=sample_tissue[match(GR_merge_dt$Sample,Sample)]$Tissue
GR_merge_dNME_only=GR_merge_dt[dNME_pval<=pval_cutoff&dMML_pval>pval_cutoff]#28863
GR_merge_dMML_only=GR_merge_dt[dMML_pval<=pval_cutoff&dNME_pval>pval_cutoff]#5989
GR_merge_both=GR_merge_dt[dNME_pval<=pval_cutoff&dMML_pval<=pval_cutoff]#818
GO_NME_only=lapply(unique(GR_merge_dNME_only$tissue),function(x) 
            GO_run(unique(unlist(GR_merge_dNME_only[tissue==x]$genes_promoter)),unique(unlist(GR_merge$genes_promoter)),cluster=x,mapping="org.Hs.eg.db"))
GO_MML_only=lapply(unique(GR_merge_dMML_only$tissue),function(x) 
            GO_run(unique(unlist(GR_merge_dMML_only[tissue==x]$genes_promoter)),unique(unlist(GR_merge$genes_promoter)),cluster=x,mapping="org.Hs.eg.db"))
GO_both=lapply(unique(GR_merge_both$tissue),function(x) 
            GO_run(unique(unlist(GR_merge_both[tissue==x]$genes_promoter)),unique(unlist(GR_merge$genes_promoter)),cluster=x,mapping="org.Hs.eg.db"))

write.csv(do.call(rbind,lapply(GO_NME_only,function(x) x[FDR<=0.8,list(Term,FDR,FC,genes,cluster)])),'../downstream/output/human_analysis/QC/dNME_only_GO.csv')
write.csv(do.call(rbind,lapply(GO_MML_only,function(x) x[FDR<=0.8,list(Term,FDR,FC,genes,cluster)])),'../downstream/output/human_analysis/QC/dMML_only_GO.csv')
write.csv(do.call(rbind,lapply(GO_both,function(x) x[FDR<=0.8,list(Term,FDR,FC,genes,cluster)])),'../downstream/output/human_analysis/QC/both_GO.csv')

write.csv(GR_merge_dNME_only[,list(dMML,dMML_pval,Sample,dNME,dNME_pval,genes_promoter,genes_body,TSS,region)],"../downstream/output/human_analysis/dMML_dNME_GO/dNME_only.csv")
write.csv(GR_merge_dMML_only[,list(dMML,dMML_pval,Sample,dNME,dNME_pval,genes_promoter,genes_body,TSS,region)],"../downstream/output/human_analysis/dMML_dNME_GO/dMML_only.csv")
write.csv(GR_merge_both[,list(dMML,dMML_pval,Sample,dNME,dNME_pval,genes_promoter,genes_body,TSS,region)],"../downstream/output/human_analysis/dMML_dNME_GO/both.csv")

theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.text.x=element_text(size=7),
                                 axis.text.y=element_text(size=7))
pdf('../downstream/output/human_analysis/QC/dNME_vs_dMML_bin2d.pdf',width=5,height=5)
#Plotting dNME dMML
print(ggplot(GR_merge_dt,
                aes(x=dMML, y=dNME))+
                xlim(c(0,1))+ylim(c(0,1))+ggtitle("dMML and dNME relationship")+
                #geom_smooth(method="loess",se=FALSE)+
                xlab("dMML")+ ylab("dNME")+geom_bin2d(bins=100)+theme_glob+
                #geom_smooth(method="lm",color="red")+
                scale_fill_gradient(name = "count", trans = "log10",high="#132B43",low="#56B1F7"))
dev.off()

pdf('../downstream/output/human_analysis/QC/dNME_vs_dMML_bin2d_sig.pdf',width=5,height=5)
#Plotting dNME dMML
print(ggplot(GR_merge_dt[GR_merge$dNME_pval<=pval_cutoff|GR_merge$dMML_pval<=pval_cutoff],
                aes(x=dMML, y=dNME))+
                xlim(c(0,1))+ylim(c(0,1))+ggtitle("dMML and dNME relationship")+
                #geom_smooth(method="loess",se=FALSE)+
                xlab("dMML")+ ylab("dNME")+geom_bin2d(bins=100)+theme_glob+
                geom_smooth(method="lm",color="red")+
                scale_fill_gradient(name = "count", trans = "log10",high="#132B43",low="#56B1F7"))
dev.off()
#dNME vs N
pdf("../downstream/output/human_analysis/QC/dNME_N.pdf")
dNME_N=data.table(dNME=GR_merge$dNME,pval=GR_merge$dNME_pval,N=as.character(GR_merge$N))
dNME_N[GR_merge$N>=10]$N=">=10"
dNME_N$N=factor(dNME_N$N,levels=c(1:10,">=10"))
print(ggplot(dNME_N,aes(x=N,y=dNME))+geom_boxplot(outlier.shape=NA))
dev.off()


