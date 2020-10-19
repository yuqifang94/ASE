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
genomic_features=readRDS(genomic_features_file)
hyper_var_all=readRDS('../downstream/output/allele_agnostic_var.rds')#cor=0.211722 
#Figure 3A
pdf('../downstream/output/graphs/Figure3/Figure3A_hypervar_dist_NME.pdf',width=3.5,height=3.5)
dist_plot_run(hyper_var_all$NME_hypervar_calc,theme_glob,ylab="NME")
dev.off()
pdf('../downstream/output/graphs/Figure3/Figure3As_hypervar_dist_MML.pdf',width=3.5,height=3.5)#0.189768 
dist_plot_run(hyper_var_all$MML_hypervar_calc,theme_glob,ylab="MML")
dev.off()
pdf('../downstream/output/graphs/Figure3/Figure3As_mean_dist_NME.pdf',width=3.5,height=3.5)#-0.1006436 
dist_plot_run(hyper_var_all$NME_meanvar_calc,theme_glob,ylab="NME")
dev.off()
pdf('../downstream/output/graphs/Figure3/Figure3As_mean_dist_MML.pdf',width=3.5,height=3.5)#-0.1101482 
dist_plot_run(hyper_var_all$MML_meanvar_calc,theme_glob,ylab="MML")
dev.off()

#Figure3B
NME_hypervar_calc=hyper_var_all$NME_hypervar_calc[abs(hyper_var_all$NME_hypervar_calc$dist)<=500]
NME_hypervar_calc=data.table(score=NME_hypervar_calc$score,hypervarquant=NME_hypervar_calc$hypervarquant001,Sample=NME_hypervar_calc$Sample)
NME_hypervar_calc=NME_hypervar_calc[,list(median_score=median(score)),by=list(hypervarquant,Sample)]
NME_hypervar_calc=NME_hypervar_calc[,list(hypervarquant=hypervarquant,NME=median_score,cor=cor(median_score,hypervarquant)),by=list(Sample)]
NME_hypervar_calc=NME_hypervar_calc[order(NME_hypervar_calc$cor,decreasing = F),]
NME_hypervar_calc$Sample=factor(NME_hypervar_calc$Sample,levels = unique(NME_hypervar_calc$Sample))
pdf('../downstream/output/graphs/Figure3/Figure3B_hypervaribility_NME.pdf',width=7,height=7)
ggplot(NME_hypervar_calc,aes(hypervarquant,Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
  xlab('Hypervaribility quantile')+ylab('Sample')+theme_glob+theme(legend.position = 'bottom')
dev.off()


#Corelation between dMML and hypervaribility
hyper_var_all_NME=hyper_var_all$NME_hypervar_calc
print(cor.test(hyper_var_all_NME$dMML[abs(hyper_var_all_NME$dist)<=500],hyper_var_all_NME$exp_stat[abs(hyper_var_all_NME$dist)<=500]))
NME_dMML_dt=data.table(NME=hyper_var_all_NME$score_original,dMML_pval=hyper_var_all_NME$dMML_pval)
NME_dMML_dt$MML_ASM="non-dMML_ASM"
NME_dMML_dt$MML_ASM[NME_dMML_dt$dMML_pval<=pval_cutoff]="dMML_ASM"
pdf('../downstream/output/graphs/FigureS4/NME_dMML.pdf',width=3.5,height=3.5)
ggplot(data=NME_dMML_dt,aes(x=NME,group=MML_ASM,color=MML_ASM))+
  geom_density(size=1)+theme_glob+theme(legend.position="bottom",legend.title = element_blank())
dev.off()
