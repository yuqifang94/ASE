rm(list=ls())
source("mainFunctions_sub.R")
library(gridExtra)
#Define ggplot theme

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
hyper_var_all=readRDS(paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV.rds'))
NME_hypervar=hyper_var_all$NME_hypervar_calc
NME_hypervar_dt=as.data.table(mcols(NME_hypervar))
NME_hypervar_dt_exp=unique(NME_hypervar_dt[!is.na(mean)&!is.na(var)&!is.na(hypervar_var)&!is.na(hypervar_logvar),list(mean,var,hypervar_var,hypervar_logvar,Sample)])
pdf('../downstream/output/human_analysis/NME_MAV/mean_var_cor.pdf')
for (sp in unique(NME_hypervar_dt_exp$Sample)){
mean_var=ggplot(NME_hypervar_dt_exp[Sample==sp],aes(x=log(mean),y=log(var)))+geom_point(alpha=0.1)+geom_smooth()+
    xlab("mean expression")+ylab("expression variance")+ggtitle(paste0("mean variance relationship\n",sp))
mean_MAV=ggplot(NME_hypervar_dt_exp[Sample==sp],aes(x=log(mean),y=hypervar_logvar))+geom_point(alpha=0.1)+geom_smooth()+
    xlab("mean expression")+ylab("MAV")+ggtitle(paste0("mean MAV relationship\n",sp))
print(grid.arrange(mean_var,mean_MAV,nrow=1))
}
dev.off()
#
pdf('../downstream/output/human_analysis/NME_MAV/mean_MAV_cor.pdf')
dev.off()