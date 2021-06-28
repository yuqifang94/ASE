rm(list=ls())
source("mainFunctions_sub.R")
#Define ggplot theme

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
# add hypervar to allele-agnostic data ------------------------------------


NME_in=readRDS(NME_agnostic_file)
MML_in=readRDS(MML_agnostic_file)
genomic_features=readRDS(genomic_features_file)
GR_calc=data.frame()
scRNA_result=data.frame()


#Rename the columns to the same as hypervar
for(fn in dir(scRNA_dir,pattern='.rds')){
   exp_hypervar_in=readRDS(paste0(scRNA_dir,fn))
  colnames(exp_hypervar_in)=c('gene_name','mean','var','fit','hypervar_var','p.value','p.adj')
  exp_hypervar_in$hypervar_logvar =log(exp_hypervar_in)
  rownames(exp_hypervar_in)=exp_hypervar_in$gene_name
  saveRDS(exp_hypervar_in,paste0(scRNA_dir,'processed/',fn))
}
MML_hypervar_calc=GRanges()
NME_hypervar_calc=GRanges()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  hyper_var_file=gsub(agnostic_dir,paste0(scRNA_dir,'processed/'),hyper_var_file)
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){
    
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    #scRNA_result=rbind(scRNA_result,sp_hyper_var)
    #Add hypervaribility inforamtion
    #Add hypervaribility inforamtion
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    
    
  }else{cat("file not exist for:",sp,'\n')}
}

#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),'../downstream/output/human_analysis/NME_MAV/allele_agnostic_var_homogeneous2_rogue_entropy.rds')
# Find number of overlapped regions ---------------------------------------
GR_merge=readRDS(GR_merge_file)
#Only use merged data for H1
GR_merge=GR_merge[!(GR_merge$Sample%in%c("rep1 - H1","rep2 - H1"))]
genomic_features=readRDS(genomic_features_file)
hyper_var_all=readRDS('../downstream/output/human_analysis/NME_MAV/allele_agnostic_var_homogeneous2_rogue_entropy.rds')#cor=0.211722 
hyper_var_all=lapply(hyper_var_all,function(x) x[x$N>=2])
#Figure 3C and D in different context
#0.1744904 
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir='../downstream/output/human_analysis/NME_MAV/expression_NME/')
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="mean",dir='../downstream/output/human_analysis/NME_MAV/expression_NME/')
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="var",dir='../downstream/output/human_analysis/NME_MAV/expression_NME/')
hyper_var_all$NME_hypervar_calc$CV=sqrt(hyper_var_all$NME_hypervar_calc$var)/hyper_var_all$NME_hypervar_calc$mean
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="CV")

dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="hypervar_logvar")
dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="mean")
dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="var")
hyper_var_all$MML_hypervar_calc$CV=sqrt(hyper_var_all$MML_hypervar_calc$var)/hyper_var_all$MML_hypervar_calc$mean
dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="CV")

