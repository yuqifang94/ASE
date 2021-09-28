source('mainFunctions_sub.R')
scRNA_dir_rogue="../downstream/input/human_analysis/expression_NME_rogue/"
NME_MAV_human_out_dir='../downstream/output/human_analysis/NME_MAV/allele_agnostic_var_homogeneous2_MAV.rds'
#Rename the columns to the same as hypervar
for(fn in dir(scRNA_dir_rogue,pattern='.rds')){
   exp_hypervar_in=readRDS(paste0(scRNA_dir_rogue,fn))
  colnames(exp_hypervar_in)=c('gene_name','mean','var','fit','hypervar_var','p.value','p.adj')
  exp_hypervar_in$hypervar_logvar =log(exp_hypervar_in)
  rownames(exp_hypervar_in)=exp_hypervar_in$gene_name
  saveRDS(exp_hypervar_in,paste0(scRNA_dir_rogue,'processed/',fn))
}

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



MML_hypervar_calc=GRanges()
NME_hypervar_calc=GRanges()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  hyper_var_file=gsub(scRNA_dir,paste0(scRNA_dir_rogue,'processed/'),hyper_var_file)
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){
    
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    
    
  }else{cat("file not exist for:",sp,'\n')}
}

#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),NME_MAV_human_out_dir)