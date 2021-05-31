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
agnostic_dir="../downstream/input/scRNA/"
scRNA_dir="../downstream/input/human_analysis/expression_NME/"
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


# # NME vs VMR currently not in use--------------------------------------------------------------
# NME_in=readRDS(NME_agnostic_file)
# #Brain
# load("../downstream/input/vmrs_hg19_brain.rda")
# vmr_HC2=vmrs_hg19$HC2
# vmr_HC1=vmrs_hg19$HC1
# names(vmr_HC2)=NULL
# names(vmr_HC1)=NULL
# #Do HC2
# vmr=do.call(c,vmr_HC2)
# saveRDS(vmr,'../downstream/output/vmr_HC2.rds')
# vmr=readRDS('../downstream/output/vmr_HC2.rds')
# NME_in_brain=NME_in[NME_in$Sample%in%c('Brain_Hippocampus_middle_paired - 149','Brain_Hippocampus_middle_paired - 150')]
# OR_quant=data.frame()
# for(percent in unique(NME_in_brain$quant_score)){
#   OR=OR_VMR(NME_in_brain,vmr,percent,NME_quant='quant_score')
#   OR_quant=rbind(OR_quant,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
#   
# }
# write.csv(OR_quant,'../downstream/output/brain_VMR.csv')
# 
# #chromHMM state get enhancer
# ah = AnnotationHub()
# ENCODE_name=ENCODE_to_sample(unique(NME_in_brain$Sample))
# ah_num=names(query(ah, c("chromhmmSegmentations", unique(ENCODE_name$ENCODE))))
# chromHMM=ah[[ah_num]]
# chromHMM_enc=chromHMM[chromHMM$abbr%in%c("7_Enh","6_EnhG")]
# NME_in_brain_enc=subsetByOverlaps(NME_in_brain,chromHMM_enc)
# vmr_enc=subsetByOverlaps(vmr,chromHMM_enc)
# #Find nearest genes
# genomic_features=readRDS(genomic_features_file)
# genomic_features=genomic_features$TSS
# NME_in_brain_enc=dist_calc(NME_in_brain_enc,genomic_features)
# vmr_enc=dist_calc(vmr_enc,genomic_features)
# #Find overlap between nearest genes
# NME_in_brain_enc_highNME=NME_in_brain_enc[NME_in_brain_enc$NME>=quantile(NME_in_brain$NME,prob=0.99)]#0.87 cutoff
# NME_in_brain_enc_highNME=as.data.table(mcols(NME_in_brain_enc_highNME))
# NME_in_brain_enc_highNME=NME_in_brain_enc_highNME[,list(NME=mean(NME),dist=mean(abs(dist))),by=list(gene)]
# NME_in_brain_enc_highNME=NME_in_brain_enc_highNME[order(NME,decreasing = T)]$gene
# bg=unique(c(vmr_enc$gene,NME_in_brain_enc$gene))
# NME_VMR=sum((bg %in%NME_in_brain_enc_highNME)&(bg %in%vmr_enc$gene))
# NME_nonVMR=sum((bg %in%NME_in_brain_enc_highNME)&!(bg %in%vmr_enc$gene))
# VMR_nonNME=sum(!(bg %in%NME_in_brain_enc_highNME)&(bg %in%vmr_enc$gene))
# nonVMR_nonNME=sum(!(bg %in%NME_in_brain_enc_highNME)&!(bg %in%vmr_enc$gene))
# 
# fisher.test(matrix(c(NME_VMR,NME_nonVMR,VMR_nonNME,nonVMR_nonNME),nrow=2))
# 
# write(unique(NME_in_brain_enc_highNME[NME_in_brain_enc_highNME%in%vmr_enc$gene]),'../downstream/output/VMR_NME.txt')
