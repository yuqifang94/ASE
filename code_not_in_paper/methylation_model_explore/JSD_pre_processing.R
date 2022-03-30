rm(list=ls())
source('mainFunctions_sub.R')


#read in JSD for MDS
in_dir='../downstream/data/mouse_JSD/'
JSD_in=GRanges()
JSD_in_ls=mclapply(dir(in_dir,pattern = '.*jsd.bedGraph'),function(x){JSD_in=read.agnostic.mouse.uc(paste(in_dir,x,sep=''))
JSD_in$JSD=JSD_in$score
return(JSD_in)},mc.cores=24)
JSD_in=fastDoCall('c',JSD_in_ls)
JSD_in$tissue=sub('_.*','',JSD_in$Sample)
JSD_in=convert_GR(JSD_in,"DT")
JSD_in_dt=dcast.data.table(JSD_in,region+N~Sample,value.var="JSD")
saveRDS(JSD_in_dt,JSD_in_fn)#74% regiOn have all data

#Filter out the region with N >=18
UC_in=readRDS(UC_in_matrix_ls_file)
JSD_in_gr=unique(granges(JSD_in))
UC_in_sub=mclapply(UC_in,function(x) {return(subsetByOverlaps(x,JSD_in,type='equal'))})
saveRDS(UC_in_sub,UC_in_QC_fn)
# EFP : 0.9786391
# Lung : 0.9792112
# NT : 0.9789671
# forebrain : 0.9785839
# heart : 0.9785712
# hindbrain : 0.978633
# intestine : 0.9788505
# kidney : 0.9789867
# limb : 0.9786628
# liver : 0.978591
# midbrain : 0.9786526
# stomach : 0.9788745



JSD_in=readRDS(JSD_in_fn)
names(JSD_in)=unlist(lapply(JSD_in,function(x) sub('_.*','',colnames(mcols(x))[1])))
for(ts in names(JSD_in)){
  colnames(mcols(JSD_in[[ts]]))=gsub('_','-',sub(paste0('-',ts),'',gsub("_all","",colnames(mcols(JSD_in[[ts]])))))
  mcols(JSD_in[[ts]])=mcols(JSD_in[[ts]])[,which(!grepl("P0",colnames(mcols(JSD_in[[ts]]))))]
  JSD_in[[ts]]=JSD_in[[ts]][which(rowSums(is.na(mcols(JSD_in[[ts]])))==0)]

}
saveRDS(JSD_in,'JSD_agnostic_mouse_matrix_dedup_N2_all_merged_ls_ft.rds')
gff_in=readGFFAsGRanges('../mm10_allele_agnostic_analysis.gff')
gff_in_sub=gff_in[as.numeric(gff_in$N)<=17]#98.3%
JSD_in_sub=mclapply(JSD_in,function(x) {return(subsetByOverlaps(x,gff_in_sub,type='equal'))})
for(ts in names(JSD_in_sub)){
  
  
  cat("top 10% quantile for",ts,'is',quantile(as.vector(as.matrix(mcols(JSD_in_sub[[ts]]))),probs=0.9),'\n')
  cat("Proportion regions letft for",ts,"is",length(JSD_in_sub[[ts]])/length(JSD_in[[ts]]),'\n')
}
# top 10% quantile for EFP is 0.2912969
# Proportion regions letft for EFP is 0.9874118
# top 10% quantile for Lung is 0.2613509
# Proportion regions letft for Lung is 0.9849098
# top 10% quantile for NT is 0.2768564
# Proportion regions letft for NT is 0.9861501
# top 10% quantile for forebrain is 0.2898555
# Proportion regions letft for forebrain is 0.9855473
# top 10% quantile for heart is 0.2799574
# Proportion regions letft for heart is 0.9858041
# top 10% quantile for hindbrain is 0.2813578
# Proportion regions letft for hindbrain is 0.9860615
# top 10% quantile for intestine is 0.2584561
# Proportion regions letft for intestine is 0.9845821
# top 10% quantile for kidney is 0.2568008
# Proportion regions letft for kidney is 0.985014
# top 10% quantile for limb is 0.2953809
# Proportion regions letft for limb is 0.9870893
# top 10% quantile for liver is 0.3264372
# Proportion regions letft for liver is 0.9849181
# top 10% quantile for midbrain is 0.2777929
# Proportion regions letft for midbrain is 0.985423
# top 10% quantile for stomach is 0.2534778
# Proportion regions letft for stomach is 0.9843143

saveRDS(JSD_in_sub,'JSD_agnostic_mouse_matrix_dedup_N2_all_merged_ls_less_equal_17CG.rds')