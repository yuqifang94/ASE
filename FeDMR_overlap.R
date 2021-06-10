rm(list=ls())
source('mainFunctions_sub.R')
#downloaded from:https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2119-x/MediaObjects/41586_2020_2119_MOESM5_ESM.zip
#Unzip and `rename craniofacial EFP *`
UC_01=readRDS(paste0(dir_cluster_in_01,'uc_0.1_1.rds'))
UC_all=readRDS(UC_in_matrix_ls_file)
feDMR_dir='../downstream/input/mouse_analysis/FeDMR/'
FeDMR_olap=data.table()
for(ts in names(UC_01)){
  feDMR_ts=GRanges()
  #read in feDMR regions
  for(tsv_in in dir(feDMR_dir,pattern=ts)){
   
    feDMR_ts=c(feDMR_ts, makeGRangesFromDataFrame(fread(paste0(feDMR_dir,tsv_in)),keep.extra.columns=T))
    
  }
  
  
}