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

  olap_all=findOverlaps(UC_all[[ts]],feDMR_ts)
  olap_UC_01=findOverlaps(convert_GR(names(UC_01[[ts]])),feDMR_ts)
  FeDMR_olap=rbind(FeDMR_olap,
                   data.table(tissue=ts,
                              total_UC_analyzed=length(UC_all[[ts]]),
                              total_UC_01=length(UC_01[[ts]]),
                              overlap_total_UC_feDMR=length(unique(queryHits(olap_all))),
                              overlap_UC_01_feDMR=length(unique(queryHits(olap_UC_01))),
                              total_FeDMR=length(feDMR_ts),
                              FeDMR_covered_all_UC=length(unique(subjectHits(olap_all))),
                              FeDMR_covered_uc_01=length(unique(subjectHits(olap_UC_01)))
                              
                              ))
  
}
FeDMR_olap$proportion_overlap_uc_01=FeDMR_olap$FeDMR_covered_uc_01/FeDMR_olap$FeDMR_covered_all_UC
FeDMR_olap$proportion_overlap_uc_all=FeDMR_olap$FeDMR_covered_all_UC/FeDMR_olap$total_FeDMR
saveRDS(FeDMR_olap,'../downstream/output/mouse_analysis/Ecker_comparison/FeDMR_olap.rds')