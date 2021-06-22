rm(list=ls())
source('mainFunctions_sub.R')
#downloaded from:https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2119-x/MediaObjects/41586_2020_2119_MOESM5_ESM.zip
#Unzip and `rename craniofacial EFP *`

UC_all=readRDS(UC_in_matrix_ls_file)
feDMR_dir='../downstream/input/mouse_analysis/FeDMR/'
FeDMR_olap=data.table()
for(ts in names(UC_all)){
  feDMR_ts=GRanges()
  #read in feDMR regions
  for(tsv_in in dir(feDMR_dir,pattern=ts)){
   
    feDMR_ts=c(feDMR_ts, makeGRangesFromDataFrame(fread(paste0(feDMR_dir,tsv_in)),keep.extra.columns=T))
    
  }
  UC_01=fread(paste0(dir_out_cluster01,ts,'.csv'))
  olap_all=findOverlaps(UC_all[[ts]],feDMR_ts)
  olap_UC_01=findOverlaps(convert_GR(UC_01$regions),feDMR_ts)
  olap_UC_01_MML=findOverlaps(convert_GR(UC_01[region_type=='MML only']$regions),feDMR_ts)
  FeDMR_olap=rbind(FeDMR_olap,
                   data.table(tissue=ts,
                              total_UC_analyzed=length(UC_all[[ts]]),
                              total_UC_01=nrow(UC_01),
                              total_UC_01_MML=sum(UC_01$region_type=='MML only',na.rm = T),
                              olap_UC_01_MML_FeDMR=length(unique(queryHits(olap_UC_01_MML))),
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
write.csv(FeDMR_olap,'../downstream/output/mouse_analysis/Ecker_comparison/FeDMR_olap.csv',row.names = F)