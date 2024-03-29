source('mainFunctions_sub.R')
UC_in_MDS=readRDS(UC_in_MDS_all_file)
analyzedRegion=UC_in_MDS$region
rm(UC_in_MDS)
analyzedRegion_gr=convert_GR(analyzedRegion,direction="GR")
mm10_CpG=getCpgSitesmm10()
analyzedRegion_gr$N=countOverlaps(analyzedRegion_gr,mm10_CpG)
mm10_CpG_covered=subsetByOverlaps(mm10_CpG,analyzedRegion_gr,minoverlap=2)
clustered=readRDS(tissue_out_filtered_fn)
