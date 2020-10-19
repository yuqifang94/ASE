# ###########Stat checking#####################################
# hetCpG_gff=readRDS(hetCpG_gff_file)
# GR_merge=readRDS(GR_merge_file)
# ###Checking regions with statistics
# covered_region=data.frame(sample=NULL,covered=NULL,gff=NULL)
# for(sp in unique(GR_merge$Sample)){
#   merge_sp=GR_merge[GR_merge$Sample==sp]
#   gff_sp=hetCpG_gff[[unique(merge_sp$Subject)]]
#   covered_sample=data.frame(sample=sp,covered=length(subsetByOverlaps(merge_sp,gff_sp)),gff=length(gff_sp))
#   covered_region=rbind(covered_region,covered_sample)
# }
# covered_region$`coverage(%)`=covered_region$covered/covered_region$gff*100