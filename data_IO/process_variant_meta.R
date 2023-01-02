rm(list=ls())
source("mainFunctions_sub.R")
#Give each SNP an ASM information for each subject
variant_meta<- function(subj,variant_in,GR_in){ #variant_in for each subject, GR_in for each subject
  cat('Processing',subj,'\n')
  GR_in_subj=GR_in[GR_in$Subject==subj]
  variant_in_subj=variant_in[[subj]]
  sp=unique(GR_in_subj$Sample)
  gr_out_sp = GRanges()
  for (sps in sp){
    cat('Processing',sps,'\n')
    gr_out_sp=c(gr_out_sp,variant_meta_sp(variant_in_subj,GR_in_subj[GR_in_subj$Sample==sps]))
    
  }
  return(gr_out_sp)
  
  
}

#For each sample, assign NME to SNP within each region
variant_meta_sp <-function(variant_subj,GR_sp){
  olap=findOverlaps(variant_subj,GR_sp,maxgap =0,type='within')
  gr_out=variant_subj[queryHits(olap)]
  #Find GR_merge olap
  olap=subjectHits(olap)
  mcols(gr_out)=cbind(mcols(gr_out),mcols(GR_sp)[olap,])
  gr_out$refNME=NA
  gr_out$altNME=NA
  gr_out$refMML=NA
  gr_out$altMML=NA
  #NME
  gr_out$refNME[gr_out$GT %in% c("0/1","0|1")]= gr_out$NME1[gr_out$GT %in% c("0/1","0|1")]
  gr_out$altNME[gr_out$GT %in% c("0/1","0|1")]= gr_out$NME2[gr_out$GT %in% c("0/1","0|1")]
  gr_out$refNME[gr_out$GT %in% c("1/0","1|0")]= gr_out$NME2[gr_out$GT %in% c("1/0","1|0")]
  gr_out$altNME[gr_out$GT %in% c("1/0","1|0")]= gr_out$NME1[gr_out$GT %in% c("1/0","1|0")]
  
  #MML
  gr_out$refMML[gr_out$GT %in% c("0/1","0|1")]= gr_out$MML1[gr_out$GT %in% c("0/1","0|1")]
  gr_out$altMML[gr_out$GT %in% c("0/1","0|1")]= gr_out$MML2[gr_out$GT %in% c("0/1","0|1")]
  gr_out$refMML[gr_out$GT %in% c("1/0","1|0")]= gr_out$MML2[gr_out$GT %in% c("1/0","1|0")]
  gr_out$altMML[gr_out$GT %in% c("1/0","1|0")]= gr_out$MML1[gr_out$GT %in% c("1/0","1|0")]
  
  gr_out$Statistic=GR_sp$Statistic[olap]
  
  # gr_out$HetCpG=FALSE
  # #This may be different from previous result?
  # gr_out$HetCpG=((gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & !(gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG')) |
  #   (!(gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & (gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG'))
  return(gr_out)
}
variant_HetCpG=readRDS(variant_HetCpG_file)
GR_merge=readRDS(GR_merge_file)
# Processing variant based result -----------------------------------------
variant_HetCpG_meta=fastDoCall('c',lapply(names(variant_HetCpG),variant_meta,variant_in=variant_HetCpG,GR_in=GR_merge))
#Trinucleotide analysis
#variant_HetCpG_meta$mask_tri=unlist(lapply(variant_HetCpG_meta$REF_tri,mask_tri))
saveRDS(variant_HetCpG_meta,variant_HetCpG_meta_file)