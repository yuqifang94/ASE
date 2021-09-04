rm(list=ls())
source('mainFunctions_sub.R')


region_in=readRDS(tissue_out_filtered_fn)
#Getting Bin enhancer
enhancer=readRDS(bin_enhancer_rds)
UC_merge=readRDS(UC_merge_max_loc_file)
#Getting target genes
region_in_enhancer=lapply(region_in,function(x) {
  olap=findOverlaps(convert_GR(x$region),enhancer)
  x=x[queryHits(olap)]
  x$gene=enhancer[subjectHits(olap)]$`Target Gene`
  return(x)
})

region_in_enhancer=lapply(region_in_enhancer,function(x){
  x=x[,list(region,tissue,cluster,gene,region_type,dMML_cor,dNME_cor)]
  UC_merge_ts=UC_merge[[unique(x$tissue)]]
  print(head(UC_merge_ts[x$region,'dNME_max_UC_pair_adj']))
  x$dNME_max_pair=UC_merge_ts[x$region,'dNME_max_UC_pair_adj']
  x$dMML_max_pair =UC_merge_ts[x$region,'dMML_max_UC_pair_adj']
  x$UC_max_pair=UC_merge_ts[x$region,'UC_max_UC_pair_adj']
  x$UC_max_time =UC_merge_ts[x$region,'UC_max_time_adj']
  x$UC_max_time=gsub(paste0(unique(x$tissue),'-|-all'),'',x$UC_max_time)
  return(x)
  
})


GO_out_all=readRDS(GO_01_enhancer_fn)

tissue_sel=names(region_in_enhancer)
for(ts in tissue_sel){
  select_top_GO_out=select_top_GO(GO_out_all$all,ts,ptcount=0,FDR_cutoff=0.1,FC_cutoff=1.5)
  ts_sel_genes=unique(unlist(strsplit(select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5]$genes,';')))
  enhancer_regions_ts=region_in_enhancer[[ts]][gene %in% ts_sel_genes][order(dNME_max_pair,decreasing=T)]
  enhancer_regions_ts$GO_terms=unlist(lapply(enhancer_regions_ts$gene,function(x) paste(select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5][grepl(x,genes)]$Term,collapse = ';')))
  enhancer_regions_ts$GO_ID=unlist(lapply(enhancer_regions_ts$gene,function(x) paste(gsub('GO:','',select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5][grepl(x,genes)]$`GO.ID`),collapse = ';')))
  region_in_enhancer[[ts]]=enhancer_regions_ts
  write.csv(enhancer_regions_ts,paste0(gene_example_dir,ts,'.csv'),row.names = F)
}
#EFP:before 2941, after 2140
#This file contains all regions we analzyed and which motif is their binding site
motif_locus_ken=readRDS(tissue_region_motif_all_regions_fn)

motif_enhancer_dNME=list()
motif_enhancer_dMML=list()
for(ts in names(region_in_enhancer)){
  enhancer_regions_ts=region_in_enhancer[[ts]]
  dMML_motif=motif_sig_Ken(ts,"dMML",motif_locus_ken,enhancer_regions_ts,dir_in_Ken=Ken_motif_folder)
  dNME_motif=motif_sig_Ken(ts,"dNME",motif_locus_ken,enhancer_regions_ts,dir_in_Ken=Ken_motif_folder)
 #motif_region relationship
  motif_enhancer_dNME[[ts]]=dNME_motif$motif_locus
  motif_enhancer_dMML[[ts]]=dMML_motif$motif_locus
  #Region-motif relationship 
  enhancer_regions_ts$dMML_motif=dMML_motif$motif_locus_dt[match(enhancer_regions_ts$region,region)]$motif
  enhancer_regions_ts$dNME_motif=dNME_motif$motif_locus_dt[match(enhancer_regions_ts$region,region)]$motif

  region_in_enhancer[[ts]]=enhancer_regions_ts
  write.csv(enhancer_regions_ts,paste0(region_motif_dir,ts,'.csv'))
}
saveRDS(region_in_enhancer,enhancer_region_fn)
