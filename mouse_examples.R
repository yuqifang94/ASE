rm(list=ls())
source('mainFunctions_sub.R')
#get UC_merge_max_loc_cluster01.rds
UC_merge_max_loc_01=readRDS(UC_merge_max_loc_01_file)
#Read in selected GO regions
region_in_fn='../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_all_region.rds'
enhancer_fn='../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds'
enhancer_bed_out_fn='../downstream/output/mouse_analysis/enhancers/bin_enhancer.bed'
UC_merge_fn='../downstream/input/mouse_analysis/UC_merge_max_loc_cluster01.rds'
enhancer_region_fn='../downstream/output/mouse_analysis/GO_analysis/all_regions_enchancer.rds'
GO_out_fn='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_01_enhancer.rds'
motif_Ken_fn='../downstream/output/mouse_analysis/motif_analysis/tissue_region_motif_all_regions.rds'
Ken_motif_folder='../downstream/input/mouse_analysis/motif_analysis/mouse_motif_enrichment_0526/'
region_motif_dir='../downstream/output/mouse_analysis/motif_analysis/region_motif/'
gene_example_dir='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/gene_examples/'
ChiP_motif_dir='../downstream/output/mouse_analysis/motif_analysis/motif_Chip_rds/'
region_in=readRDS(region_in_fn)
#Getting Bin enhancer
enhancer=readRDS(enhancer_fn)
region_in_enhancer=lapply(region_in,function(x) {
  olap=findOverlaps(convert_GR(x$region),enhancer)
  x=x[queryHits(olap)]
  x$gene=enhancer[subjectHits(olap)]$`Target Gene`
  return(x)
})
enhancer_bed_out=data.table(chr=seqnames(enhancer),start=start(enhancer),end=end(enhancer),genes=enhancer$`Target Gene`)
write.table(enhancer_bed_out,enhancer_bed_out_fn,row.names = F,col.names =F,quote=F )
UC_merge=readRDS(UC_merge_fn)

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
saveRDS(region_in_enhancer,enhancer_region_fn)
# #Mouse gene example
GO_out_all=readRDS(GO_out_fn)
enhancer_regions=readRDS(enhancer_region_fn)
enhancer_regions_all=list()
tissue_sel=names(enhancer_regions)
for(ts in tissue_sel){
  select_top_GO_out=select_top_GO(GO_out_all$all,ts,ptcount=0,FDR_cutoff=0.1,FC_cutoff=1.5)
  ts_sel_genes=unique(unlist(strsplit(select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5]$genes,';')))
  enhancer_regions_ts=enhancer_regions[[ts]][gene %in% ts_sel_genes][order(dNME_max_pair,decreasing=T)]
  enhancer_regions_ts$GO_terms=unlist(lapply(enhancer_regions_ts$gene,function(x) paste(select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5][grepl(x,genes)]$Term,collapse = ';')))
  enhancer_regions_ts$GO_ID=unlist(lapply(enhancer_regions_ts$gene,function(x) paste(gsub('GO:','',select_top_GO_out$tissue_all_merged[FDR<=0.1&FC>=1.5][grepl(x,genes)]$`GO.ID`),collapse = ';')))
  enhancer_regions_all[[ts]]=enhancer_regions_ts
  write.csv(enhancer_regions_ts,paste0(gene_example_dir,ts,'.csv'),row.names = F)
}
# saveRDS(enhancer_regions_all,enhancer_region_fn)
# mouse_motif_overlap -----------------------------------------------------
# tissue_sel=c("forebrain","heart",'limb','midbrain','hindbrain','liver','EFP')
# GO_out_all=readRDS('../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run/GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_enhancer.rds')
# GO_out_all=lapply(GO_out_all$all,function(x) do.call(rbind,lapply(x,function(y) 
#   cbind(y$csv_in_ts_clu[,gsub('UC-','',colnames(y$csv_in_ts_clu)) %in% paste0('E',10:20,'.5-E',11:21,'.5'),with=FALSE],
#         y$csv_in_ts_clu[,gsub('dMML-','',colnames(y$csv_in_ts_clu)) %in% paste0('E',10:20,'.5-E',11:21,'.5'),with=FALSE],
#         y$csv_in_ts_clu[,gsub('dNME-','',colnames(y$csv_in_ts_clu)) %in% paste0('E',10:20,'.5-E',11:21,'.5'),with=FALSE],
#    #y$csv_in_ts_clu[,grepl("UC-|dNME-|dMML-",colnames(y$csv_in_ts_clu)),with=FALSE],
#     y$csv_in_ts_clu[,list(cluster,region,region_type)]))))
# cor_in=readRDS('../downstream/output/mouse_analysis/correlation/tissue_out_N17_kmeans_10run_filtered_ref11.rds')
# GO_out_all=lapply(GO_out_all,function(x) {
#   uc_dt=  x[,grepl("UC-",colnames(x)),with=FALSE]
#   dNME_dt=  x[,grepl("dNME-",colnames(x)),with=FALSE]
#   dMML_dt=  x[,grepl("dMML-",colnames(x)),with=FALSE]
#   uc_max=apply(uc_dt,1,which.max)
#   x$UC_max_time=gsub('UC-','',colnames(uc_dt)[uc_max])
#   x$dNME_max_UC_pair=as.data.frame(dNME_dt)[cbind(seq_along(uc_max), uc_max)]
#   x$UC_max_UC_pair=as.data.frame(uc_dt)[cbind(seq_along(uc_max), uc_max)]
#   x$dMML_max_UC_pair=as.data.frame(dMML_dt)[cbind(seq_along(uc_max), uc_max)]
#   return(x)
# })
# enhancer_regions=readRDS('../downstream/output/mouse_analysis/enhancers/GO_regions_enchancer.rds')

# CpG_mm10=CpG_mm10=getCpgSitesmm10()
# end(CpG_mm10)=start(CpG_mm10)+1
# motif_locus_ken_CG=lapply(motif_locus_ken,function(motifs) {
#   for(mt in names(motifs)){
#     if(length(motifs[[mt]])>0){
#     motifs[[mt]]$mouse_region=motifs[[mt]]$region
#     motifs[[mt]]$region=NULL
#     motifs[[mt]]$motif=mt
#     motifs[[mt]]$NCG=countOverlaps(motifs[[mt]],CpG_mm10,minoverlap = 2)
#     }
#     
#   }
#   return(do.call(c,motifs))
#   
# })
# motif_locus_ken_CG=lapply(motif_locus_ken_CG,function(x){
#   
#   names(x)=NULL
#   return(do.call('c',x))
#   
# }
#   
#   )
# motif_locus_ken_CG_olap=lapply(motif_locus_ken_CG,function(x) x[x$NCG>0])
# motif_locus_ken_CG_olap=lapply(motif_locus_ken_CG_olap,convert_GR,dir="DT")
# 
# 
# saveRDS(motif_locus_ken_CG,'../downstream/output/mouse_analysis/motif_analysis/tissue_region_motif_all_locus_CG.rds')
motif_locus_ken=readRDS(motif_Ken_fn)
enhancer_regions=readRDS(enhancer_region_fn)
motif_enhancer_dNME=list()
motif_enhancer_dMML=list()
for(ts in names(enhancer_regions)){
  
  enhancer_regions_ts=enhancer_regions[[ts]]
  # enhancer_regions_ts$dNME_max_UC_pair =GO_out_all[[ts]][match(enhancer_regions_ts$region,region)]$dNME_max_UC_pair
  # enhancer_regions_ts$dMML_max_UC_pair =GO_out_all[[ts]][match(enhancer_regions_ts$region,region)]$dMML_max_UC_pair
  # enhancer_regions_ts$dNME_UC_cor=cor_in[[ts]][match(enhancer_regions_ts$region,region)]$dNME_cor
  # enhancer_regions_ts$dMML_UC_cor=cor_in[[ts]][match(enhancer_regions_ts$region,region)]$dMML_cor

  dMML_motif=motif_sig_Ken(ts,"dMML",motif_locus_ken,enhancer_regions_ts,dir_in_Ken=Ken_motif_folder)
  dNME_motif=motif_sig_Ken(ts,"dNME",motif_locus_ken,enhancer_regions_ts,dir_in_Ken=Ken_motif_folder)

 #motif_region relationship
  
  motif_enhancer_dNME[[ts]]=dNME_motif$motif_locus
  motif_enhancer_dMML[[ts]]=dMML_motif$motif_locus
  #Region-motif relationship 
  enhancer_regions_ts$dMML_motif=dMML_motif$motif_locus_dt[match(enhancer_regions_ts$region,region)]$motif
  enhancer_regions_ts$dNME_motif=dNME_motif$motif_locus_dt[match(enhancer_regions_ts$region,region)]$motif

  enhancer_regions[[ts]]=enhancer_regions_ts
  write.csv(enhancer_regions_ts,paste0(region_motif_dir,ts,'.csv'))
  


}
saveRDS(enhancer_regions,enhancer_region_fn)




#create pubMed searching queue based on tissue

enhancer_regions=readRDS(enhancer_region_fn)
#table for key items
#In format of  "Dlx3[Title/Abstract] AND craniofacial[Title/Abstract]AND (Human OR Mouse)"
#"forebrain" "heart"     "limb"      "midbrain"  "hindbrain" "liver"     "EFP"

tissue_sel=names(enhancer_regions)
key_terms=c(
  "((craniofacial[Title/Abstract]) OR (face[Title/Abstract]) OR (head[Title/Abstract]))",
  "(forebrain[Title/Abstract])",
  "(heart[Title/Abstract])",
  "(hindbrain[Title/Abstract])",
  "(limb[Title/Abstract])",
  "(liver[Title/Abstract])",
  "(midbrain[Title/Abstract])"
 
  
  
)
names(key_terms)=tissue_sel
enhancer_regions_motif_dNME_all=list()
enhancer_regions_motif_dMML_all=list()
for(ts in tissue_sel){
  cat('Processing:',ts, 'with key:',key_terms[ts],'\n')
  tt1=proc.time()[[3]]
  GO_in=fread(paste0(gene_example_dir,unique(ts),'.csv'))
  enhancer_regions[[ts]]$GO_ID=GO_in[match(enhancer_regions[[ts]]$region,region)]$GO_ID
  enhancer_regions_motif=enhancer_regions[[ts]][(!is.na(dNME_motif)|!is.na(dMML_motif))&(!is.na(GO_ID))]
  enhancer_regions_motif=enhancer_regions_motif[order(dNME_max_pair,decreasing=T)]
  selected_genes=unique(enhancer_regions_motif$gene)
  gene_pub_med=do.call(rbind,lapply(selected_genes,pubmed_rec,keys=key_terms[ts]))
  if(!is.null(gene_pub_med)){
    #maximum countof 300
    PMID_unique=unique(gene_pub_med$PMID)
    
    #if(length(PMID_unique)>300){
      #Counting number of citations
    # breaks=ceiling(length(PMID_unique)/300)
    # PMID_cut=cut(1:length(PMID_unique),breaks,label=FALSE)
    # #Sys.sleep(300)
    # pmc_ref_count=c()
    # for(idx in 1:breaks){
    # 
    #   pmc_ref_count=c(pmc_ref_count,do.call(c,lapply(entrez_summary(db="pubmed", id=PMID_unique[PMID_cut==idx]),function(x) x$pmcrefcount)))
    #   Sys.sleep(120)
    # }
    # }else( pmc_ref_count=do.call(c,lapply(entrez_summary(db="pubmed", id=PMID_unique),function(x) x$pmcrefcount)))
    # 
    # gene_pub_med$pmc_count=pmc_ref_count[gene_pub_med$PMID]
    #gene_pub_med_collapse=gene_pub_med[,list(PMID=paste(PMID,collapse=";"),pmc_count=paste(pmc_count,collapse = ';'),title=paste(title,collapse=' ; ')),by=gene]
    gene_pub_med_collapse=gene_pub_med[,list(PMID=paste(PMID,collapse=";"),title=paste(title,collapse=' ; ')),by=gene]
    #Only genes with Pubmed literature
    enhancer_regions_motif_dNME=enhancer_regions_motif[!is.na(dNME_motif)]
    enhancer_regions_motif_dMML=enhancer_regions_motif[!is.na(dMML_motif)]
    enhancer_regions_motif_dNME=cbind(enhancer_regions_motif_dNME,gene_pub_med_collapse[match(enhancer_regions_motif_dNME$gene,gene)])
    enhancer_regions_motif_dMML=cbind(enhancer_regions_motif_dMML,gene_pub_med_collapse[match(enhancer_regions_motif_dMML$gene,gene)])
  }
  else{
    enhancer_regions_motif_dNME=enhancer_regions_motif[!is.na(dNME_motif)]
  enhancer_regions_motif_dMML=enhancer_regions_motif[!is.na(dMML_motif)]
  enhancer_regions_motif_dMML$PMID=NA
  enhancer_regions_motif_dNME$PMID=NA
  enhancer_regions_motif_dMML$title=NA
  enhancer_regions_motif_dNME$title=NA
  }

    cat("Writing result\n")
    write.csv(enhancer_regions_motif_dNME[,list(region,tissue,dNME_motif,cluster,gene,region_type,
                                                dNME_max_pair ,dMML_max_pair,UC_max_time,
                                               UC_max_pair,dNME_cor,
                                               dMML_cor,PMID,title,GO_ID)],
              paste0(region_motif_dir,ts,'_dNME_Pubmed_annotated.csv'),
              row.names = F,quote=F)
    write.csv(enhancer_regions_motif_dMML[,list(region,tissue,dNME_motif,cluster,gene,region_type,
                                                dNME_max_pair ,dMML_max_pair,UC_max_time,
                                                UC_max_pair,dNME_cor,
                                                dMML_cor,PMID,title,GO_ID)],
              paste0(region_motif_dir,ts,'_dMML_Pubmed_annotated.csv'),
              row.names = F,quote=F)
    enhancer_regions_motif_dNME_all[[ts]]=enhancer_regions_motif_dNME
    enhancer_regions_motif_dMML_all[[ts]]=enhancer_regions_motif_dMML
  
 
}

saveRDS(enhancer_regions_motif_dNME_all,'../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dNME_all.rds')
saveRDS(enhancer_regions_motif_dMML_all,'../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dMML_all.rds')

enhancer_regions_motif_dNME_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dNME_all.rds')
enhancer_regions_motif_dMML_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dMML_all.rds')

enhancer_regions_motif_dNME_all_dNME=do.call(rbind,enhancer_regions_motif_dNME_all)
enhancer_regions_motif_dNME_all_dNME=enhancer_regions_motif_dNME_all_dNME[region_type=='NME only',list(tissue,dNME_motif,cluster,gene,dNME_max_pair,dNME_cor,dMML_max_pair,dMML_cor,UC_max_time,PMID)]
colnames(enhancer_regions_motif_dNME_all_dNME)=c('Tissue','Motif','Cluster','Gene','dNME','dNME-UC correlation','dMML','dMML-UC correlation','Stage','PMID')
enhancer_regions_motif_dMML_all_dMML=do.call(rbind,enhancer_regions_motif_dMML_all)
enhancer_regions_motif_dMML_all_dMML=enhancer_regions_motif_dMML_all_dMML[region_type=='MML only',list(tissue,dMML_motif,cluster,gene,dNME_max_pair,dNME_cor,dMML_max_pair,dMML_cor,UC_max_time,PMID)]
colnames(enhancer_regions_motif_dMML_all_dMML)=c('Tissue','Motif','Cluster','Gene','dNME','dNME-UC correlation','dMML','dMML-UC correlation','Stage','PMID')
write.csv(enhancer_regions_motif_dNME_all_dNME,'../downstream/output/mouse_analysis/motif_analysis/merged_motif_dNME.csv',row.names = F)
write.csv(enhancer_regions_motif_dMML_all_dMML,'../downstream/output/mouse_analysis/motif_analysis/merged_motif_dMML.csv',row.names = F)

# Writing output ----------------------------------------------------------
# 
# 
# for(ts in names(enhancer_regions_motif_dNME_all)){
#   write.csv(enhancer_regions_motif_dNME_all[[ts]][,list(region,tissue,dNME_motif,cluster,gene,region_type,
#                                               dNME_max_UC_pair ,dMML_max_UC_pair,UC_max_time,
#                                               UC_max_time,UC_max_pair,dNME_UC_cor,
#                                               dMML_UC_cor,dNME_UC_cor,PMID,title,GO_ID)],
#             paste0(region_motif_dir,ts,'_dNME_Pubmed_annotated.csv'),
#             row.names = F,quote=F)
#   
#   
# }
# write.csv(do.call(rbind,lapply(enhancer_regions_motif_dNME_all,function(x) 
#   x[region_type=="dNME_only",list(tissue,dNME_motif,gene,UC_max_time,
#                                   dNME_max_UC_pair,dNME_UC_cor,
#                                   dMML_max_UC_pair,dMML_UC_cor,
#                                   UC_max_pair,
#                                   PMID,GO_ID)])),
#   '../downstream/output/mouse_analysis/motif_analysis/region_motif/all_regions_dNME_only.csv',row.names = F,quote=F)
# write.csv(do.call(rbind,lapply(enhancer_regions_motif_dMML_all,function(x) 
#   x[region_type=="dMML_only",list(tissue,dMML_motif,gene,UC_max_time,
#                                   dNME_max_UC_pair,dNME_UC_cor,
#                                   dMML_max_UC_pair,dMML_UC_cor,
#                                   UC_max_pair,
#                                   PMID,GO_ID)])),
#   '../downstream/output/mouse_analysis/motif_analysis/region_motif/all_regions_dMML_only.csv',row.names = F,quote=F)
# #Read in and rank the dNME motif and dMML motif separately
# dir_motif_gene='../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/'
# #dNME only motif, need to have annotated GO terms
# motif_dNME_out=data.table()
# for(fn in dir(dir_motif_gene,pattern='dNME')){
#   csv_in=fread(paste0(dir_motif_gene,fn))
#   # need to have annotated GO terms
#   GO_in=fread(paste0('../downstream/output/mouse_analysis/GO_analysis/gene_examples/',unique(csv_in$tissue),'.csv'))
#   csv_in=csv_in[region %in% GO_in$region&region_type=="dNME_only"]
#   csv_in$GO_terms=GO_in[match(csv_in$region,region)]$GO_terms
#   #Change to dNME accordingly
#   csv_in$rank1=csv_in[,list(rank=rank(dNME_UC_cor)),by=motif]$rank
#   csv_in$rank2=csv_in[,list(rank=rank(dNME_max_UC_pair)),by=motif]$rank
#   motif_dNME_out=rbind(motif_dNME_out,csv_in[,.SD[order(rank1+rank2,decreasing=T)],by=motif])
# 
# }
# write.csv(motif_dNME_out[,list(motif,region,gene,tissue,
#                                 UC_max_time ,
#                                 dNME_max_UC_pair,UC_max_pair,dMML_max_UC_pair,
#                                 dNME_UC_cor,dMML_UC_cor,
#                                 cluster,GO_terms)],'../downstream/output/mouse_analysis/motif_analysis/mouse_motif_example_GO_dNME.csv',
#           row.names = F)
# 
# 
# #Check the distribution in forebrain
# 
# forebrain_dMML=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/Forebrain_dMML_temp.csv')
# forebrain_dNME=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/Forebrain_dNME_temp.csv')
# 

# Reformat motif from Ken -------------------------------------------------

dir_Ken='../downstream/input/mouse_analysis/motif_analysis/mouse_motif_enrichment_0526/'
dMML_motifs=data.table()
for(fn in dir(dir_Ken,pattern='dMML')){
  ts=gsub('_.*','',fn)
  csv_in=fread(paste0(dir_Ken,fn))
  csv_in$tissue=ts
  csv_in$motif=gsub('.*_','',csv_in$motif)
  dMML_motifs=rbind(dMML_motifs,csv_in[,list(tissue,motif,human_high_NME,residual,log_OR_dMML,normalized_log_OR_dMML,FDR_dMML)])
}
write.csv(dMML_motifs,'../downstream/output/mouse_analysis/motif_analysis/Ken_dMML_all.csv')
dNME_motifs=data.table()
for(fn in dir(dir_Ken,pattern='dNME')){
  ts=gsub('_.*','',fn)
  csv_in=fread(paste0(dir_Ken,fn))
  csv_in$tissue=ts
  csv_in$motif=gsub('.*_','',csv_in$motif)
  dNME_motifs=rbind(dNME_motifs,csv_in[,list(tissue,motif,human_high_NME,residual,log_OR_dNME,normalized_log_OR_dNME,FDR_dNME)])
}
write.csv(dNME_motifs,'../downstream/output/mouse_analysis/motif_analysis/Ken_dNME_all.csv')
# # archived ----------------------------------------------------------------
# 
# 
# #dNME region
# region_sample=list()
# set.seed(123)
# for(ts in names(GO_out_all)){
# 
#   motif_dNME=Ken_motif_merge(readRDS(paste0('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_locus/',ts,'_motif_site_dNME.rds')))
#   if(!is.null(motif_dNME)){
#     GO_out_all_ts=GO_out_all[[ts]]
#     GO_out_all_ts_dNME=GO_out_all_ts[region_type=="dNME_only"]
# 
#     GO_out_all_ts_dNME$motif_dNME=convert_GR(motif_dNME,direction="DT")[match(GO_out_all_ts_dNME$region,region)]$motif
#     GO_out_all_ts_dNME$tissue=ts
#     GO_out_all_ts_dNME_motif=GO_out_all_ts_dNME[!is.na(motif_dNME)]
#     region_sample[[ts]]=GO_out_all_ts_dNME_motif[sample(1:nrow(GO_out_all_ts_dNME_motif),20)]
#   }
# }
# #dMML region
# region_sample=list()
# set.seed(123)
# for(ts in names(GO_out_all)){
# 
#   motif_dMML=Ken_motif_merge(readRDS(paste0('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_locus/',ts,'_motif_site_dMML.rds')))
#   if(!is.null(motif_dMML)){
#     GO_out_all_ts=GO_out_all[[ts]]
#     GO_out_all_ts_dMML=GO_out_all_ts[region_type=="dMML_only"]
# 
#     GO_out_all_ts_dMML$motif_dMML=convert_GR(motif_dMML,direction="DT")[match(GO_out_all_ts_dMML$region,region)]$motif
#     GO_out_all_ts_dMML$tissue=ts
#     GO_out_all_ts_dMML_motif=GO_out_all_ts_dMML[!is.na(motif_dMML)]
#     if(nrow(GO_out_all_ts_dMML_motif)>20){
#     region_sample[[ts]]=GO_out_all_ts_dMML_motif[sample(1:nrow(GO_out_all_ts_dMML_motif),20)]
#     }else(region_sample[[ts]]=GO_out_all_ts_dMML_motif)
#   }
# }
# write.csv(region_sample$heart[,list(region,UC_max_time,UC_max_UC_pair,dMML_max_UC_pair,dNME_max_UC_pair,motif_dMML,tissue)],
#           '../downstream/output/mouse_analysis/examples/dMML_heart_GO_motif.csv')
# Make Bedfile for Ken's motif binding site -------------------------------
Ken_motif_locus=readRDS(motif_Ken_fn)
for(stat_in in c("dMML","dNME")){
    for(tissue in names(Ken_motif_locus)){
      motif_sig=fread(paste0(Ken_motif_folder,tissue,'_OR_residual_',stat_in,'.csv'))
      bed_file_path=paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/',stat_in,'/',tissue,'/')
      motif_locus_in=Ken_motif_locus[[tissue]][motif_sig$motif]
      ifelse(!dir.exists(file.path(bed_file_path)), dir.create(file.path(bed_file_path)), FALSE)
      lapply(names(motif_locus_in),function(x) {export.bed(motif_locus_in[[x]],paste0(bed_file_path,
                                                                                      gsub("::","_",gsub('.*_','',x)),'.bed'))})
    
    }
}
# Selecct motif in factor book --------------------------------------------
#Get all motif binding site for motifs in Ken's list and factor book and find examples
#Need to have few brown SNPs 
#Need to have large examples
#Looking for genes related to those tissue
#Start with heart and forebrain
#Factorbook
# #CHip_atlas: GATA4 mouse embyro
# factor_in=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_heart_GATA4.bed',skip = 1,sep='\t')
# colnames(factor_in)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
# factor_in=factor_in[,list(seqnames,start,end,metadata,log10qval)]
# factor_in$metadata=gsub('%20','_',factor_in$metadata)
# factor_in$metadata=gsub('%3','_',factor_in$metadata)
# #gata4 in heart
# factor_in_GATA4=factor_in[grepl("GATA",factor_in$metadata)]
# factor_in_GATA4_gr=makeGRangesFromDataFrame(factor_in_GATA4)




# Check mosue ChiP-seq result ---------------------------------------------
enhancer_regions_motif_dNME_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dNME_all_V4.rds')
enhancer_regions_motif_dMML_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dMML_all_V4.rds')

#embyro Heart
factor_in_heart_embyro=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_embyro_heart_allTF.bed',skip = 1,sep='\t')
colnames(factor_in_heart_embyro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_heart_embyro=factor_in_heart_embyro[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_heart_embyro=factor_olap('heart',factor_in_heart_embyro,Ken_motif_folder,stage="embyro_heart")
ChIP_olap_heart_embyro=ChIP_olap(factor_olap_heart_embyro$factor_in_dNME,factor_olap_heart_embyro$factor_in_dMML,
                                 enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(enhancer_regions_heart_embyro$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_consistent_dNME.rds'))


#Adult heart
factor_in_heart=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_heart_all_TF.bed',skip = 1,sep='\t')
colnames(factor_in_heart)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_heart=factor_in_heart[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_heart_adult=factor_olap('heart',factor_in_heart,Ken_motif_folder,stage="adult")
ChIP_olap_heart_adult=ChIP_olap(factor_olap_heart_adult$factor_in_dNME,factor_olap_heart_adult$factor_in_dMML,
                                 enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(ChIP_olap_heart_adult$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_adult_consistent_dNME.rds'))

#embyro all:
factor_in_embyro=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_embyro_all_TF.bed',skip = 1,sep='\t')
colnames(factor_in_embyro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_embyro=factor_in_embyro[,list(seqnames,start,end,metadata,log10qval)]

#Heart embyro all
factor_olap_heart_embyro_all=factor_olap('heart',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_heart_embyro_all=ChIP_olap(factor_olap_heart_embyro_all$factor_in_dNME,factor_olap_heart_embyro_all$factor_in_dMML,
                                enhancer_regions_motif_dNME_all$heart,enhancer_regions_motif_dMML_all$heart)
saveRDS(ChIP_olap_heart_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_heart_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_heart_embyro_all_consistent_dMML.rds'))
#Forebrain embyro all
factor_olap_forebrain_embyro_all=factor_olap('forebrain',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_forebrain_embyro_all=ChIP_olap(factor_olap_forebrain_embyro_all$factor_in_dNME,factor_olap_forebrain_embyro_all$factor_in_dMML,
                                     enhancer_regions_motif_dNME_all$forebrain,enhancer_regions_motif_dMML_all$forebrain)
saveRDS(ChIP_olap_forebrain_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_forebrain_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_embyro_all_consistent_dMML.rds'))
#EFP embyro all
factor_olap_EFP_embyro_all=factor_olap('EFP',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_EFP_embyro_all=ChIP_olap(factor_olap_EFP_embyro_all$factor_in_dNME,factor_olap_EFP_embyro_all$factor_in_dMML,
                                         enhancer_regions_motif_dNME_all$EFP,enhancer_regions_motif_dMML_all$EFP)
saveRDS(ChIP_olap_EFP_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_EFP_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_EFP_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_EFP_embyro_all_consistent_dMML.rds'))
#limb embyro all
factor_olap_limb_embyro_all=factor_olap('limb',factor_in_embyro,Ken_motif_folder,stage="embyro_all")
ChIP_olap_limb_embyro_all=ChIP_olap(factor_olap_limb_embyro_all$factor_in_dNME,factor_olap_limb_embyro_all$factor_in_dMML,
                                   enhancer_regions_motif_dNME_all$limb,enhancer_regions_motif_dMML_all$limb)
saveRDS(ChIP_olap_limb_embyro_all$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_limb_embyro_all_consistent_dNME.rds'))
saveRDS(ChIP_olap_limb_embyro_all$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_limb_embyro_all_consistent_dMML.rds'))
#Adult neuro
factor_in_neuro=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_neuro_all_TF.bed',skip = 1,sep='\t')
colnames(factor_in_neuro)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
factor_in_neuro=factor_in_neuro[,list(seqnames,start,end,metadata,log10qval)]

factor_olap_forebrain_adult=factor_olap('forebrain',factor_in_embyro,Ken_motif_folder,stage="adult")
ChIP_olap_forebrain_adult=ChIP_olap(factor_olap_forebrain_adult$factor_in_dNME,factor_olap_forebrain_adult$factor_in_dMML,
                                         enhancer_regions_motif_dNME_all$forebrain,enhancer_regions_motif_dMML_all$forebrain)
saveRDS(ChIP_olap_forebrain_adult$dNME_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_adult_consistent_dNME.rds'))
saveRDS(ChIP_olap_forebrain_adult$dMML_region,paste0(ChiP_motif_dir,'enhancer_regions_forebrain_adult_consistent_dMML.rds'))

# #Liver
# #Embyro
# factor_in_liver=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_embyro_liver.bed',skip = 1,sep='\t')
# colnames(factor_in_liver)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
# factor_in_liver=factor_in_liver[,list(seqnames,start,end,metadata,log10qval)]
# factor_in_liver$metadata=gsub('%20','_',factor_in_liver$metadata)
# factor_in_liver$metadata=gsub('%3','_',factor_in_liver$metadata)
# forebrain_ken_dNME=fread(paste0(Ken_motif_folder,'liver_OR_residual_dNME.csv'))
# forebrain_ken_dMML = fread(paste0(Ken_motif_folder,'liver_OR_residual_dMML.csv'))
# factor_in_liver_dNME=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dNME$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# factor_in_liver_dMML=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dMML$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# 
# for(motif in gsub('.*_','',liver_ken_dNME$motif)){
#   motif_ChIP=factor_in_liver_dNME[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_liver_embyro_dNME_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }
# for(motif in gsub('.*_','',liver_ken_dMML$motif)){
#   motif_ChIP=factor_in_liver_dMML[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dMML/mm10_liver_embyro_dMML_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }
# #Adult
# factor_in_liver=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_liver.bed',skip = 1,sep='\t')
# colnames(factor_in_liver)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
# factor_in_liver=factor_in_liver[,list(seqnames,start,end,metadata,log10qval)]
# factor_in_liver$metadata=gsub('%20','_',factor_in_liver$metadata)
# factor_in_liver$metadata=gsub('%3','_',factor_in_liver$metadata)
# liver_ken_dNME=fread('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_result_20210406/liver_OR_residual_dNME.csv')
# liver_ken_dMML = fread('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_result_20210406/liver_OR_residual_dMML.csv')
# factor_in_liver_dNME=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dNME$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# factor_in_liver_dMML=factor_in_liver[grepl(paste(gsub('.*_','',liver_ken_dMML$motif),collapse = "|"),factor_in_liver$metadata,ignore.case=T)]
# 
# for(motif in gsub('.*_','',liver_ken_dNME$motif)){
#   motif_ChIP=factor_in_liver_dNME[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_liver_adult_dNME_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }
# for(motif in gsub('.*_','',liver_ken_dMML$motif)){
#   motif_ChIP=factor_in_liver_dMML[grepl(motif,metadata,ignore.case = T)]
#   if(nrow(motif_ChIP)>0){
#     write.table(motif_ChIP,
#                 paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dMML/mm10_forebrain_adult_dMML_ChiPatlas_',gsub('::','_',motif),'.bed'),
#                 col.names = F,row.names = F,quote=F)
#   }
#   
# }


# Find examples by looking at overlap between motif and ChIP data -----------
enhancer_regions_motif_dNME_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dNME_all.rds')
enhancer_regions_motif_dMML_all=readRDS('../downstream/output/mouse_analysis/motif_analysis/enhancer_regions_motif_dMML_all.rds')
#Heart 
factor_in_embyro_heart_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_embyro_heart_dNME.rds')))
factor_in_heart_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_heart_dNME.rds')))
heart_dNME_pubmed=enhancer_regions_motif_dNME_all$heart
heart_dNME_pubmed=heart_dNME_pubmed[,list(region,tissue,cluster,gene,UC_max_time,region_type,PMID,dNME_max_pair,dMML_max_pair,UC_max_pair,dMML_motif,dNME_motif)]
olap=findOverlaps(convert_GR(heart_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_embyro_heart_dNME))
heart_dNME_pubmed_motif=heart_dNME_pubmed[queryHits(olap)][order(dNME_max_pair,decreasing=T)]

olap=findOverlaps(convert_GR(heart_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_heart_dNME))
heart_dNME_pubmed_motif=heart_dNME_pubmed[queryHits(olap)][order(dNME_max_pair,decreasing=T)]
#Forebrain
factor_in_embyro_forebrain_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_embyro_forebrain_dNME.rds')))
factor_in_forebrain_dNME=readRDS(paste0(paste0(ChiP_motif_dir,'factor_in_adult_forebrain_dNME.rds')))
forebrain_dNME_pubmed=enhancer_regions_motif_dNME_all$forebrain
forebrain_dNME_pubmed=forebrain_dNME_pubmed[,list(region,tissue,cluster,gene,UC_max_time,region_type,PMID,dNME_max_pair,dMML_max_pair,UC_max_pair,dMML_motif,dNME_motif)]
olap=findOverlaps(convert_GR(forebrain_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_embyro_forebrain_dNME))
forebrain_dNME_pubmed_motif=forebrain_dNME_pubmed[queryHits(olap)]
forebrain_dNME_pubmed_motif$motif_ChIP=factor_in_embyro_forebrain_dNME[subjectHits(olap)]$metadata
forebrain_dNME_pubmed_motif$motif_in_ChIP=forebrain_dNME_pubmed_motif[,list(motif_in_ChIP=grepl(gsub(';','|',dNME_motif),motif_ChIP,ignore.case = T)),by=list(region,motif_ChIP)]$motif_in_ChIP

olap=findOverlaps(convert_GR(forebrain_dNME_pubmed$region),makeGRangesFromDataFrame(factor_in_forebrain_dNME))
forebrain_dNME_pubmed_motif=forebrain_dNME_pubmed[queryHits(olap)]
forebrain_dNME_pubmed_motif$motif_ChIP=factor_in_forebrain_dNME[subjectHits(olap)]$metadata
forebrain_dNME_pubmed_motif$motif_in_ChIP=forebrain_dNME_pubmed_motif[,list(motif_in_ChIP=grepl(gsub(';','|',dNME_motif),motif_ChIP,ignore.case = T)),by=list(region,motif_ChIP)]$motif_in_ChIP
# archived ----------------------------------------------------------------



liver_regions=fread('../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dNME.csv')
olap=findOverlaps(convert_GR(liver_regions$region),GRanges(seqname=factor_in_liver_dNME$seqnames,IRanges(start=factor_in_liver_dNME$start,end=factor_in_liver_dNME$end)))

write.csv(liver_regions[queryHits(olap)],'../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dNME_ChiP.csv')

liver_regions=fread('../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dMML.csv')
olap=findOverlaps(convert_GR(liver_regions$region),GRanges(seqname=factor_in_liver_dMML$seqnames,IRanges(start=factor_in_liver_dMML$start,end=factor_in_liver_dMML$end)))

write.csv(liver_regions[queryHits(olap)],'../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_liver_dMML_ChiP.csv')


#Subset the motif
#Ken motif
heart_Ken_binding_site=readRDS('../downstream/input/mouse_analysis/motif_analysis/Ken_motif_binding_site/tissue_region_motif_all.rds')
heart_enhancer_GO_motif=

olap_motif=findOverlaps(convert_GR(region_in_enhancer_GO$heart$region),factor_in_GATA4_gr)

unique(region_in_enhancer$heart[queryHits(olap_motif)])[order(dNME_max_pair,decreasing = T)]
unique(region_in_enhancer_GO$heart[queryHits(olap_motif)][grepl("heart|cardi|aort|ventricu|vascula|muscle|actin",GO_result),
                                   list(region,cluster,dNME_max_pair,dNME_max_time,UC_max_pair,dMML_max_pair,gene),with=T])[order(dNME_max_pair,decreasing=T)]




#Examples for correlation
both_example = data.frame(value=as.numeric(UC_merge$heart["chr4:97739493-97739742",!grepl("max",colnames(UC_merge$heart))]),
                          stage_stat=colnames(UC_merge$heart["chr4:97739493-97739742",!grepl("max",colnames(UC_merge$heart))]))

both_example$stat=gsub('-.*','',both_example$stage_stat)
both_example$stage=sub('dMML-|dNME-|UC-','',both_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(both_example[both_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()

dNME_only_example = data.frame(value=as.numeric(UC_merge$heart["chr10:80117344-80117593",!grepl("max",colnames(UC_merge$heart))]),
                               stage_stat=colnames(UC_merge$heart["chr10:80117344-80117593",!grepl("max",colnames(UC_merge$heart))]))

dNME_only_example$stat=gsub('-.*','',dNME_only_example$stage_stat)
dNME_only_example$stage=sub('dMML-|dNME-|UC-','',dNME_only_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(dNME_only_example[dNME_only_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()

dMML_only_example = data.frame(value=as.numeric(UC_merge$heart["chr10:62313679-62313928",!grepl("max",colnames(UC_merge$heart))]),
                               stage_stat=colnames(UC_merge$heart["chr10:62313679-62313928",!grepl("max",colnames(UC_merge$heart))]))

dMML_only_example$stat=gsub('-.*','',dMML_only_example$stage_stat)
dMML_only_example$stage=sub('dMML-|dNME-|UC-','',dMML_only_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(dMML_only_example[dMML_only_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()


neither_example = data.frame(value=as.numeric(UC_merge$heart["chr1:151764174-151764423",!grepl("max",colnames(UC_merge$heart))]),
                             stage_stat=colnames(UC_merge$heart["chr1:151764174-151764423",!grepl("max",colnames(UC_merge$heart))]))

neither_example$stat=gsub('-.*','',neither_example$stage_stat)
neither_example$stage=sub('dMML-|dNME-|UC-','',neither_example$stage_stat)
selected_stages=paste0("E",10:15,".5-E",11:16,'.5')
ggplot(neither_example[neither_example$stage%in%selected_stages,],aes(x=stage,y=value,group=stat,color=stat))+geom_line()+geom_point()
#IGV examples

neuro_motif_ChIP=fread('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_forebrain_adult_selected_dNME_ChiPatlas.bed')
forebrain_motif=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/forebrain_dNME_Pubmed_annotated_region_only.csv')
olap=findOverlaps(convert_GR(forebrain_motif$region),makeGRangesFromDataFrame(neuro_motif_ChIP,seqnames.field = "V1",start.field = "V2",end.field = "V3"))
forebrain_motif=forebrain_motif[unique(queryHits(olap))][order(dNME_max_UC_pair,decreasing = T)]
write.csv(forebrain_motif,'../downstream/output/mouse_analysis/motif_analysis/region_motif/forebrain_dNME_Pubmed_annotated_region_only_ChIp_olap.csv')

heart_motif_ChIP=fread('../downstream/input/mouse_analysis/motif_analysis/chipatlas/mm10_heart_GATA4.bed',skip = 1,sep='\t')
heart_motif_ChIP=fread('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/dNME/mm10_heart_embyro_selected_dNME_ChiPatlas.bed',sep=' ')
heart_motif=fread('../downstream/output/mouse_analysis/motif_analysis/region_motif/heart_motif_region_only.csv')
olap=findOverlaps(convert_GR(heart_motif$region),makeGRangesFromDataFrame(heart_motif_ChIP,seqnames.field = "V1",start.field = "V2",end.field = "V3"))
heart_motif=heart_motif[unique(queryHits(olap))][order(dNME_max_UC_pair,decreasing = T)]
write.csv(heart_motif,'../downstream/output/mouse_analysis/motif_analysis/region_motif/heart_dNME_Pubmed_annotated_region_only_ChIp_olap2.csv')


