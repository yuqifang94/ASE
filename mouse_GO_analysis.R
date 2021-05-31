rm(list=ls())
source('mainFunctions_sub.R')
# GO analysis -------------------------------------------------------------
#Define input parameters
dir_out_cluster='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_1/'
folder_out=paste0(dir_out_cluster,'region_assigned/')
dir_out_GO='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/'
UC_merge=readRDS('../downstream/input/mouse_analysis/UC_only_all_regions.rds')#Define all analyzed regions, were using UC_merge_max_loc_cluster01.rds,4626
cutoff_fn='01'
#Runnning
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
#prepare enhancer background gene list
uc_gr=lapply(UC_merge,function(x) rownames(x))
uc_gr=Reduce(intersect,uc_gr)
uc_gr=convert_GR(uc_gr)
enhancer=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)
#Prepare promoter background gene
tss=get_mm10_tss()

bg_promoter=names(subsetByOverlaps(tss,uc_gr,maxgap = 2000))

# #filter by repeats
# 
# dir_in='../downstream/input/mouse_analysis/clustering/tissue_specific/currently_in_use/ts_cluster_0_1_non_repeats/region_assigned/'
# re_web=readRDS('../downstream/output/mouse_analysis/repeats/re_web.rds')
# cor_dt_filtered=readRDS('../downstream/output/mouse_analysis/correlation/correlation_dt_N17_kmeans_10run_filtered_all_regions_non_repeats.rds')
# names(cor_dt_filtered)=NULL
# olap_enhancer=findOverlaps(enhancer,re_web,minoverlap = mean(width(convert_GR(do.call('rbind',cor_dt_filtered)$region))))
# enhancer_bg=subsetByOverlaps(enhancer[unique(queryHits(olap_enhancer))],uc_gr)
# bg_enhancer=unique(enhancer_bg$`Target Gene`)
# olap_tss=findOverlaps(enhancer,re_web,maxgap = 2000)
# bg_promoter=names(subsetByOverlaps(tss[unique(queryHits(olap_tss))],uc_gr,maxgap = 2000))
#GO run

for(enc_type in c("enhancer","promoter")){
  enc_type="promoter"
  GO_out_all=list()
  for(region_type in c("all","NME only","Neither","Both","MML only")){
    GO_out_all[[region_type]]=list()
    for(ts in tissue_all){
      if(enc_type=="enhancer"){
        bg=bg_enhancer
      }else 
        if(enc_type=="promoter"){
          bg=bg_promoter
          
        }
      GO_out_all[[region_type]][[ts]]=GO_run_tissue(ts,folder_out,enc_type=enc_type,region_type_sel=region_type,bg=bg,DNase=F)
      GO_out_all[[region_type]][[ts]]=lapply(GO_out_all[[region_type]][[ts]],function(x){
        return(list(GO_out_cluster_all=x$GO_out_cluster_all,
                    csv_in_ts_clu=cbind(x$csv_in_ts_clu,as.data.table(UC_merge[[ts]][x$csv_in_ts_clu$region,!grepl('max',colnames(UC_merge[[ts]]))]))))
        
        
      })
      
    }
  }
  # UC_maxUC is wrong
  saveRDS(GO_out_all,paste0(dir_out_GO,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_',enc_type,'.rds'))
}

# Plot heatmaps -----------------------------------------------------------

enc_type='enhancer'
GO_out_all=readRDS(paste0(dir_out_GO,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_',enc_type,'.rds'))
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
#Plot all terms in a single plot
for(region_type in names(GO_out_all)){
  plot_GO_heatmap_all(tissue_all,GO_out_all[[region_type]],region_type=region_type,enc_type="enhancer",ptcount=0,FDR_cutoff=0.2,
                      dir_plot=paste0('../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/UC_',cutoff_fn,'/'))
  
}
enc_type='promoter'
cutoff_fn='01'
GO_out_all=readRDS(paste0(dir_out_GO,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_',enc_type,'.rds'))
for(region_type in names(GO_out_all)){
  plot_GO_heatmap_all(tissue_all,GO_out_all[[region_type]],region_type=region_type,enc_type="promoter",ptcount=0,FDR_cutoff=0.2,
                      dir_plot=paste0('../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/UC_',cutoff_fn,'/'))
  
}


# Write csv output --------------------------------------------------------

chrs <- names(Mmusculus)[1:21]#2276796
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
GO_sheets(GO_out_all,"enhancer",dMML_cor=dMML_cor,dNME_cor=dNME_cor,mm10_CpG=cpgr,FDR_cutoff = 0.2,out_dir='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/GO_sheets/')
GO_sheets(GO_out_all,"promoter",dMML_cor=dMML_cor,dNME_cor=dNME_cor,mm10_CpG=cpgr,FDR_cutoff = 0.2,out_dir='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/GO_sheets/')

# # Write motif analysis result ---------------------------------------------
# motif_dir='../downstream/input/mouse_analysis/motif_analysis/mouse_motif_enrichment_0510/'
# motif_all=data.table()
# for(fn in dir(motif_dir,pattern='csv')){
#   stat=gsub('.csv','',gsub('.*_','',fn))
#   csv_in=fread(paste0(motif_dir,fn))
#   csv_in$motif=gsub('.*_','',csv_in$motif)
#   csv_out=csv_in[,grepl("motif|FDR",colnames(csv_in)),with=F]
#   colnames(csv_out)=c('motif','FDR')
#   csv_out$stat=stat
#   csv_out$tissue=gsub('_.*','',fn)
#   motif_all=rbind(motif_all,csv_out)
#   }
# write.csv(motif_all[FDR<=0.1,list(tissue,motif,FDR,stat)],'../downstream/output/mouse_analysis/motif_analysis/motif_all_Ken.csv')
