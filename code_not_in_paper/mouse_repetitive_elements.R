# re_UCSC=import.bed('../downstream/input/repetitive_element/mm10_repeat_masker_UCSC.bed')
# re_UCSC=re_UCSC[seqnames(re_UCSC)%in%paste0("chr",1:19)]
re_web_in=fread('../downstream/input/mouse_analysis/repetitive_element/mm10_repeat_masker_web.out')
re_web=re_web_in[V5 %in%paste0("chr",1:19),list(V5,V6,V7,V10,V11)]
re_web=makeGRangesFromDataFrame(re_web,keep.extra.columns=T,seqnames.field = "V5",start.field = "V6",
                                end.field = "V7")
colnames(mcols(re_web))=c("repeat_type","repeat_class")
#sum(width(unique(re_web))): 38.5% of mouse genome (2730871774)
#filter out RNA related repats
re_web=re_web[!grepl("RNA",re_web$repeat_class)]
re_web=re_web[!re_web$repeat_class%in%c("Other","Unknown")]
#16k Simplerepeats
re_web$repeat_type[re_web$repeat_class=="Simple_repeat"]="Simple_repeat"
re_web_class_type=as.data.table(unique(mcols(re_web)))
saveRDS(re_web,'../downstream/output/mouse_analysis/repeats/re_web.rds')
#Looking for enrichment of repetitive element
#Among all the regions we analyzed, how many are in UC sig list, how many are in repetitive
uc=readRDS('../downstream/output/uc_matrix_DNase.rds')
#All regions
all_regions=unique(do.call(c,lapply(uc,rownames)))
all_regins_gr=convert_GR(all_regions)
folder_input="../downstream/input/ts_cluster_0_1/"
enhancer=readRDS('../downstream/output/bin_enhancer.rds')
repeats_output=list()
for(csv_in in dir(folder_input,pattern="csv")){
  cat("Processing",csv_in,'\n')
  tt1=proc.time()[[3]]
  repeats_output_ts=data.table()
  ts_in=fread(paste0(folder_input,csv_in))
  ts_region=convert_GR(ts_in$region)
  #Looking at enhancers
  #ts_region=subsetByOverlaps(ts_region,enhancer)
  
  non_ts_region=convert_GR(all_regions[!(all_regions%in%convert_GR(ts_region,direction="DT"))])
  ts_region=ts_region[seqnames(ts_region)%in%paste0("chr",1:19)]
  non_ts_region=non_ts_region[seqnames(non_ts_region)%in%paste0("chr",1:19)]
  tissue= gsub(".csv","",csv_in)
  repeats_output[[tissue]]=repeat_olap_individual(ts_region,non_ts_region,re_web,re_web_class_type,tissue)
  gc()
  cat("Finish processing in ",proc.time()[[3]]-tt1,'\n')
}
saveRDS(repeats_output,'../downstream/output/mouse_analysis/repeats/repeats_output.rds')

repeats_output=readRDS('../downstream/output/mouse_analysis/repeats/repeats_output.rds')
repeats_output=lapply(repeats_output,function(x) {
  x$percent_overlap_UC=x$observed/(x$observed+x$ts_region_non_repeat)
  return(x)}
  
  )
repeats_output_sig_percent=lapply(repeats_output,function(x) {
  sum(x[FDR<=0.1]$observed )/(x[1]$observed+x[1]$ts_region_non_repeat )
  
})
mean(unlist(repeats_output_sig_percent),na.rm=T)#0.006 = 0.6%
write.csv(do.call(rbind,repeats_output)[order(FDR)],'../downstream/output/mouse_analysis/repeats/all_repeat_overlap.csv')

#Average size of the repeats
re_web$repeat_size=width(re_web)
re_web_dt=as.data.table(mcols(re_web))
re_web_dt_size=re_web_dt[,list(mean_size=mean(repeat_size),sd_size=sd(repeat_size),repeat_class=unique(repeat_class)),by=list(repeat_type)]

repeats_output=lapply(repeats_output,function(x){
  x$mean_size=re_web_dt_size[match(x$repeat_type,repeat_type)]$mean_size
  x$sd_size=re_web_dt_size[match(x$repeat_type,repeat_type)]$sd_size
  return(x)
  
})
write.csv(do.call(rbind,lapply(repeats_output,function(x) x[order(FDR,decreasing = F)])),'../downstream/output/repeat_overlap.csv')
#Do it for each tissue, all repeats
repeats_output_ts_allrep=data.table()

for(csv_in in dir(folder_input,pattern="csv")){
  cat("Processing",csv_in,'\n')
  tt1=proc.time()[[3]]
  repeats_output_ts=data.table()
  ts_in=fread(paste0(folder_input,csv_in))
  ts_region=convert_GR(ts_in$region)
  #Looking at enhancers
  #ts_region=subsetByOverlaps(ts_region,enhancer)
  
  non_ts_region=convert_GR(all_regions[!(all_regions%in%convert_GR(ts_region,direction="DT"))])
  ts_region=ts_region[seqnames(ts_region)%in%paste0("chr",1:19)]
  non_ts_region=non_ts_region[seqnames(non_ts_region)%in%paste0("chr",1:19)]
  
  
  repeats_output_ts_allrep=rbind(repeats_output_ts_allrep,repeat_olap_all(ts_region,non_ts_region,re_web))
  gc()
  cat("Finish processing in ",proc.time()[[3]]-tt1,'\n')
}
repeats_output_ts_allrep$FC=(repeats_output_ts_allrep$ts_repeat/repeats_output_ts_allrep$ts_non_repeat)/
  (repeats_output_ts_allrep$non_ts_repeat/repeats_output_ts_allrep$non_ts_non_repeat)

# Human analysis ----------------------------------------------------------

re_web_in=fread('../downstream/input/human_analysis//Repeats/hg19.fa.out')
re_web=re_web_in[V5 %in%paste0("chr",1:19),list(V5,V6,V7,V10,V11)]
re_web=makeGRangesFromDataFrame(re_web,keep.extra.columns=T,seqnames.field = "V5",start.field = "V6",
                                end.field = "V7")
colnames(mcols(re_web))=c("repeat_type","repeat_class")
#sum(width(unique(re_web))): 38.5% of mouse genome (2730871774)
#filter out RNA related repats
re_web=re_web[!grepl("RNA",re_web$repeat_class)]
re_web=re_web[!re_web$repeat_class%in%c("Other","Unknown")]
#16k Simplerepeats
re_web$repeat_type[re_web$repeat_class=="Simple_repeat"]="Simple_repeat"
re_web_class_type=as.data.table(unique(mcols(re_web)))
GR_merge=readRDS(GR_merge_file)

#dMML region
repeats_output_dMML=repeat_olap_individual(GR_merge[GR_merge$dMML_pval<=pval_cutoff],GR_merge[GR_merge$dMML_pval>pval_cutoff],
                                           re_web,re_web_class_type,'dMML',ncores=10)
saveRDS(repeats_output_dMML,'../downstream/output/human_analysis/repeats/repeats_output_dMML.rds')
#dNME region
repeats_output_dNME=repeat_olap_individual(GR_merge[GR_merge$dNME_pval<=pval_cutoff],GR_merge[GR_merge$dNME_pval>pval_cutoff],
                                           re_web,re_web_class_type,'dNME',ncores=10)
saveRDS(repeats_output_dNME,'../downstream/output/human_analysis/repeats/repeats_output_dNME.rds')
repeats_output_dMML_all=repeat_olap_all(GR_merge[GR_merge$dMML_pval<=pval_cutoff],GR_merge[GR_merge$dMML_pval>pval_cutoff],re_web)
repeats_output_dNME_all=repeat_olap_all(GR_merge[GR_merge$dNME_pval<=pval_cutoff],GR_merge[GR_merge$dNME_pval>pval_cutoff],re_web)
repeats_output_dNME_all$tissue=NULL
repeats_output_dNME_all$proportion=repeats_output_dNME_all$ts_repeat/repeats_output_dNME_all$total_ts_reion
repeats_output_dNME_all$FC=(repeats_output_dNME_all$ts_repeat/repeats_output_dNME_all$ts_non_repeat)/
  (repeats_output_dNME_all$non_ts_repeat/repeats_output_dNME_all$non_ts_non_repeat)
repeats_output_dMML_all$tissue=NULL
repeats_output_dMML_all$proportion=repeats_output_dMML_all$ts_repeat/repeats_output_dMML_all$total_ts_reion
repeats_output_dMML_all$FC=(repeats_output_dMML_all$ts_repeat/repeats_output_dMML_all$ts_non_repeat)/
  (repeats_output_dMML_all$non_ts_repeat/repeats_output_dMML_all$non_ts_non_repeat)

repeats_output_dNME=readRDS('../downstream/output/human_analysis/repeats/repeats_output_dNME.rds')
write.csv(repeats_output_dNME,'../downstream/output/human_analysis/repeats/repeats_output_dNME_minolap.csv')
repeats_output_dMML=readRDS('../downstream/output/human_analysis/repeats/repeats_output_dMML.rds')
write.csv(repeats_output_dMML,'../downstream/output/human_analysis/repeats/repeats_output_dMML_min_olap.csv')


# Human analysis ----------------------------------------------------------


