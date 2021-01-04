source('mainFunctions_sub.R')
chromHMM_only=readRDS('../downstream/input/cluster/chromHMM_enhancer.rds')
all=readRDS('../downstream/input/cluster/jsd.rds')
dir_in='../downstream/input/mm10_cluster_all/'

dir_out='../downstream/input/mm10_cluster_chromHMM/'
conf_matrix_all=list()
for(ts in names(chromHMM_only)){
  chromHMM_only_in=chromHMM_only[[ts]]
  all_in=all[[ts]]
  conf_matrix=data.table(region=names(chromHMM_only_in),cluster_chromHMM=chromHMM_only_in,
                         cluster_all=all_in[names(chromHMM_only_in)])
  conf_matrix_all[[ts]]=dcast.data.table(conf_matrix,cluster_chromHMM~cluster_all)
}

saveRDS(conf_matrix_all,'../downstream/output/conf_matrix_all.rds')
for(ts in names(chromHMM_only)){
  csv_in_ts=fread(paste0(dir_in,ts,".csv"))
  chromHMM_only_in=chromHMM_only[[ts]]
  csv_in_ts=csv_in_ts[chromHMM_enhancer ==TRUE]
  csv_in_ts$cluster_all=csv_in_ts$cluster
  csv_in_ts$cluster=chromHMM_only_in[csv_in_ts$region]
  write.csv(csv_in_ts,paste0(dir_out,ts,'.csv'))
}



# csv_in_ts[region%in%names(chromHMM_only$forebrain)]
# chromHMM=readRDS('../downstream/output/chromHMM_enhancer.rds')
# forebrain_enhancer=convert_GR(names(chromHMM_only$forebrain))
# subsetByOverlaps(forebrain_enhancer,chromHMM[chromHMM$tissue=='forebrain'])
# UC_raw=readRDS('../downstream/output/uc_matrix_DNase.rds')
# UC_raw=UC_raw$forebrain
# timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
# UC_raw=UC_raw[,colnames(UC_raw)%in%timeorder]
# chromHMM_only_forebrain=sort(chromHMM_only$forebrain)
# chromHMM_only_forebrain=chromHMM_only_forebrain[names(chromHMM_only_forebrain) %in% rownames(UC_raw)]
# UC=UC_raw[names(chromHMM_only_forebrain),]
# scalematrix <- function(data) {
#   cm <- rowMeans(data)
#   csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
#   (data - cm) / csd
# }
# pheatmap(scalematrix(UC),cluster_rows = F,cluster_cols = F,scale = "none",show_colnames = T,show_rownames = F,file='test.pdf')
