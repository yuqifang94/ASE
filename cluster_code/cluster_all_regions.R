source('mainFunctions_sub.R')
UC_in_matrix_ls=readRDS(UC_in_matrix_ls_file)
UC_in_matrix_ls_mt=lapply(UC_in_matrix_ls,convert_GR,direction='matrix')
all_analyzed_regions=Reduce(intersect,lapply(UC_in_matrix_ls,function(x) rownames(x)))

UC_in_matrix_ls_red=mclapply(UC_in_matrix_ls_mt,
        function(x) {
            
            x[all_analyzed_regions,]#?for first occurance
        }
        
        
,mc.cores=7)
UC_in_matrix_ls_red=do.call(cbind,UC_in_matrix_ls_red)#4989824
UC_in_matrix_ls_red_ft=UC_in_matrix_ls_red[rowSums(is.na(UC_in_matrix_ls_red))==0,]#4876367
UC_cutoff=0.1
UC_in_matrix_ls_red_ft=UC_in_matrix_ls_red_ft[rowSums(UC_in_matrix_ls_red_ft>UC_cutoff)>0,]#2085884
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
UC_in_matrix_ls_red_ft=UC_in_matrix_ls_red_ft[,sub('-all','',sub('.*?-','',colnames(UC_in_matrix_ls_red))) %in% timeorder]
saveRDS(UC_in_matrix_ls_red_ft,'../downstream/output/mouse_analysis/clustering/UC_in_matrix_ls_red_ft.rds')
for(i in 1:10){
    set.seed(i)
    cat("Processing ",i, 'literation\n')
    #https://stackoverflow.com/questions/21382681/kmeans-quick-transfer-stage-steps-exceeded-maximum
    clu <- kmeans(UC_in_matrix_ls_red_ft,70,iter.max = 10000,algorithm="MacQueen")$cluster
    saveRDS(clu,paste0('../downstream/output/mouse_analysis/clustering/UC_cluster_all_time_',i,'.rds'))
    #2 runs remains
}

library(RColorBrewer)
library(pheatmap)
library(gplots)

  #clu=readRDS(cluster_region_out_fn)
  set.seed(123)
  clu <- kmeans(UC_in_matrix_ls_red_ft,70,iter.max = 10000,algorithm="MacQueen")$cluster
  saveRDS(clu,paste0('../downstream/output/mouse_analysis/clustering/UC_cluster_all_time_',123,'.rds'))
  clu=readRDS(paste0('../downstream/output/mouse_analysis/clustering/UC_cluster_all_time_',123,'.rds'))
  clu=sort(clu)
  UC_in_matrix_ls_red_ft=readRDS('../downstream/output/mouse_analysis/clustering/UC_in_matrix_ls_red_ft.rds')
  mat_out=UC_in_matrix_ls_red_ft[names(clu),]
   rowann_out <- data.frame(cluster=clu,stringsAsFactors = F)
   row_gap=cumsum(table(clu) ) 
  #Refine plotting parameters
  colann <- data.frame(time=sub('-all','',sub('.*?-','',colnames(mat_out))),tissue=sub('-.*','',colnames(mat_out)),stringsAsFactors = F)
  rownames(colann) <- colnames(mat_out)
  c1 <- mouse_color()
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  c2 <- sample(color, 70)
  names(c2) <- 1:70
  c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
  names(c4) <- sort(unique(colann[,1]))
  #remove row with all NA 
  
 # jpeg(paste0(figure_path,'clustering_all_regions_01.jpg'))
  pheatmap(scalematrix(mat_out),cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
           annotation_col = colann,show_colnames = F,show_rownames = F,
           gaps_row = row_gap,gaps_col = cumsum(rle(colann[,2])$lengths),
           annotation_colors = list(tissue=c1,cluster=c2,time=c4
                                    #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)
           ),
           filename=paste0(figure_path,'clustering_all_regions_01.jpeg'))
 # dev.off()