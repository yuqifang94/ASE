source('mainFunctions_sub.R')
UC_in_matrix_ls_red_ft_fn='../downstream/output/mouse_analysis/clustering/UC_in_matrix_ls_red_ft.rds'
order_clu<-function(UC_in_mt,clu){
  clum <- rowsum(UC_in_mt,clu)/as.vector(table(clu))
  n <- names(clu)
  maxp <- apply(clum,1,function(i) {
    i[i < 0] <- 0
    i <- i-min(i)
    i <- i/max(i)
    sum(i*c(1:length(i)))/sum(i)
  })
  clu <- rank(maxp)[clu]
  names(clu) <- n
  return(clu)
}


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
saveRDS(UC_in_matrix_ls_red_ft,UC_in_matrix_ls_red_ft_fn)

UC_in_matrix_ls_red_ft=readRDS(UC_in_matrix_ls_red_ft_fn)
for(i in 1:10){
    set.seed(i)
    cat("Processing ",i, 'literation\n')
    #https://stackoverflow.com/questions/21382681/kmeans-quick-transfer-stage-steps-exceeded-maximum
    clu <- order_clu(UC_in_matrix_ls_red_ft,kmeans(UC_in_matrix_ls_red_ft,70,iter.max = 10000)$cluster)
    saveRDS(clu,paste0('../downstream/output/mouse_analysis/clustering/UC_cluster_all_time_',i,'.rds'))
    #2 runs remains
}

library(RColorBrewer)
library(pheatmap)
library(gplots)

  #clu=readRDS(cluster_region_out_fn)
  set.seed(123)
  clu <- kmeans(UC_in_matrix_ls_red_ft,70,iter.max = 10000)$cluster


  


  #This is reorder clu to the order of maximum UC in column

  saveRDS(clu,paste0('../downstream/output/mouse_analysis/clustering/UC_cluster_all_time_',123,'.rds'))
  
  #Local run
  UC_in_matrix_ls_red_ft=readRDS('../downstream/output/mouse_analysis/clustering/UC_in_matrix_ls_red_ft.rds')
  clu=readRDS(paste0('../downstream/output/mouse_analysis/clustering/UC_cluster_all_time_',123,'.rds'))
  clu=order_clu(UC_in_matrix_ls_red_ft,clu)
  #sub_sample=sample(1:length(clu),round(length(clu)/100))
  #clu=sort(clu[sub_sample])
  clu=sort(clu)
  #GO analysis
  #Background
  bg=convert_GR(names(clu))
  enhancer=readRDS(bin_enhancer_rds)
  bg=enhancer$`Target Gene`
  cluster_GO=list()
  for(i in 1:70){
  clu_i_GR= convert_GR(names(clu[clu==i]))
  enhancer_clu_i=subsetByOverlaps(enhancer,clu_i_GR)
  
  cluster_GO[[i]]=GO_run(unique(enhancer_clu_i$`Target Gene`),bg,cluster=i)
  }
  saveRDS(cluster_GO,'../downstream/output/mouse_analysis/GO_analysis/cluster_all_regions_GO.rds')
  cluster_GO=readRDS('../downstream/output/mouse_analysis/GO_analysis/cluster_all_regions_GO.rds')
#Plotting
  mat_out=UC_in_matrix_ls_red_ft[names(clu),]
  

   rowann_out <- data.frame(cluster=clu,stringsAsFactors = F)
   #row_gap=cumsum(rle(as.numeric(clu))$lengths) 
  #Refine plotting parameters
  colann <- data.frame(time=sub('-all','',sub('.*?-','',colnames(mat_out))),tissue=sub('-.*','',colnames(mat_out)),stringsAsFactors = F)
  rownames(colann) <- colnames(mat_out)
  c1 <- mouse_color()
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
library(gplots)

  c2 <- col2hex(sample(color, 70))
  names(c2) <- 1:70
  c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
  names(c4) <- sort(unique(colann[,1]))
  #remove row with all NA 
  mat_out_sc=scalematrix(mat_out)
 # jpeg(paste0(figure_path,'clustering_all_regions_01.jpg'))
  pheatmap(mat_out_sc,cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
           annotation_col = colann,show_colnames = F,show_rownames = F,
           gaps_col = cumsum(rle(colann[,2])$lengths),
           annotation_colors = list(tissue=c1,cluster=c2,time=c4
                                    #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)
           ),
           filename=paste0(figure_path,'clustering_all_regions_01.jpeg'))
 # dev.off()
  for(i in 1:70){
    pheatmap(mat_out_sc[names(clu[clu==i]),],cluster_rows=F,cluster_cols=F,show_rownames = F,
             filename=paste0(figure_path,'All_region_70_cluster/clustering_all_regions_01_C',i,'.jpeg'),
             annotation_col = colann, gaps_col = cumsum(rle(colann[,2])$lengths), annotation_colors = list(tissue=c1,time=c4))
  }