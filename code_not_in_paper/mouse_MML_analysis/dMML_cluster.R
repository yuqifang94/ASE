source('mainFunctions_sub.R')
cut=0.134375#From tissue-specific region counting
UC_in=readRDS(UC_merge_file)
d=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=20)
# lapply(d,nrow)
# $EFP
# [1] 4876367
# $forebrain
# [1] 4876367
# $heart
# [1] 4876367
# $hindbrain
# [1] 4876367
# $limb
# [1] 4876367

# $liver
# [1] 4876367
# $midbrain
# [1] 4876367
aid <- sapply(names(d),function(i) {
    names(which(rowSums(d[[i]] > cut) > 0))
  })  
#lapply(aid,length)
# $EFP
# [1] 935740
# $forebrain
# [1] 1165702
# $heart
# [1] 994001
# $hindbrain
# [1] 1132357
# $limb
# [1] 970433
# $liver
# [1] 1093959
# $midbrain
# [1] 1113027
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
for (seed in 1:10) {
  cat('Processing:',seed,'\n')
  cluster_d <- sapply(names(d),function(i) {
    #This is for one and only one
    #sid <- setdiff(aid[[i]],unlist(aid[names(aid)!=i]))
    sid=aid[[i]]
    i <- d[[i]]
    i <- i[sid,]
    colnames(i)=gsub('dMML-','',colnames(i))
    i <- i[,colnames(i) %in% timeorder]
    i <- scalematrix(i)
    i <- i[complete.cases(i),]
    set.seed(seed)
    clu <- kmeans(i,10,iter.max = 10000)$cluster
    n <- names(clu)
    clum <- rowsum(i,clu)/as.vector(table(clu))
    maxp <- apply(clum,1,function(i) {
      i[i < 0] <- 0
      i <- i-min(i)
      i <- i/max(i)
      sum(i*c(1:length(i)))/sum(i)
    })
    clu <- rank(maxp)[clu]
    names(clu) <- n
    clu
  })
  dir_uc=paste0(dir_cluster_in_non_ts,'dmml_',sub('\\.','',as.character(cut)),'/')
    ifelse(!dir.exists(file.path(dir_uc)), dir.create(file.path(dir_uc)), FALSE)
  saveRDS(cluster_d,file=paste0(dir_uc,'dmml_',cut,'_',seed,'.rds'))
}

#Plotting non-tissue-specific heatmap
figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_non_ts_dMML.jpeg')
 
  # Plot heatmap ------------------------------------------------------------
  cat('Plotting heatmap\n')
  library(RColorBrewer)
  library(pheatmap)
  library(gplots)

  # dmml <-readRDS(dmml_cor_file)
  # dnme <-readRDS(dnme_cor_file)
  tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
  timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
  seed=1
  clu=readRDS(paste0(dir_uc,'dmml_',cut,'_',seed,'.rds'))
  d=d[tissue_all]
  d=lapply(d,function(x) {
    
    colnames(x)=gsub(paste0('dMML-|-all|',paste(tissue_all,'-',sep='',collapse = '|')),'',colnames(x))
    return(x)
  })
   dMML_all=d
  d <- sapply(d,function(i) {
    
    i <- i[rowSums(i) > 0,]
    i <- i[,colnames(i) %in% timeorder]
    i <- i[,order(match(colnames(i),timeorder))]
    
    #i <- scalematrix(i)
    i <- i[complete.cases(i),]
  })
  
  mat_out=matrix(ncol=39,nrow=0)
  rowann_out=data.frame()
  row_gap=c(0)
  for (n in names(d)) {
    cl <- clu[[n]]
    cl <- sort(cl)
    mat <- do.call(cbind,sapply(tissue_all,function(i) {
      tmp <- matrix(NA,nrow=length(cl),ncol=ncol(d[[i]]),dimnames = list(names(cl),colnames(d[[i]])))
      
      rn <- intersect(names(cl),rownames(d[[i]]))
      tmp[rn,] <- as.matrix(d[[i]][rn,])
      
      colnames(tmp) <- paste0(i,':',colnames(tmp))
      
      tmp
    }))
    na_ma=-which(rowSums(is.na(mat))>0)
    if(length(na_ma)>0){
      mat= mat[na_ma,]
      rowann <- data.frame(tissue_r=n,cluster=sub(':.*','',cl),
                           #dMMLJSDcor=dmml[[n]][rownames(mat)],
                           #dNMEJSDcor=dnme[[n]][rownames(mat)],
                           stringsAsFactors = F)
      rowann=rowann[na_ma,]
    }else{
      
      rowann <- data.frame(tissue_r=n,cluster=sub(':.*','',cl),
                           #dMMLJSDcor=dmml[[n]][rownames(mat)],
                           #dNMEJSDcor=dnme[[n]][rownames(mat)],
                           stringsAsFactors = F)
    }
    #Note for non-tissue specific UC 01, add tissue name to mat
    rownames(mat)=paste0(n,'-',rownames(mat))
    mat_out=rbind(mat_out,mat)
    
    rownames(rowann) <- rownames(mat)
    rowann <- rowann[,ncol(rowann):1]
    rowann_out=rbind(rowann_out,rowann)
    #row_gap=c(row_gap,row_gap[length(row_gap)]+cumsum(rle(sub(':.*','',cl))$lengths))
    row_gap=c(row_gap,row_gap[length(row_gap)]+nrow(mat))
  }
  
  #Refine plotting parameters
  colann <- data.frame(time=sub('.*:','',colnames(mat_out)),tissue=sub(':.*','',colnames(mat_out)),stringsAsFactors = F)
  rownames(colann) <- colnames(mat_out)
  c1 <- mouse_color()
  c2 <- brewer.pal(10,'Set3')
  names(c2) <- 1:10
  c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
  names(c4) <- sort(unique(colann[,1]))
  #remove row with all NA 

  sub_sp=sort(sample(1:nrow(mat_out),round(nrow(mat_out))))
    mat_out_sc=scalematrix(mat_out[sub_sp,])
  jpeg(figure_name,width=4000,height=30000,res=500)
  pheatmap(mat_out_sc,cluster_rows = F,
           annotation_row = rowann_out[sub_sp,],
           cluster_cols = F,
           annotation_col = colann,show_colnames = F,show_rownames = F,
           #gaps_row = row_gap[-1],
           gaps_col = cumsum(rle(colann[,2])$lengths),
           annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4)
                                    #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10))
                                    
           #filename=figure_name
  )
  dev.off()
  


