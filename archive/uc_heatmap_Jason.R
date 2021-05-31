rm(list=ls())
library(RColorBrewer)
library(pheatmap)
library(gplots)
source('mainFunctions_sub.R')
d <- readRDS('../downstream/output/uc_matrix_DNase.rds')
dmml <- readRDS('../downstream/input/dmmlcor.rds')
dnme <- readRDS('../downstream/input/dNMEcor.rds')
#dmml <- sapply(dmml,abs)
#dnme <- sapply(dnme,abs)
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
#clu <- readRDS('../downstream/input/Jason_UC_cluster/uc_0.1.rds')
clu=readRDS('../downstream/input/hcluster_result/uc_0.1.rds')
d=d[tissue_all]
d <- sapply(d,function(i) {
  i <- i[rowSums(i) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  i <- i[,order(match(colnames(i),timeorder))]

  #i <- scalematrix(i)
  i <- i[complete.cases(i),]
})

mat_out=matrix(ncol=39)
rowann_out=data.frame()
row_gap=c(0)
for (n in names(d)) {
  cl <- clu[[n]]
  cl <- sort(cl)
  mat <- do.call(cbind,sapply(tissue_all,function(i) {
    tmp <- matrix(NA,nrow=length(cl),ncol=ncol(d[[i]]),dimnames = list(names(cl),colnames(d[[i]])))
    rn <- intersect(names(cl),rownames(d[[i]]))
    tmp[rn,] <- d[[i]][rn,]
    colnames(tmp) <- paste0(i,':',colnames(tmp))
    tmp
  }))
  mat_out=rbind(mat_out,mat)
  rowann <- data.frame(tissue_r=n,cluster=sub(':.*','',cl),dMMLJSDcor=dmml[[n]][rownames(mat)],
                       dNMEJSDcor=dnme[[n]][rownames(mat)],stringsAsFactors = F)
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
# tiff(paste0('../downstream/output/heatmap_acrosstissue/',n,'.tiff'),width=3000,height=3000,res=300)
# #png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
# pheatmap(mat,cluster_rows = F,annotation_row = rowann,cluster_cols = F,
#          annotation_col = colann,show_colnames = F,show_rownames = F,
#          gaps_row = row_gap,gaps_col = cumsum(rle(colann[,1])$lengths),
#          annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4,dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)))
# dev.off()

tiff('../downstream/output/heatmap_acrosstissue/all_sc_N17_ft_hclust.tiff',width=5000,height=5000,res=300)
#png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
pheatmap(scalematrix(mat_out),cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
         annotation_col = colann,show_colnames = F,show_rownames = F,
         gaps_row = row_gap,gaps_col = cumsum(rle(colann[,2])$lengths),
         annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4,dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)))
dev.off()
