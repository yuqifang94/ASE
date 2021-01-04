library(RColorBrewer)
library(pheatmap)
d <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/proc/dnase.rds')
mml <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/proc/mml.rds')
nme <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/proc/nme.rds')

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}

timeorder <- sapply(1:20,function(i) paste0('day',i,'_5-day',i+1,'_5'))
clu <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/cluster/cluster.rds')
d <- sapply(d,function(i) {
  i <- i[rowSums(i > 0.1) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  i <- i[,order(match(colnames(i),timeorder))]
  scalematrix <- function(data) {
    cm <- rowMeans(data)
    csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
    (data - cm) / csd
  }
  i <- scalematrix(i)
  i <- i[complete.cases(i),]
})

for (n in names(d)) {
  cl <- clu[[n]]
  cl <- sort(cl)
  mat <- do.call(cbind,sapply(c(n,setdiff(names(d),n)),function(i) {
    tmp <- matrix(NA,nrow=length(cl),ncol=ncol(d[[i]]),dimnames = list(names(cl),colnames(d[[i]])))
    rn <- intersect(names(cl),rownames(d[[i]]))
    tmp[rn,] <- d[[i]][rn,]
    colnames(tmp) <- paste0(i,':',colnames(tmp))
    tmp
  }))
  dmml <- sapply(colnames(d[[n]]),function(i) {
    time <- strsplit(i,'-')[[1]]
    mml[rownames(d[[n]]),paste0(n,'-',time[1],'-all')]-mml[rownames(d[[n]]),paste0(n,'-',time[2],'-all')]
  })
  dnme <- sapply(colnames(d[[n]]),function(i) {
    time <- strsplit(i,'-')[[1]]
    nme[rownames(d[[n]]),paste0(n,'-',time[1],'-all')]-nme[rownames(d[[n]]),paste0(n,'-',time[2],'-all')]
  })
  dmmlcor <- corfunc(dmml,d[[n]])
  dnmecor <- corfunc(dnme,d[[n]])
  colann <- data.frame(tissue=sub(':.*','',colnames(mat)),time=sub('.*:','',colnames(mat)),stringsAsFactors = F)
  rownames(colann) <- colnames(mat)
  rowann <- data.frame(cluster=as.character(cl),dMML=dmmlcor,dNME=dnmecor,stringsAsFactors = F)
  rownames(rowann) <- names(cl)
  
  c1 <- brewer.pal(length(c(n,setdiff(names(d),n))),'Set1')
  names(c1) <- c(n,setdiff(names(d),n))
  c2 <- brewer.pal(max(clu[[1]]),'Set3')
  names(c2) <- 1:max(clu[[1]])
  
  png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/',n,'.png'),width = 800,height=800)
  pheatmap(mat,cluster_rows = F,annotation_row = rowann,cluster_cols = F,annotation_col = colann,show_colnames = F,show_rownames = F,
           gaps_row = cumsum(table(cl)),gaps_col = cumsum(rle(colann[,1])$lengths),annotation_colors = list(tissue=c1,cluster=c2))
  dev.off()
}
