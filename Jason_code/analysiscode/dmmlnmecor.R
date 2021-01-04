library(RColorBrewer)
library(pheatmap)
library(ggplot2)
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
  dmmlcortype <- ifelse(dmmlcor > 0.5,'dMMLpositive',ifelse(dmmlcor < -0.5,'dMMLnegative','dMMLzero'))
  dnmecortype <- ifelse(dnmecor > 0.5,'dNMEpositive',ifelse(dnmecor < -0.5,'dNMEnegative','dNMEzero'))
  dtype <- paste0(dmmlcortype,'_',dnmecortype)
  pd <- data.frame(type=dtype,cluster=as.character(cl),stringsAsFactors = F)
  pd$cluster <- factor(pd$cluster,levels=1:10)
  pdf(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/dmmlnmecor/',n,'.pdf'),width = 8,height=6)
  print(ggplot(pd[!grepl('NA',pd[,1]),],aes(type)) + geom_bar() + theme_classic() + facet_wrap(~cluster,scales = 'free_x',nrow=2) + coord_flip())
  dev.off()
}
