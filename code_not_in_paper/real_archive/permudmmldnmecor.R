library(RColorBrewer)
library(pheatmap)
library(ggplot2)
d <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/proc/uc.rds')
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

timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
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

dmmlcor <- dnmecor <- list()
for (n in names(d)) {
  dmml <- sapply(colnames(d[[n]]),function(i) {
    time <- strsplit(i,'-')[[1]]
    abs(mml[rownames(d[[n]]),paste0(n,'-',time[1])]-mml[rownames(d[[n]]),paste0(n,'-',time[2])])
  })
  dnme <- sapply(colnames(d[[n]]),function(i) {
    time <- strsplit(i,'-')[[1]]
    abs(nme[rownames(d[[n]]),paste0(n,'-',time[1])]-nme[rownames(d[[n]]),paste0(n,'-',time[2])])
  })
  mat <- do.call(rbind,combinat::permn(1:ncol(d[[n]])))
  iden <- apply(mat,1,function(i) {
    mean(i==sort(i))
  })
  mat <- mat[iden!=1,]
  set.seed(12345)
  mat <- mat[sample(1:nrow(mat),20),]
  sampid <- t(mat)
  print(sum(apply(sampid,2,function(i) {
    mean(i==sort(i))
  })==1))
  dmmlcor[[n]] <- sapply(1:20,function(i) corfunc(dmml,d[[n]][,sampid[,i]]))
  dnmecor[[n]] <- sapply(1:20,function(i) corfunc(dnme,d[[n]][,sampid[,i]]))
}

saveRDS(dmmlcor,file='/dcl01/hongkai/data/zji4/ase/mouse/res/dmmldnmecor/permudmmlcor.rds')
saveRDS(dnmecor,file='/dcl01/hongkai/data/zji4/ase/mouse/res/dmmldnmecor/permudnmecor.rds')



