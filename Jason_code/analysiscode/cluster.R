d <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/proc/dnase.rds')
timeorder <- sapply(1:20,function(i) paste0('day',i,'_5-day',i+1,'_5'))

d <- sapply(d,function(i) {
  i <- i[rowSums(i > 0.1) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  scalematrix <- function(data) {
    cm <- rowMeans(data)
    csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
    (data - cm) / csd
  }
  i <- scalematrix(i)
  i <- i[complete.cases(i),]
  set.seed(12345)
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
saveRDS(d,file='/dcl01/hongkai/data/zji4/ase/mouse/cluster/cluster.rds')

