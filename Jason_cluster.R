suppressMessages(library(GenomicRanges))

d <- readRDS('allele_agnostic_uc_complement/UC_agnostic_mouse_matrix_dedup_N2_all_non_MDS.rds')

d <- sapply(d,function(gr) {
  grv <- paste0(as.character(seqnames(gr)),':',start(gr),'-',end(gr))
  am <- as.matrix(mcols(gr))
  rownames(am) <- grv
  am <- am[,!grepl('P0',colnames(am))]
  am <- am[complete.cases(am),]
})
d <- d[sapply(d,ncol) > 6]
names(d) <- sapply(d,function(i) sub('-.*','',colnames(i)[1]))
for (i in 1:length(d)) colnames(d[[i]]) <- sub('-all','',sub(paste0(names(d)[i],'-'),'',colnames(d[[i]])))
saveRDS('fulluc.rds')

#### clustering
library(RColorBrewer)
library(pheatmap)


for (seed in 1:10) {
  cut <-0.01
  seed=1
  d <- readRDS('fulluc.rds')
  timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
  
  aid <- sapply(names(d),function(i) {
    names(which(rowSums(d[[i]] > cut) > 0))
  })  
  
  d2 <- sapply(names(d),function(i) {
    ###### this is the line to ensure it's > cut in only one tissue
    sid <- setdiff(aid[[i]],unlist(aid[names(aid)!=i]))
    i <- d[[i]]
    i <- i[sid,]
    i <- i[,colnames(i) %in% timeorder]
    scalematrix <- function(data) {
      cm <- rowMeans(data)
      csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
      (data - cm) / csd
    }
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
  saveRDS(d2,file=paste0('uc_',cut,'_',seed,'.rds'))
}
