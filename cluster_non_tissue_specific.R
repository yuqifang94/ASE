source('mainFunctions_sub.R')
cut <- 0.1
d=readRDS(UC_in_matrix_cluster_file)

# $EFP
# [1] 5008585
# $forebrain
# [1] 5005549
# $heart
# [1] 5006020
# $hindbrain
# [1] 5006806
# $limb
# [1] 5006777
# $liver
# [1] 4995202
# $midbrain
# [1] 5004404
aid <- sapply(names(d),function(i) {
    names(which(rowSums(d[[i]] > cut) > 0))
  })  
  timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
#   $EFP
# [1] 951833
# $forebrain
# [1] 1184078
# $heart
# [1] 1007666
# $hindbrain
# [1] 1149669
# $limb
# [1] 988370
# $liver
# [1] 1111038
# $midbrain
# [1] 1129165
for (seed in 1:10) {
  cat('Processing:',seed,'\n')
  cluster_d <- sapply(names(d),function(i) {
    #This is for one and only one
    #sid <- setdiff(aid[[i]],unlist(aid[names(aid)!=i]))
    sid=aid[[i]]
    i <- d[[i]]
    colnames(i)=sub('-all','',sub('.*?-','',colnames(i)))
    i <- i[sid,]
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
  dir_uc=paste0(dir_cluster_in,'uc_',sub('\\.','',as.character(cut)),'/')
    ifelse(!dir.exists(file.path(dir_uc)), dir.create(file.path(dir_uc)), FALSE)
  saveRDS(cluster_d,file=paste0(dir_uc,cut,'_',seed,'.rds'))
}




