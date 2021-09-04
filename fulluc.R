source('mainFunctions_sub.R')
d <- readRDS(UC_in_matrix_ls_file)

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
k <- table(unlist(sapply(d,rownames)))
id <- names(k)[k==length(d)]
d <- sapply(d,function(i) i[id,],simplify = F)
saveRDS(d,file=UC_in_matrix_cluster_file)


