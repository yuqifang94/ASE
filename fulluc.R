suppressMessages(library(GenomicRanges))

d <- readRDS('/home-4/zji4@jhu.edu/scratch/andy_ASE/data/raw/UC_agnostic_mouse_matrix_dedup_N2_all_non_MDS.rds')

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
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/andy_ASE/data/proc/fulluc.rds')


