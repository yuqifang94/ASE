library(pheatmap)
library(ggplot2)
suppressMessages(library(GenomicRanges))

d <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/raw/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')
r <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/raw/mm10_DNase_N2.rds')
grr <- paste0(as.character(seqnames(r)),':',start(r),'-',end(r))
dnase <- grr[mcols(r)$region_type=='DNase']

d <- sapply(d,function(gr) {
  grv <- paste0(as.character(seqnames(gr)),':',start(gr),'-',end(gr))
  am <- as.matrix(mcols(gr))
  rownames(am) <- grv
  am <- am[rownames(am) %in% dnase,]
  am <- am[,!grepl('day0',colnames(am))]
  am <- am[complete.cases(am),]
})
d <- d[sapply(d,ncol) > 6]
names(d) <- sapply(d,function(i) sub('-.*','',colnames(i)[1]))
for (i in 1:length(d)) colnames(d[[i]]) <- sub('-all','',sub(paste0(names(d)[i],'-'),'',colnames(d[[i]])))
saveRDS(d,file='/dcl01/hongkai/data/zji4/ase/mouse/data/proc/dnase.rds')

