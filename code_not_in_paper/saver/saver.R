library(SAVER)
af <- list.files('/home-4/zji4@jhu.edu/scratch/scdata/HCL/proc/count')
rc <- unlist(sapply(list.files('/home-4/zji4@jhu.edu/scratch/scdata/HCL/proc/rc/',full.names = T),readRDS,USE.NAMES = F))
f <- af[as.numeric(commandArgs(trailingOnly = T))]
if (!file.exists(paste0('/home-4/zji4@jhu.edu/scratch/scdata/HCL/proc/saver/',f))) {
  d <- as.matrix(readRDS(paste0('/home-4/zji4@jhu.edu/scratch/scdata/HCL/proc/count/',f)))
  lib <- rc[colnames(d)]#raw count of all samples genes that in this sample, cols: cell, rows: gene
  lib <- lib/median(rc)#rc is the total cell count per cell, so lib here is a total cell count devidied by median cell count for all cells
  d <- t(t(d)/lib)
  res <- saver(d,ncores=5,size.factor=1)$estimate
  saveRDS(res,file=paste0('/home-4/zji4@jhu.edu/scratch/scdata/HCL/proc/saver/',f))
}
