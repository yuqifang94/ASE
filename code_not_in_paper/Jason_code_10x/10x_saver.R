library(SAVER)
m <- readRDS('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/proc/10x.rds')
m <- 2^m-1
p <- sub(':.*','',colnames(m))
i <- as.numeric(commandArgs(trailingOnly = T))
tp <- unique(p)[i]
m <- m[,p==tp]
m <- m[rowSums(m) > 0,]
d <- saver(m,ncores=10,size.factor=1)$estimate  
saveRDS(d,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/saver/10x/',tp,'.rds'))

