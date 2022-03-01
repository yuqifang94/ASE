library(SAVER)
m <- readRDS('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/proc/c1.rds')
m <- 2^m-1
p <- sub(':.*','',colnames(m))
for (tp in unique(p)) {
  d <- saver(m[,p==tp],ncores=10,size.factor=1)$estimate  
  saveRDS(d,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/saver/c1/',tp,'.rds'))
}

