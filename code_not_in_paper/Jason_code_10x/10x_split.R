m <- readRDS('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/proc/10x.rds')
p <- sub(':.*','',colnames(m))
for (tp in unique(p)) {
  d <- m[,p==tp]
  saveRDS(d,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/split/10x/',tp,'.rds'))
}
