af <- list.files('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/saver/10x')
for (f in af) {
  expr <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/saver/10x/',f))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
  expr <- expr[!grepl('^Rpl|^Rps',rownames(expr)),]
  cv <- apply(expr,1,sd)/rowMeans(expr)
  saveRDS(cv,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/cv/10x/',f))
}

