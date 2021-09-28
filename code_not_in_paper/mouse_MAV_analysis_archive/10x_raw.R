af <- list.files('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/split/10x')
for (f in af) {
  expr <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/split/10x/',f))
  expr <- expr[rowMeans(expr > 0) > 0.1,]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.frame(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  saveRDS(res,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/res/10x/',f))
}

