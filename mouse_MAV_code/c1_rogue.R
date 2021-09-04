suppressMessages(library(ROGUE))
suppressMessages(library(tidyverse))
af <- list.files('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/split/c1')
for (f in af) {
  expr <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/split/c1/',f))
  expr <- 2^(expr)-1
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
  expr <- expr[!grepl('^Rpl|^Rps',rownames(expr)),]
  res <- data.frame(SE_fun(expr))
  saveRDS(res,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/rogue/c1/',f))
}


