#This code is adapted from Jason only chaning path
af <- list.files('../downstream/data/mouse_c1_saver/')
for (f in af) {
  expr <- readRDS(paste0('../downstream/data/mouse_c1_saver/',f))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
  expr <- expr[!grepl('^Rpl|^Rps',rownames(expr)),]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.frame(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  saveRDS(res,file=paste0(dir_scRNA_mouse,f))
}

