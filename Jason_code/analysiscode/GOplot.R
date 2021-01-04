library(gridExtra)
allres <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/GO/GO.rds')

for (type in unique(allres[,1])) {
  d <- NULL
  for (sc in 1:10) {
    res <- allres[allres[,1]==type & allres[,2]==sc,]
    res <- res[res$FDR < 0.05 & res$Annotated >=10,]
    res <- res[order(res$FDR,-res$FC),]
    d <- rbind(d,data.frame(cluster=sc,term=res[1:3,4]))
  }
  pdf(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/GO/',type,'.pdf'),height=9)
  grid.table(d,rows=NULL)
  dev.off()
}

