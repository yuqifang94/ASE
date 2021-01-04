library(data.table)
suppressMessages(library(GenomicRanges))
suppressMessages(library(topGO))
gtf <- fread('/dcl01/hongkai/data/zji4/resource/gtf/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
names(gr) <- gn
pro <- promoters(gr,upstream=2000,downstream=1000)
d <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/cluster/cluster.rds')
allres <- NULL
for (n in names(d)) {
  print(n)
  clu <- d[[n]]
  clugr <- names(clu)
  clugr <- GRanges(seqnames=sub(':.*','',clugr),IRanges(start=as.numeric(sub('.*:','',sub('-.*','',clugr))),end=as.numeric(sub('.*-','',clugr))))
  back <- unique(names(pro)[as.matrix(findOverlaps(pro,clugr))[,1]])
  for (sc in 1:max(clu)) {
    clugr <- names(clu)[clu==sc]
    clugr <- GRanges(seqnames=sub(':.*','',clugr),IRanges(start=as.numeric(sub('.*:','',sub('-.*','',clugr))),end=as.numeric(sub('.*-','',clugr))))
    gl <- unique(names(pro)[as.matrix(findOverlaps(pro,clugr))[,1]])
    geneList <- factor(as.integer(back %in% gl))
    names(geneList) <- back
    
    suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)})
    sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
    sigres <- sigres[sigres$Annotated >= 10,]
    sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
    sigres <- sigres[as.numeric(sigres$FDR) <= 0.25,]
    fc <- (sigres[,"Significant"]/sum(GOdata@allScores[GOdata@feasible]==1))/(sigres[,"Annotated"]/sum(GOdata@feasible))
    sigres <- data.frame(sigres,FC=fc)
    sigres <- sigres[order(sigres$FDR,-sigres$FC),]
    if (nrow(sigres) > 0) {
      allres <- rbind(allres,data.frame(tissue=n,cluster=sc,sigres,stringsAsFactors = F))
    }
  }
}
saveRDS(allres,file='/dcl01/hongkai/data/zji4/ase/mouse/GO/GO.rds')
