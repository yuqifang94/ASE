suppressMessages(library(rtracklayer))
library(data.table)
suppressMessages(library(GenomicRanges))
#for (pattern in c('run2','runSL')) {
gtf <- fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))

# var <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/mousereprogram/res/fit/hypervar/',pattern,'.rds'))[[2]]
# ct <- sub('_.*','',colnames(var))
# var <- sapply(unique(ct),function(i) rowMeans(var[,which(ct==i)]))
# var <- var[intersect(row.names(var),gn),]

library(splines)

expr <- expr[intersect(row.names(expr),gn),]
ct <- sub('_.*','',colnames(expr))
m <- sapply(unique(ct),function(i) {
  rowMeans(expr[,which(ct==i)])
})

var <- sapply(unique(ct),function(i) {
  apply(expr[,which(ct==i)],1,var)
})

hypervar <- sapply(unique(ct),function(i) {
  m <- rowMeans(expr[,which(ct==i)])
  lm <- log(m)
  var <- apply(expr[,which(ct==i)],1,var)
  resid <- resid(lm(var~bs(lm)))
})

#start <- ifelse(gtf[,7]=='+',gtf[,4]-10000,gtf[,5]+5000)
#end <- ifelse(gtf[,7]=='+',gtf[,4]-5000,gtf[,5]+10000)
#progr <- GRanges(seqnames=gtf[,1],IRanges(start=start,end=end),strand=gtf[,7])
#names(progr) <- gn

gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- gn
for (reg in c(1000,2500,5000)) {
  progr <- promoters(gene,upstream=reg*2,downstream=reg)
  progr <- progr[row.names(var)]
  
  nme <- sapply(colnames(var),function(i) {
    d <- import(paste0('/home-4/zji4@jhu.edu/scratch/mousereprogram/data/wgbs/mm10/NME-',sub('-2$','-1',i),'.bw'))
    nmv <- unlist(mcols(d))
    o <- as.matrix(findOverlaps(progr,d))
    o <- data.frame(names(progr)[o[,1]],nmv[o[,2]])
    vv <- tapply(o[,2],list(o[,1]),mean)
    tar <- rep(NA,length(names(progr)))
    names(tar) <- names(progr)
    tar[names(vv)] <- vv
    tar
  })
  mcorv <- sapply(intersect(row.names(nme),row.names(m)),function(i) cor(m[i,],nme[i,]))
  varcorv <- sapply(intersect(row.names(nme),row.names(m)),function(i) cor(var[i,],nme[i,]))
  hypervarcorv <- sapply(intersect(row.names(nme),row.names(m)),function(i) cor(hypervar[i,],nme[i,]))
  saveRDS(mcorv,file=paste0('/home-4/zji4@jhu.edu/scratch/mousereprogram/res/wgbs/varcor/acrosstime/',reg,'_mean.rds'))
  saveRDS(varcorv,file=paste0('/home-4/zji4@jhu.edu/scratch/mousereprogram/res/wgbs/varcor/acrosstime/',reg,'_var.rds'))
  saveRDS(hypervarcorv,file=paste0('/home-4/zji4@jhu.edu/scratch/mousereprogram/res/wgbs/varcor/acrosstime/',reg,'_hypervar.rds'))
}

