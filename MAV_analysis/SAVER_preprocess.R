#From Jason's code
library(data.table)
rm(list=ls())
d <- fread('../downstream/data/mouseLimb/exprMatrix.tsv',data.table=F)
m <- fread('../downstream/data/mouseLimb/meta.tsv',data.table = F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
m$merged_celltype=trimws(gsub(paste(1:10,collapse="|"),"",m$cell_type))
m$cellType_stage=paste0(m$merged_celltype,"-E",m$stage)
print(identical(m[,1],colnames(d)))
colnames(d) <- paste0(m$cellType_stage,':cell',1:ncol(d))
saveRDS(d,'../downstream/data/mouseLimb/limb_10x.rds')