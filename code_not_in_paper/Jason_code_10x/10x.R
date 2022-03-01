library(data.table)
d <- fread('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/raw/10x/exprMatrix.tsv',data.table=F)
m <- fread('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/raw/10x/meta.tsv',data.table = F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
print(identical(m[,1],colnames(d)))
colnames(d) <- paste0(m$stage,':cell',1:ncol(d))
saveRDS(d,file='/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/proc/10x.rds')

