#Uuse SAVER in MARCC, conda_R with Jason's code at /scratch/users/zji4@jhu.edu/andy_ASE/hypervar/encmouse/data/code
#MARCC R/4.0.2
library(SAVER)
library(data.table)
#All takes too long, need to separate by stage
m=readRDS('../../downstream/data/mouseLimb/limb_10x.rds')
stages=gsub('.*-|:.*','',colnames(m))
i <- as.numeric(commandArgs(trailingOnly = T))#This is nth stage
stages_selected <- unique(stages)[i]
cat("Processing stage:",stages_selected,'\n')
m <- 2^m-1#dim: 43346 90637
m <- m[,stages==stages_selected]
m <- m[rowSums(m) > 0,]#26384 90637
d <- saver(m[1:200,],ncores=10,size.factor=1)$estimate  
saveRDS(d,file=paste0('../../downstream/data/mouseLimb/10xLimbSaver_',stages_selected,'.rds'))