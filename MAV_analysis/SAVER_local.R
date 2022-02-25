#It seems SAVER can only be on conda_R
library(SAVER)
library(data.table)
rawCount=readRDS("../downstream/data/mouseLimb/mouse_limb_scRNA_10x_rawCount.rds")

annotation=fread('../downstream/data/mouseLimb/meta.tsv')
#Normalizing across all smples
totalCount=colSums(rawCount)
totalCount=totalCount/median(totalCount)
#Normalzing the total count of each cell by the median total count
res <- saver(rawCount,ncores=20,size.factor=1)$estimate
