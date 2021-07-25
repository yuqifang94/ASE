source('mainFunctions_sub.R')
library("org.Mm.eg.db")

enhancer_Bin=readRDS(bin_enhancer_rds)
#Reading in RNA data
xx <- as.list(org.Mm.egENSEMBLTRANS)
RNA_mouse_dir='../downstream/data/RNA_tsv/'
for(fn in ,dir(RNA_mouse_dir,pattern='.tsv')){
  RNA_in=fread(paste0(RNA_mouse_dir,fn))



}