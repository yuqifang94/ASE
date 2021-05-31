library(data.table)
dir_in='../downstream/input/enhancer_atlas/'
dat_in=data.table()
for(fn in dir(dir_in)){
dat_in=rbind(dat_in,fread(paste0('../downstream/input/enhancer_atlas/',fn)))
}
enhancer_in=data.table(targe_region=gsub('_.*','',dat_in$V1),
                       target_gene=lapply(strsplit(sub('.*_','',dat_in$V1),'\\$'),function(x) x[[2]]))
enhancer_gr=convert_GR(enhancer_in$targe_region)
enhancer_gr$gene=enhancer_in$target_gene#total 41493
bin_enhancer=readRDS('../downstream/output/bin_enhancer.rds')#432 overlap, total 21141
library(liftOver)
#https://hgdownload-test.gi.ucsc.edu/goldenPath/mm9/liftOver/
ch = import.chain('../downstream/input/mm9ToMm10.over.chain')
seqlevelsStyle(cur) = "UCSC"  # necessary
enhancer_gr_mm10 = unlist(liftOver(enhancer_gr, ch))
subsetByOverlaps(bin_enhancer,resize(enhancer_gr_mm10,2000,fix="center"))#7403 of 21141 overlap
olap=findOverlaps(bin_enhancer,enhancer_gr_mm10)
bin_enhancer$enhancer_atlas=NA
bin_enhancer$enhancer_atlas[queryHits(olap)]=enhancer_gr_mm10[subjectHits(olap)]$gene#160 same gene?
