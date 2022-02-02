library(bsseq)
coverage=read.table('coverage.txt',col.names=c('chr','start','end','cov'))
coverage_gr=makeGRangesFromDataFrame(coverage,keep.extra.columns = TRUE)
file="HUES64_ectoderm_21_ectoderm_paired_cat_dup_rh.CpG_report.txt"
##Ignore this, you need to separate positive and negative strain then calculate coverage
dat_report <- read.table(file, row.names = NULL,
                  col.names = c("chr", "pos", "strand", "M", "U", "context", "context2"),
                  colClasses = c("character", "numeric", "character","integer", "integer", "character", "character"))
gr.cg <- GRanges(dat_report $chr, IRanges(dat_report$pos, dat_report$pos))
gr.cg$score=dat_report$M+dat_report$U
rtracklayer::export.bed(gr.cg,"test_CpG.bed")