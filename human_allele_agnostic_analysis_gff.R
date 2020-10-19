rm(list=ls())
source('mainFunctions_sub.R')
#gff size =250bp
extend_size=20000
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- GenomicFeatures::genes(txdb)
TSS<-promoters(genes,upstream=100,downstream=100)
TSS$gene_symbol=AnnotationDbi::select(Homo.sapiens,key=as.character(TSS$gene_id),
                                      keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL

#extend 50k
TSS=TSS[seqnames(TSS) %in% paste0('chr',1:22)]
start(TSS)=start(TSS)-extend_size
end(TSS)=end(TSS)+extend_size+1
TSS=trim(TSS)
TSS_break=subdivideGRanges(TSS,250)
chrs <- names(Hsapiens)[1:22]#2276796
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 1)))) #use first location
olap=findOverlaps(TSS_break,cpgr)
iG_df=data.frame(TSS_hit=queryHits(olap),CpG_start=start(cpgr)[subjectHits(olap)])
CpG_df_agg=aggregate(CpG_start~TSS_hit,CpG_df,function(x) paste(x,collapse=', '))
CpG_df_agg$N=unlist(lapply(CpG_df_agg$CpG_start,function(x) length(strsplit(x,', ')[[1]])))
CpG_df_agg$CpG_start=paste("[",CpG_df_agg[,2],"]",sep='')
TSS_break_out=TSS_break[CpG_df_agg[,1]]
strand(TSS_break_out)='*'
TSS_break_out$N=CpG_df_agg$N
TSS_break_out$CpGs=CpG_df_agg$CpG_start
#hg19:https://www.encodeproject.org/files/ENCFF001TDO/
blacklist_region=import.bed('../downstream/input/hg19_blacklist.bed.gz')
#TSS_break_out=TSS_break_out[TSS_break_out$N>1]
olap=findOverlaps(TSS_break_out,blacklist_region)
TSS_break_out=TSS_break_out[-queryHits(olap)]
TSS_break_out_gff=granges(TSS_break_out)
mcols(TSS_break_out_gff)=mcols(TSS_break_out)[,c('N','CpGs')]
export.gff3(sort(TSS_break_out_gff),'../downstream/output/allele_agnostic_gff/human_20kb_allele_agnostic_250bp.gff')

#hg19
#about 4000 regions are affected, 1600200
# sed -i 's/%2c/,/g' hg19_FANTOM_TSS_3kb_250bp.gff
# sed -i 's/];/]/g' hg19_FANTOM_TSS_3kb_250bp.gff
# sed -i 's/chr//g' hg19_FANTOM_TSS_3kb_250bp.gff
# subj=(H9 HUES64 skin03 STL001 STL002 STL003 STL011 H1 HuFGM02 112 149 150)
# for i in "${subj[@]}"; do cp hg19_FANTOM_TSS_3kb_250bp.gff ~/work/shared/CpelAsm/data/${i}/cpelasm/${i}_allele_agnostic_analysis.gff; done