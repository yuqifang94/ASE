rm(list=ls())
source("mainFunctions_sub.R")

#DNase vs non-DNase mouse
  
  DNAase=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp.rds')
  control=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp_control.rds')
  PRC_binding=readRDS('../downstream/output/mm10_PRC.rds')
  PRC_binding=granges(PRC_binding)
  # enhancers=readRDS('../downstream/input/enhancers_HACER.rds')
  rds_save_file="../downstream/input/mm10_agnostic.rds"
  
  chrs <- names(Mmusculus)[1:21]#2276796
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
  cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- GenomicFeatures::genes(txdb)
  TSS<-promoters(genes,upstream=0,downstream=0)
  TSS$gene_name=AnnotationDbi::select(Mus.musculus,key=as.character(TSS$gene_id),
                                      keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
  
  PRC_binding$region_type="PRC"
  DNAase$region_type='DNase'
  control$region_type='control'
  TSS_break=c(DNAase,control,PRC_binding)
  TSS_break=dist_calc(TSS_break,TSS)

  #https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/
  blacklist_region=import.bed('../downstream/input/mm10.blacklist.bed.gz')
  gff_gen(TSS_break,cpgr,blacklist_region,'../downstream/output/mm10_allele_agnostic_analysis.rds','../downstream/output/mm10_allele_agnostic_analysis.gff')
# Human agnostic analysis -------------------------------------------------
  #Get all CpG location
  chrs <- names(Hsapiens)[1:24]
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 1)))) #use first location
  #https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/
  hg19_bl=import.bed('../downstream/input/hg19_blacklist.bed.gz')
  
  # 20kb from hg19 TSS ------------------------------------------------------
  extend_size=20000
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- GenomicFeatures::genes(txdb)
  TSS<-promoters(genes,upstream=100,downstream=100)
  TSS$gene_symbol=AnnotationDbi::select(Homo.sapiens,key=as.character(TSS$gene_id),
                                        keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
  TSS=TSS[seqnames(TSS) %in% paste0('chr',1:22)]
  start(TSS)=start(TSS)-extend_size
  end(TSS)=end(TSS)+extend_size+1
  TSS=trim(TSS)
  TSS_break=subdivideGRanges(TSS,250)
  gff_gen(TSS_break,cpgr,hg19_bl,'../downstream/output/hg19_TSS_20kb_250bp.rds','../downstream/output/hg19_TSS_20kb_250bp.gff')

  # DNase vs non-DNase ------------------------------------------------------
  DNase=readRDS('../downstream/input/DNase_hg19_250bp.rds')
  control=readRDS('../downstream/input/DNase_hg19_250bp_control.rds')
  TSS_break=c(DNase,control)
  gff_gen(TSS_break,cpgr,hg19_bl,'../downstream/output/hg19_DNase_250bp.rds','../downstream/output/hg19_DNase_250bp.gff')

  # All region have SNP -----------------------------------------------------
  GR_merge=readRDS(GR_merge_file)
  gff_gen(granges(unique(GR_merge)),cpgr,hg19_bl,'../downstream/output/human_ASM_region_allele_agnostic_250bp.gff','../downstream/output/human_ASM_region_allele_agnostic_250bp.gff')


#hg19
#about 4000 regions are affected, 1600200
# sed -i 's/%2c/,/g' hg19_TSS_20kb_250bp.gff
# sed -i 's/];/]/g' hg19_TSS_20kb_250bp.gff
# sed -i 's/chr//g' hg19_TSS_20kb_250bp.gff
# subj=(H9 HUES64 skin03 STL001 STL002 STL003 STL011 H1 HuFGM02 112 149 150)
# for i in "${subj[@]}"; do cp hg19_TSS_20kb_250bp.gff ~/work/shared/CpelAsm/data/${i}/cpelasm/${i}_allele_agnostic_analysis.gff; done
#mm10
# sed -i 's/%2c/,/g' mm10_allele_agnostic_analysis.gff
# sed -i 's/];$/]/g' mm10_allele_agnostic_analysis.gff
  # sed -i 's/rtracklayer/\./g' mm10_allele_agnostic_analysis.gff
  # sed -i 's/sequence_feature/\./g' mm10_allele_agnostic_analysis.gff

#filter regions


# generating sbatching arguments ------------------------------------------

  
  