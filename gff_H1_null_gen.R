rm(list=ls())
source("mainFunctions_sub.R")
library(exomeCopy)
library(Mus.musculus)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
species='mm10_DNase'
#######Making gff file that are within 5k of TSS
# enhancers=read.table('../downstream/input/enhancer_HACER.txt',sep='\t',stringsAsFactors = F)
# colnames(enhancers)=c('Enhancer_ID','chr','start','end','center','FANTOM5',
#                       'asssociated_gene_FANTOM5','associated_gene_50kb','associated_gene-4DGenome',
#                       'Cell/Tissue','Detection_method','PMID','closest_gene','distance','Technique','CellType',
#                       'Genome','source','Normalized-count','density','VISTA','Ensembl','ENCODE','chromHMM')
# #enhancers=as.data.frame(enhancers,stringsAsFactors=F)
# enhancers=makeGRangesFromDataFrame(enhancers,keep.extra.columns = T)

if(species=='mm10'){
FANTOM=import.bed('../downstream/input/FANTOM_mm10.bed.gz')
# enhancers=readRDS('../downstream/input/enhancers_HACER.rds')
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- GenomicFeatures::genes(txdb)
TSS<-promoters(genes,upstream=0,downstream=0)
TSS$gene_symbol=AnnotationDbi::select(Mus.musculus,key=as.character(TSS$gene_id),
                                      keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
out_name="../downstream/output/mm10_FANTOM_TSS_3kb_250bp.gff"
rds_save_file="../downstream/output/mm10_FANTOM_TSS_3kb_250bp.rds"
chrs <- names(Mmusculus)[1:19]#2276796
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr <- do.call(c, lapply(1:19, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
}else if(species=='hg19'){
  
     FANTOM=import.bed('../downstream/input/FANTOM_hg19.bed.gz')
    # enhancers=readRDS('../downstream/input/enhancers_HACER.rds')
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    genes <- GenomicFeatures::genes(txdb)
    TSS<-promoters(genes,upstream=0,downstream=0)
    TSS$gene_symbol=AnnotationDbi::select(Homo.sapiens,key=as.character(TSS$gene_id),
                                          keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
    out_name="../downstream/output/hg19_FANTOM_TSS_3kb_250bp.gff"
    rds_save_file="../downstream/output/hg19_FANTOM_TSS_3kb_250bp.rds"
    
    chrs <- names(Hsapiens)[1:22]#2276796
    cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
    cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 1)))) #use first location
}else if(species=='mm10_DNase'){
  
  DNAase=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp.rds')
  control=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp_control.rds')
  
  # enhancers=readRDS('../downstream/input/enhancers_HACER.rds')
   rds_save_file="../downstream/input/mm10_DNase.rds"
  
  chrs <- names(Mmusculus)[1:21]#2276796
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
  cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location

  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- GenomicFeatures::genes(txdb)
  TSS<-promoters(genes,upstream=0,downstream=0)
  TSS$gene_name=AnnotationDbi::select(Mus.musculus,key=as.character(TSS$gene_id),
                                        keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
  
  DNAase$region_type='DNase'
  control$region_type='control'
  TSS_break=c(DNAase,control)
  TSS_break=dist_calc(TSS_break,TSS)
  out_name="../downstream/output/mm10_DNase_3kb_250bp.gff"

}else if(species=='mm10_PRC'){
  #PRC2 regions

  read_chromHMM_bed<-function(bed_dir,rep){
    bed_out=GRanges()
    for(fn in dir(bed_dir,pattern='.bed.gz')){
      #get sample name etc
      fn_split=strsplit(fn,'_')[[1]]
      stage=gsub('e','day',fn_split[1])
      stage=gsub('P','day',stage)
      stage=gsub('\\.','\\_',stage)
      tissue=gsub('facial-prominence','EFP',fn_split[2])
      tissue=gsub('neural-tube','NT',tissue)
      bed_in=read.table(paste(bed_dir,fn,sep=''))
      colnames(bed_in)=c('chr','start','end','chrom_num','chrom_state')
      bed_in=makeGRangesFromDataFrame(bed_in,keep.extra.columns = T)
      #bed_in=reduce(bed_in[bed_in$chrom_state%in%c('Hc-P','Pr-B')])
      bed_in$stage=stage
      bed_in$tissue=tissue
      bed_in$rep=rep
      bed_in$Sample=paste(stage,tissue,sep='-')
      bed_out=c(bed_out,bed_in)
    }
    return(bed_out)
  }
  pooled_PRC= read_chromHMM_bed('../downstream/input/chromHMM_mm10/pooled/','pooled')
  
  #rep1_PRC= read_chromHMM_bed('../downstream/input/chromHMM_mm10/rep1/','rep1')
  #rep2_PRC= read_chromHMM_bed('../downstream/input/chromHMM_mm10/rep2/','rep2')
  TSS_break=subdivideGRanges(pooled_PRC,250)
  # enhancers=readRDS('../downstream/input/enhancers_HACER.rds')
  rds_save_file="../downstream/output/mm10_PRC.rds"
  
  chrs <- names(Mmusculus)[1:21]#2276796
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
  cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- GenomicFeatures::genes(txdb)
  TSS<-promoters(genes,upstream=0,downstream=0)
  TSS$gene_symbol=AnnotationDbi::select(Mus.musculus,key=as.character(TSS$gene_id),
                                        keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
  TSS_break=dist_calc(TSS_break,TSS)
 
  out_name="../downstream/output/mm10_PRC_250bp.gff"
}

olap=findOverlaps(TSS_break,cpgr)
CpG_df=data.frame(TSS_hit=queryHits(olap),CpG_start=start(cpgr)[subjectHits(olap)])
CpG_df_agg=aggregate(CpG_start~TSS_hit,CpG_df,function(x) paste(x,collapse=', '))
CpG_df_agg$N=unlist(lapply(CpG_df_agg$CpG_start,function(x) length(strsplit(x,', ')[[1]])))
CpG_df_agg$CpG_start=paste("[",CpG_df_agg[,2],"]",sep='')
TSS_break_out=TSS_break[CpG_df_agg[,1]]
strand(TSS_break_out)='*'
TSS_break_out$N=CpG_df_agg$N
TSS_break_out$CpGs=CpG_df_agg$CpG_start
#filter out blacklist region and region with N=1
blacklist_region=import.bed('../downstream/input/mm10.blacklist.bed.gz')
#TSS_break_out=TSS_break_out[TSS_break_out$N>1]
olap=findOverlaps(TSS_break_out,blacklist_region)
TSS_break_out=TSS_break_out[-queryHits(olap)]
TSS_break_out_gff=granges(TSS_break_out)
mcols(TSS_break_out_gff)=mcols(TSS_break_out)[,c('N','CpGs')]
export.gff3(sort(TSS_break_out_gff),out_name)

all_chromatin= read_chromHMM_bed('../downstream/input/chromHMM_mm10/pooled/','pooled')
#all_chromatin_marker=all_chromatin[all_chromatin$chrom_state%in%c('En-Sd','En-Sp','En-W')]
all_chromatin_marker=all_chromatin[all_chromatin$chrom_state%in%c('Pr-A','Pr-W')]
TSS_break_out=readRDS('../downstream/output/mm10_DNase.rds')
#create a annotation file of which Sample it is.
for (sp in unique(all_chromatin_marker$Sample)){
  mcols(TSS_break_out)[,sp]=NA
  sp_olap=findOverlaps(TSS_break_out,all_chromatin_marker[all_chromatin_marker$Sample==sp])
  mcols(TSS_break_out)[unique(queryHits(sp_olap)),sp]=TRUE
}
#Annotate chromatin states

sp_olap=findOverlaps(TSS_break_out,all_chromatin_marker)
sp_df=data.table(qt=queryHits(sp_olap),sample=all_chromatin_marker$Sample[subjectHits(sp_olap)],
                 state_chrom =all_chromatin_marker$chrom_state[subjectHits(sp_olap)])
paste_unique<-function(x) {paste(unique(x),collapse=',')}
sp_df_dc=dcast.data.table(sp_df,qt~sample,value.var = "state_chrom",fun.aggregate =paste,collapse=',',fill = NA)
mcols(TSS_break_out)=cbind(mcols(TSS_break_out),sp_df_dc[,-1])

saveRDS(TSS_break_out,'../downstream/output/mm10_DNase_all_chromHMM.rds')
saveRDS(all_chromatin_marker,'../downstream/output/mm10_chromHMM_promoter.rds')
# saveRDS(TSS_break_out,'../downstream/output/mm10_DNase_enhancer.rds')
# saveRDS(all_chromatin_marker,'../downstream/output/mm10_chromHMM_enhancer.rds')
saveRDS(TSS_break_out,rds_save_file)
#hg19
#about 4000 regions are affected, 1600200
# sed -i 's/%2c/,/g' hg19_FANTOM_TSS_3kb_250bp.gff
# sed -i 's/];/]/g' hg19_FANTOM_TSS_3kb_250bp.gff
# sed -i 's/chr//g' hg19_FANTOM_TSS_3kb_250bp.gff
# subj=(H9 HUES64 skin03 STL001 STL002 STL003 STL011 H1 HuFGM02 112 149 150)
# for i in "${subj[@]}"; do cp hg19_FANTOM_TSS_3kb_250bp.gff ~/work/shared/CpelAsm/data/${i}/cpelasm/${i}_allele_agnostic_analysis.gff; done
#mm10
# sed -i 's/%2c/,/g' mm10_allele_agnostic_analysis.gff
# sed -i 's/];$/]/g' mm10_allele_agnostic_analysis.gff
#mm10 DNase
# sed -i 's/%2c/,/g' mm10_DNase_3kb_250bp.gff
# sed -i 's/];/]/g' mm10_DNase_3kb_250bp.gff
# sed -i 's/rtracklayer/\./g' mm10_DNase_3kb_250bp.gff
# sed -i 's/sequence_feature/\./g' mm10_DNase_3kb_250bp.gff
#mm10 PRC
# sed -i 's/%2c/,/g' mm10_PRC_250bp.gff
# sed -i 's/];/]/g' mm10_PRC_250bp.gff
# sed -i 's/rtracklayer/\./g' mm10_PRC_250bp.gff
# sed -i 's/sequence_feature/\./g' mm10_PRC_250bp.gff

#filter regions
