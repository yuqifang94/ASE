#######Making gff file that are within 5k of TSS
# enhancers=read.table('../downstream/input/enhancer_HACER.txt',sep='\t',stringsAsFactors = F)
# colnames(enhancers)=c('Enhancer_ID','chr','start','end','center','FANTOM5',
#                       'asssociated_gene_FANTOM5','associated_gene_50kb','associated_gene-4DGenome',
#                       'Cell/Tissue','Detection_method','PMID','closest_gene','distance','Technique','CellType',
#                       'Genome','source','Normalized-count','density','VISTA','Ensembl','ENCODE','chromHMM')
# #enhancers=as.data.frame(enhancers,stringsAsFactors=F)
# enhancers=makeGRangesFromDataFrame(enhancers,keep.extra.columns = T)
# else if(species=='mm10_PRC'){
#   #PRC2 regions
# 
#   read_chromHMM_bed<-function(bed_dir,rep){
#     bed_out=GRanges()
#     for(fn in dir(bed_dir,pattern='.bed.gz')){
#       #get sample name etc
#       fn_split=strsplit(fn,'_')[[1]]
#       stage=gsub('e','day',fn_split[1])
#       stage=gsub('P','day',stage)
#       stage=gsub('\\.','\\_',stage)
#       tissue=gsub('facial-prominence','EFP',fn_split[2])
#       tissue=gsub('neural-tube','NT',tissue)
#       bed_in=read.table(paste(bed_dir,fn,sep=''))
#       colnames(bed_in)=c('chr','start','end','chrom_num','chrom_state')
#       bed_in=makeGRangesFromDataFrame(bed_in,keep.extra.columns = T)
#       #bed_in=reduce(bed_in[bed_in$chrom_state%in%c('Hc-P','Pr-B')])
#       bed_in$stage=stage
#       bed_in$tissue=tissue
#       bed_in$rep=rep
#       bed_in$Sample=paste(stage,tissue,sep='-')
#       bed_out=c(bed_out,bed_in)
#     }
#     return(bed_out)
#   }
#   pooled_PRC= read_chromHMM_bed('../downstream/input/chromHMM_mm10/pooled/','pooled')
#   
#   #rep1_PRC= read_chromHMM_bed('../downstream/input/chromHMM_mm10/rep1/','rep1')
#   #rep2_PRC= read_chromHMM_bed('../downstream/input/chromHMM_mm10/rep2/','rep2')
#   TSS_break=subdivideGRanges(pooled_PRC,250)
#   # enhancers=readRDS('../downstream/input/enhancers_HACER.rds')
#   rds_save_file="../downstream/output/mm10_PRC.rds"
#   
#   chrs <- names(Mmusculus)[1:21]#2276796
#   cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
#   cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
#   
#   txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#   genes <- GenomicFeatures::genes(txdb)
#   TSS<-promoters(genes,upstream=0,downstream=0)
#   TSS$gene_symbol=AnnotationDbi::select(Mus.musculus,key=as.character(TSS$gene_id),
#                                         keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
#   TSS_break=dist_calc(TSS_break,TSS)
#  
#   out_name="../downstream/output/mm10_PRC_250bp.gff"
# }


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