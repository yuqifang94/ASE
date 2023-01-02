rm(list=ls())
source("mainFunctions_sub.R")
#Get features
getGeneralFeats_CpG <- function(CpGdir,enhancerDir='',chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Features included
  featureNickNames <- c("genome-wide","CpG island","CpG shore","CpG shelf","CpG open sea",
                        "gene body","exon","intron","intergenic")
  
  # Define list of feature GRs
  outGR <- GRangesList()
  GRtemp <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  
  outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
  #Redefining CpG islands using hidden Markov models 
  CpG_all <- readRDS(paste(CpGdir,"CpG_hg19.rds",sep=""))
  CpG_all<-setGenomeLengths(CpG_all)
  #Could also use UCSC genome browser CpG file,this file is from Jordi Abante
  cpg_islands <- readRDS(paste(CpGdir,"cpg_islands_hg19.rds",sep=""))
  cpg_islands<-subsetByOverlaps(CpG_all,cpg_islands)
  outGR[["CpG island"]] <- setGenomeLengths(cpg_islands)
  
  # extract the shore defined by 2000 bp upstream and downstream of cpg islands
  shore1 <- flank(cpg_islands, 2000)
  shore2 <- flank(cpg_islands,2000,FALSE)
  shore1_2 <- reduce(c(shore1,shore2))
  
  # extract the features (ranges) that are present in shores only and not in
  # cpg_islands (ie., shores not overlapping islands)
  cpgi_shores <- setdiff(shore1_2, cpg_islands)
  olap=findOverlaps(CpG_all,cpgi_shores)
  cpgi_shores<-subsetByOverlaps(CpG_all,cpgi_shores)
  outGR[["CpG shore"]] <- setGenomeLengths(cpgi_shores)
  
  # extract the shore defined by 4000 bp upstream and downstream of cpg islands
  shelves1 <- flank(cpg_islands, 4000)
  shelves2 <- flank(cpg_islands,4000,FALSE)
  shelves1_2 <- reduce(c(shelves1,shelves2))
  
  # create a set of ranges consisting CpG Islands, Shores
  island_shores <- c(cpg_islands,cpgi_shores)
  
  # extract the features (ranges) that are present in shelves only
  # and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
  cpgi_shelves <- setdiff(shelves1_2, island_shores)
  cpgi_shelves<-subsetByOverlaps(CpG_all,cpgi_shelves)
  outGR[["CpG shelf"]] <- setGenomeLengths(cpgi_shelves)
  
  # Open sea
  open_sea <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]],outGR[["CpG shelf"]]))
  open_sea<-subsetByOverlaps(CpG_all,open_sea)
  outGR[["CpG open sea"]] <- setGenomeLengths(open_sea)
  
  # Enhancers 
  #enhancers <- import.bed(paste(enhancerDir,"enhancers.bed",sep=""))[,c()]
  
  #outGR[["enhancer"]] <- setGenomeLengths(enhancers)
  
  # Other generic features
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- GenomicFeatures::genes(txdb)
  outGR[["gene body"]] <- setGenomeLengths(genes)
  exons <- GenomicFeatures::exons(txdb)
  outGR[["exon"]] <- setGenomeLengths(exons[,c()])
  introns <- GenomicFeatures::intronicParts(txdb)
  outGR[["intron"]] <- setGenomeLengths(introns[,c()])
  intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
  outGR[["intergenic"]] <- setGenomeLengths(intergenic)
  #Use annotation hub for TSS, promoter have something to do with strand
  proms <- promoters(genes,upstream=2000,downstream=1000)#ask Michael about this
  outGR[["promoter"]] <- setGenomeLengths(proms)
  TSS<-promoters(genes,upstream=0,downstream=0)
  outGR[["TSS"]] <- setGenomeLengths(TSS)
  # Gene name mapping
  geneBodyNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["gene body"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["gene body"]]$gene_name <- geneBodyNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["promoter"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["promoter"]]$gene_name <- promNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["TSS"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["TSS"]]$gene_name <- promNameMap$SYMBOL
  # Return
  return(outGR)
  
}
CpG_hg19=getCpgSitesHg19()
saveRDS(CpG_hg19,"../downstream/input/human_analysis/CpG_hg19.rds")
saveRDS(getGeneralFeats_CpG("../downstream/input/human_analysis/"),genomic_features_file)