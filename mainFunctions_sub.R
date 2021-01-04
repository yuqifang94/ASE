# Dependencies in use -----------------------------------------------------
if (!requireNamespace("rtracklayer", quietly = TRUE))
{BiocManager::install("rtracklayer")}
library(rtracklayer)
if (!requireNamespace("GenomicRanges", quietly = TRUE))
{BiocManager::install("GenomicRanges")}
library(GenomicRanges)
if (!requireNamespace("data.table", quietly = TRUE))
{BiocManager::install("data.table")}
library(data.table)
if (!requireNamespace("Gmisc", quietly = TRUE))
{BiocManager::install("Gmisc")}
library(Gmisc)
if (!requireNamespace("matrixStats", quietly = TRUE))
{BiocManager::install("matrixStats")}
library(matrixStats)
if (!requireNamespace("parallel", quietly = TRUE))
{BiocManager::install("parallel")}
library(parallel)
if (!requireNamespace("readr", quietly = TRUE))
{install.packages("readr")}
library(readr)
if (!requireNamespace("readxl", quietly = TRUE))
{install.packages("readxl")}
library(readxl)
if (!requireNamespace("reshape2", quietly = TRUE))
{BiocManager::install("reshape2")}
library(reshape2)

#plot 
if (!requireNamespace("pheatmap", quietly = TRUE))
{BiocManager::install("pheatmap")}
library(pheatmap)
if (!requireNamespace("ggfortify", quietly = TRUE))
{BiocManager::install("ggfortify")}
library(ggfortify)
if (!requireNamespace("RColorBrewer", quietly = TRUE))
{BiocManager::install("RColorBrewer")}
library(RColorBrewer)
if (!requireNamespace("ggplot2", quietly = TRUE))
{BiocManager::install("ggplot2")}
library(ggplot2)
if (!requireNamespace("ggpubr", quietly = TRUE))
{BiocManager::install("ggpubr")}
library(ggpubr)
if (!requireNamespace("shiny", quietly = TRUE))
{BiocManager::install("shiny",version="devel")}
library(shiny)
if (!requireNamespace("lattice", quietly = TRUE))
{install.packages("lattice")}
library(lattice)
if (!requireNamespace("gridExtra", quietly = TRUE))
{BiocManager::install("gridExtra")}
library(gridExtra)
if (!requireNamespace("grid", quietly = TRUE))
{BiocManager::install("grid")}
library(grid)
if (!requireNamespace("pander", quietly = TRUE))
{install.packages("pander")}
library(pander)
if (!requireNamespace("preprocessCore", quietly = TRUE))
{install.packages("preprocessCore")}
library(preprocessCore)
#hg19 genome 
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))
{BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")}
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
if (!requireNamespace("Homo.sapiens", quietly = TRUE))
{BiocManager::install("Homo.sapiens")}
library(Homo.sapiens)
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
{BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
library(BSgenome.Hsapiens.UCSC.hg19)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
{BiocManager::install("org.Hs.eg.db")}
library(org.Hs.eg.db)

#hg38 genome
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38", quietly = TRUE)){
  BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")}
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)

#mm10 genome related
if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)){
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
if (!requireNamespace("Mus.musculus", quietly = TRUE)){
  BiocManager::install("Mus.musculus")}
library(Mus.musculus)
if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)){
  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")}
library(BSgenome.Mmusculus.UCSC.mm10)

#genome tools
if (!requireNamespace("liftOver", quietly = TRUE)){
  BiocManager::install("liftOver")}
library(liftOver)

#Accessibility analysis
if (!requireNamespace("limma", quietly = TRUE))
{BiocManager::install("limma")}
library(limma)
if (!requireNamespace("edgeR", quietly = TRUE))
{BiocManager::install("edgeR")}
library(edgeR)
if (!requireNamespace("tximport", quietly = TRUE))
{BiocManager::install("tximport")}
library(tximport)
if (!requireNamespace("DESeq2", quietly = TRUE))
{BiocManager::install("DESeq2")}
library(DESeq2)

#Annotations
if (!requireNamespace("annotatr", quietly = TRUE))
{BiocManager::install("annotatr",INSTALL_opts = c('--no-lock'))}
library(annotatr)
if (!requireNamespace("AnnotationHub", quietly = TRUE))
{BiocManager::install("AnnotationHub")}
library(AnnotationHub)
if (!requireNamespace("GenomicFeatures", quietly = TRUE))
{BiocManager::install("GenomicFeatures")}
library(GenomicFeatures)
if (!requireNamespace("VariantAnnotation", quietly = TRUE))
{BiocManager::install("VariantAnnotation")}
library(VariantAnnotation)
if (!requireNamespace("topGO", quietly = TRUE))
{BiocManager::install("topGO")}
library(topGO)



# Not in use --------------------------------------------------------------
# 
# if (!requireNamespace("Repitools", quietly = TRUE))
# {BiocManager::install("Repitools")}
# library(Repitools)
# if(!require(psych)){install.packages("psych")}
# library('psych')
# if(!require(vcd)){install.packages("vcd")}
# if(!require(DescTools)){install.packages("DescTools")}
# if(!require(rcompanion)){install.packages("rcompanion")}
# if(!require(ggpubr)){install.packages("ggpubr")}
# library("ggpubr")
# if(!require(readxl)){install.packages("readxl")}
# library(readxl) 
# if(!require(ggpointdensity)){install.packages("ggpointdensity")}
# library(ggpointdensity)
# if (!requireNamespace("rethinking", quietly = TRUE))
# {install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
#   library(devtools)
#   devtools::install_github("rmcelreath/rethinking")}
# library(rethinking)
# if (!requireNamespace("gwasrapidd", quietly = TRUE)){
#   remotes::install_github("ramiromagno/gwasrapidd")}
# library(gwasrapidd)
# library(qusage)
# if (!requireNamespace("gwascat", quietly = TRUE)){
# BiocManager::install("gwascat")}
# library(gwascat)


#global pval cutoff file names
pval_cutoff=0.1 #Onuchic use 0.1
gff_in_file='../downstream/input/gff_in.rds'
variant_HetCpG_file='../downstream/input/variant_HetCpG_new.rds'
GR_file='../downstream/output/GRs_final1.rds'
GR_allele_file='../downstream/output/GRs_allele_final1.rds'
hetCpG_gff_file='../downstream/input/hetCpG_gff_final1.rds'
GR_merge_file="../downstream/output/GR_merge_final12_ls.rds"
variant_HetCpG_meta_file='../downstream/output/variant_HetCpG_meta_final1_ls.rds'
genomic_features_file="../downstream/input/genomic_features2020.rds"
NME_agnostic_file="../downstream/input/NME_allele_agnostic_merge_20k_homogeneous_excluding_dMML2.rds"
MML_agnostic_file="../downstream/input/MML_allele_agnostic_merge_20k_homogeneous2.rds"
motif_gene_file='../downstream/output/motif_all_JASPAR_default.rds' #For all SNP
#Checking chromosome name
chr_check<-function(gr_in){
  #Get plus one location
  chr_check=grepl('chr',seqlevels(gr_in))#checking if seqlevels contain chr
  if(any(!chr_check)){seqlevels(gr_in)[which(!chr_check)]=paste("chr",seqlevels(gr_in)[which(!chr_check)],sep='')}
  return(gr_in)
}
#Read in gff file
readAllGff <- function(inDir,subjects,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Loop over all GFF files
  outGR <- GRanges()
  #subjects <- c("H9","HUES64","skin03","HuFGM02","STL001","STL002","STL003","STL011")
  for (sub in subjects) {
    # Import GFF file
    cat('importing',sub,'\n')
    tmpGR <- readSubGff(inDir,sub,chrsOfInterest)
    start(tmpGR)=start(tmpGR)
    # Append to output GR
    outGR <- c(outGR,tmpGR)
  }
  
  # Return GR with all haplotypes
  return(outGR)
  
}
#Read in gff file for each subject
readSubGff <- function(inDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Import GFF
  outGR <- import.gff(paste(inDir,sub,"_het.cpelasm.gff",sep=""))
  
  # Retain a the required subset
  outGR <- outGR[,c("N","CpGs","hetCpGg1","hetCpGg2")]
  
  # Add subject metadata column
  outGR$Subject <- sub
  
  # Make N numeric and filter
  outGR$N <- as.numeric(outGR$N)
  
  # Add genome info 
  outGR <- setGenomeLengths(outGR,chrsOfInterest=chrsOfInterest)
  
  # Return GR with all haplotypes
  return(outGR)
  
}
###Here most current version of Het CpG analysis

#From vcf file, extract het CpG information
extractHetCpG<-function(vcfDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  cat('Processing subject:', sub,'\n')
  tt1=proc.time()[[3]]
  vcf <- readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19")
  
  gt <- as.vector(geno(vcf)$GT)
  vcf <- rowRanges(vcf)
  vcf$GT <- gt
  vcf$snpId <- paste(sub,seq(1:length(vcf)),sep="-")
  # Keep only relevant variables
  vcf <- vcf[,c("REF","ALT","GT","snpId")]
  vcf$REF <- as.character(vcf$REF)
  vcf$ALT <- as.character(unlist(vcf$ALT))
  names(vcf)=NULL
  # Delete labels
  vcf=het_CpG_df(vcf)
  vcf$sub=sub
  cat('Done processing',sub,'using:', proc.time()[[3]]-tt1,'\n')
  return(vcf)
  
}
#ID heterozygous CpG from vcf file
het_CpG_df<-function(var_gr){
  var_gr=chr_check(var_gr)
  plus_loc=as.character(Views(Hsapiens,GenomicRanges::shift(var_gr,1)))
  minus_loc=as.character(Views(Hsapiens,GenomicRanges::shift(var_gr,-1)))
  #get dinucleotide for ref, alt, plus and minus, find some examples region randomly: check if those match
  var_gr$REF_plus=paste0(var_gr$REF,plus_loc)
  var_gr$REF_minus=paste0(minus_loc,var_gr$REF)
  var_gr$ALT_plus=paste0(var_gr$ALT,plus_loc)
  var_gr$ALT_minus=paste0(minus_loc,var_gr$ALT)
  #get trinucleotide
  var_gr$REF_tri=paste0(minus_loc,var_gr$REF,plus_loc)
  var_gr$ALT_tri=paste0(minus_loc,var_gr$ALT,plus_loc)
  #currently not in use
  # #check if heterogygouze: note rowSum =2 have trinucleotide form CGG with ref =G alt =C
  # var_gr$npmCG=rowSums(as.data.table(mcols(var_gr))
  #                      [,.(str_count(REF_tri,pattern="CG"),str_count(ALT_tri,pattern="CG"))])
  # var_gr$HetCpG=var_gr$npmCG>0
  #Add CG content for genome1 and genome2 based on GT
  var_gr$genome1_plus=NA
  var_gr$genome1_minus=NA
  var_gr$genome1_tri=NA
  var_gr$genome2_plus=NA
  var_gr$genome2_minus=NA
  var_gr$genome2_tri=NA
  ####Fix the issue with genome file order does not equal to the ref/alt order, calculate the CG or het CG in genome 1 based on actual ref/alt order
  #Genome1
  var_gr$genome1_plus[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_plus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome1_minus[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_minus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome1_tri[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_tri[var_gr$GT %in% c("0/1","0|1")]
  
  var_gr$genome1_plus[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_plus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome1_minus[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_minus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome1_tri[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_tri[var_gr$GT %in% c("1/0","1|0")]
  #Genome2
  var_gr$genome2_plus[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_plus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome2_minus[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_minus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome2_tri[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_tri[var_gr$GT %in% c("0/1","0|1")]
  
  var_gr$genome2_plus[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_plus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome2_minus[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_minus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome2_tri[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_tri[var_gr$GT %in% c("1/0","1|0")]
  return(var_gr)
}

# Function to set genome and chromosome lengths to a GR object
setGenomeLengths <- function(GR,chrsOfInterest=paste("chr",1:22,sep=""),genome_in="hg19"){
  # Get genome info
  GR=chr_check(GR)
  genome.seqinfo <- seqinfo(getBSgenome(genome_in))
  genome.seqinfo <- genome.seqinfo[chrsOfInterest]
  GR <- GR[seqnames(GR) %in% chrsOfInterest]
  genome(GR) <- genome(genome.seqinfo)
  seqlevels(GR) <- seqlevels(genome.seqinfo)
  seqlengths(GR) <- seqlengths(genome.seqinfo)
  
  return(GR)
}
#read in all sample tissue diff
import.subject<-function(inDir,calc='diff'){
  #for calc: diff -> dMML etc, allele -> NME etc
  bed_in=dir(inDir,pattern="bedGraph")
  sample_in=unique(sub('_phased.*','',bed_in))
  GRs=GRanges()
  for (sp in sample_in) {
    subjects=sub('_.*','',sp)
    tissues=sub(paste0(subjects,"_"),'',sp)
    # Print sample being loaded
    cat("Loading:",subjects,'-',tissues,'\n')
    if (calc=='diff'){
      GR.in=read.diffGR(subjects,tissues,inDir)
    }else if(calc=='allele'){
      GR.in=read.alleleGR(subjects,tissues,inDir)
    }else {cat('Wrong calc \n')}
    GRs=c(GRs,GR.in)
  }
  return(GRs)
  
}
#Function to read in single GR object:
read.diffGR<-function(subjects,tissues,inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  #Initialization
  GRs=GRanges()
  #Make sure the inputs are unique
  # dmml
  filename=paste(inDir,subjects,"_",tissues,"_phased_tmml_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dMML')
  
  GRs <- c(GRs,GR.in)
  # dnme
  filename=paste(inDir,subjects,"_",tissues,"_phased_tnme_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dNME')
  GRs <- c(GRs,GR.in)
  # uc
  filename=paste(inDir,subjects,"_",tissues,"_phased_tpdm_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'UC')
  GRs <- c(GRs,GR.in)
  #Check if pvalue available
  GRs <- GRs[!is.na(GRs$pvalue)]
  # Add sample field
  GRs$Sample <- paste(tissues,"-",GRs$Subject)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  return(GRs)
}
#file_ends can be 
#c('dmml_pvals,dnme_pvals,uc_pvals')
#c('mml1,mml2,nme1,nme2')
read.alleleGR<-function(subjects,tissues,inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  #MML1
  GRs=GRanges()
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'MML')
  GR.in$Genome="1"
  GRs <- c(GRs, GR.in)
  #MML2
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'MML')
  GR.in$Genome="2"
  GRs <- c(GRs, GR.in)
  #NME1
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="1"
  GRs <- c(GRs, GR.in)
  #NME2
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="2"
  GRs <- c(GRs, GR.in)
  # Add sample field
  GRs$Sample <- paste(tissues,"-",subjects)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  GRs$K=GRs$'NA.1'
  GRs$'NA.1'=NULL
  return(GRs)
}
#Read in each bed file, for new method, no need to resize
import.ASMbed<-function(subjects,tissues,filename,pvalue=TRUE,Statistic,chrsOfInterest=paste("chr",1:22,sep="")){
  GR <- import.bedGraph(filename)
  #fit  bedGraph reads, import.bedGraph will remove 1 from start
  #Check if all files are 0 based or 1 based? Check on genome browser, UCSC: check SNP location (0 based)
  start(GR)=start(GR)-1
  GR_out=chr_check(GR)
  GR_out$ASM=NULL
  GR_out$Data=NULL
  GR_out$Subject <- subjects
  GR_out$tissue<-tissues
  GR_out$Statistic <- Statistic
  GR_out$Value<- GR$score
  #for differential analysis, make sure the GR_out are unique
  if(pvalue){
    GR_out$pvalue <- as.numeric(GR$NA.)
    GR_out <- GR_out[!duplicated(GR[,c()])]
    }
  else{GR_out$N <- GR$NA.} #not diff for new samples
  return(GR_out)
}
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
  #Could also use UCSC genome browser CpG file
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
  proms <- promoters(genes,upstream=2000,downstream=1000)
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

#Put hetCpG count into each sample in gr_allele

##redo instead of adding allele information to allele CpG, add it to GR_merge
#This is mainly modifed to fit the output of CPEL ASM

stat_merge<-function(gr_in,allele_in,vcf_in,CpG){
  #Check merge behavior, 
  dMML=gr_in[gr_in$Statistic=="dMML"]
  dNME=gr_in[gr_in$Statistic=="dNME"]
  UC=gr_in[gr_in$Statistic=="UC"]
  gr=unique(granges(gr_in))
  olap_dMML=findOverlaps(gr,dMML,type="equal")
  gr$dMML[queryHits(olap_dMML)]=dMML$Value[subjectHits(olap_dMML)]
  gr$dMML_pval[queryHits(olap_dMML)]=dMML$pvalue[subjectHits(olap_dMML)]
  gr$Sample=unique(gr_in$Sample)
  
  olap_dNME=findOverlaps(gr,dNME,type="equal")
  gr$dNME[queryHits(olap_dNME)]=dNME$Value[subjectHits(olap_dNME)]
  gr$dNME_pval[queryHits(olap_dNME)]=dNME$pvalue[subjectHits(olap_dNME)]
  
  olap_UC=findOverlaps(gr,UC,type="equal")
  gr$UC[queryHits(olap_UC)]=UC$Value[subjectHits(olap_UC)]
  gr$UC_pval[queryHits(olap_UC)]=UC$pvalue[subjectHits(olap_UC)]
  
  gr$NME1=gr$MML1=gr$NME2=gr$MML2=gr$N=NA
  #genome 1 NME
  olap_NME1=findOverlaps(gr,allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME'],type="equal")
  gr$NME1[queryHits(olap_NME1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME']$Value[subjectHits(olap_NME1)]
  gr$N[queryHits(olap_NME1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME']$N[subjectHits(olap_NME1)]
  #genome 2 NME
  olap_NME2=findOverlaps(gr,allele_in[allele_in$Genome==2 & allele_in$Statistic=='NME'],type="equal")
  gr$NME2[queryHits(olap_NME2)]=allele_in[allele_in$Genome==2 & allele_in$Statistic=='NME']$Value[subjectHits(olap_NME2)]
  #genome 1 MML
  olap_MML1=findOverlaps(gr,allele_in[allele_in$Genome==1 & allele_in$Statistic=='MML'],type="equal")
  gr$MML1[queryHits(olap_MML1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='MML']$Value[subjectHits(olap_MML1)]
  #genome 1 MML
  olap_MML2=findOverlaps(gr,allele_in[allele_in$Genome==2 & allele_in$Statistic=='MML'],type="equal")
  gr$MML2[queryHits(olap_MML2)]=allele_in[allele_in$Genome==2 & allele_in$Statistic=='MML']$Value[subjectHits(olap_MML2)]
  gr$Subject=unique(gr_in$Subject)
  gr$tissue=unique(gr_in$tissue)
  #Add g1cg etc information to the gr
  gr=add_hetCPG(gr,vcf_in,CpG)
  return(gr)
}

add_hetCPG<-function(gr,vcf_in,CpG){
  olap=findOverlaps(vcf_in,gr,type='within',select='all')
  #Count number of CG, here g1CG=genome1 and g2CG=genome2, however, we need to split based on the GT, also ref
  vcf_in$g1CG=as.numeric(grepl("CG",vcf_in$genome1_tri))
  vcf_in$g2CG=as.numeric(grepl("CG",vcf_in$genome2_tri))
  vcf_in$refCG=as.numeric(grepl("CG",vcf_in$REF_tri))
  vcf_in$altCG=as.numeric(grepl("CG",vcf_in$ALT_tri))
  #Count number of Het CpG here *Check 
  df_sub=data.table(subjHits=subjectHits(olap),
                    g1CG=vcf_in$g1CG[queryHits(olap)],g2CG=vcf_in$g2CG[queryHits(olap)],
                    refCG=vcf_in$refCG[queryHits(olap)],altCG=vcf_in$altCG[queryHits(olap)])
  
  agg_sub=df_sub[,.(g1CG=sum(g1CG),g2CG=sum(g2CG),refCG=sum(refCG),altCG=sum(altCG)),by=subjHits]
  gr$g1CG[agg_sub$subjHits]=agg_sub$g1CG#agg_sub$subjHits is the unique subject hits
  gr$g2CG[agg_sub$subjHits]=agg_sub$g2CG
  gr$refCG[agg_sub$subjHits]=agg_sub$refCG#agg_sub$subjHits is the unique subject hits
  gr$altCG[agg_sub$subjHits]=agg_sub$altCG
  gr$N_hg19=countOverlaps(gr,CpG)
  #count number of CG lost from reference
  gr$N_nonhet=gr$N_hg19-countOverlaps(gr,vcf_in[(vcf_in$refCG-vcf_in$altCG)>0])
  
  return(gr)
}


#add genes to GR
add_gene_GR<-function(GR_merge_in,feature_in,feature_names){
  mcols(GR_merge_in)[[feature_names]]=NA
  olap_promoter=findOverlaps(GR_merge_in,feature_in)
  df_idx=data.table(qt=queryHits(olap_promoter),
                    genes=as.character(feature_in$gene_name[subjectHits(olap_promoter)]),
                    stringsAsFactors = F)
  df_idx=df_idx[,list(genes=list(genes)),by=list(df_idx$qt)]
  colnames(df_idx)=c('qt','genes')
  df_idx$genes=lapply(df_idx$genes,as.character)
  mcols(GR_merge_in)[[feature_names]][df_idx$qt]=df_idx$genes
  return(GR_merge_in)
  
}

#Count Number of Het CpG at extened regions, each extend 500 bp
hetCGallele_merged<-function(sub,gr_merge,CpG,vcf_in,gene_size=500){
  cat('Analyzing',sub,'\n')
  t1=proc.time()[[3]]
  #Import vcf file
  sub_vcf=vcf_in[[sub]]
  sub_allele=gr_merge[gr_merge$Subject==sub]
  #gr_allele got resized with start and end +1, use +2 to include equal start & end, no longer needed for new output?
  #sub_het=resize(sub_het, width(sub_het) + 4, fix="center")
  sub_allele$gff_size=width(sub_allele)
  gr_out=GR_resize_merged(sub_allele,CpG,sub_vcf,gene_size=gene_size,sub)#Check for calculating density
  cat('Finish analyzing',sub,proc.time()[[3]]-t1,'\n')
  return(gr_out)
}
#Generate density from GR allele, at least need 200 bp for density
#Count number of hetCG at each allele GR from result, add gff size, and N
GR_resize_merged<-function(GR.in,CpG_sites,hetCpG,sub,gene_size=500){
  ##From definitino of CpG island, use 200 bp regions
  GR.extend=resize(GR.in,width=width(GR.in)+gene_size,fix='center')
  GR.in$CG_hg19_extend=countOverlaps(GR.extend,CpG_sites)
  #any SNP contain CG in ref genome
  GR.in$CG_nonhet_extend=GR.in$CG_hg19_extend-countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$REF_tri)])
  gr_seq=getSeq(Hsapiens,GR.extend,as.character=T)
  GR.in$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
  #Count CpG in genome 1
  GR.in$CG_allele_extend_g1=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$genome1_tri)])
  GR.in$CG_allele_extend_g2=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$genome2_tri)])
  GR.in$gff_size_extend=width(GR.extend)
  return(GR.in)
}

#Count expted CG numbers and return data.frame
countCGOR<-function(x){ #x=input seq
  #calculate Odds ratio for expected CG vs actual CG
  #Expected CG number C * number G/total length
  # Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes". Journal of Molecular Biology.
  #Wiki, actual: ((number of C + Number of G)/2)^2/length of genomics Normalized CpG content, whole genome ~25%
  NC=countPattern('C',x)
  NG=countPattern('G',x)
  #NCG=countPattern('CG',x)
  CG_exp=NC*NG/nchar(x) #PC*PG*length
  CG_exp_norm=((NC+NG)/2)^2/nchar(x) #assuming PC=PG = (NC+NG)/2/length
  return(CG_exp_norm)
}


#Give each SNP an ASM information for each subject
variant_meta<- function(subj,variant_in,GR_in){ #variant_in for each subject, GR_in for each subject
  cat('Processing',subj,'\n')
  GR_in_subj=GR_in[GR_in$Subject==subj]
  variant_in_subj=variant_in[[subj]]
  sp=unique(GR_in_subj$Sample)
  gr_out_sp = GRanges()
  for (sps in sp){
    cat('Processing',sps,'\n')
    gr_out_sp=c(gr_out_sp,variant_meta_sp(variant_in_subj,GR_in_subj[GR_in_subj$Sample==sps]))
    
  }
  return(gr_out_sp)
  
  
}

#For each sample, assign NME to SNP within each region
variant_meta_sp <-function(variant_subj,GR_sp){
  olap=findOverlaps(variant_subj,GR_sp,maxgap =0,type='within')
  gr_out=variant_subj[queryHits(olap)]
  #Find GR_merge olap
  olap=subjectHits(olap)
  mcols(gr_out)=cbind(mcols(gr_out),mcols(GR_sp)[olap,])
  gr_out$refNME=NA
  gr_out$altNME=NA
  gr_out$refMML=NA
  gr_out$altMML=NA
  #NME
  gr_out$refNME[gr_out$GT %in% c("0/1","0|1")]= gr_out$NME1[gr_out$GT %in% c("0/1","0|1")]
  gr_out$altNME[gr_out$GT %in% c("0/1","0|1")]= gr_out$NME2[gr_out$GT %in% c("0/1","0|1")]
  gr_out$refNME[gr_out$GT %in% c("1/0","1|0")]= gr_out$NME2[gr_out$GT %in% c("1/0","1|0")]
  gr_out$altNME[gr_out$GT %in% c("1/0","1|0")]= gr_out$NME1[gr_out$GT %in% c("1/0","1|0")]
  
  #MML
  gr_out$refMML[gr_out$GT %in% c("0/1","0|1")]= gr_out$MML1[gr_out$GT %in% c("0/1","0|1")]
  gr_out$altMML[gr_out$GT %in% c("0/1","0|1")]= gr_out$MML2[gr_out$GT %in% c("0/1","0|1")]
  gr_out$refMML[gr_out$GT %in% c("1/0","1|0")]= gr_out$MML2[gr_out$GT %in% c("1/0","1|0")]
  gr_out$altMML[gr_out$GT %in% c("1/0","1|0")]= gr_out$MML1[gr_out$GT %in% c("1/0","1|0")]
  
  gr_out$Statistic=GR_sp$Statistic[olap]
  
  # gr_out$HetCpG=FALSE
  # #This may be different from previous result?
  # gr_out$HetCpG=((gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & !(gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG')) |
  #   (!(gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & (gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG'))
  return(gr_out)
}

#Read in allele-agnositc model
read.agnostic<-function(file_in,GR_merge_in=NULL,allele_include=T,olap_type="any",all_regions=NA,sample_in=NA,hyper_var_file=NA){
  stat=toupper(sub('.*_','',sub('.bedGraph','',file_in)))
  informME_in=read.bedGraph.informME(file_in)
  if(length(GR_merge_in)>0){
    GR_merge_in=GR_merge_in[GR_merge_in$Sample==sample_in]
  #Find overlapped region
    olap=findOverlaps(informME_in,GR_merge_in,type=olap_type)
    if(length(olap)>0) {informME_in=informME_in[-queryHits(olap)]}
  }
  #add GR_merge data
  if(allele_include){
  
    cat('Percent overlap with dNME region:',length(unique(queryHits(olap)))/length(informME_in)*100,'%\n')
    informME_in$score_original=informME_in$score
 
    #replace value instead of remove regions
    olap=findOverlaps(all_regions,GR_merge_in)
    asm_replacement=data.table(idx=queryHits(olap), 
                               score= rowMeans(as.matrix(elementMetadata(GR_merge_in)[paste(stat,c('1','2'),sep='')]))[subjectHits(olap)],
                               dMML=GR_merge_in$dMML[subjectHits(olap)],
                               dMML_pval=GR_merge_in$dMML_pval[subjectHits(olap)],
                               N=GR_merge_in$N[subjectHits(olap)])
    asm_replacement=asm_replacement[,list(score=mean(score),dMML=mean(dMML),dMML_pval=mean(dMML_pval)),by=list(idx)]
    all_regions=all_regions[asm_replacement$idx]
    all_regions$score=asm_replacement$score
    all_regions$K=NA
    all_regions$dMML=asm_replacement$dMML
    all_regions$dMML_pval=asm_replacement$dMML_pval
    informME_in=c(informME_in,all_regions)
  }
 informME_in=informME_in[!is.infinite(informME_in$score)]
  informME_in$Sample=sample_in
  informME_in$hyper_var_fn=hyper_var_file
  return(informME_in)

}
#read in hypervaribility from file
read_hypervar<-function(hyper_var_in){
  hyper_var=data.table()
  for (fn in hyper_var_in){
    fn_in=readRDS(fn)
    fn_in$gene_name=rownames(fn_in)
    fn_in=as.data.table(fn_in)
    rownames(fn_in)=NULL
    hyper_var=rbind(hyper_var,fn_in)
  }
  
  #Find the dataset that have most regions overlap with NME regions, check correlation?
  hyper_var=hyper_var[,list(mean=mean(mean,na.rm=T),var=mean(var,na.rm=T),
                            hypervar_var=mean(hypervar_var,na.rm=T),hypervar_logvar=mean(hypervar_logvar,na.rm=T)),
                      by=list(hyper_var$gene_name)]
  colnames(hyper_var)[1]='gene_name'
  return(hyper_var)
}
#Calculating NME vs hypervaribility
dist_plot_calc<-function(informME_in,genes_hypervar,genomic_features,enhancer=FALSE){
  if(enhancer){informME_in_dist=informME_in}else{informME_in_dist=dist_calc(informME_in,genomic_features$TSS)}
  mcols(informME_in_dist)=cbind(mcols(informME_in_dist),genes_hypervar[,-1][match(informME_in_dist$gene,genes_hypervar$gene_name)])
  
  return(informME_in_dist)
}
#calculate distance to the genomic features
dist_calc<-function(gr_in,gr_feature){
  
  dist_nearest=nearest(resize(gr_in,1,fix="center"),resize(gr_feature,1,fix="start"))
  
  gr_non_na=which(!is.na(dist_nearest))
  gr_feature=gr_feature[dist_nearest[gr_non_na]]
  gr_in=gr_in[gr_non_na]
  gr_in$dist=NA
  gr_in$gene=NA
  
  sgn <- as.integer(ifelse(strand(gr_feature)=="+",1,-1))
  gr_in$dist=sgn*(start(resize(gr_in,1,fix="center"))-start(resize(gr_feature,1,fix="start")))
  gr_in$gene=gr_feature$gene_name
  
  return(gr_in)
}
read.bedGraph.informME<-function(file_in){
  
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    return(informME_in)
  }
}
#read in mouse enhancer
read_chromHMM_bed<-function(bed_dir,rep){
  bed_out=GRanges()
  for(fn in dir(bed_dir,pattern='.bed.gz')){
    #get sample name etc
    fn_split=strsplit(fn,'_')[[1]]
    stage=gsub('e','E',fn_split[1])
    tissue=gsub('facial-prominence','EFP',fn_split[2])
    tissue=gsub('neural-tube','NT',tissue)
    bed_in=read.table(paste(bed_dir,fn,sep=''))
    colnames(bed_in)=c('chr','start','end','chrom_num','chrom_state')
    bed_in=makeGRangesFromDataFrame(bed_in,keep.extra.columns = T)
    bed_in$stage=stage
    bed_in$tissue=tissue
    bed_in$rep=rep
    bed_in$Sample=paste(stage,tissue,sep='-')
    bed_out=c(bed_out,bed_in)
  }
  return(bed_out)
}
#Enrichment analysis for genomic features in dMML and dNME
testEnrichmentFeature_stat<-function(dataGR,featureGR,maxgap=0,output_ft=1){
  # Find ranges overlapping with feature
  olaps <- findOverlaps(dataGR,featureGR,type="any",select="all",maxgap = maxgap)
  
  indFeature <- queryHits(olaps)
  featurestatistic <- dataGR[indFeature]
  complementarystatistic <- dataGR[-indFeature]
  
  # Enrichment of in feature
  #featurestatistic <- featureData[featureData$Statistic==statistic]
  #complementarystatistic <- complementaryData[complementaryData$Statistic==statistic]
  contTablestatistic <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTablestatistic) <- c("Feature","Complementary")
  contTablestatistic[1,]$ASM <- sum(featurestatistic$ASM=="Yes")
  contTablestatistic[1,]$nonASM <- sum(featurestatistic$ASM=="No")
  contTablestatistic[2,]$ASM <- sum(complementarystatistic$ASM=="Yes")
  contTablestatistic[2,]$nonASM <- sum(complementarystatistic$ASM=="No")
  #print(contTablestatistic)
  #Return overlap >=output_ft
  if(contTablestatistic[1,1]>=output_ft){
  return(list(contTablestatistic,fisher.test(contTablestatistic)))
  }
}
#MAE enrich
MAE_enrich<-function(GR_merge,pval_cutoff,genes='genes_promoter',stat='dMML_pval',MAE=MAE){
  #GR_merge=GR_merge[!is.na(GR_merge$genes_promoter)]
  GR_merge=elementMetadata(GR_merge)
  stat_gene=GR_merge[[genes]][GR_merge[[stat]]<=pval_cutoff]
  non_stat_gene=GR_merge[[genes]][GR_merge[[stat]]>pval_cutoff]
  stat_MAE=sum(unlist(lapply(stat_gene,function(x) any(x %in% MAE))))
  nonstat_MAE=sum(unlist(lapply(non_stat_gene,function(x) any(x %in% MAE)))) 
  nonstat_nonMAE=sum(!unlist(lapply(non_stat_gene,function(x) any(x %in% MAE))))
  stat_nonMAE=sum(!unlist(lapply(stat_gene,function(x) any(x %in% MAE))))
  print(matrix(c(stat_MAE,stat_nonMAE,nonstat_MAE,nonstat_nonMAE),nrow=2))
  ft=fisher.test(matrix(c(stat_MAE,stat_nonMAE,nonstat_MAE,nonstat_nonMAE),nrow=2))
  return(data.frame(OR=ft$estimate,pvalue=ft$p.value,lowerCI= ft$conf.int[1],upperCI=ft$conf.int[2]))
}


#ChromHMM* check
ENCODE_to_sample<-function(sample_in){
  #Assign code to each sample
  ENCODE_number_to_sample=data.frame(sample=sample_in,stringsAsFactors = F)
  ENCODE_number_to_sample$subject=unlist(lapply(ENCODE_number_to_sample$sample,function(x) strsplit(x,' - ')[[1]][2]))
  ENCODE_number_to_sample$tissue=unlist(lapply(ENCODE_number_to_sample$sample,function(x) gsub(c('_single','_paired'),'',strsplit(x,' - ')[[1]][1])))
  #Assign known code to given sample
  ENCODE_number_to_sample$ENCODE=NA
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$subject=='H1']='E003'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='stem_27_undifferentiated_paired - HUES64']='E016'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='42_embryonic_stem_cell_single - H9']='E008'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='foreskin_melanocyte_paired - skin03']='E061'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='foreskin_keratinocyte_paired - skin03']='E058'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Small_Intestine']='E109'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Gastric']='E094'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Left_Ventricle']='E095'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Lung']='E096'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Psoas_Muscle']='E100'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Right_Ventricle']='E105'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Sigmoid_Colon']='E106'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Spleen']='E113'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Thymus']='E112'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue %in% c('Adipose','Adipose_Tissue')]='E063'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Aorta']='E065'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Esophagus']='E079'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Adrenal_Gland']='E080'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Ovary']='E097'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Pancreas']='E087'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Liver']='E066'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Right_Atrium']='E014'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='ectoderm_paired - HUES64']='E012'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='endoerm_27_paired - HUES64']='E011'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='mesoderm_23_paired - HUES64']='E013'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='brain_germinal_matrix_tissue_paired - HuFGM02']='E070'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='Brain_substantia_nigra_paired - 112']='E074'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='Brain_Hippocampus_middle_paired - 149']='E071'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='Brain_Hippocampus_middle_paired - 150']='E071'
  
  return(ENCODE_number_to_sample)
}
##Use  CMH test for chromHMM analysis
chromHMM_OR<-function(GR_merge,chromHMM,sample_name,pval_cutoff=0.1,stat="dNME_pval"){
  
  GR_merge_sp=GR_merge[GR_merge$Sample%in%sample_name]
  GR_merge_sp$ASM="No"
  GR_merge_sp$ASM[elementMetadata(GR_merge_sp)[,stat]<=pval_cutoff]="Yes"
  out_df=data.frame()
  count_table=list()
  count_table_N=data.frame()
  for(states in unique(chromHMM$name)){
    OR=testEnrichmentFeature_stat(GR_merge_sp,chromHMM[chromHMM$name==states])
    #Get contengency table
    count_table[[states]]=OR[[1]]
    OR=OR[[2]]
    #result directly from fisher.test for each state
    out_df=rbind(out_df,
                 data.frame(state=states,OR=OR$estimate,p_value=OR$p.value,
                            lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
    # }
  }
  return(list(out_df,count_table))
}
#Combine contengency table
chromHMM_combine<-function(chromHMM_in){
  cont_table_all=list()
  #ChromHMM, list object, each is a sample, within sample ,each is a state
  for (states in names(chromHMM_in[[1]][[2]])){
    #extract 2x2 table for each states for each sample, return a list of sample with its contengency table
    chromHMM_in_cont=lapply(chromHMM_in, function(x,states) x[[2]][[states]],states=states)
    #For each state, construct a CMH table, 4 columns
    cont_table_all[[states]]=do.call(rbind,lapply(seq_along(chromHMM_in_cont),melt_cont,cont_in=chromHMM_in_cont))
    print(states)
    
  }
  cont_table_all_CMH=lapply(cont_table_all,CMH_test)
  cont_table_all_CMH_df=data.frame(states=names(cont_table_all_CMH), 
                                   OR=unlist(lapply(cont_table_all_CMH, function(x) x$estimate[[1]])), 
                                   p_value=unlist(lapply(cont_table_all_CMH, function(x) x$p.value)),
                                   lower_CI=unlist(lapply(cont_table_all_CMH,function(x) x$conf.int[1])),
                                   upper_CI=unlist(lapply(cont_table_all_CMH,function(x) x$conf.int[2])))
  
  return(cont_table_all_CMH_df)
}

melt_cont<-function(i,cont_in){
  sp=names(cont_in)[i]
  # cont_out=cont_in[[i]]
  # cont_out$subject=factor(paste(sp,cont_out$N,sep='_'))
  # return(cont_out)
  count_out=as.numeric(unlist(cont_in[[i]]))
  if(length(count_out)>0){
  #print(count_out)
  #cutoff of smallest number in the table
  if (all(count_out>0)){
    return(data.frame(subject=factor(rep(sp,4)),
                      ASM=factor(c('ASM','ASM','Non-ASM','Non-ASM')),
                      feature=factor(c('Feature','Non_feature','Feature','Non_feature')),
                      count=count_out
    ))
  }
  }
}

CMH_test<-function(df_in,CMH_eqn=count~ASM+feature+subject){
  if (nrow(df_in)>0){
    CMH_table=xtabs(CMH_eqn,data=df_in)
    
    if (length(dim(CMH_table)[3]!=0)){
      
      if(dim(CMH_table)[3]>1){
        fs=mantelhaen.test(CMH_table)
      }
      else if(dim(CMH_table)[3]==1){
        #print(as.matrix(CMH_table[,,1]))
        #print(count_all)
        fs=fisher.test(as.matrix(CMH_table[,,1]))
        
      }
      return(fs)
    }
    
  }
}
#SNP OR
OR_calc<-function(tb_in,SNP,SNP_name,pval_cutoff=NA){
  
  larger_NME_SNP=sum(tb_in[[SNP_name]]==SNP&(tb_in$refNME<tb_in$altNME))
  larger_NME_nonSNP=sum(tb_in[[SNP_name]]!=SNP&(tb_in$refNME<tb_in$altNME))
  lower_NME_SNP=sum(tb_in[[SNP_name]]==SNP&(tb_in$refNME>tb_in$altNME))
  lower_NME_nonSNP=sum(tb_in[[SNP_name]]!=SNP&(tb_in$refNME>tb_in$altNME))
  cont_table=matrix(c(larger_NME_SNP,larger_NME_nonSNP,lower_NME_SNP,lower_NME_nonSNP),nrow=2)
  OR=fisher.test(cont_table)
  #print(cont_table)
  return(data.table(OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.int[1],upperCI=OR$conf.int[2],SNP=SNP))
}

dist_plot_run<-function(informME_in_dist,theme_glob,ylab,stat_in,cutoff=pval_cutoff,dir='../downstream/output/graphs/Figure3/'){
  
  #informME_in_dist=informME_in_dist[-which(informME_in_dist$dMML_pval<=cutoff)]
  
  informME_in_dist$exp_stat=informME_in_dist[[stat_in]]
  informME_in_dist=informME_in_dist[!is.na(exp_stat)]
  plot_informME_dat=informME_in_dist[,list(score=score,stat_in=stat_in,dist=dist,exp_stat=exp_stat,gene=gene,
                                           quant=findInterval(exp_stat,quantile(unique(data.table(gene=gene,exp_stat=exp_stat))$exp_stat,prob=c(0,0.25,0.5,0.75),na.rm=T)),
                                           hypervarquant=findInterval(exp_stat,quantile(unique(data.table(gene=gene,exp_stat=exp_stat))$exp_stat,prob=seq(0.01,1,0.01),na.rm=T))/100),
                                     #scorequant001=findInterval(exp_stat,quantile(score,prob=c(0.01,1,0.01),na.rm=T))/100
                                     by=list(Sample)]
  rm(informME_in_dist)
  quant=c("0-25%","25%-50%","50%-75%","75%-100%")
  plot_informME_dat$quant=quant[plot_informME_dat$quant]
  plot_informME_dat$quant=as.factor(plot_informME_dat$quant)
  dist_plot=list()
  #print(unique(plot_informME_dat[,c(3,4)]))
  pdf(paste0(dir,'Figure3A_',ylab,'_',stat_in,'_dist.pdf'),width=3.5,height=3.5)
  for(sp in unique(plot_informME_dat$Sample)){
  print(ggplot(plot_informME_dat[abs(dist)<=3000&Sample==sp],aes(x=dist,y=score,color=quant))+theme_glob+
    geom_smooth(size=1,na.rm=TRUE,se=TRUE)+ggtitle(sp)+
    theme(legend.position="bottom")+ labs(color = "quantile")+
    xlab("Distance to TSS")+ylab(ylab)+guides(color=guide_legend(nrow=2,byrow=TRUE)))
  }
 dev.off()
  #find where min NME is 
  #plot_informME_dat_median=plot_informME_dat[abs(dist)<=3000,list(median_score=median(score)),by=list(dist,quant)]
  #min_dist=plot_informME_dat_median[which.min(median_score)]$dist
  #plot_informME_dat_ft=plot_informME_dat[dist<=min_dist+250&dist>=min_dist-250]
  plot_informME_dat_ft=plot_informME_dat[dist>=0&dist<=500]
  #print(cor.test(plot_informME_dat_cor$score,plot_informME_dat_cor$exp_stat))
  #Use correlation between hypervaribility quantile and NME
  # plot_informME_dat_cor=plot_informME_dat_ft[,list(median_score=median(score)),
  #                                         by=list(exp_stat,gene,hypervarquant,Sample)]
  cor_mean=plot_informME_dat_ft[,list(cor=cor(score,hypervarquant),pvalue=cor.test(score,hypervarquant)$p.value),by=list(Sample)]
  write.csv(cor_mean,paste0(dir,'Figure3A_',ylab,'_',stat_in,'_500bp_tb.csv'),row.names = F)
  print(mean(cor_mean$cor))
  # print(cor.test(plot_informME_dat_cor$median_score,plot_informME_dat_cor$hypervarquant))
  heatmap_informME_dat=plot_informME_dat_ft[,list(median_score=median(score)),by=list(hypervarquant,Sample)]
  heatmap_informME_dat=heatmap_informME_dat[,list(hypervarquant=hypervarquant,NME=median_score,cor=cor(median_score,hypervarquant)),by=list(Sample)]
  heatmap_informME_dat=heatmap_informME_dat[order(heatmap_informME_dat$cor,decreasing = F),]
  heatmap_informME_dat$Sample=factor(heatmap_informME_dat$Sample,levels = unique(heatmap_informME_dat$Sample))
  heatmap_plot=ggplot(heatmap_informME_dat,aes(hypervarquant,Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
    xlab('quantile')+ylab('Sample')+theme_glob+theme(legend.position = 'bottom')
  ggsave(paste0(dir,'Figure3A_',ylab,'_',stat_in,'_heatmap.pdf'),heatmap_plot,device='pdf',width=7,height=7,units="in")
  # cor_plot=ggplot(plot_informME_dat,aes(x=exp_stat,y=dat))+
  #   geom_smooth(size=1,na.rm=TRUE,se=TRUE)+theme_glob+
  #   theme(legend.position="bottom")+ggtitle('Correlation to hypervaribility within 500 bp upstream from TSS')+geom_point(alpha=0.1)
  # 
  # print(cor_plot)
  return(plot_informME_dat)
}
direction_calc_enriched_subj<-function(motif_loc,variant_all,gene_in,pval_cutoff=0.1,stat="NME"){
  
  #print(motif_sig_df)
  # motif_sig_df=motif_sig_df[which(motif_sig_df$qval<=0.2&motif_sig_df$OR>1),]
  motif_direction_out=mclapply(gene_in,direction_enriched_sample,
                               variant_gene=variant_all,motif_gene_subj=motif_loc,pval_cutoff=pval_cutoff,stat=stat,mc.cores =1)
  
  
  
  return(do.call(rbind,motif_direction_out))
}

direction_enriched_sample<-function(tf,variant_gene,motif_gene_subj,pval_cutoff,nperm=0,stat="NME"){
  motif_gene_subj=motif_gene_subj[motif_gene_subj$geneSymbol==tf]
  variant_gene=variant_gene[mcols(variant_gene)[[paste0('d',stat,'_pval')]]<=pval_cutoff]
  olap=findOverlaps(variant_gene,motif_gene_subj)
  variant_gene=variant_gene[queryHits(olap)]
  variant_gene$alleleDiff=motif_gene_subj$alleleDiff[subjectHits(olap)]
  variant_gene=variant_gene[!is.na(variant_gene$alleleDiff)]
  
  #alleleDiff is calculated use ref - alt, prefer low ent ones
  variant_gene$stat_diff=mcols(variant_gene)[[paste0('alt',stat)]]-mcols(variant_gene)[[paste0('ref',stat)]]
  same_dir=sum(sign(variant_gene$alleleDiff)== sign(variant_gene$stat_diff),na.rm = TRUE)
  opposite_dir=sum(sign(variant_gene$alleleDiff)!= sign(variant_gene$stat_diff),na.rm = TRUE)
  # same_dir=sum(sign(variant_gene$alleleDiff)== sign(variant_gene$MML_diff),na.rm = TRUE)
  # opposite_dir=sum(sign(variant_gene$alleleDiff)!= sign(variant_gene$MML_diff),na.rm = TRUE)
  total_data=same_dir+opposite_dir
  
  variant_gene_df=data.frame(alleleDiff=sign(variant_gene$alleleDiff),stat_diff=sign(variant_gene$stat_diff))
  len_x=nrow(variant_gene_df)
  if(nperm>0){
    same_dir_perm=replicate(nperm,
                            sum(sample(variant_gene_df$alleleDiff,len_x,replace = F)==sample(variant_gene_df$stat_diff,len_x,replace = F)))
    same_dir_perm_prob=same_dir_perm/total_data
  }else{same_dir_perm_prob=-1}
  if(same_dir >0 &opposite_dir>0){
    
    binom=binom.test(same_dir,(same_dir+opposite_dir),0.5)
    #print(binom)
    # binom=summary(lm(abs(variant_gene$alleleDiff)~abs(variant_gene$dMML)))
    #binom$p.value=pf(binom$fstatistic[1],df1 = binom$fstatistic[2],df2 = binom$fstatistic[3],lower.tail = F)
    # return(data.frame(TF=unique(motif_gene_subj$geneSymbol),total_data=same_dir+opposite_dir,same_dir=same_dir,opposite_dir=opposite_dir,
    #                   binom.pval=binom$p.value,prob=binom$estimate[[1]],NSNP=length(variant_gene),stringsAsFactors = F))
    prob_binom=binom$estimate[[1]]
    if(nperm>0){binom.pval=sum(abs(same_dir_perm_prob-0.5)>=abs(prob_binom-0.5))/nperm}else{binom.pval=NA}
    # if(prob_binom>0.5){
    #   binom.pval=sum(same_dir_perm_prob>=(same_dir/total_data)|same_dir_perm_prob<=(1-(same_dir/total_data)))/nperm
    #   
    # }else if(prob_binom<0.5){binom.pval=sum(same_dir_perm_prob<=(same_dir/total_data)|same_dir_perm_prob>=(1-(same_dir/total_data)))/nperm}
    #cat(tf,':',binom.pval,'\n')
    return(data.table(TF=unique(motif_gene_subj$geneSymbol),total_data=total_data,same_dir=same_dir,opposite_dir=opposite_dir,
                      binom.pval_perm=binom.pval,binom.pval=binom$p.value,prob=prob_binom,NSNP=length(variant_gene),stringsAsFactors = F))
  }
}
OMIM_annotation<-function(motif_in,OMIM){
  motif_in=motif_in[order(motif_in$Proportion,decreasing=T),]
  motif_in$OMIM=NA
  for(tf in unique(motif_in$TF)){
    tf_in=gsub('\\(var.2\\)','',tf)
    tf_in=gsub('\\(var.3\\)','',tf_in)
    
    tf_in=unlist(strsplit(tf_in,'::'))
    #print(tf_in)
    OMIM_disease=OMIM$Phenotypes[which(unlist(lapply(OMIM$`Gene Symbols`,function(x) any(x%in% tf_in))))]
    if(length(OMIM_disease)>0){
      motif_in$OMIM[motif_in$TF==tf]=paste(OMIM_disease,collapse=';')}
    
  }
  
  return(motif_in)
}
NME_dNME_ken<-function(motif_in,GR_in,stat_in){
  tt1=proc.time()[[3]]
  olap=findOverlaps(motif_in,GR_in,select='all')
  if(length(GR_in$NME1)>0){GR_in$NME=(GR_in$NME1+GR_in$NME2)/2}
  subj_olap=subjectHits(olap)
  stat_in_df=data.table(qt=queryHits(olap),stat_in=elementMetadata(GR_in)[[stat_in]][subj_olap],
                        sample=GR_in$Sample[subj_olap],stringsAsFactors = T)
  sample_all=unique(GR_in$Sample)
  gc()
  stat_in_df_stat=dcast.data.table(stat_in_df,qt~sample,value.var = "stat_in",fun.aggregate = mean)
  
  for(sp in sample_all){
    elementMetadata(motif_in)[[sp]]=NA
    if(!is.null(stat_in_df_stat[[sp]])){
      elementMetadata(motif_in)[[sp]][stat_in_df_stat$qt]=stat_in_df_stat[[sp]]
    }
    
  }
  #add column of NA if motif don't have all samples
  print(proc.time()[[3]]-tt1)
  gc()
  return(motif_in)
}

read.agnostic.mouse<-function(in_dir,tissue,stage,stat_type,replicate){
  file_in=paste(in_dir,'mm10_',tissue,'_',stage,'_',replicate,'_allele_agnostic_',stat_type,'.bedGraph',sep='')
  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    informME_in$tissue=tissue
    stage=gsub('_5','.5',stage)
    stage=gsub('day','E',stage)
    stage=gsub('E0','P0',stage)
    informME_in$stage=stage
    informME_in$bioreplicate=replicate
    informME_in$Sample=paste(tissue,stage,replicate,sep='-')
    
    return(informME_in)
  }
}
#UC
read.agnostic.mouse.uc<-function(file_in,matrix=FALSE,fileter_N=1,gff_in=NA){
  
  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    #process file name
    file_in=strsplit(file_in,'\\/')[[1]]
    file_in=file_in[length(file_in)]
    file_in=gsub('_uc.bedGraph','',file_in)
    file_in=gsub('_jsd.bedGraph','',file_in)
    comp= strsplit(file_in,'-vs-')[[1]]
    comp= strsplit(file_in,'-vs-')[[1]]
    strain=unlist(lapply(strsplit(comp,'_'),function(x) x[1]))
    #if  contain BL6DBA, use ref is BL6DBA
    strain=ifelse('BL6DBA'%in%strain,'BL6DBA','mm10')
    comp=unlist(lapply(strsplit(comp,'_'),function(x) paste(x[-1],collapse = '_')))
    comp=comp[comp!='']
    comp_stage=unlist(lapply(comp,function(x) {x_split=strsplit(x,'_')[[1]]
    x_split=x_split[-length(x_split)][-1]
    x_split=paste(x_split,collapse = '_')
    return(x_split)}))
    tissue1=strsplit(comp[1],'_')[[1]][1]
    tissue2=strsplit(comp[2],'_')[[1]][1]
    #if BL6DBA, the 1st comp_stage is empty
    
    comp_stage=gsub('_5','.5',comp_stage)
    comp_stage=gsub('day','E',comp_stage)
    comp_stage=gsub('E0','P0',comp_stage)
    replicate=strsplit(comp[1],'_')[[1]][length(strsplit(comp[1],'_')[[1]])]
    replicate=gsub('merged','',replicate)
    informME_in$Sample=paste0(tissue1,'-',comp_stage[1],'-',comp_stage[2],'-',replicate)
    informME_in=informME_in[informME_in$N>=fileter_N]
    cat('Minimum N:',min(informME_in$N),'\n')
    #informME_in$Ref=strain
    if(matrix){
      informME_in$Sample=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
      informME_in_dt=as.data.table(mcols(informME_in))[,c("score","Sample")]
      colnames(informME_in_dt)=c("UC","Sample")
      informME_in_dt$UC=as.numeric(informME_in_dt$UC)
      informME_in_dt$region=paste0(seqnames(informME_in),":",start(informME_in),"-",end(informME_in))
      informME_in_dt=informME_in_dt[match(gff_in,region),"UC"]
      colnames(informME_in_dt)=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
      return(informME_in_dt)
    }
    else{return(informME_in)}
  }
}

agnostic_matrix_conversion<-function(gr_in,stat='NME'){
  gr_out=granges(unique(gr_in))
  olap=findOverlaps(gr_in,gr_out,type='equal')
  stat_in_df=elementMetadata(gr_in[queryHits(olap)])[c(stat,'Sample')]
  stat_in_df$idx=NA
  stat_in_df$idx[queryHits(olap)]=subjectHits(olap)
  stat_in_df=as.data.table(stat_in_df)
  
  stat_in_df_stat=dcast.data.table(data=stat_in_df,formula=idx~Sample,value.var = stat,fun.aggregate=mean)#remove agg.fun for new run
  gr_out=gr_out[stat_in_df_stat$idx]
  mcols(gr_out)=stat_in_df_stat[,-1]
  return(gr_out)
  
}

cor_dMML_dNME_enrich<-function(all_region_ts,prob_cutoff,FeDMR){
  all_region_ts$change=FALSE
  cutoff_cor=quantile(all_region_ts$cor,prob=prob_cutoff)
  all_region_ts$change[all_region_ts$cor>=cutoff_cor]=TRUE
  #Find regions that have FeDMR in at least 1 tissue
  
  olap=findOverlaps(all_region_ts,FeDMR)
  all_region_ts$FeDMR=FALSE
  all_region_ts$FeDMR[queryHits(olap)]=TRUE
  cont_table=matrix(c(sum(all_region_ts$change&all_region_ts$FeDMR),
                      sum(all_region_ts$change&(!all_region_ts$FeDMR)),
                      sum((!all_region_ts$change)&(all_region_ts$FeDMR)),
                      sum((!all_region_ts$change)&(!all_region_ts$FeDMR))),nrow=2)
  OR=fisher.test(cont_table)
  print(cont_table)
  return(data.table(OR=OR$estimate,pvalue=OR$p.value,cutoff_cor=cutoff_cor,
                    lowerCI=OR$conf.int[1],upperCI=OR$conf.int[2]))
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

#GO annotation
GO_run<-function(gl,back,cluster,ptcount=0){
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                                  annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")})
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  fc <- ((sigres[,"Significant"]+ptcount)/(sum(GOdata@allScores[GOdata@feasible]==1)+ptcount))/((sigres[,"Annotated"]+ptcount)/(sum(GOdata@feasible)+ptcount))
  sigres <- data.frame(sigres,FC=fc)
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
  sigres=as.data.table(sigres)
  siggene_forID=lapply(sigres$GO.ID,function(x,GOdata){
    gene=sigGenes(GOdata)[sigGenes(GOdata)%in%unlist(genesInTerm(GOdata, x))]
    gl_dt=data.table(rank=1:length(gl),gene=gl)
    mt=match(gl_dt$gene,gene)
    mt=mt[!is.na(mt)]
    highest_rank=NA
    return(list(paste(gene[mt],collapse =";"),highest_rank))
    
    
  },GOdata=GOdata)
  siggene=unlist(lapply(siggene_forID,function(x) x[[1]]))
  max_rank=unlist(lapply(siggene_forID,function(x) x[[2]]))
  if(nrow(sigres)>0){ 
    sigres$genes=siggene
    sigres$higest_ranks=max_rank
    
  }
  sigres$cluster=cluster
  return(sigres)
}
GO_run_tissue<-function(tissue,dir_in,nme_cor,mml_cor,enc_type){
  GO_out_all=list()
  csv_files=dir(paste0('../downstream/input/',dir_in),pattern="csv")
  for (ts in tissue){
    cat("Processing:",ts,'\n')
    fn=paste0(ts,'.csv')
    #read in csv file for given tissue
    csv_in_ts=fread(paste0('../downstream/input/',dir_in,'/',fn))
    csv_in_ts=csv_in_ts[order(dNME_maxJSD_rank,decreasing = F)]
    # Getting enhancer
    if(enc_type=="chromHMM_enhancer"){csv_in_ts=csv_in_ts[csv_in_ts$chromHMM_enhancer]}else
      if(enc_type=="non_chromHMM_enhancer"){csv_in_ts=csv_in_ts[!csv_in_ts$chromHMM_enhancer]}else 
        if(enc_type=="all_regions"){csv_in_ts=csv_in_ts}
    #GO annotation
    if(nrow(csv_in_ts)>1){
      #GO annotation for each cluster
      csv_out=lapply(1:10,function(clu){
        sp=paste0(ts,'-',clu)
        csv_in_ts_clu=csv_in_ts[cluster==clu]
        csv_in_ts_clu=csv_in_ts_clu[order(dNME_maxJSD_rank,decreasing=F)]
        #Add NME and mml cor
        csv_in_ts_clu$nme_cor=nme_cor[[ts]][match(csv_in_ts_clu$region,names(nme_cor[[ts]]))]
        csv_in_ts_clu$mml_cor=mml_cor[[ts]][match(csv_in_ts_clu$region,names(mml_cor[[ts]]))]
        
        if(nrow(csv_in_ts_clu)>1){
          
          #GO annotation for chromHMM 
          GO_out_cluster=GO_run(csv_in_ts_clu$gene,unique(csv_in_ts$gene),cluster=clu)
          write.csv(GO_out_cluster[FC>=1.5],
                    file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GO/',dir_in,'/',sp,'_cluster_GO.csv'),row.names = F,quote = T)
          csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
          #GO annotation for other types of regions
          GO_out_cluster_NME=GO_run(csv_in_ts_clu[nme_cor>=0.7]$gene,unique(csv_in_ts$gene),cluster=clu)
          GO_out_cluster_NME_only=GO_run(csv_in_ts_clu[nme_cor>=0.7&mml_cor<0.7]$gene,unique(csv_in_ts$gene),cluster=clu)
          GO_out_cluster_MML=GO_run(csv_in_ts_clu[mml_cor>=0.7]$gene,unique(csv_in_ts$gene),cluster=clu)
          GO_out_cluster_MML_only=GO_run(csv_in_ts_clu[mml_cor>=0.7&nme_cor<0.7]$gene,unique(csv_in_ts$gene),cluster=clu)
          GO_out_cluster_non_MML=GO_run(csv_in_ts_clu[mml_cor<0.7]$gene,unique(csv_in_ts$gene),cluster=clu)
          GO_out_cluster_NME_MML=GO_run(csv_in_ts_clu[nme_cor>=0.7&mml_cor>=0.7]$gene,unique(csv_in_ts$gene),cluster=clu)
          GO_out_cluster_non_NME_non_MML=GO_run(csv_in_ts_clu[nme_cor<0.7&mml_cor<0.7]$gene,unique(csv_in_ts$gene),cluster=clu)
          
          #return(list(csv_in_ts_clu=csv_in_ts_clu,GO_out_cluster=GO_out_cluster))
          return(list(GO_out_cluster_all=GO_out_cluster,GO_out_cluster_NME=GO_out_cluster_NME,GO_out_cluster_non_MML=GO_out_cluster_non_MML,csv_in_ts_clu=csv_in_ts_clu,
                      GO_out_cluster_NME_only=GO_out_cluster_NME_only,GO_out_cluster_MML=GO_out_cluster_MML,GO_out_cluster_MML_only=GO_out_cluster_MML_only,
                      GO_out_cluster_non_NME_non_MML=GO_out_cluster_non_NME_non_MML,GO_out_cluster_NME_MML=GO_out_cluster_NME_MML))
        }
        
      })
      GO_out_all[[ts]]=csv_out
      write.csv(fastDoCall('rbind',lapply(csv_out,function(x) x$csv_in_ts_clu))[order(dNME_maxJSD,decreasing=T)],
                file=paste0('../downstream/output/mm10_result/',enc_type,'/enhancer_gene_list/',dir_in,'/',ts,'_all.csv'))
    }
    
  }
  return(GO_out_all)
}


dcast_matrix<-function(dt_in,value_in,colnames_order=colnames_order){
  dt_in=dcast.data.table(dt_in,Term~cluster,value.var  = value_in)
  dt_in_mt=as.matrix(dt_in[,-1])
  rownames(dt_in_mt)=dt_in$Term
  dt_in_mt=dt_in_mt[order(max.col(dt_in_mt,)),]
  dt_in_mt=dt_in_mt[,colnames_order]
  return(dt_in_mt)
}
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
library(pheatmap)
plot_GO_heatmap<-function(selected_tissue,GO_anno,GO_out){
  GO_tissue=list()
  for(tissue in selected_tissue){
    GO_in=GO_out[[tissue]]
    GO_in=fastDoCall('rbind',lapply(GO_in,function(x) {
      x=x[[GO_anno]]
      x$sig_num=sum(x$FDR<=0.1)
      return(x)
    }))
    
    
    GO_in=GO_in[,.(GO.ID,Term,classicFisher,FDR,FC,cluster,sig_num)]
    GO_in_top=do.call(c,lapply(1:10,function(x) {
      
      return(GO_in[cluster==x][FC>=1.5&FDR<=0.1][order(FDR,-FC,decreasing=F)][1:5]$GO.ID)
    }))
    GO_in=GO_in[GO.ID%in%GO_in_top]
    GO_in$log10FDR=-log10(GO_in$FDR)
    GO_in_main= dcast_matrix(GO_in,"FC")
    GO_in_main=GO_in_main[order(max.col(GO_in_main),decreasing = F),]
    GO_in_FDR= dcast_matrix(GO_in,"FDR")
    #GO_in_FDR_log10= dcast_matrix(GO_in,"log10FDR")
    GO_in_FDR[GO_in_FDR<=0.1]="*"
    GO_in_FDR[GO_in_FDR>0.1]=""
    GO_in_FDR=GO_in_FDR[rownames(GO_in_main),]
    col_label=unique(GO_in[,.(cluster,sig_num)])
    col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
    c2 <- brewer.pal(10,'Set3')
    names(c2) <- 1:10
    breaksList = seq(-1, 1, by = 0.01)
    colann= data.frame(cluster=as.character(1:10))
    #pdf(paste0('../downstream/output/graphs/Figure6/GO_', tissue,'_',GO_anno,'_FC.pdf'),width=25,height=14)
    pheatmap(scalematrix(GO_in_main),cluster_rows =F,cluster_cols = F,
             show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR,border_color = NA,
             color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),
             filename=paste0('../downstream/output/graphs/Figure6/all_regions_chromHMM/',GO_anno,'/GO_', tissue,'_',GO_anno,'_FC_chromHMM.pdf'),
             cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
             fontsize=30,legend = F,labels_col=col_label$sig_num)
    #dev.off()
    GO_tissue[[tissue]]=GO_in_main
    #,breaks=breaksList,color=col(length(breaksList))
  }
  
}
matrix_conv<-function(dt_in,value.var){
  out_dc=dcast.data.table(dt_in,region~stage,value.var=value.var)
  rn=out_dc$region
  out_dc=as.matrix(out_dc[,-1])
  rownames(out_dc)=rn
  return(out_dc)
}
# Currently Not in use ----------------------------------------------------


#collapase variants
variants_collapase<-function(varsDiff){
  variants <- paste(as.character(varsDiff$REF),as.character(unlist(varsDiff$ALT)),sep="-")
  # Combine same variants (ALT/REF -> REF/ALT)
  variants[variants %in% c("A-C","C-A")] <- "A-C"
  variants[variants %in% c("A-G","G-A")] <- "A-G"
  variants[variants %in% c("A-T","T-A")] <- "A-T"
  variants[variants %in% c("C-G","G-C")] <- "C-G"
  variants[variants %in% c("C-T","T-C")] <- "C-T"
  variants[variants %in% c("G-T","T-G")] <- "G-T"
  return(variants)
}

#Merge allele agnositc data 
allele_agnostic_merge<-function(GR_in,nme_in,mml_in,pval_cutoff=0.1){
  NME_in=import.bedGraph(nme_in)
  if(all(seqlevels(NME_in)==gsub('chr','',seqlevels(NME_in)))){seqlevels(NME_in)=paste('chr',seqlevels(NME_in),sep='')}
  MML_in=import.bedGraph(mml_in)
  if(all(seqlevels(MML_in)==gsub('chr','',seqlevels(MML_in)))){seqlevels(MML_in)=paste('chr',seqlevels(MML_in),sep='')}
  olap=findOverlaps(NME_in,GR_in)
  agnostic_diff_df=data.frame(NME=NME_in$score[queryHits(olap)],dNME=GR_in$dNME[subjectHits(olap)],
                              dNME_pval=GR_in$dNME_pval[subjectHits(olap)],dMML=GR_in$dMML[subjectHits(olap)],
                              dMML_pval=GR_in$dMML_pval[subjectHits(olap)],MML=MML_in$score[queryHits(olap)])
  agnostic_diff_df$dNME_ASM=agnostic_diff_df$dNME_pval<=pval_cutoff
  agnostic_diff_df$dMML_ASM=agnostic_diff_df$dMML_pval<=pval_cutoff
  agnostic_diff_df$dNME_ASM_non_dMML=agnostic_diff_df$dMML_pval>pval_cutoff & agnostic_diff_df$dNME_pval<=pval_cutoff
  agnostic_diff_df$sample=unique(GR_in$Sample)
  return(agnostic_diff_df)
  
}

#Processing RNA-seq data
RNA_seq_process<-function(dir_in,fds,condition_rep=3){
  files=paste(dir_in,fds,"/t_data.ctab",sep="")
  tmp=read_tsv(files[1])
  tx2gene <- tmp[, c("t_name", "gene_name")]
  txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)
  txi[1:3]=lapply(txi[1:3],function(x) {colnames(x)=fds 
  return(x)})
  sampleTable <- data.frame(condition = factor(rep(c("genome1", "genome2"), each =condition_rep)))
  rownames(sampleTable) <- colnames(txi$counts)
  ####Perfrom differential RNA analysis
  dds_RNA<- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
  if(ncol(as.data.frame(assay(dds_RNA)))==2){
    res_RNA<-as.data.frame(cpm(assay(dds_RNA)))
    res_RNA$gene_name=rownames(res_RNA)
    res_RNA=res_RNA[(res_RNA[,1]!=0&res_RNA[,2]!=0),]
    res_RNA=res_RNA[(res_RNA[,1]>10&res_RNA[,2]>10),]
    res_RNA$log2FoldChange=log2((res_RNA[,2])/(res_RNA[,1]))
  }else{
    dds_RNA<-DESeq(dds_RNA)
    res_RNA<-results(dds_RNA,name= "condition_genome2_vs_genome1")
  }
  return(res_RNA)
}

#Tissue to germlayer
tissue_to_germlayer<-function(GR_input){
  GR_input$germlayer=NA
  tissue_ectoderm=c("foreskin_keratinocyte_paired",
                    "foreskin_melanocyte_paired",
                    "ectoderm_paired",
                    "brain_cerebellum_tissue_paired",
                    "brain_germinal_matrix_tissue_paired",
                    "Brain_substantia_nigra_paired",
                    "Brain_Hippocampus_middle_paired" )
  tissue_mesoderm=c("mesoderm_23_paired","Adipose_single",
                    "Left_Ventricle_single","Psoas_Muscle_single" ,
                    "Right_Ventricle_single","Right_Atrium_single","Spleen_single",
                    "Adrenal_Gland_single","Aorta_single","Ovary_single")
  tissue_endoderm=c("Small_Intestine_single","Lung_single","endoerm_27_paired",
                    "Bladder_single" ,"Gastric_single", "Sigmoid_Colon_single",
                    "Thymus_single","Esophagus_single", "Pancreas_single" ,"Liver_single")
  tissue_ESC=c("rep1","rep2","merged","42_embryonic_stem_cell_single" , "stem_27_undifferentiated_paired")
  GR_input$germlayer[GR_input$tissue %in% tissue_ectoderm]='ectoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_mesoderm]='mesoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_endoderm]='endoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_ESC]='ESC'
  return(GR_input)
}








OR_VMR<-function(NME_dat,vmr,percent,NME_quant='quant_score'){
  NME_dat$VMR=FALSE
  olap=findOverlaps(NME_dat,vmr)
  NME_dat$VMR[unique(queryHits(olap))]=TRUE
  NME_VMR=sum((elementMetadata(NME_dat)[[NME_quant]]%in%percent)&NME_dat$VMR)
  nonNME_VMR=sum((!elementMetadata(NME_dat)[[NME_quant]]%in%percent)&NME_dat$VMR)
  nonNME_nonVMR=sum((!elementMetadata(NME_dat)[[NME_quant]]%in%percent)&!NME_dat$VMR)
  NME_nonVMR=sum((elementMetadata(NME_dat)[[NME_quant]]%in%percent)&!NME_dat$VMR)
  #print(matrix(c(NME_VMR,nonNME_VMR,NME_nonVMR,nonNME_nonVMR),nrow = 2))
  fisher.test(matrix(c(NME_VMR,nonNME_VMR,NME_nonVMR,nonNME_nonVMR),nrow = 2))
}

get_traits_GWAS<-function(variant_in,trait_gr_in,pval_cutoff=0.1,count_cutoff=3,stat='dNME_pval',CMH=FALSE,maxgap=500,ncores=15){
  traits_ls=mclapply(trait_gr_in,get_traits_GWAS_all_trait,
                   variant_in=variant_in,pval_cutoff=pval_cutoff,count_cutoff=count_cutoff,stat=stat,CMH=CMH,maxgap=maxgap,
                   mc.cores=ncores)
  return(traits_ls)
  #list(do.call(rbind,lapply(traits_ls,function(x)x[[1]])),do.call(rbind,lapply(traits_ls,function(x)x[[2]])))
}
get_traits_GWAS_all_trait<-function(trait_gr,variant_in,pval_cutoff,count_cutoff,stat,CMH,maxgap=500){
  #trait_gr=trait_gr[trait_gr$trait%in%trait]
  trait=unique(trait_gr$`DISEASE/TRAIT`)
  OR_output=data.frame()
  CMH_df=data.frame()
  dNME_sig=variant_in[elementMetadata(variant_in)[,stat]<=pval_cutoff]
  dNME_non_sig=variant_in[elementMetadata(variant_in)[,stat]>pval_cutoff]
  dNME_traits_gr=findOverlaps(dNME_sig,trait_gr,maxgap = maxgap)
  # dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
  # non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
  # dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
  # non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
  dNME_trait=length(unique(queryHits(dNME_traits_gr)))
  dNME_non_trait=length(dNME_sig)-dNME_trait
  non_dNME_trait=length(subsetByOverlaps(dNME_non_sig,trait_gr,maxgap = maxgap))
  non_dNME_non_trait=length(dNME_non_sig)-non_dNME_trait
  trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
  journal=paste(unique(trait_gr$JOURNAL),collapse = ',')
  trait_gr_out=trait_gr[subjectHits(dNME_traits_gr)]
  mcols(trait_gr_out)=mcols(trait_gr_out)[c('MAPPED_GENE','STRONGEST SNP-RISK ALLELE','DISEASE/TRAIT')]
  names(mcols(trait_gr_out))=c('genes','risk allele','traits')
  rm(dNME_sig)
  rm(dNME_non_sig)
  rm(trait_gr)
  gc()
  if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
    
    CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
                      count=trait_count,subject='All')
    cat('processing trait',trait,'\n')
    
    if(length(CMH_df)>0){
      if (CMH){return(CMH_df)}
      else{
        
        OR=CMH_test(CMH_df)
        OR_output=cbind(t(CMH_df$count),data.frame(trait=trait,OR=OR$estimate,
                                                   p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
        colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
        OR_output$journal=journal
        
        return(list(OR_output,trait_gr_out))
      }
    }
  }
}
#Modify this to adapt trait analysis
get_traits_GWAS_all_trait_single<-function(variant_in,trait_gr,pval_cutoff,count_cutoff,stat,CMH,maxgap){
 
  OR_output=data.frame()
  CMH_df=data.frame()
  dNME_sig=variant_in[elementMetadata(variant_in)[,stat]<=pval_cutoff]
  dNME_non_sig=variant_in[elementMetadata(variant_in)[,stat]>pval_cutoff]
  
  # dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
  # non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
  # dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
  # non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
  dNME_trait=length(subsetByOverlaps(dNME_sig,trait_gr,maxgap = maxgap))
  dNME_non_trait=length(dNME_sig)-dNME_trait
  non_dNME_trait=length(subsetByOverlaps(dNME_non_sig,trait_gr,maxgap = maxgap))
  non_dNME_non_trait=length(dNME_non_sig)-non_dNME_trait
  trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
  rm(dNME_sig)
  rm(dNME_non_sig)
  rm(trait_gr)
  gc()
  if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
    
    CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
                      count=trait_count,subject='All')

    
    if(length(CMH_df)>0){
      if (CMH){return(CMH_df)}
      else{
        
        OR=CMH_test(CMH_df)
        OR_output=cbind(t(CMH_df$count),data.frame(OR=OR$estimate,
                                                   p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
        colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
        OR_output$journal=paste(unique(trait_gr$JOURNAL),collapse = ',')
        return(OR_output)
      }
    }
  }
}



get_traits_GWAS_trait<-function(trait,variant_in,pval_cutoff,count_cutoff,stat,CMH){
  traits_sp_ls=lapply(unique(variant_in$germlayer),get_traits_GWAS_sp_trait,trait=trait,
                      variant_in=variant_in,pval_cutoff=pval_cutoff,count_cutoff=count_cutoff,stat=stat,CMH=CMH)
  traits_sp=do.call(rbind,traits_sp_ls)
  if(CMH&length(traits_sp)>0){
    OR=CMH_test(traits_sp)
    OR_output=data.frame(trait=trait,OR=OR$estimate,p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2])
    return(OR_output)
  }else if(!CMH){
    return(traits_sp)
  }


}




get_traits_GWAS_sp_trait<-function(sp,trait,variant_in,pval_cutoff,count_cutoff,stat,CMH){
 
      OR_output=data.frame()
      CMH_df=data.frame()
      variant_in_sp=variant_in[variant_in$germlayer==sp&!is.na(variant_in$trait),]
 
      dNME_sig=variant_in_sp[,stat]<=pval_cutoff
      dNME_non_sig=variant_in_sp[,stat]>pval_cutoff
    
      dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
      non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
      dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
      non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
      trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
      if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
       
        CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
                                       count=trait_count,subject=sp)
        cat('processing',sp,'with trait',trait,'\n')

      if(length(CMH_df)>0){
        if (CMH){return(CMH_df)}
        else{
          OR=CMH_test(CMH_df)
          OR_output=cbind(t(CMH_df$count),data.frame(subject=sp,trait=trait,OR=OR$estimate,
                                                     p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
          colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
          return(OR_output)
        }
      }
     }
}

# Currently not in use ----------------------------------------------------


ASM_het_enrichment<-function(gr_in){
  ASM_het=sum(gr_in$ASM=='Yes' & gr_in$HetCpG)
  ASM_non_het=sum(gr_in$ASM=='Yes' & !gr_in$HetCpG)
  non_ASM_non_het=sum(gr_in$ASM=='No' & !gr_in$HetCpG)
  non_ASM_het=sum(gr_in$ASM=='No' & gr_in$HetCpG)
  cont_table=matrix(c(ASM_het,ASM_non_het,non_ASM_het,non_ASM_non_het),nrow=2)
  print(cont_table)
  fisher.test(cont_table)
}
testEnrichmentFeature <- function(dataGR,featureGR){
  
  # Find ranges overlapping with feature
  olaps <- findOverlaps(dataGR,featureGR,type="any",select="all")
  indFeature <- queryHits(olaps)
  featureData <- dataGR[indFeature]
  complementaryData <- dataGR[-indFeature]
  
  # Enrichment of dMML-HASM in feature
  featureDmml <- featureData[featureData$Statistic=="dMML"]
  complementaryDmml <- complementaryData[complementaryData$Statistic=="dMML"]
  contTableDmml <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableDmml) <- c("Feature","Complementary")
  contTableDmml[1,]$ASM <- sum(featureDmml$ASM=="Yes")
  contTableDmml[1,]$nonASM <- sum(featureDmml$ASM=="No")
  contTableDmml[2,]$ASM <- sum(complementaryDmml$ASM=="Yes")
  contTableDmml[2,]$nonASM <- sum(complementaryDmml$ASM=="No")
  dmmlFisher <- fisher.test(contTableDmml)
  
  # Enrichment of dNME-HASM in feature
  featureDnme <- featureData[featureData$Statistic=="dNME"]
  complementaryDnme <- complementaryData[complementaryData$Statistic=="dNME"]
  contTableDnme <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableDnme) <- c("Feature","Complementary")
  contTableDnme[1,]$ASM <- sum(featureDnme$ASM=="Yes")
  contTableDnme[1,]$nonASM <- sum(featureDnme$ASM=="No")
  contTableDnme[2,]$ASM <- sum(complementaryDnme$ASM=="Yes")
  contTableDnme[2,]$nonASM <- sum(complementaryDnme$ASM=="No")
  dnmeFisher <- fisher.test(contTableDnme)
  
  # Enrichment of UC-HASM in feature
  featureUc <- featureData[featureData$Statistic=="UC"]
  complementaryUc <- complementaryData[complementaryData$Statistic=="UC"]
  contTableUc <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableUc) <- c("Feature","Complementary")
  contTableUc[1,]$ASM <- sum(featureUc$ASM=="Yes")
  contTableUc[1,]$nonASM <- sum(featureUc$ASM=="No")
  contTableUc[2,]$ASM <- sum(complementaryUc$ASM=="Yes")
  contTableUc[2,]$nonASM <- sum(complementaryUc$ASM=="No")
  ucFisher <- fisher.test(contTableUc)
  
  # Return list of Fisher's test
  return(list(dmmlFisher,dnmeFisher,ucFisher))
  
}

getGeneralFeats_mm9 <- function(CpGdir,enhancerDir='',chrsOfInterest=paste("chr",1:21,sep="")){
  
  # Features included
  featureNickNames <- c("genome-wide","CpG island","CpG shore","CpG shelf","CpG open sea",
                        "gene body","exon","intron","intergenic")
  
  # Define list of feature GRs
  outGR <- GRangesList()
  GRtemp <- unlist(tileGenome(seqinfo(Mus),ntile=1))
  
  outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
  #Redefining CpG islands using hidden Markov models 
  CpG_all <- readRDS(paste(CpGdir,"CpG_hg19.rds",sep=""))
  CpG_all<-setGenomeLengths(CpG_all)
  #Could also use UCSC genome browser CpG file
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
  genes <- genes(txdb)
  outGR[["gene body"]] <- setGenomeLengths(genes)
  exons <- exons(txdb)
  outGR[["exon"]] <- setGenomeLengths(exons[,c()])
  introns <- intronicParts(txdb)
  outGR[["intron"]] <- setGenomeLengths(introns[,c()])
  intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
  outGR[["intergenic"]] <- setGenomeLengths(intergenic)
  #Use annotation hub for TSS
  proms <- promoters(genes,upstream=2000,downstream=2000)
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
#Get CpG sites from hg19
getCpgSitesH19 <- function(chrsOfInterest=paste("chr",1:22,sep="")){
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x) GRanges(x,IRanges(start(matchPattern("CG", Hsapiens[[x]])),width=2)))
  
  # Set genome and seqlengths
  cgs <- setGenomeLengths(do.call('c',cgs))
  # Return CpG site GR
  return(cgs)
}

#Get CpG density for each chromosome
getCpgdensH19 <- function(chrsOfInterest=paste("chr",1:22,sep="")){
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cgs_dist=lapply(cgs,function(x) x[2:length(x)]-x[1:length(x)-1])
  cgs_df=data.frame(CG_number=unlist(lapply(cgs,length)),
                    CG_dist=unlist(lapply(cgs_dist,mean)),
                    total_length=unlist(lapply(chrsOfInterest,function(x) length(Hsapiens[[x]]))))
  cgs_df$CG_density=cgs_df$CG_number/cgs_df$total_length
  return(list(cgs_df,cgs_dist))
}



# Function to plot a GR using Gviz
plotGR <- function(CpGdir,enhancerDir,GR,startHight,highSize=500,reverseStrand=FALSE,chr="chr11",lim=c(2010000,2022500)){
  
  # Get genome
  gen <- "hg19"

  # GRanges to intersect with and keep the relevant data
  windowGR <- GRanges(seqnames=chr,ranges=IRanges(start=lim[1],end=lim[2]),strand="*")
  
  # Subset GR
  GR <- subsetByOverlaps(GR,windowGR)

  # Create data tracks
  dmmlTrack <- DataTrack(GR[GR$Statistic=="dMML",c("Value")],name="dMML")
  dnmeTrack <- DataTrack(GR[GR$Statistic=="dNME",c("Value")],name="dNME")
  ucTrack <- DataTrack(GR[GR$Statistic=="UC",c("Value")],name="UC")
  
  # Gene track 1
  bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="grch37.ensembl.org",
                path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
  biomTrack <- BiomartGeneRegionTrack(genome=gen,chromosome=chr,start=lim[1],end=lim[2],name="ENSEMBL",biomart=bm)

  # Add CpG island annotation track
  genomicFeatures <- getGeneralFeats(CpGdir,enhancerDir)
  cpgIslands <- genomicFeatures[["CpG island"]]
  cpgIslands <- subsetByOverlaps(cpgIslands,windowGR)
  islandTrack <- AnnotationTrack(cpgIslands,name="CpG islands")
  
  # Chromosome information tracks
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome=gen,chromosome=chr)
  
  # Highlight
  ht <- HighlightTrack(trackList=list(dmmlTrack,dnmeTrack,ucTrack),start=startHight,width=highSize,chromosome=chr)
  
  # Return plot
  plotTracks(list(itrack,gtrack,biomTrack,islandTrack,ht),from=lim[1],to=lim[2],
             transcriptAnnotation="symbol",type=c("gradient"),stacking="squish",reverseStrand=reverseStrand,
             collapseTranscripts = "meta")
  
}
# Function to assign genome
assignGenomeFile <- function(row) {
  
  genome = NA
  if((row[1]=="REF")&(row[2] %in% c("0|1","0/1"))){
    genome = 1
  } else if((row[1]=="ALT")&(row[2] %in% c("0|1","0/1"))){
    genome = 2
  } else if((row[1]=="REF")&(row[2] %in% c("1|0","1/0"))) {
    genome = 2
  } else if((row[1]=="ALT")&(row[2] %in% c("1|0","1/0"))){
    genome = 1
  } else {
    # nothing
  }
  
  return(genome)
}

# Function that returns bool if variant results in a CpG site
hetCpgSite <- function(row) {
  
  # Check if we have C or G, otherwise return false
  if(row[3] %in% c("C","G")){
    # continue
  } else {
    return(FALSE)
  }
  
  # Initialize output
  hetCpg <- FALSE
  context <- toString(getSeq(Hsapiens,row[1],start=as.numeric(row[2])-1,end=as.numeric(row[2])+1,strand="+"))
  if((row[3]=="C")&(substr(context,3,3)=="G")){
    hetCpg <- TRUE
  } else if((row[3]=="G")&(substr(context,1,1)=="C")){
    hetCpg <- TRUE
  } else {
    # nothing
  }
  
  # Return binary vector
  return(hetCpg)
  
}

#Find min width, no longer needed for new output
gff_min<-function(op_more,olap_sub,qt,gffwid){
  gff_idx=olap_sub[which(qt==op_more)]
  return(gff_idx[which.min(gffwid[gff_idx])])
}

###End of calculation



#generate and export bed file

ASM_bed_gen_sp<-function(gr_ASM,gr_allele,sp,outdir){
  gr_ASM=gr_ASM[gr_ASM$CpGdiff!=0]
  gr_ASM=gr_ASM[order(gr_ASM$diff,decreasing = TRUE)]
  #do it for each sample
  #sp='Gastric - STL001'
  sp_ASM=gr_ASM[gr_ASM$Sample==sp]#stats for this subject
  sp_allele=gr_allele[gr_allele$Sample==sp & gr_allele$HetCpG]
  sp_bed=granges(sp_ASM)
  sp_bed$diff=sp_ASM$diff
  olap1=findOverlaps(sp_allele[sp_allele$Genome==1],sp_bed,type='equal')
  olap2=findOverlaps(sp_allele[sp_allele$Genome==2],sp_bed,type='equal')
  sp_bed$A1[subjectHits(olap1)]=sp_allele[sp_allele$Genome==1]$Value[queryHits(olap1)]
  sp_bed$A2[subjectHits(olap2)]=sp_allele[sp_allele$Genome==2]$Value[queryHits(olap2)]
  #Export bed file
  export_bed(sp_bed,'diff',paste(outdir,sp,'_hetASM_diff.bedGraph',sep=''))
  export_bed(sp_bed,'A1',paste(outdir,sp,'_hetASM_A1.bedGraph',sep=''))
  export_bed(sp_bed,'A2',paste(outdir,sp,'_hetASM_A2.bedGraph',sep=''))
}
export_bed<-function(gr_in,dat,out_name){
  bed_out=granges(gr_in)
  bed_out$score=elementMetadata(gr_in)[,dat]
  export(bed_out,out_name,format='bedGraph')
}




#Generate motif from variant file

#From vcf file, extract het CpG information
extractmotif<-function(vcfDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  cat('Processing subject:', sub,'\n')
  tt1=proc.time()[[3]]
  #genomeGr <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  #genomeGr <- setGenomeLengths(genomeGr)
  vcf <- readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19")
  gt <- as.vector(geno(vcf)$GT)
  vcf <- rowRanges(vcf)
  vcf$GT <- gt
  vcf$snpId <- paste(sub,seq(1:length(vcf)),sep="-")
  # Keep only relevant variables
  vcf <- vcf[,c("REF","ALT","GT","snpId")]
  vcf$REF <- as.character(vcf$REF)
  vcf$ALT <- as.character(unlist(vcf$ALT))
  names(vcf)=NULL
  # Delete labels
  vcf=motif_df(vcf)
  vcf$sub=sub
  cat('Done processing',sub,'using:', proc.time()[[3]]-tt1,'\n')
  return(vcf)
  
}
#ID motif from vcf file
motif_df<-function(var_gr){
  #Get 11 nucleotide location, 3 nucleotide 4-6, 5 nulceotide 3-7, 7 nucleotide 2-8
  nucleo_11_alt<-nucleo_11_X<-nucleo_11_ref<-as.matrix(getSeq(Hsapiens,resize(var_gr,11,fix='center')))
  nucleo_11_X[,6]='N'
  nucleo_11_alt[,6]=var_gr$ALT
  #Running too long. using matrix operation
  #cat matrix into 1 string
  nucleo_11_ref=apply(nucleo_11_ref,1,function(x) paste(x,collapse=''))
  nucleo_11_X=apply(nucleo_11_X,1,function(x) paste(x,collapse=''))
  nucleo_11_alt=apply(nucleo_11_alt,1,function(x) paste(x,collapse=''))
  #Get 3,5,7,9 nucleotide location
  nucleo_9_ref=substr(nucleo_11_ref,start=2,stop=10)
  nucleo_7_ref=substr(nucleo_11_ref,start=3,stop=9)
  nucleo_5_ref=substr(nucleo_11_ref,start=4,stop=8)
  nucleo_3_ref=substr(nucleo_11_ref,start=5,stop=7)
  
  nucleo_9_alt=substr(nucleo_11_alt,start=2,stop=10)
  nucleo_7_alt=substr(nucleo_11_alt,start=3,stop=9)
  nucleo_5_alt=substr(nucleo_11_alt,start=4,stop=8)
  nucleo_3_alt=substr(nucleo_11_alt,start=5,stop=7)
  
  nucleo_9_X=substr(nucleo_11_X,start=2,stop=10)
  nucleo_7_X=substr(nucleo_11_X,start=3,stop=9)
  nucleo_5_X=substr(nucleo_11_X,start=4,stop=8)
  nucleo_3_X=substr(nucleo_11_X,start=5,stop=7)
  
  #Assign value to var_gr
  var_df=data.frame(nucleo_11_ref=nucleo_11_ref,
  nucleo_9_ref=nucleo_9_ref,
  nucleo_7_ref=nucleo_7_ref,
  nucleo_5_ref=nucleo_5_ref,
  nucleo_3_ref= nucleo_3_ref,
  
  nucleo_11_alt=nucleo_11_alt,
  nucleo_9_alt=nucleo_9_alt,
  nucleo_7_alt=nucleo_7_alt,
  nucleo_5_alt=nucleo_5_alt,
  nucleo_3_alt=nucleo_3_alt,
  
  nucleo_11_X=nucleo_11_X,
  nucleo_9_X=nucleo_9_X,
  nucleo_7_X=nucleo_7_X,
  nucleo_5_X=nucleo_5_X,
  nucleo_3_X=nucleo_3_X,stringsAsFactors = FALSE)
  var_gr=makeGRangesFromDataFrame(cbind(var_gr,var_df),keep.extra.columns = TRUE)
  return(var_gr)
}

#For each sample, add dMML and dNME information
extract_diff_values<-function(sp,diff,variant){
  #Get subject information for sp
  subj= strsplit(sp,' - ')[[1]][2]
  variant=variant[[subj]]
  #For this sample, extract the loci
  #dMML
  dMML=diff[diff$Statistic=='dMML']
  olap_dMML=findOverlaps(variant,dMML,type='within')
  outGR=variant[queryHits(olap_dMML)]
  outGR$dMML=dMML$Value[subjectHits(olap_dMML)]
  outGR$dMML_pval=dMML$pvalue[subjectHits(olap_dMML)]
  #dNME
  dNME=diff[diff$Statistic=='dNME']
  olap_dNME=findOverlaps(outGR,dNME,type='within')
  outGR$dNME=NA
  outGR$dNME_pval=NA
  outGR[queryHits(olap_dNME)]$dNME=dNME$Value[subjectHits(olap_dNME)]
  outGR[queryHits(olap_dNME)]$dNME_pval=dNME$pvalue[subjectHits(olap_dNME)]
  outGR$Sample=sp
  return(outGR)
}


#Reshape each motif and keep extra column
reshape_sample_variant<-function(variant){
  variant=as.data.frame(variant)
  melt_id_var=colnames(variant)
  melt_id_var=melt_id_var[-which(melt_id_var=="REF" | melt_id_var=="ALT")]
 
  # Melt: origin = ref/alt, variant= ATCG
  variant <- melt(data=variant,id.vars=melt_id_var,
                measure.vars=c("REF","ALT"),
                variable.name="Origin",value.name="Variant")
  #Assign genome
  variant$Genome <- apply(variant[,c("Origin","GT")],1,assignGenomeFile)
  #Reshape ref columns
  ref_id=c('nucleo_11_ref','nucleo_9_ref','nucleo_7_ref','nucleo_5_ref','nucleo_3_ref')
  alt_id=c('nucleo_11_alt','nucleo_9_alt','nucleo_7_alt','nucleo_5_alt','nucleo_3_alt')
  new_id=c('nucleo_11','nucleo_9','nucleo_7','nucleo_5','nucleo_3')
  #reassign ID
  variant[,new_id]=NA
  variant[variant$Genome==1,new_id]=variant[variant$Genome==1,ref_id]
  variant[variant$Genome==2,new_id]=variant[variant$Genome==2,alt_id]
  variant=variant[,-which(colnames(variant)%in%c(ref_id,alt_id))]
  #Make Granges
  variant=makeGRangesFromDataFrame(variant,keep.extra.columns = TRUE)
  return(variant)
}
assignGenomeFile <- function(row) {
  
  genome = NA
  if((row[1]=="REF")&(row[2] %in% c("0|1","0/1"))){
    genome = 1
  } else if((row[1]=="ALT")&(row[2] %in% c("0|1","0/1"))){
    genome = 2
  } else if((row[1]=="REF")&(row[2] %in% c("1|0","1/0"))) {
    genome = 2
  } else if((row[1]=="ALT")&(row[2] %in% c("1|0","1/0"))){
    genome = 1
  } else {
    # nothing
  }
  
  return(genome)
}
#Extract allele values for each allele
extract_allele_value<-function(outGR,cpelAllele){
  # Cross resulting GR with MML of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with MML of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML2[outGR$Genome=="1"] <- NA
  
  # Consolidate MML1 and MML2 columns into single column
  outGR$MML <- rowSums(cbind(outGR$MML1,outGR$MML2), na.rm=T)
  outGR$MML1 <- outGR$MML2 <- NULL
  
  # Cross resulting GR with NME of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with NME of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME2[outGR$Genome=="1"] <- NA
  
  # Consolidate NME1 and NME2 columns into single column
  outGR$NME <- rowSums(cbind(outGR$NME1,outGR$NME2), na.rm=T)
  outGR$NME1 <- outGR$NME2 <- NULL

  # Return
  return(outGR)
}
#Enrichment test
motif_enrichment<-function(GR_allele,pval_cutoff=0.1,p_stat,motif,motif_type){#pstat either dNME_pval or dMML_pval
  GR_allele=as.data.frame(GR_allele)
  GR_allele$ASM=NA
  GR_allele$ASM[GR_allele[,p_stat]<=pval_cutoff] =  TRUE
  GR_allele$ASM[GR_allele[,p_stat]>pval_cutoff] = FALSE
  ASM_motif=sum(GR_allele$ASM & GR_allele[,motif_type]==motif)
  ASM_not_motif=sum(GR_allele$ASM & !GR_allele[,motif_type]==motif)
  notASM_motif=sum(!GR_allele$ASM & GR_allele[,motif_type]==motif)
  notASM_notmotif=sum(!GR_allele$ASM & !GR_allele[,motif_type]==motif)
  cont_table=matrix(c(ASM_motif,ASM_not_motif,notASM_motif,notASM_notmotif),nrow=2,byrow=TRUE)
  return(fisher.test(cont_table))
}

# #for given sample
# variant_meta_sp<-function(variant_subj,GR_st){
#   st=unique(GR_st$Statistic)
#   gr_out_st = GRanges()
#   for (statistics in st){
#     gr_out_st=c(gr_out_st,variant_meta_sp_st(variant_subj,GR_st[GR_st$Statistic==statistics]))
#     
#   }
#   return(gr_out_st)
# }

#Enrichment of variants
variant_enrich<-function(variant_in,cutoff=0.1){
  variant_in$ASM=FALSE
  variant_in$ASM[variant_in$dNME_pval<=cutoff]=TRUE
  cont_table=matrix(c(sum(variant_in$variant & variant_in$ASM),
                      sum(variant_in$variant & !variant_in$ASM),
                      sum(!variant_in$variant & variant_in$ASM),
                      sum(!variant_in$variant & !variant_in$ASM)),
                    nrow=2,byrow=TRUE)
  colnames(cont_table)=c('ASM','Not ASM')
  rownames(cont_table)=c('Vairant','Not variant')
  print(cont_table)
  return(fisher.test(cont_table))
}
#Subset for SNP-containing ranges
SNP_conmtaining_hap<-function(gr_in,variant_in){
  SNP_sub=GRanges()
  SNP_not=GRanges()
  for (subj in unique(gr_in$Subject)){
    olap=findOverlaps(gr_in[gr_in$Subject==subj],variant_in[[subj]])
    SNP_sub=c(SNP_sub,subsetByOverlaps(gr_in[gr_in$Subject==subj], variant_in[[subj]]))
    SNP_not=c(SNP_not,gr_in[gr_in$Subject==subj][-queryHits(olap)])
  }
return(list(SNP_containing=SNP_sub,Non_SNP_containing=SNP_not))
}
#ASM_het_CpG_enrichment
ASM_het_enrich<-function(gr_in,title){
OR_df=data.frame(sp=NULL,OR=NULL,lower_CI=NULL,upper_CI=NULL)
for (sp in unique(gr_in$Sample)){
  OR=ASM_het_enrichment(gr_in[gr_in$Sample==sp])
  OR_df=rbind(OR_df,data.frame(sp=sp,subjects=strsplit(sp,' - ')[[1]][2],
            OR=OR$estimate,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
}
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
ggplot(OR_df,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 13)+
  ggtitle(title)+xlab('Sample name')+ylab('Odds Ratio')+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+  theme_bar

}

#Plot global distribution of density vs CG type
plot_density<-function(CG_df,ylab,title,ylim,xlab='CG type'){
ggplot(CG_df,aes(x=CG_type,y=density,fill=CG_type))+#scale_fill_manual(values = c("blue","red"))+
  geom_boxplot(outlier.shape = NA)+xlab(xlab)+ylab(ylab)+ylim(ylim)+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
  ggtitle(title)+theme(legend.title = element_blank())
}
#CpG feature enrichment plot
genome_feature_plot<-function(gr,feature,Stats,title,ylim=c(0,2),pval_cutoff=0.1){
  groups=unique(gr$Group)
  gr$ASM=NA
  gr$ASM[which(elementMetadata(gr)[,paste(Stats,"pval",sep="_")]<=pval_cutoff)]="Yes"
  gr$ASM[which(elementMetadata(gr)[,paste(Stats,"pval",sep="_")]>pval_cutoff)]="No"
  OR_df=data.frame(group=groups,OR=0,lower_CI=0,upper_CI=0)
  for(gp in groups){
    OR_out=testEnrichmentFeature_stat(gr[gr$Group==gp],feature)[[2]]
    OR_df$OR[OR_df$group==gp]=OR_out$estimate
    OR_df$lower_CI[OR_df$group==gp]=OR_out$conf.int[1]
    OR_df$upper_CI[OR_df$group==gp]=OR_out$conf.int[2]
  }  
  print(OR_df)
  theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
  ggplot(OR_df,aes(x=group,y=OR,fill=group)) + ylim(ylim)+
    geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
    geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                  position=position_dodge(.9))+theme_bar
    
}

#calculating distance

multiple_reduce<-function(qt_idx,hits){

  hits_red=hits[queryHits(hits)==qt_idx]

  idx_min=which.min(elementMetadata(hits_red)$distance)
  return(c(subjectHits(hits_red)[idx_min],elementMetadata(hits_red)$distance[idx_min]))
  
}
compare_dist<-function(x){
  if(any(is.na(x))){
    return(x[!is.na(x)])
  }else{
    return(x[which.min(abs(x))])
  }
  
}
compare_dist_gene<-function(x){
 
  if(any(is.na(x[1:2]))){
    return(x[which(!is.na(x[1:2]))+2])
  }else{
    return(x[which.min(abs(as.numeric(x[1:2])))+2])
  }
  
  
}
#Distance to given granges
gr_distance<-function(gr_in,gr_feature,xlab,main,ylim){
  gr_dist=GRanges()
  for (subj in unique(gr_in$Sample)){
    gr_dist=c(gr_dist,dist_calc(gr_in[gr_in$Sample == subj],gr_feature))
              }
  gr_dist$dist_round=gr_dist$dist_round/1000
  #gr_all_close=gr_dist[abs(gr_dist$dist_round)<=50]
  gr_all_close=gr_dist[abs(gr_dist$dist_round)<=15]
  gr_count=table(gr_all_close$dist_round)
  gr_plot_df=data.frame(dist=as.numeric(names(gr_count)),percent_ASM=gr_count/length(gr_all_close))
  plot(gr_plot_df$dist,gr_plot_df$percent_ASM.Freq,pch=1,cex=0.8,ylab='Proportion of ASM',xlab=xlab,main=main,ylim=ylim,xlim=c(-15,15))
  lines(gr_plot_df$dist,gr_plot_df$percent_ASM.Freq,lwd=1.5)
  abline(h=mean(gr_plot_df$percent_ASM.Freq),lty=2,lwd=4)
  return(gr_all_close)
}

readEnhancer <- function(enhancerDir){
  
  #Function to convert colnames to granges
  load(paste(enhancerDir,"enhancers_intersect.RData",sep=""))
  enhancer_gr_all <- do.call('c',lapply(rownames(max_states),rownames2Granges))
  
  # Return
  return(enhancer_gr_all)
}
#Find the overlap event in ASM
olap_ASM<-function(GR_in){
  dMML=GR_in[GR_in$Statistic=='dMML']
  dNME=GR_in[GR_in$Statistic=='dNME']
  UC=GR_in[GR_in$Statistic=='UC']
  dMML_dNME=length(subsetByOverlaps(dMML,dNME))
  dMML_UC=length(subsetByOverlaps(dMML,UC))
  dNME_UC=length(subsetByOverlaps(dMML,dNME))
  olap3=length(subsetByOverlaps(subsetByOverlaps(UC,dNME),dMML))
  return(data.frame(dMML_dNME=dMML_dNME-olap3,dMML_UC=dMML_UC-olap3,
                    dNME_UC=dNME_UC-olap3,olap3=olap3,dMML=length(dMML),
                    dNME=length(dNME),UC=length(UC),
                    dMML_nonolap=length(dMML)-dMML_dNME-dMML_UC-olap3,
                    dNME_nonolap=length(dNME)-dMML_dNME-dNME_UC-olap3,
                    UC_nonolap=length(UC)-dNME_UC-dMML_UC-olap3,
                    sample=unique(GR_in$Sample),subject=unique(GR_in$Subject)))
}

#Calculate enrichment of each variants in ASM
variants_OR<-function(variant_gr,variant,statistics,cutoff=0.1){
  variant_gr$pvalue=elementMetadata(variant_gr)[,paste(statistics,'pval',sep="_")]
  invariant=variant_gr[variant_gr$variants==variant]
  nonvariant=variant_gr[variant_gr$variants!=variant]
  variant_ASM=sum(invariant$pvalue<=cutoff)
  variant_nonASM=sum(invariant$pvalue>cutoff)
  nonvariant_ASM=sum(nonvariant$pvalue<=cutoff)
  nonvariant_nonASM=sum(nonvariant$pvalue>cutoff)
  cont_table=matrix(c(variant_ASM,variant_nonASM,nonvariant_ASM,nonvariant_nonASM),nrow=2)
  return(fisher.test(cont_table))
}
#calculate odds ratio of trinucleotide
tri_nucleo_OR<-function(gr,tri,stat_type,cutoff=0.05){
  gr$pvalue=elementMetadata(gr)[,paste(stat_type,'pval',sep='_')]
  tri_gr=gr[gr$mask_tri==tri]
  nontri=gr[gr$mask_tri!=tri]
  tri_ASM=sum(tri_gr$pvalue<=cutoff)
  tri_nonASM=sum(tri_gr$pvalue>cutoff)
  nontri_ASM=sum(nontri$pvalue<=cutoff)
  nontri_nonASM=sum(nontri$pvalue>cutoff)
  cont_table=matrix(c(tri_ASM,tri_nonASM,nontri_ASM,nontri_nonASM),byrow = T,nrow=2)
  return(fisher.test(cont_table))
}
#Merge 3 stats for each sample


density_plot_hyper<-function(GR_in,genes_hypervaribility_in,genomic_features,quant_dat="hyper_var_promoter"){
  #$hyper_var is used to calculating quantile
  GR_in=dist_calc(GR_in,genomic_features$TSS,k_round=100)
  plot_out=data.frame(NME=GR_in$mean_NME,dNME=GR_in$dNME,dist=GR_in$dist,dNME_pval=GR_in$dNME_pval,
                      MML=GR_in$mean_MML,dMML=GR_in$dMML,dMML_pval=GR_in$dMML_pval,
                      hyper_var=elementMetadata(GR_in)[,quant_dat])
  plot_out$quant=findInterval(plot_out$hyper_var,quantile(plot_out$hyper_var,prob=c(0,0.25,0.5,0.75),na.rm=T))
  quant=c("0-25%","25%-50%","50%-75%","75%-100%")
  plot_out$quant=quant[plot_out$quant]
  plot_out=plot_out[!is.na(plot_out$quant),]
  print(table(plot_out$quant))
  #within 1k of promoter
  NME=ggplot(plot_out,aes(x=NME,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
    ggtitle("NME distribution")
  dNME=ggplot(plot_out,aes(x=dNME,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
    ggtitle("dNME distribution")
  print(grid.arrange(NME,dNME,nrow=1))
  MML=ggplot(plot_out,aes(x=MML,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
    ggtitle("MML distribution")
  dMML=ggplot(plot_out,aes(x=dMML,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
    ggtitle("dMML distribution")
  print(grid.arrange(MML,dMML,nrow=1))
  dNME_dot=ggplot(plot_out,aes(x=hyper_var,y=dNME))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
    ggtitle("dNME vs hypervaribility")
  NME_dot=ggplot(plot_out,aes(x=hyper_var,y=NME))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
    ggtitle("NME vs hypervaribility")
  print(grid.arrange(NME_dot,dNME_dot,nrow=1))
  dMML_dot=ggplot(plot_out,aes(x=hyper_var,y=dMML))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
    ggtitle("dMML vs hypervaribility")
  MML_dot=ggplot(plot_out,aes(x=hyper_var,y=MML))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
    ggtitle("MML vs hypervaribility")
  print(grid.arrange(MML_dot,dMML_dot,nrow=1))
  #print(table(NME_plot_out[abs(NME_plot_out$dist)<=1000,"quant"]))
  dist_NME=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=NME,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("NME distribution vs distance to TSS")
  dist_dNME=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=dNME,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("dNME distribution vs distance to TSS")
  print(grid.arrange(dist_NME,dist_dNME,nrow=1))
  dist_MML=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=MML,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("MML distribution vs distance to TSS")
  dist_dMML=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=dMML,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("dMML distribution vs distance to TSS")
  print(grid.arrange(dist_MML,dist_dMML,nrow=1))
}  


#Look for gene being bind by those TF
motif_family_enrich<-function(hits,background,family_anno,pse_count=25){
  family_anno$Family=family_anno$Family
  bg_non_hit=background[!background %in% hits]
  df_out=data.frame(OR=NULL,pvalue=NULL,lowerCI=NULL,upperCI=NULL,family=NULL,hits_in=NULL)
  for(fam in unique(family_anno$Family)){
    
    family=family_anno$Name[family_anno$Family==fam]
    
    hit_in_family<-sum(hits %in% family)+pse_count
    hit_not_family<-sum(!hits %in% family)+pse_count
    bg_in_family<-sum(bg_non_hit %in% family)+pse_count
    bg_not_family<-sum(!bg_non_hit %in% family)+pse_count

    cont_table=matrix(c(hit_in_family,hit_not_family,bg_in_family,bg_not_family),nrow=2)
    fs_pse=fisher.test(cont_table)
    fs=fisher.test(cont_table-(pse_count-1))
    df_out=rbind(df_out,data.frame(OR=fs_pse$estimate[[1]],pvalue=fs$p.value,
                                   lowerCI=fs$conf.int[[1]],upperCI=fs$conf.int[[2]],family=fam,
                                   hits_in=hit_in_family-pse_count,hit_not_family=hit_not_family-pse_count,
                                   bg_in_family=bg_in_family-pse_count,bg_not_family=bg_not_family-pse_count,
                                   variance=fs_pse$conf.int[[2]]-fs_pse$conf.int[[1]]))
  }
  df_out=df_out[df_out$hits_in>1,]
  df_out$qvalue=p.adjust(df_out$pvalue,method='BH')
  return(df_out[order(df_out$pvalue,decreasing=F),])
}

split_ranges<-function(TSS_in,size=350){
  n=width(TSS_in)
  n=n-size
  start_TSS_in=start(TSS_in)
  strand_in=strand(TSS_in)
  chr_in=seqnames(TSS_in)
  #gene_name_TSS=TSS_in$gene_symbol
  #i=1
  for (num in seq(size,n,size)){
    start_TSS_in=c(start_TSS_in,start(TSS_in)+num)
    strand_in=c(strand_in,strand(TSS_in))
    chr_in=c(chr_in,seqnames(TSS_in))
    #gene_name_TSS=c(gene_name_TSS,paste(TSS_in$gene_symbol,i,sep='-'))
     #i=i+1
  }
  
  TSS_in_break=data.frame(seqnames=chr_in,start=start_TSS_in,end=start_TSS_in+size-1,strand=strand_in)
  TSS_in_break=makeGRangesFromDataFrame(TSS_in_break[which(TSS_in_break$start>0),],keep.extra.columns = T)
  TSS_in_break=TSS_in_break[seqnames(TSS_in_break)%in% seqlevels(TSS_in_break)[1:24]]#select only auto+xy
  return(TSS_in_break)
}
#Need to fix it, previous setting =200
# split_ranges_larger_5k<-function(TSS_in,size=350){
#   
#   TSS_out=split_ranges_5k(TSS_in,))
#   return(TSS_out)
# }

make_gr_dnase<-function(x){
  x=strsplit(x,':')[[1]]
  range_sp=strsplit(x[2],"-")[[1]]
  return(makeGRangesFromDataFrame(data.frame(seqnames=x[1],start=range_sp[1],end=range_sp[2])))
}

stat_hyper_enrich<-function(GR_merge_sp,stat='NME'){
  stat_in_quant=quantile(c(elementMetadata(GR_merge)[,paste(stat,'1',sep='')],elementMetadata(GR_merge)[,paste(stat,'2',sep='')]),prob=0.75)
  GR_merge_sp$hyper_var_genes=GR_merge_sp$hyper_var_promoter>=GR_merge_sp$hyper_var_upper
  GR_merge_sp=GR_merge_sp[!is.na(GR_merge_sp$hyper_var_genes)]
  # (elementMetadata(GR_merge_sp)[,paste(stat,'1',sep='')]>=stat_in_quant|
  #     elementMetadata(GR_merge_sp)[,paste(stat,'2',sep='')]>=stat_in_quant)
  ASM_hypervar=sum(GR_merge_sp$dNME_pval<=0.1 &GR_merge_sp$hyper_var_genes)
  ASM_not_hypervar=sum(GR_merge_sp$dNME_pval<=0.1&!GR_merge_sp$hyper_var_genes)
  
  notASM_hypervar=sum(GR_merge_sp$dNME_pval>0.1  &GR_merge_sp$hyper_var_genes)
  notASM_not_hypervar=sum(GR_merge_sp$dNME_pval>0.1 &
                            !GR_merge_sp$hyper_var_genes)
  cont_table=matrix(c(ASM_hypervar,ASM_not_hypervar,notASM_hypervar,notASM_not_hypervar),nrow=2)
 
  table_CMH=data.frame(subject=factor(tissue),
                       ASM=factor(c("ASM","ASM","NonASM","NonASM")),
                       feature=factor(c("Hypervar","NonHypervar","Hypervar","NonHypervar")),
                       count=c(ASM_hypervar,ASM_not_hypervar,notASM_hypervar,notASM_not_hypervar))
  print(table_CMH)

    if(all(table_CMH$count>1)){
      print(fisher.test(cont_table))
      return(table_CMH)
    }
  
}


RNA_df<-function(GR_in,RNA_in){
  #Checking overlapping with H1
  cat('Number of genes covered:',sum(GR_in$genes_promoter %in% rownames(RNA_in)),'\n')
  rna_asm_GR=rownames(RNA_in)[rownames(RNA_in) %in% GR_in$genes_promoter]
  #Make a summary df for those genes, GR_merge have most dNME
  rna_asm_hyper=data.table(dNME_promo=GR_in$dNME[GR_in$genes_promoter %in% rna_asm_GR])
  rna_asm_hyper$dNME_pval=GR_in$dNME_pval[GR_in$genes_promoter %in% rna_asm_GR]
  rna_asm_hyper$genes=GR_in$genes_promoter[GR_in$genes_promoter %in% rna_asm_GR]
  rna_asm_hyper$ASE_log2FC=RNA_in$log2FoldChange[match(rna_asm_hyper$genes,rownames(RNA_in))]
  rna_asm_hyper$dMML_pval=GR_in$dMML_pval[GR_in$genes_promoter %in% rna_asm_GR]
  rna_asm_hyper$dMML=GR_in$dMML[GR_in$genes_promoter %in% rna_asm_GR]
  rna_asm_hyper$dMML_relative=(GR_in$MML2-GR_in$MML1)[GR_in$genes_promoter %in% rna_asm_GR]
  rna_asm_hyper$dNME=GR_in$dNME[GR_in$genes_promoter %in% rna_asm_GR]
  rna_asm_hyper$dNME_relative=(GR_in$NME2-GR_in$NME1)[GR_in$genes_promoter %in% rna_asm_GR]
  return(rna_asm_hyper)
}


DEVG_analysis<-function(GR_in,tissue,stat_var='hyper_var_TSS'){
  if(!is.na(tissue)){GR_in=GR_in[GR_in$tissue==tissue]}
  GR_in$hyper_var_genes=elementMetadata(GR_in)[[stat_var]]>=GR_in$hyper_var_upper
  NME_quant=quantile(c(GR_in$NME1,GR_in$NME2),prob=0.75)
  high_NME_two_hyprvarible=GR_in[which(GR_in$NME1>=NME_quant&GR_in$NME2>=NME_quant& GR_in$hyper_var_genes)]
  high_NME_one_hyprvarible=GR_in[which((GR_in$NME1>=NME_quant|GR_in$NME2>=NME_quant)
                                       &(GR_in$dNME_pval<=0.01)&
                                         GR_in$hyper_var_genes)]
  low_NME_dNME=GR_in[which(GR_in$NME1<NME_quant&GR_in$NME2<NME_quant&
                             GR_in$dNME_pval<=pval_cutoff& GR_in$hyper_var_genes)]
  low_dNME_non_dNME=GR_in[which(GR_in$NME1<NME_quant&GR_in$NME2<NME_quant&
                                  GR_in$dNME_pval>pval_cutoff& GR_in$hyper_var_genes)]
  if(!is.na(tissue)){
    return(data.frame(high_NME_two_hyprvarible=length(high_NME_two_hyprvarible),
                      high_NME_one_hyprvarible=length(high_NME_one_hyprvarible),
                      low_NME_dNME=length(low_NME_dNME),
                      low_dNME_non_dNME=length(low_dNME_non_dNME),
                      tissue=tissue))
  }else{
    
    return(list(high_NME_two_hyprvarible=high_NME_two_hyprvarible,
                high_NME_one_hyprvarible=high_NME_one_hyprvarible,
                low_NME_dNME=low_NME_dNME,low_dNME_non_dNME=low_dNME_non_dNME))
  }
}


motif_reformat<-function(motif_in){
  motif_in$TF=gsub('\\(var.2\\)','',motif_in$TF)
  motif_in$TF=gsub('\\(var.3\\)','',motif_in$TF)
  motif_in=unlist(strsplit(motif_in$TF,"::"))
  return(motif_in)
}
motif_pref<-function(variant_in,motif_gene,motif_dir_in){
  olap=findOverlaps(variant_in,motif_gene)
  variant_in$alleleDiff=NA
  variant_in$alleleDiff[queryHits(olap)]=motif_gene$alleleDiff[subjectHits(olap)]
  motif_sig=motif_gene[motif_gene$geneSymbol %in% motif_dir_in$TF[motif_dir_in$qvalue<=0.1]]
  
  variant_in=subsetByOverlaps(variant_in,motif_sig)
  #variant_in=variant_in[variant_in$dNME_pval<=pval_cutoff]
  variant_in=variant_in[sign(variant_in$NME2-variant_in$NME1) == sign(variant_in$alleleDiff)]
  return(variant_in)
}
convert_GR<-function(x){
  strand=unlist(lapply(strsplit(x,','),function(x)x[2]))
  strand=ifelse(is.na(strand),"*",strand)
  x=ifelse(is.na(strand),x,sub(paste0(',\\','-'),'',sub(paste0(',\\','+'),'',x)))

  gr=GRanges(seqnames=sub(':.*','',x),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',x))),
                     end=as.numeric(sub('.*-','',x))),strand=strand)
  return(gr)
  
}
trait_variant=function(x,traits) {
  variant=NA
  
  variant=tryCatch(get_variants(efo_id = x),error=function(e) NA)
  if(!is.na(variant)){
    if(length(variant@variants$variant_id)>0){
      cat('Processing',traits$trait[traits$efo_id==x],'\n')
      
      var_out=as.data.frame(variant@variants)
      var_out=var_out[!is.na(var_out$chromosome_position)&!is.na(var_out$chromosome_name),]
      if(length(var_out$chromosome_position)>0&length(var_out$chromosome_name)>0&length(var_out$chromosome_name)==length(var_out$chromosome_position)){
        var_gr=makeGRangesFromDataFrame(var_out,
                                        seqnames.field = 'chromosome_name',start.field = 'chromosome_position',end.field = 'chromosome_position')
        seqlevels(var_gr)=paste('chr',seqlevels(var_gr),sep='')
        var_gr$efo_id=traits$trait[traits$efo_id== x]
        var_gr$rsid=variant@variants$variant_id
        #check genomic version
        snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
        var_out=snpsById(snps,variant_id[1])
        if(length(subsetByOverlaps(var_gr[1],var_out))==0){
          
        }
    return(var_out)
      }
    }
  }
}

agnostic_matrix_conversion<-function(gr_in,stat='NME'){
  gr_out=granges(unique(gr_in))
  olap=findOverlaps(gr_in,gr_out,type='equal')
  stat_in_df=elementMetadata(gr_in[queryHits(olap)])[c(stat,'Sample')]
  stat_in_df$idx=NA
  stat_in_df$idx[queryHits(olap)]=subjectHits(olap)
  stat_in_df=as.data.table(stat_in_df)
 
  stat_in_df_stat=dcast.data.table(data=stat_in_df,formula=idx~Sample,value.var = stat,fun.aggregate=mean)#remove agg.fun for new run
  gr_out=gr_out[stat_in_df_stat$idx]
  mcols(gr_out)=stat_in_df_stat[,-1]
  return(gr_out)
 
}
PCA_df_prep<-function(UC_in_matrix_sub){
  UC_in_matrix_sub_df=as.data.frame(t(as.matrix(mcols(UC_in_matrix_sub))))
  UC_in_matrix_sub_PCA=prcomp(UC_in_matrix_sub_df, scale. = TRUE)
  UC_in_matrix_sub_df$sample=unlist(lapply(strsplit(rownames(UC_in_matrix_sub_df),'-'),function(x) paste(x[-length(x)],collapse = '-')))
 
  return(list(UC_in_matrix_sub_PCA,UC_in_matrix_sub_df))
}

col_fun <- colorRampPalette(
  c(
    ##rev(brewer.pal(8,"RdYlBu"))[1:4],
    rev(brewer.pal(8,"RdYlBu"))[1:2],
    #"white",
    rev(brewer.pal(8,"RdYlBu"))[5:8]
  )
)

#UC
read.agnostic.mouse.uc<-function(file_in,matrix=FALSE,fileter_N=1,gff_in=NA){
  
  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    #process file name
    file_in=strsplit(file_in,'\\/')[[1]]
    file_in=file_in[length(file_in)]
    comp= strsplit(strsplit(file_in,'_uc.bedGraph')[[1]],'-vs-')[[1]]
    strain=unlist(lapply(strsplit(comp,'_'),function(x) x[1]))
    #if  contain BL6DBA, use ref is BL6DBA
    strain=ifelse('BL6DBA'%in%strain,'BL6DBA','mm10')
    comp=unlist(lapply(strsplit(comp,'_'),function(x) paste(x[-1],collapse = '_')))
    comp=comp[comp!='']
    comp_stage=unlist(lapply(comp,function(x) {x_split=strsplit(x,'_')[[1]]
    x_split=x_split[-length(x_split)][-1]
    x_split=paste(x_split,collapse = '_')
    return(x_split)}))
    tissue1=strsplit(comp[1],'_')[[1]][1]
    tissue2=strsplit(comp[2],'_')[[1]][1]
    #if BL6DBA, the 1st comp_stage is empty
    
    comp_stage=gsub('_5','.5',comp_stage)
    comp_stage=gsub('day','E',comp_stage)
    comp_stage=gsub('E0','P0',comp_stage)
    replicate=strsplit(comp[1],'_')[[1]][length(strsplit(comp[1],'_')[[1]])]
    replicate=gsub('merged','',replicate)
    informME_in$Sample=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
    informME_in=informME_in[informME_in$N>=fileter_N]
    cat('Minimum N:',min(informME_in$N),'\n')
    #informME_in$Ref=strain
    if(matrix){
      informME_in_dt=as.data.table(mcols(informME_in))[,c("score","Sample")]
      colnames(informME_in_dt)=c("UC","Sample")
      informME_in_dt$UC=as.numeric(informME_in_dt$UC)
      informME_in_dt$region=paste0(seqnames(informME_in),":",start(informME_in),"-",end(informME_in))
      informME_in_dt=informME_in_dt[match(gff_in,region),"UC"]
      colnames(informME_in_dt)=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
      return(informME_in_dt)
    }
    else{return(informME_in)}
  }
}

#Get CpG sites from hg19
getCpgSitesmm10 <- function(chrsOfInterest=paste("chr",1:19,sep="")){
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x)  GRanges(x,IRanges(start(matchPattern("CG", Mmusculus[[x]])),with=2)))
  # Set genome and seqlengths
  cgs <- setGenomeLengths(do.call('c',cgs),chrsOfInterest=chrsOfInterest,genome_in="mm10")
  # Return CpG site GR
  return(cgs)
}
