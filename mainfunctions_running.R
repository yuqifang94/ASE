#Dependencies
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
if (!requireNamespace("readr", quietly = TRUE))
{install.packages("readr")}
library(readr)
if (!requireNamespace("readxl", quietly = TRUE))
{install.packages("readxl")}
library(readxl)
if (!requireNamespace("annotatr", quietly = TRUE))
{BiocManager::install("annotatr",INSTALL_opts = c('--no-lock'))}
library(annotatr)
if (!requireNamespace("AnnotationHub", quietly = TRUE))
{BiocManager::install("AnnotationHub")}
library(AnnotationHub)
if (!requireNamespace("shiny", quietly = TRUE))
{BiocManager::install("shiny",version="devel")}
library(shiny)
if (!requireNamespace("rtracklayer", quietly = TRUE))
{BiocManager::install("rtracklayer")}
library(rtracklayer)
library(GenomicRanges)
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))
{BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")}
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
if (!requireNamespace("Homo.sapiens", quietly = TRUE))
{BiocManager::install("Homo.sapiens")}
library(Homo.sapiens)
if (!requireNamespace("GenomicFeatures", quietly = TRUE))
{BiocManager::install("GenomicFeatures")}
library(GenomicFeatures)
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
{BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
library(BSgenome.Hsapiens.UCSC.hg19)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
{BiocManager::install("org.Hs.eg.db")}
library(org.Hs.eg.db)
if (!requireNamespace("VariantAnnotation", quietly = TRUE))
{BiocManager::install("VariantAnnotation")}
library(VariantAnnotation)
if (!requireNamespace("reshape2", quietly = TRUE))
{BiocManager::install("reshape2")}
library(reshape2)
if (!requireNamespace("gridExtra", quietly = TRUE))
{BiocManager::install("gridExtra")}
library(gridExtra)
if (!requireNamespace("grid", quietly = TRUE))
{BiocManager::install("grid")}
library(grid)
if (!requireNamespace("ggplot2", quietly = TRUE))
{install.packages("ggplot2")}
library(ggplot2)
if (!requireNamespace("lattice", quietly = TRUE))
{install.packages("lattice")}
library(lattice)
if (!requireNamespace("Repitools", quietly = TRUE))
{BiocManager::install("Repitools")}
library(Repitools)
if (!requireNamespace("topGO", quietly = TRUE))
{BiocManager::install("topGO")}
library(topGO)
if(!require(psych)){install.packages("psych")}
library('psych')
if(!require(vcd)){install.packages("vcd")}
if(!require(DescTools)){install.packages("DescTools")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(ggpubr)){install.packages("ggpubr")}
library("ggpubr")
if(!require(readxl)){install.packages("readxl")}
library(readxl) 
if(!require(ggpointdensity)){install.packages("ggpointdensity")}
library(ggpointdensity)
if (!requireNamespace("rethinking", quietly = TRUE))
{install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
  library(devtools)
  devtools::install_github("rmcelreath/rethinking")}
library(rethinking)
if (!requireNamespace("gwasrapidd", quietly = TRUE)){
  remotes::install_github("ramiromagno/gwasrapidd")}
library(gwasrapidd)
library(qusage)
if (!requireNamespace("gwascat", quietly = TRUE)){
  BiocManager::install("gwascat")}
library(gwascat)
if (!requireNamespace("liftOver", quietly = TRUE)){
  BiocManager::install("liftOver")}
library(liftOver)
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38", quietly = TRUE)){
  BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")}
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
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
NME_agnostic_file="../downstream/input/NME_allele_agnostic_merge_2k.rds"
MML_agnostic_file="../downstream/input/MML_allele_agnostic_merge_2k.rds"
motif_gene_file='../downstream/output/motif_all_JASPAR_default.rds' #For all SNP
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
  vcf=het_CpG_df(vcf)
  vcf$sub=sub
  cat('Done processing',sub,'using:', proc.time()[[3]]-tt1,'\n')
  return(vcf)
  
}
#ID heterozygous CpG from vcf file
het_CpG_df<-function(var_gr){
  var_gr=chr_check(var_gr)
  plus_loc=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,1)))
  minus_loc=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,-1)))
  #get dinucleotide for ref, alt, plus and minus, find some examples region randomly: check if those match
  var_gr$REF_plus=paste_nucleotide(var_gr$REF,plus_loc,'plus')
  var_gr$REF_minus=paste_nucleotide(var_gr$REF,minus_loc,'minus')
  var_gr$ALT_plus=paste_nucleotide(var_gr$ALT,plus_loc,'plus')
  var_gr$ALT_minus=paste_nucleotide(var_gr$ALT,minus_loc,'minus')
  #get trinucleotide
  var_gr$REF_tri=paste_trinucleotide(minus_loc,var_gr$REF,plus_loc)
  var_gr$ALT_tri=paste_trinucleotide(minus_loc,var_gr$ALT,plus_loc)
  #check if heterogygouze: note rowSum =2 have trinucleotide form CGG with ref =G alt =C
  var_gr$npmCG=rowSums(as.data.frame(var_gr)[,c('REF_plus','REF_minus','ALT_plus','ALT_minus')]=='CG')
  var_gr$HetCpG=var_gr$npmCG>0
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
#Get dinucleotide 
paste_nucleotide<-function(ref,seq,direction){
  df=data.frame(ref=ref,seq=unlist(seq))
  if(direction=='plus'){return(paste(df$ref,df$seq,sep=''))}
  else if(direction=='minus'){return(paste(df$seq,df$ref,sep=''))}
  else{print('wrong direction')}
}
#Get trinucleotide
paste_trinucleotide<-function(minus,ref,plus){
  df=data.frame(minus=unlist(minus),ref=ref,plus=unlist(plus))
  return(paste(df$minus,df$ref,df$plus,sep=''))
}
# Function to set genome and chromosome lengths to a GR object
setGenomeLengths <- function(GR,chrsOfInterest=paste("chr",1:22,sep="")){
  # Get genome info
  GR=chr_check(GR)
  hg19<-getBSgenome("hg19")
  genome.seqinfo <- seqinfo(hg19)
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
  # H1
  H1_subject <- rep("H1",3)
  H1_subject_labels <- rep("CL3",3)
  H1_tissues <- c(
    "rep1","rep2",
    "merged"
  )             
  H1_tissue_labels <- c(
    "Embryonic Stem Cell", "Embryonic Stem Cell", "Embryonic Stem Cell"
  )
  H1_gtex_labels <- c(
    "", "Embryonic Stem Cell", "Embryonic Stem Cell"
  )
  H1_diff_labels <- c(
    "Undifferentiated","Undifferentiated","Undifferentiated"
  )
  # GM12878
  # GM12878_subject <- rep("GM12878",1)
  # GM12878_subject_labels <- rep("CL4",1)
  # GM12878_tissues <- c(
  #   #"1","2",
  #   "merged"
  # )             
  # GM12878_tissue_labels <- c(
  #   "EBV B-cell","EBV B-cell","EBV B-cell"
  # )
  # GM12878_gtex_labels <- c(
  #   "","",""
  # )
  # GM12878_diff_labels <- c(
  #   "","",""
  # )
  #HuFGM02
  HuFGM02_subject <- rep("HuFGM02",2)
  HuFGM02_subject_labels <- rep("Brain1",2)
  HuFGM02_tissues <- c(
    "brain_germinal_matrix_tissue_paired",
    "brain_cerebellum_tissue_paired"
  )             
  HuFGM02_tissue_labels <- c(
    "brain_germinal_matrix",
    "brain_cerebellum_tissue"
  )
  HuFGM02_gtex_labels <- c(
    "",""
  )
  HuFGM02_diff_labels <- c(
    "differentiated",'differentiated'
  )
  #112
  b112_subject <- rep("112",1)
  b112_subject_labels <- rep("Brain2",1)
  b112_tissues <- c(
    "Brain_substantia_nigra_paired"
  )             
  b112_tissue_labels <- c(
    "Brain_substantia_nigra"
  )
  b112_gtex_labels <- c(
    ""
  )
  b112_diff_labels <- c(
    "differentiated"
  )
  
  #149
  b149_subject <- rep("149",1)
  b149_subject_labels <- rep("Brain3",1)
  b149_tissues <- c(
    "Brain_Hippocampus_middle_paired"
  )             
  b149_tissue_labels <- c(
    "Brain_Hippocampus_middle"
  )
  b149_gtex_labels <- c(
    ""
  )
  b149_diff_labels <- c(
    "differentiated"
  )
  #150
  b150_subject <- rep("150",1)
  b150_subject_labels <- rep("Brain4",1)
  b150_tissues <- c(
    "Brain_Hippocampus_middle_paired"
  )             
  b150_tissue_labels <- c(
    "Brain_Hippocampus_middle"
  )
  b150_gtex_labels <- c(
    ""
  )
  b150_diff_labels <- c(
    "differentiated"
  )
  # H9
  H9_subject <- rep("H9",1)
  H9_subject_labels <- rep("CL1",1)
  H9_tissues <- c(
    "42_embryonic_stem_cell_single"
  )             
  H9_tissue_labels <- c(
    "Embryonic Stem Cell"
  )
  H9_gtex_labels <- c(
    ""
  )
  H9_diff_labels <- c(
    "Undifferentiated"
  )
  
  # HUES64
  HUES64_subject <- rep("HUES64",4)
  HUES64_subject_labels <- rep("CL2",4)
  HUES64_tissues <- c(
    "stem_27_undifferentiated_paired",
    "ectoderm_paired",
    "endoerm_27_paired",
    "mesoderm_23_paired"
  )
  HUES64_tissue_labels <- c(
    "Embyonic Stem Cell",
    "Ectoderm",
    "Endoderm",
    "Mesoderm"
  )
  HUES64_gtex_labels <- c("","","","")
  HUES64_diff_labels <- c(
    "Undifferentiated",
    "Semidifferentiated",
    "Semidifferentiated",
    "Semidifferentiated"
  )
  
  # skin03
  skin03_subject <- rep("skin03",2)
  skin03_subject_labels <- rep("D11",2)
  skin03_tissues <- c(
    "foreskin_keratinocyte_paired",
    "foreskin_melanocyte_paired"
  )
  skin03_tissue_labels <- c(
    "Foreskin Keratinocyte",
    "Foreskin Melanocyte"
  )
  skin03_gtex_labels <- c("","")
  skin03_diff_labels <- rep("Differentiated",length(skin03_subject))
  
  # STL001
  stl001_subject <- rep("STL001",11)
  stl001_subject_labels <- rep("D5",11)
  stl001_tissues <- c(
    "Adipose_single",
    "Small_Intestine_single",
    "Bladder_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Lung_single",
    "Psoas_Muscle_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Spleen_single",
    "Thymus_single"
  )
  stl001_tissue_labels <- c(
    "Adipose",
    "Small Intestine",
    "Bladder",
    "Gastric",
    "Left Ventricle",
    "Lung",
    "Psoas Muscle",
    "Right Ventricle",
    "Sigmoid Colon",
    "Spleen",
    "Thymus"
  )
  stl001_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "",
    "",
    "",
    "Heart_Left_Ventricle",
    "Lung",
    "",
    "",
    "Colon_Transverse",
    "",
    ""
  )
  stl001_diff_labels <- rep("Differentiated",length(stl001_subject))
  
  # STL002
  stl002_subject <- rep("STL002",11)
  stl002_subject_labels <- rep("D6",11)
  stl002_tissues <- c(
    "Adipose_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Lung_single",
    "Ovary_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Small_Intestine_single",
    "Spleen_single"
  )
  stl002_tissue_labels <- c(
    "Adipose",
    "Adrenal Gland",
    "Aorta",
    "Esophagus",
    "Gastric",
    "Lung",
    "Ovary",
    "Pancreas",
    "Psoas Muscle",
    "Small Intestine",
    "Spleen"
  )
  stl002_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Esophagus_Mucosa",
    "",
    "Lung",
    "Uterus",
    "Pancreas",
    "",
    "",
    ""
  )
  stl002_diff_labels <- rep("Differentiated",length(stl002_subject))
  
  # STL003
  stl003_subject <- rep("STL003",13)
  stl003_subject_labels <- rep("D7",13)
  stl003_tissues <- c(
    "Adipose_Tissue_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Right_Atrium_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Small_Intestine_single",
    "Spleen_single"
  )
  stl003_tissue_labels <- c(
    "Adipose",
    "Adrenal Gland",
    "Aorta",
    "Esophagus",
    "Gastric",
    "Left Ventricle",
    "Pancreas",
    "Psoas Muscle",
    "Right Atrium",
    "Right Ventricle",
    "Sigmoid Colon",
    "Small Intestine",
    "Spleen"
  )
  stl003_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Esophagus",
    "",
    "Heart_Left_Ventricle",
    "Pancreas",
    "",
    "",
    "",
    "Colon_Transverse",
    "",
    ""
  )
  stl003_diff_labels <- rep("Differentiated",length(stl003_subject))
  
  # STL011
  stl011_subject <- rep("STL011",1)
  stl011_subject_labels <- rep("D8",1)
  stl011_tissues <- c(
    "Liver_single"
  )
  stl011_tissue_labels <- c(
    "Liver"
  )
  stl011_gtex_labels <- c(
    "Liver"
  )
  stl011_diff_labels <- rep("Differentiated",length(stl011_subject))
  
  # Create single vectors
  subjects <- c(H1_subject,H9_subject,HUES64_subject,skin03_subject,stl001_subject,
                stl002_subject,stl003_subject,stl011_subject,HuFGM02_subject,b112_subject,b150_subject,b149_subject)
  tissues <- c(H1_tissues,H9_tissues,HUES64_tissues,skin03_tissues,stl001_tissues,
               stl002_tissues,stl003_tissues,stl011_tissues,HuFGM02_tissues,b112_tissues,b150_tissues,b149_tissues)
  subject_labels <- c(H1_subject_labels,H9_subject_labels,HUES64_subject_labels,skin03_subject_labels,
                      stl001_subject_labels,stl002_subject_labels,stl003_subject_labels,
                      stl011_subject_labels,HuFGM02_subject_labels,b112_subject_labels,b150_subject_labels,b149_subject_labels)
  tissue_labels <- c(H1_tissue_labels,H9_tissue_labels,HUES64_tissue_labels,skin03_tissue_labels,
                     stl001_tissue_labels,stl002_tissue_labels,stl003_tissue_labels,
                     stl011_tissue_labels,HuFGM02_tissue_labels,b112_tissue_labels,b150_tissue_labels,b149_tissue_labels)
  gtex_labels <- c(H1_gtex_labels,H9_gtex_labels,HUES64_gtex_labels,skin03_gtex_labels,
                   stl001_gtex_labels,stl002_gtex_labels,stl003_gtex_labels,
                   stl011_gtex_labels,HuFGM02_gtex_labels,b112_gtex_labels,b150_gtex_labels,b149_gtex_labels)
  diff_labels <- c(H1_diff_labels,H9_diff_labels,HUES64_diff_labels,skin03_diff_labels,stl001_diff_labels,
                   stl002_diff_labels,stl003_diff_labels,stl011_diff_labels,HuFGM02_diff_labels,b112_diff_labels,b149_diff_labels,b150_diff_labels)
  GRs=GRanges()
  for (i in 1:length(subjects)) {
    # Print sample being loaded
    print(paste("Loading sample:",subjects[i],tissues[i]))
    if (calc=='diff'){
      GR.in=read.diffGR(subjects[i],tissues[i],inDir,cutoff=0.05)
    }else if(calc=='allele'){
      GR.in=read.alleleGR(subjects[i],tissues[i],inDir)
    }else {cat('Wrong calc \n')}
    GR.in$SubjectLabel <- subject_labels[i]
    GR.in$Tissue <- tissue_labels[i]
    GR.in$GTEx <- gtex_labels[i]
    GR.in$State <- diff_labels[i]
    GRs=append(GRs,GR.in)
  }
  return(GRs)
  
}
#Function to read in single GR object:
read.diffGR<-function(subjects,tissues,inDir,cutoff=0.05,chrsOfInterest=paste("chr",1:22,sep="")){
  #Initialization
  GRs=GRanges()
  #Make sure the inputs are unique
  # dmml
  filename=paste(inDir,subjects,"_",tissues,"_phased_tmml_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dMML')
  
  GRs <- append(GRs,GR.in)
  # dnme
  filename=paste(inDir,subjects,"_",tissues,"_phased_tnme_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dNME')
  GRs <- append(GRs,GR.in)
  # uc
  filename=paste(inDir,subjects,"_",tissues,"_phased_tpdm_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'UC')
  GRs <- append(GRs,GR.in)
  #Check ASM
  GRs <- GRs[!is.na(GRs$pvalue)]
  GRs$ASM <- NA
  GRs[GRs$pvalue<=cutoff]$ASM <- "Yes"
  GRs[GRs$pvalue>cutoff]$ASM <- "No"
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
  GRs <- append(GRs, GR.in)
  #MML2
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'MML')
  GR.in$Genome="2"
  GRs <- append(GRs, GR.in)
  #NME1
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="1"
  GRs <- append(GRs, GR.in)
  #NME2
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="2"
  GRs <- append(GRs, GR.in)
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
  GR_out$sample_name<-paste(tissues,subjects, sep=' - ')
  GR_out$Statistic <- Statistic
  GR_out$Value<- GR$score
  #for differential analysis, make sure the GR_out are unique
  if(pvalue){
    GR_out$pvalue <- as.numeric(GR$NA.)
    GR_out <- GR_out[!duplicated(GR[,c()])]
  }
  else{GR_out$N <- GR$NA.} #not diff for new samples
  
  #GR_out <- resize(GR_out, width(GR_out) + 1, fix="end")
  #GR_out <- resize(GR_out, width(GR_out) + 1, fix="start")
  GR_out <- setGenomeLengths(GR_out,chrsOfInterest=chrsOfInterest)
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
#Count number of CpG in gff, find examples to check
gff_hetCpG_count<-function(sub,gff_in,vcf_in,CpG){
  cat('Processing',sub,'\n')
  #Read in gff file
  gff_sub=gff_in[gff_in$Subject==sub]
  #Read in vcf het CpG information
  vcf_sub=vcf_in[[sub]]
  #Filter het CpG regions
  vcf_sub_het=vcf_sub[vcf_sub$HetCpG]
  #For each gff region count overlaps
  #gff_het_count=countOverlaps(gff_sub,vcf_sub_het)
  #cat(sub,':',sum(gff_het_count>1)/length(gff_het_count)*100, '% all gff region have more than 1 het CpG\n', sep='')
  #cat(sub,':',sum(gff_het_count>1)/sum(gff_het_count>0)*100, '% all het CpG region have more than 1 het CpG\n', sep='')
  #cat(sub,':',sum(gff_het_count>0)/length(gff_het_count)*100, '% gff region have het CpG\n', sep='')
  olap=findOverlaps(vcf_sub_het,gff_sub,type='within',select='all')
  #Count number of CG, here g1CG=genome1 and g2CG=genome2, however, we need to split based on the GT, also ref
  vcf_sub_het$g1CG=(vcf_sub_het$genome1_plus=='CG')+(vcf_sub_het$genome1_minus=='CG')
  vcf_sub_het$g2CG=(vcf_sub_het$genome2_plus=='CG')+(vcf_sub_het$genome2_minus=='CG')
  vcf_sub_het$refCG=(vcf_sub_het$REF_plus=='CG')+(vcf_sub_het$REF_minus=='CG')
  vcf_sub_het$altCG=(vcf_sub_het$ALT_plus=='CG')+(vcf_sub_het$ALT_minus=='CG')
  #Add Het CpG information
  gff_sub$HetCpG=FALSE
  gff_sub[subjectHits(olap)]$HetCpG=TRUE
  #g1CG, SNP have CG in ref 
  #g2CG SNP have CG in alt
  gff_sub$g1CG=0
  gff_sub$g2CG=0
  gff_sub$refCG=0
  gff_sub$altCG=0
  #Count number of Het CpG here *Check 
  df_sub=data.frame(subjHits=subjectHits(olap),
                    g1CG=vcf_sub_het$g1CG[queryHits(olap)],g2CG=vcf_sub_het$g2CG[queryHits(olap)],
                    refCG=vcf_sub_het$refCG[queryHits(olap)],altCG=vcf_sub_het$altCG[queryHits(olap)])
  
  agg_sub=aggregate(df_sub,by=list(df_sub$subjHits),FUN=sum)
  gff_sub$g1CG[agg_sub$Group.1]=agg_sub$g1CG#agg_sub$Group.1 is the unique subject hits
  gff_sub$g2CG[agg_sub$Group.1]=agg_sub$g2CG
  gff_sub$refCG[agg_sub$Group.1]=agg_sub$refCG#agg_sub$Group.1 is the unique subject hits
  gff_sub$altCG[agg_sub$Group.1]=agg_sub$altCG
  gff_sub$N_hg19=countOverlaps(gff_sub,CpG)
  gff_sub$N_nonhet=gff_sub$N_hg19-countOverlaps(gff_sub,vcf_sub_het[vcf_sub_het$refCG>0])
  
  return(gff_sub)
}
#Merge stats into single granges object
stat_merge<-function(gr_in,allele_in){
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
  
  gr$NME1=gr$MML1=gr$NME2=gr$MML2=NA
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
  if (length(UC$Value)>0){gr$UC=UC$Value}
  if (length(UC$pvalue)>0){gr$UC=UC$pvalue}
  
  return(gr)
}
#add genes to GR in the format of list
add_gene_GR<-function(GR_merge_in,feature_in,feature_names){
  elementMetadata(GR_merge_in)[[feature_names]]=NA
  olap_promoter=findOverlaps(GR_merge_in,feature_in)
  df_idx=data.frame(qt=queryHits(olap_promoter),
                    genes=as.character(feature_in$gene_name[subjectHits(olap_promoter)]),
                    stringsAsFactors = F)
  df_idx=aggregate(df_idx$genes,by=list(df_idx$qt),FUN=list)
  colnames(df_idx)=c('qt','genes')
  df_idx$genes=lapply(df_idx$genes,as.character)
  elementMetadata(GR_merge_in)[[feature_names]][df_idx$qt]=df_idx$genes
  return(GR_merge_in)
}
#Add hypervaribility to each region
add_hyper_var<-function(GR_in,hyper_var_in,sample_name,upper_cutoff=0.75,lower_cutoff=0.25){
  GR_sp=GR_in[GR_in$Sample %in%sample_name]
  hyper_var=read_hypervar(hyper_var_in)
  
  #check output
  GR_in$hyper_var_promoter[GR_in$Sample %in%sample_name]=unlist(lapply(GR_in$genes_promoter[GR_in$Sample %in%sample_name],
                                                                       function(x) mean(hyper_var$hypervar_logvar[hyper_var$gene_name %in% x])))
  
  GR_in$hyper_var_TSS[GR_in$Sample %in%sample_name]=unlist(lapply(GR_in$TSS[GR_in$Sample %in%sample_name],
                                                                  function(x) mean(hyper_var$hypervar_logvar[hyper_var$gene_name %in% x])))
  GR_in$hyper_var_body[GR_in$Sample %in%sample_name]=unlist(lapply(GR_in$genes_body[GR_in$Sample %in%sample_name],
                                                                   function(x) mean(hyper_var$hypervar_logvar[hyper_var$gene_name %in% x])))
  GR_in$hyper_var_lower[GR_in$Sample %in%sample_name]=quantile(hyper_var$hypervar_logvar,prob=lower_cutoff)
  GR_in$hyper_var_upper[GR_in$Sample %in%sample_name]=quantile(hyper_var$hypervar_logvar,prob=upper_cutoff)
  GR_in$hyper_var_median[GR_in$Sample %in%sample_name]=median(hyper_var$hypervar_logvar)
  GR_in$hyper_var_fn[GR_in$Sample %in%sample_name]=rep(list(hyper_var_in),length(GR_in$hyper_var_fn[GR_in$Sample %in%sample_name]))
  return(GR_in)
}
read_hypervar<-function(hyper_var_in){
  hyper_var=data.frame()
  for (fn in hyper_var_in){
    fn_in=readRDS(fn)
    fn_in$gene_name=rownames(fn_in)
    rownames(fn_in)=NULL
    hyper_var=rbind(hyper_var,fn_in)
  }
  
  #Find the dataset that have most regions overlap with NME regions, check correlation?
  hyper_var=aggregate(hyper_var[,c('mean','var','hypervar_var','hypervar_logvar')],by=list(hyper_var$gene_name),FUN=mean,na.rm=T)
  colnames(hyper_var)[1]='gene_name'
  return(hyper_var)
}
#Count number of Heterozygous CpG at each gff region, notice the boundary condition, keep in mind gff file is shorterned by 1 due to bedGraph adjustment
hetCGallele_merged<-function(sub,gr_merge,gff,CpG,vcf_in,gene_size=500){
  cat('Analyzing',sub,'\n')
  t1=proc.time()[[3]]
  #Import vcf file
  sub_vcf=vcf_in[[sub]]
  het_vcf=sub_vcf[sub_vcf$HetCpG]
  sub_allele=gr_merge[gr_merge$Subject==sub]
  sub_het=gff[[sub]]
  #gr_allele got resized with start and end +1, use +2 to include equal start & end, no longer needed for new output?
  #sub_het=resize(sub_het, width(sub_het) + 4, fix="center")
  gr_merge_olap=findOverlaps(sub_allele,sub_het,type='equal')
  sub_allele=sub_allele[queryHits(gr_merge_olap)]
  sub_allele$g1CG=sub_het$g1CG[subjectHits(gr_merge_olap)]
  sub_allele$g2CG=sub_het$g2CG[subjectHits(gr_merge_olap)]
  sub_allele$N_hg19=NULL
  sub_allele$N_hg19=sub_het$N_hg19[subjectHits(gr_merge_olap)]
  sub_allele$N_nonhet=NULL
  sub_allele$N_nonhet=sub_het$N_nonhet[subjectHits(gr_merge_olap)]
  sub_allele$gff_size=width(sub_allele)
  gr_out=GR_resize_merged(sub_allele,CpG,het_vcf,gene_size=gene_size,sub)#Check for calculating density
  cat('Finish analyzing',sub,proc.time()[[3]]-t1,'\n')
  return(gr_out)
}
#Generate density from GR allele, at least need 200 bp for density
#Count number of hetCG at each allele GR from result, add gff size, and N
GR_resize_merged<-function(GR.in,CpG_sites,hetCpG,sub,gene_size=500){
  ##From definitino of CpG island, use 200 bp regions
  GR.extend=resize(GR.in,width=width(GR.in)+gene_size,fix='center')
  GR.in$CG_hg19_extend=countOverlaps(GR.extend,CpG_sites)
  #Change here for ref
  GR.in$CG_nonhet_extend=GR.in$CG_hg19_extend-countOverlaps(GR.extend,hetCpG[hetCpG$genome1_plus=='CG'|hetCpG$genome1_minus=='CG'])
  gr_seq=getSeq(Hsapiens,GR.extend,as.character=T)
  GR.in$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
  #Count CpG in genome 1
  GR.in$CG_allele_extend_g1=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[hetCpG$genome1_plus=='CG'|hetCpG$genome1_minus=='CG'])
  GR.in$CG_allele_extend_g2=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[hetCpG$genome2_plus=='CG'|hetCpG$genome2_minus=='CG'])
  GR.in$gff_size_extend=width(GR.extend)
  return(GR.in)
}
#Give each SNP an ASM information for each subject
variant_meta<- function(subj,variant_in,GR_in){ #variant_in for each subject, GR_in for each subject
  cat('Processing',subj,'\n')
  GR_in_subj=GR_in[GR_in$Subject==subj]
  variant_in_subj=variant_in[[subj]]
  sp=unique(GR_in_subj$Sample)
  gr_out_sp = GRanges()
  for (sps in sp){
    cat('Processing',sp,'\n')
    gr_out_sp=c(gr_out_sp,variant_meta_sp(variant_in_subj,GR_in_subj[GR_in_subj$Sample==sps]))
    
  }
  return(gr_out_sp)
  
  
}

#for given stat: need to debug
variant_meta_sp <-function(variant_subj,GR_sp){
  olap=findOverlaps(variant_subj,GR_sp,maxgap =0,type='within')
  gr_out=variant_subj[queryHits(olap)]
  #Find the overlap region
  #gr_out=variant_subj[queryHits(olap)]
  olap=subjectHits(olap)
  gr_out$Sample=GR_sp$Sample[olap]
  #gr_out$pvalue=GR_sp$pvalue[olap]
  gr_out$dNME_pval=GR_sp$dNME_pval[olap]
  gr_out$dMML_pval=GR_sp$dMML_pval[olap]
  gr_out$UC_pval=GR_sp$UC_pval[olap]
  gr_out$dMML=GR_sp$dMML[olap]
  gr_out$dNME=GR_sp$dNME[olap]
  gr_out$UC=GR_sp$UC[olap]
  gr_out$NME1=GR_sp$NME1[olap]
  gr_out$NME2=GR_sp$NME2[olap]
  gr_out$MML1=GR_sp$MML1[olap]
  gr_out$MML2=GR_sp$MML2[olap]
  gr_out$g1CG=GR_sp$g1CG[olap]
  gr_out$g2CG=GR_sp$g2CG[olap]
  gr_out$tissue=GR_sp$tissue[olap]
  gr_out$N=GR_sp$N[olap]
  gr_out$genes_promoter=GR_sp$genes_promoter[olap]
  gr_out$genes_body=GR_sp$genes_body[olap]
  gr_out$TSS=GR_sp$TSS[olap]
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
  
  gr_out$HetCpG=FALSE
  #This may be different from previous result?
  gr_out$HetCpG=((gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & !(gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG')) |
    (!(gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & (gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG'))
  return(gr_out)
}
#Read in allele-agnositc model for hg19 analysis
read.agnostic<-function(file_in,GR_merge_in,stat='NME',allele_include=T){
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    olap=findOverlaps(GR_merge_in,informME_in)
    print(length(unique(subjectHits(olap)))/length(informME_in))
    #add GR_merge data
    if(allele_include){
      #replace value instead of remove regions
      informME_in$score[subjectHits(olap)]=
        rowMeans(as.matrix(elementMetadata(GR_merge_in)[paste(stat,c('1','2'),sep='')]))[queryHits(olap)]
    }else( informME_in=informME_in[-subjectHits(olap)])
    
    quant=c("0-25%","25%-50%","50%-75%","75%-100%")
    informME_in$quant_score=findInterval(informME_in$score,quantile(informME_in$score,prob=c(0,0.25,0.5,0.75),na.rm=T))
    informME_in$quant_score=quant[informME_in$quant_score]
    informME_in$Sample=unique(GR_merge_in$Sample)
    return(informME_in)
  }
}
read.agnostic.mouse<-function(in_dir,tissue,stage,stat_type,replicate){
  file_in=paste(in_dir,'mm10_',tissue,'_',stage,'_merged',replicate,'_allele_agnostic_',stat_type,'.bedGraph',sep='')
  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    informME_in$tissue=tissue
    informME_in$stage=stage
    informME_in$bioreplicate=replicate
    informME_in$Sample=paste(tissue,stage,replicate,sep='-')
    
    return(informME_in)
  }
}
agnostic_matrix_conversion<-function(gr_in,stat='NME'){
  gr_out=granges(unique(gr_in))
  olap=findOverlaps(gr_in,gr_out,type='equal')
  stat_in_df=elementMetadata(gr_in[queryHits(olap)])[c(stat,'Sample')]
  stat_in_df$idx=NA
  stat_in_df$idx[queryHits(olap)]=subjectHits(olap)
  stat_in_df=as.data.frame(stat_in_df)
  
  stat_in_df_stat=dcast(data=stat_in_df,formula=idx~Sample,value.var = stat,fun.aggregate=mean)#remove agg.fun for new run
  
  for(sp in unique(colnames(stat_in_df_stat))[-1]){
    elementMetadata(gr_out)[[sp]]=NA
    elementMetadata(gr_out)[[sp]][stat_in_df_stat$idx]=stat_in_df_stat[,sp]
    
  }
  return(gr_out)
  
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
