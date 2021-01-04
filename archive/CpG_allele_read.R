# Genomics
rm(list=ls())
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(VariantAnnotation)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(Repitools)
# Source main functions
#setwd("~/code/HASM-MetaAnalysis/")
source("mainFunctions_sub.R")
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003","STL011",'GM12878','H1')
#Read in gff files
cat('Processing gff file\n')
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
saveRDS(gff_in,'../downstream/input/gff_in_new.rds')
#Extract Heterozygous CpGs from vcf files
cat('processing variant file\n')
variant_HetCpG=lapply(subjects,function(x) extractHetCpG('../downstream/data/vcfFiles/',x))
names(variant_HetCpG)=subjects
saveRDS(variant_HetCpG,'../downstream/input/variant_HetCpG_new.rds')
#For each region in the gff, Cont number of Het CPG in it
cat('Extracting Het CpG\n')
hetCpG_gff=lapply(subjects,gff_hetCpG_count,gff_in=gff_in,vcf_in=variant_HetCpG,CpG=CpG_hg19)
names(hetCpG_gff)=subjects
saveRDS(hetCpG_gff,'../downstream/input/hetCpG_gff_new.rds')
#Read in CpG location and read in data for each allele
gr_allele=readRDS('../downstream/output/GR.all.allele.H1.GM12878.rds')

#For each region in output, add heterogyzous CpG number in it
cat('Processing allelic file')
gr_allele_CpG=GRanges()
for (subj in subjects){gr_allele_CpG=c(gr_allele_CpG,hetCGallele_sub(subj,gr_allele,hetCpG_gff,CpG_hg19,variant_HetCpG))}

saveRDS(gr_allele_CpG,'../downstream/output/gr_allele_CpG_new.rds')

gr_allele_CpG_resize=GRanges()
for (subj in subjects){gr_allele_CpG_resize=c(gr_allele_CpG_resize,GR_resize_sub(subj,gr_allele_CpG,CpG_hg19,variant_HetCpG))}


saveRDS(gr_allele_CpG_resize,'../downstream/output/gr_allele_CpG_new2.rds')
