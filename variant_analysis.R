# Genomics
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(VariantAnnotation)
library(reshape2)
# Source main functions
#setwd("~/code/HASM-MetaAnalysis/")
source("mainFunctions.R")
# All Directories
# args <- commandArgs(trailingOnly = TRUE)
# inDir <- "~/data/jabante1/ASM/Onuchic-2018/CpelAsmOutput/"
# vcfDir <- "/scratch/groups/afeinbe2/shared/JuliASM/data/vcfFiles/"
# outDir <- "/data/afeinbe2/jabante1/ASM/Onuchic-2018/VariantAnalysis/"
# sub <- args[1]
# # Load allele results
# cpelAllele <- resultsCpelAllele(inDir)
# #readRDS
# # Get subjecst and tissues
# subjectsTissues <- getSubjectsTissues(inDir)#a list of subject and Tissue
# #get all tissues from that subject
# whichTissues <- subjectsTissues$Subject == sub
# tissues <- subjectsTissues$Tissues[whichTissues]
# Loop over tissues

GRs=readRDS('../downstream/output/GRs_allele.rds')
outDir='../downstream/output/'
vcfDir='../downstream/data/vcfFiles/'
inDir='../downstream/data/gff_file/'
subjects=unique(GRs$Subject)
for (sub in subjects){
  varSub <- GRangesList()
  cpelAllele=GRs[GRs$Subject==sub]
  tissues=unique(cpelAllele$Tissue)
  cat(paste('Processing subject:',sub,'\n'))
  varSub=mclapply(tissues,function(x) crossVcfEtGffEtCpelAsm(vcfDir,inDir,cpelAllele,sub,x),mc.cores=6)
  # for(tissue in tissues){#do it for each tissue and sample combination
  #   varSub[[tissue]] <- crossVcfEtGffEtCpelAsm(vcfDir,inDir,cpelAllele,sub,tissue)
  # }
  
  # Save RDS
  saveRDS(varSub,file=paste(outDir,sub,"-variants_unique",".rds",sep=""))
  #cat(paste('Finish processing subject:',sub,'using:',proc.time()[3]-t1,'sec','\n'))
}
var_all=mclapply(subjects,function(x) extractHetCpG(vcfDir,x),mc.cores=6)
#make table
# sub='H9'
# tissue='Embryonic Stem Cell'
# cpelAllele=readRDS('../downstream/output/GRs_allele.rds')
# vcfDir='../downstream/data/vcfFiles/'
# inDir='../downstream/data/gff_file/'
# varSub <- crossVcfEtGffEtCpelAsm(vcfDir,inDir,cpelAllele,sub,tissue)

GR=readRDS('../downstream/output/GRs_no_brain.rds')
gr_allele=readRDS('../downstream/output/GRs_allele.rds')
RDS_dir='../downstream/output/Het_cpg/RDS_in/'
dat_check_NME=GRanges()
dat_check_MML=GRanges()
GR_sub_ASM_NME=GRanges()
GR_sub_ASM_MML=GRanges()
dat_check_NME_all=GRanges()
dat_check_MML_all=GRanges()
dat_all=GRanges()
for(var_dat in dir(RDS_dir)){
  sub=gsub('-variants_unique.rds','',var_dat)
  cat('Checking subject:',sub,'\n')
  dat=readRDS(paste(RDS_dir,var_dat,sep=''))
  #dat=GR_allele_correct_sub(dat,gr_allele,sub)
  GR_sub=GR[GR$Subject==sub]
  dat_check_NME=c(dat_check_NME,hetCpG_check_subject(dat,GR_sub,'NME',ASM=TRUE,sub))
  dat_check_MML=c(dat_check_MML,hetCpG_check_subject(dat,GR_sub,'MML',ASM=TRUE,sub))
  dat_check_NME_all=c(dat_check_NME_all,hetCpG_check_subject(dat,GR_sub,'NME',ASM=FALSE,sub))
  dat_check_MML_all=c(dat_check_MML_all,hetCpG_check_subject(dat,GR_sub,'MML',ASM=FALSE,sub))
  GR_sub_ASM_NME=c(GR_sub_ASM_NME,subsetByOverlaps(GR_sub[GR_sub$ASM=='Yes'],dat_check_NME))
  GR_sub_ASM_MML=c(GR_sub_ASM_MML,subsetByOverlaps(GR_sub[GR_sub$ASM=='Yes'],dat_check_MML))
  dat_all=c(dat_all,do.call('c',dat))
}
saveRDS(dat_check_NME,'../downstream/output/Het_cpg/ASM_NME_het_hom_difference_unique.rds')
saveRDS(dat_check_MML,'../downstream/output/Het_cpg/ASM_MML_het_hom_difference_unique.rds')
saveRDS(GR_sub_ASM_NME,'../downstream/output/Het_cpg/ASM_NME_GR_unique.rds')
saveRDS(GR_sub_ASM_MML,'../downstream/output/Het_cpg/ASM_MML_GR_unique.rds')
saveRDS(dat_all,'../downstream/output/Het_cpg/ASM_dat_unique.rds')

dat_check_NME=readRDS('../downstream/output/Het_cpg/ASM_NME_het_hom_difference_unique.rds')
dat_check_MML=readRDS('../downstream/output/Het_cpg/ASM_MML_het_hom_difference_unique.rds')
GR_sub_ASM_MML=readRDS('../downstream/output/Het_cpg/ASM_MML_GR_unique.rds')
GR_sub_ASM_NME=readRDS('../downstream/output/Het_cpg/ASM_NME_GR_unique.rds')

NME_allele=gr_allele_CpG[gr_allele_CpG$Statistic=='NME']
NME_allele=add_ASM(NME_allele,GR[GR$Statistic=='dNME'])
NME_allele_diff_all=allele_diff(NME_allele)
NME_allele_diff_ASM_calc=allele_diff(NME_allele[which(NME_allele$ASM=='Yes')])
NME_allele_diff_ASM=NME_allele_diff_ASM_calc[[1]]
NME_allele_diff_ASM$type='Non Het CpG'
NME_allele_diff_ASM$type[NME_allele_diff_ASM$CpGdiff!=0]='Het CpG'

#olap=findOverlaps(NME_allele_diff_ASM,dat_check_NME)
#NME_allele_diff_ASM_nohet=NME_allele_diff_ASM[-queryHits(olap)]

#NME_allele_diff_ASM_nohet$type='Non het CpG'
#dat_check_NME$type='Het CpG'

#Use density plot
#NME_df=as.data.frame(rbind(elementMetadata(NME_allele_diff_ASM_nohet)[,c('diff','type')],elementMetadata(dat_check_NME)[,c('diff','type')]))
NME_df=data.frame(NME_allele_diff_ASM$diff,NME_allele_diff_ASM$type)
colnames(NME_df)=c('Allele_difference','Allele_type')
ggplot(NME_df,aes(x=Allele_difference,fill=Allele_type))+
  geom_density(alpha=0.6)+xlab('Allele Difference')+
  ggtitle('NME Allele difference at different type of allele at ASM')

#More CpG genome2-genome1
NME_allele_gr_ASM=NME_allele_diff_ASM_calc[[2]]
NME_het_df=rbind(data.frame(NME=NME_allele_gr_ASM$Value[NME_allele_gr_ASM$CpGstat=='More'],type='more CpG'),
                 data.frame(NME=NME_allele_gr_ASM$Value[NME_allele_gr_ASM$CpGstat=='Less'],type='Less CpG'))


ggplot(NME_het_df,aes(x=NME,y=..scaled..,fill=type))+
  geom_density(alpha=0.6)+xlab('NME')+ggtitle('NME at different alleles at ASM')+
  theme(legend.position="bottom")


#Plot specific allele type
ggplot(NME_df[NME_df$Allele_type=='Het CpG',],aes(x=Allele_difference,fill=Allele_type))+
  geom_density(alpha=0.6)+xlab('Allele Difference')+theme(legend.position="none")+
  ggtitle('Allelic dNME difference at heterozygous CpG')

ggplot(NME_df[NME_df$Allele_type=='Non het CpG',],aes(x=Allele_difference,fill=Allele_type))+
  geom_density(alpha=0.6)+xlab('Allele Difference')+theme(legend.position="none")+
  ggtitle('Allelic dNME at non-heterozygous CpG')

#Only look at those allele with CpG vs nonCpG
#NME_het_df=rbind(data.frame(NME=dat_check_NME$cpg,type='CpG'),data.frame(NME=dat_check_NME$nonCpG,type='Non CpG'))


NME_het_df_all=rbind(data.frame(NME=dat_check_NME_all$cpg,type='CpG'),data.frame(NME=dat_check_NME_all$nonCpG,type='Non CpG'))
ggplot(NME_het_df_all,aes(x=NME,y=..scaled..,fill=type))+
  geom_density(alpha=0.6)+xlab('NME')+ggtitle('NME at different alleles at all loci')+
  theme(legend.position="bottom")

MML_allele=gr_allele[gr_allele$Statistic=='MML']
MML_allele=add_ASM(MML_allele,GR[GR$Statistic=='dMML'])
MML_allele_diff_ASM=allele_diff(MML_allele[which(MML_allele$ASM=='Yes')])
olap=findOverlaps(MML_allele_diff_ASM,dat_check_MML)
MML_allele_diff_ASM_nohet=MML_allele_diff_ASM[-queryHits(olap)]
MML_allele_diff_ASM_nohet$type='Non het CpG'
dat_check_MML$type='Het CpG'
MML_df=as.data.frame(rbind(elementMetadata(MML_allele_diff_ASM_nohet)[,c('diff','type')],elementMetadata(dat_check_MML)[,c('diff','type')]))
colnames(MML_df)=c('Allele_difference','Allele_type')
ggplot(MML_df,aes(x=Allele_difference,fill=Allele_type))+
  geom_density(alpha=0.6)+xlab('Allele Difference')+
  ggtitle('MML Allele difference at different type of allele at ASM')


#Plot specific allele type
ggplot(MML_df[MML_df$Allele_type=='Het CpG',],aes(x=Allele_difference,y=..scaled..,fill=Allele_type))+
  geom_density(alpha=0.6)+xlab('Allele Difference')+theme(legend.position="none")+
  ggtitle('Allele MML difference at heterozygous CpG')

ggplot(MML_df[MML_df$Allele_type=='Non het CpG',],aes(x=Allele_difference,y=..scaled..,fill=Allele_type))+
  geom_density(alpha=0.6)+xlab('Allele Difference')+theme(legend.position="none")+
  ggtitle('Allele MML difference at non-heterozygous CpG')


#Only look at those allele with CpG vs nonCpG
MML_het_df=rbind(data.frame(MML=dat_check_MML$cpg,type='CpG'),data.frame(MML=dat_check_MML$nonCpG,type='Non CpG'))
ggplot(MML_het_df,aes(x=MML,y=..scaled..,fill=type))+
  geom_density(alpha=0.6,)+xlab('MML')+ggtitle('MML at different alleles at ASM')+
  theme(legend.position="bottom")

MML_het_df_all=rbind(data.frame(MML=dat_check_MML_all$cpg,type='CpG'),data.frame(MML=dat_check_MML_all$nonCpG,type='Non CpG'))
ggplot(MML_het_df_all,aes(x=MML,y=..scaled..,fill=type))+
  geom_density(alpha=0.6,)+xlab('MML')+ggtitle('MML at different alleles at all loci')+
  theme(legend.position="bottom")

theme_bar=theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                legend.position="bottom")

NME_mean_diff=diff_mean_sample(dat_check_NME)
NME_mean_diff$subject=do.call('c',lapply(strsplit(as.character(NME_mean_diff$sample),split='-'),function(x) x[1]))
NME_mean_diff$tissue=do.call('c',lapply(strsplit(as.character(NME_mean_diff$sample),split='-'),function(x) x[2]))
ggplot(NME_mean_diff,aes(x=sample,y=diff,fill=subject))+geom_bar(stat='identity')+theme_bar+ylab('Difference')+
  ggtitle('NME difference between Non-CpG and CpG')
#ggplot(NME_mean_diff,aes(x=sample,y=diff,fill=tissue))+geom_bar(stat='identity')+theme_bar+ylab('NME difference at between non-CpG and CpG')

MML_mean_diff=diff_mean_sample(dat_check_MML)
MML_mean_diff$subject=do.call('c',lapply(strsplit(as.character(MML_mean_diff$sample),split='-'),function(x) x[1]))
ggplot(MML_mean_diff,aes(x=sample,y=diff,fill=subject))+geom_bar(stat='identity')+theme_bar+ylab('Difference')+
  ggtitle('MML difference between Non-CpG and CpG')


#Enrichment analysis
#All CpG analyzed extract from gff file
gff_in= readAllGff('../downstream/data/gff_file/')
gff_CpG=gff_to_CpG_loop(gff_in)
gff_CpG=resize(gff_CpG,2)
saveRDS(gff_CpG,'../downstream/input/gff_CpG.rds')
saveRDS(gff_in,'../downstream/input/gff_in.rds')
genomic_features <- readRDS("../downstream/input/genomic_features.rds")
#MML
MML=GR[GR$Statistic=='dMML']
MML_ASM=MML[MML$ASM=='Yes']
MML_sea_OR=calc_OR(genomic_features[['CpG open sea']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_island_OR=calc_OR(genomic_features[['CpG island']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_shore_OR=calc_OR(genomic_features[['CpG shore']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_shelf_OR=calc_OR(genomic_features[['CpG shelf']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_gene_OR=calc_OR(genomic_features[['gene body']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_exon_OR=calc_OR(genomic_features[['exon']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_intron_OR=calc_OR(genomic_features[['intron']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_promoter_OR=calc_OR(genomic_features[['promoter']],MML_ASM,GR_sub_ASM_MML,gff_CpG)
MML_intergenic_OR=calc_OR(genomic_features[['intergenic']],MML_ASM,GR_sub_ASM_MML,gff_CpG)

#NME
NME=GR[GR$Statistic=='dNME']
NME_ASM=NME[NME$ASM=='Yes']
NME_sea_OR=calc_OR(genomic_features[['CpG open sea']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_island_OR=calc_OR(genomic_features[['CpG island']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_shore_OR=calc_OR(genomic_features[['CpG shore']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_shelf_OR=calc_OR(genomic_features[['CpG shelf']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_gene_OR=calc_OR(genomic_features[['gene body']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_exon_OR=calc_OR(genomic_features[['exon']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_intron_OR=calc_OR(genomic_features[['intron']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_promoter_OR=calc_OR(genomic_features[['promoter']],NME_ASM,GR_sub_ASM_NME,gff_CpG)
NME_intergenic_OR=calc_OR(genomic_features[['intergenic']],NME_ASM,GR_sub_ASM_NME,gff_CpG)


#GO analysis
het_NME_gene=subsetByOverlaps(genomic_features[['promoter']],dat_check_NME)
het_NME_gene_all=subsetByOverlaps(genomic_features[['promoter']],NME_ASM)
het_NME_GO=GO_anno(het_NME_gene$gene_name,het_NME_gene_all$gene_name)
het_NME_GO_table=het_NME_GO[[1]]
het_NME_GO_table[het_NME_GO_table$Significant>=10&het_NME_GO_table$weightF<0.1,]

het_MML_gene=subsetByOverlaps(genomic_features[['promoter']],dat_check_MML)
het_MML_gene_all=subsetByOverlaps(genomic_features[['promoter']],MML_ASM)
het_MML_GO=GO_anno(het_MML_gene$gene_name,het_MML_gene_all$gene_name)
het_MML_GO_table=het_MML_GO[[1]]
het_MML_GO_table[het_MML_GO_table$weightF<0.2,]


table(scoresInTerm(het_NME_GO[[2]], 'GO:0030307',use.names=TRUE)) #use score =2
#CHIP-seq analysis


#Find number of het CpG in Gff regions
gff_in=readRDS('../downstream/input/gff_in.rds')
vcf_in=readRDS('../downstream/input/var_all.rds')
gr_allele=readRDS('../downstream/output/GRs_allele.rds')
gff_out=list()
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003","STL011" )
hetCpG_vcf=mclapply(subjects,gff_hetCpG_count,gff_in=gff_in,vcf_in=vcf_in,mc.cores = 7)
hetCpG_Vcf=readRDS('../downstream/input/hetCpG_vcf.rds')
gr_allele_CpG=GRanges()
for (subj in subjects){gr_allele_CpG=c(gr_allele_CpG,hetCGallele_sub(subj,gr_allele,hetCpG_Vcf))}
