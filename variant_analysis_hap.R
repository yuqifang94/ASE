#Dependencies
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
# Genomics
# Source main functions
#setwd("~/code/HASM-MetaAnalysis/")
source("mainFunctions_sub.R")
#global pval cutoff,remove GM12878 in summary analysis
pval_cutoff=0.1 #Onuchic use 0.1
##############################due to the output of CpelASM, it would be better to use gff ranges and assign all GR parameters to gff###############
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
gff_in=readRDS('../downstream/input/gff_in_new.rds')
#Note the first line cannot be NA
GR=import.subject('../downstream/data/Julia_bedfile/clean_file/',gff_in)
GR_allele=import.subject('../downstream/data/bedGraph_allele/',gff_in,calc='allele')
##############################Variant analysis: Counting enrichment of each variant from all variant ##############################################
variant_HetCpG_meta=do.call('c',lapply(names(variant_HetCpG),variant_meta_subj,variant_in=variant_HetCpG,GR_in=GR))
saveRDS(variant_HetCpG_meta,'../downstream/output/ASM_enrich_meta.rds')
variant_HetCpG_meta=readRDS('../downstream/output/ASM_enrich_meta.rds')
variant_HetCpG_meta$ASM=NA
variant_HetCpG_meta$ASM[variant_HetCpG_meta$pvalue<=pval_cutoff]='Yes'
variant_HetCpG_meta$ASM[variant_HetCpG_meta$pvalue>pval_cutoff]='No'
variant_HetCpG_meta=variant_HetCpG_meta[!is.na(variant_HetCpG_meta$ASM)]
variant_HetCpG_meta=variant_HetCpG_meta[variant_HetCpG_meta$Sample %in% c( unique(variant_HetCpG_meta$Sample)[1:43])]
ASM_het_enrichment(variant_HetCpG_meta)
ASM_enrich_all=list()

for (st in unique(variant_HetCpG_meta$Statistic)){
  ASM_enrich=data.frame(sp=NULL,OR=NULL)
  for (sp in unique(variant_HetCpG_meta$Sample)){
    ASM_enrich=rbind(ASM_enrich,data.frame(sp=sp,subjects=strsplit(sp,' - ')[[1]][2],
                                           OR=ASM_het_enrichment(variant_HetCpG_meta[variant_HetCpG_meta$Sample==sp & variant_HetCpG_meta$Statistic==st])$estimate))
  }
  ASM_enrich_all[[st]]=ASM_enrich
}

#Box plot of all types of statistics
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
ggplot(ASM_enrich_all$dMML,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 8)+
  ggtitle('dMML Het CpG enrichment')+xlab('Sample name')+ylab('Odds Ratio')+theme_bar

ggplot(ASM_enrich_all$dNME,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 8)+
  ggtitle('dNME Het CpG enrichment')+xlab('Sample name')+ylab('Odds Ratio')+theme_bar

ggplot(ASM_enrich_all$UC,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 8) +
  ggtitle('UC Het CpG enrichment')+xlab('Sample name')+ylab('Odds Ratio')+theme_bar
##############################Het Cpg Analysis: Calculating density and het CpG information #######################################################
########Find number of het CpG in Gff regions
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003","STL011","GM12878","H1")
subjects=c("GM12878","H1")
subjects='STL011'
#Read in gff files
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
gff_in=readRDS('../downstream/input/gff_in_new.rds')
#Extract Heterozygous CpGs from vcf files
variant_HetCpG=lapply(subjects,function(x) extractHetCpG('../downstream/vcfFiles/',x))
variant_HetCpG=readRDS('../downstream/input/variant_HetCpG_new.rds')
names(variant_HetCpG)=subjects
#For each region in the gff, Cont number of Het CPG in it
hetCpG_gff=lapply(subjects,gff_hetCpG_count,gff_in=gff_in,vcf_in=variant_HetCpG,CpG=CpG_hg19)
names(hetCpG_gff)=subjects
hetCpG_gff=readRDS('../downstream/input/hetCpG_gff.rds')
#Read in CpG location and read in data for each allele
gr_allele=readRDS('../downstream/output/GRs_allele.rds')
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
#For each region in output, add heterogyzous CpG number in it
gr_allele_CpG=GRanges()
for (subj in subjects){gr_allele_CpG=c(gr_allele_CpG,hetCGallele_sub(subj,gr_allele,hetCpG_gff,CpG_hg19,variant_HetCpG))}
#Resize and calculate number of CGs
gr_allele_CpG_resize=GRanges()
for (subj in subjects){gr_allele_CpG_resize=c(gr_allele_CpG_resize,GR_resize_sub(subj,gr_allele_CpG,CpG_hg19,variant_HetCpG))}
#This one should not have density in it
saveRDS(gr_allele_CpG,'../downstream/input/gr_allele_CpG_500.rds')

#Reading in GR
GR=readRDS('../downstream/output/GR.all.diff.H1.GM12878.rds')#Use ASM cutoff=0.05
gr_allele_CpG=readRDS('../downstream/output/gr_allele_CpG_new2.rds')#dropping one CpG region in future
gr_allele_CpG$TpA_CpG=gr_allele_CpG$TpA_count_extend/gr_allele_CpG$CG_hg19_extend
#Find ways to calculate CpG density, here use number of CG/size of the region
#Plot density distribution
#hg19 CG density
cpg_hg19_density=getCpgdensH19()
hist(log10(unlist(cpg_hg19_density[[2]])),xlab='log10(distance)',main='Distance between CpG site across hg19')
##############################Het Cpg Analysis: Het CpG distribution vs non Het CpG #######################################################

############################## NME analysis #######################################################
NME_allele=gr_allele_CpG[gr_allele_CpG$Statistic=='NME']
#Assign ASM information and find the ASM region
NME_allele=add_ASM(NME_allele,GR[GR$Statistic=='dNME'])
NME_allele=NME_allele[width(NME_allele)>2]
NME_allele_ASM_calcvariant_HetCpG=readRDS('../downstream/input/variant_HetCpG_new.rds')
#Subset by SNP-containing haplotype?
NME_allele_SNP=SNP_conmtaining_hap(NME_allele,variant_HetCpG)
NME_allele_ASM=NME_allele[which(NME_allele$ASM=='Yes')]
NME_allele_ASM=NME_allele_ASM[width(NME_allele_ASM)>2]
#Looking for Het CpG effect over ASM different regions 
NME_allele_ASM_calc=allele_calc_plot(NME_allele_ASM,'NME','../downstream/output/dNME/','ASM Region compare SNP containing')
NME_allele_ASM_calc=allele_diff(NME_allele[which(NME_allele$ASM=='Yes')])
NME_allele_calc=allele_diff(NME_allele)
allele_plot(NME_allele_ASM_calc,'NME','../downstream/output/dNME/','ASM Region')
NME_allele_all_calc=allele_calc_plot(NME_allele,'NME','../downstream/output/dNME/','all')
saveRDS(NME_allele_ASM_calc,'../downstream/output/dNME/NME_ASM_het_calc.rds')
saveRDS(NME_allele_all_calc,'../downstream/output/dNME/NME_all_het_calc.rds')
#Looking for Het CpG effect over differnet size of regions or number of CpG
#Generate bed files for genome browser
outdir='../downstream/output/dNME/'
NME_ASM=NME_allele_ASM_calc[[1]]
NME_all=NME_allele_calc[[1]]
NME_ASM=readRDS('../downstream/output/dNME/NME_ASM_het_calc.rds')
NME_ASM=NME_ASM[[1]]
NME_ASM_m1=NME_ASM[width(NME_ASM)>2]
plot(ecdf(abs(NME_ASM$diff[NME_ASM$CpGdiff!=0])),main='cumulative distribution of dNME at ASM',xlab='dNME',xlim=c(0.3,1),lwd=3)
lines(ecdf(abs(NME_ASM$diff[NME_ASM$CpGdiff==0])),col='red',lwd=3)
plot(ecdf(abs(NME_all$diff[NME_all$CpGdiff!=0])),main='cumulative distribution of dNME',xlab='dNME',xlim=c(0,1),lwd=3)
lines(ecdf(abs(NME_all$diff[NME_all$CpGdiff==0])),col='red',lwd=3)

for (sp in unique(NME_ASM_het$Sample)){ASM_bed_gen_sp(NME_ASM_het,gr_allele_CpG,sp,'../downstream/output/dNME/bedallele/')}
##############################Science paper: find average NME at 2 alleles at dMML region #################################################
NME_all=NME_allele_calc[[2]]
NME_all=NME_all[width(NME_all)>2]
NME_mean_dMML=c()
#Looking at dMML ASM
for (sp in unique(GR$Sample)){NME_mean_dMML=c(NME_mean_dMML,
                                         subsetByOverlaps(NME_all[NME_all$Sample==sp],GR[GR$Sample==sp & GR$Statistic=='dMML' &GR$pvalue<=pval_cutoff])$Value)
}
NME_ASM=data.frame(NME_mean=NME_mean,NME_type='dMML ASM',line_type='solid',stringsAsFactors = F)
#Looking at imprinted region dMML >0.9
NME_mean_imp=c()
for (sp in unique(GR$Sample)){NME_mean_imp=c(NME_mean_imp,subsetByOverlaps(NME_all[NME_all$Sample==sp],
                                                           GR[GR$Sample==sp & GR$Statistic=='dMML' &GR$pvalue<=pval_cutoff & GR$Value>=0.9])$Value)
}
NME_ASM_imp=data.frame(NME_mean=NME_mean_imp,NME_type='dMML ASM bimodal',line_type='solid',stringsAsFactors = F)
NME_mean_non_ASM=c()
#Looking at dMML nonASM
for (sp in unique(GR$Sample)){NME_mean_non_ASM=c(NME_mean_non_ASM,subsetByOverlaps(NME_all[NME_all$Sample==sp],
                                                                                   GR[GR$Sample==sp & GR$Statistic=='dMML' &GR$pvalue>pval_cutoff])$Value)}
NME_non_ASM=data.frame(NME_mean=NME_mean_non_ASM,NME_type='Non dMML ASM',line_type='dashed',stringsAsFactors = F)
dNME_ASM=data.frame(NME_mean=NME_allele_ASM_calc[[2]]$Value,NME_type='dNME ASM',line_type='soild')
NME_mean_df=rbind(NME_ASM,NME_ASM_imp,NME_non_ASM,dNME_ASM)

ggplot(NME_mean_df,aes(x=NME_mean,color=NME_type))+
  geom_density(aes(linetype=line_type),alpha=0.6,size=1)+xlab('Mean NME')+ggtitle('Average NME at 2 alleles')+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+ylim(0,4)+scale_color_manual(values=c("blue","red","black",'purple'))+
  geom_hline(yintercept=0, colour="white", size=1)


ggplot(NME_mean_df,aes(x=NME_type,y=NME_mean,fill=NME_type))+
  geom_boxplot()+xlab('Type of ASM')+ggtitle('NME at 2 alleles')+ylab('NME')+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+ylim(0,1)
############################## MML analysis #######################################################
MML_allele=gr_allele_CpG[gr_allele_CpG$Statistic=='MML']
MML_allele=add_ASM(MML_allele,GR[GR$Statistic=='dMML'])
MML_allele_ASM=MML_allele[which(MML_allele$ASM=='Yes')]
MML_allele_calc=allele_diff(MML_allele)
MML_allele_ASM_calc=allele_calc_plot(MML_allele_ASM,'MML','../downstream/output/dMML/','ASM Region')
allele_plot(MML_allele_ASM_calc,'MML','../downstream/output/dMML/','ASM Region')
MML_allele_all_calc=allele_calc_plot(MML_allele,'MML','../downstream/output/dMML/','all')
saveRDS(MML_allele_ASM_calc,'../downstream/output/dMML/MML_ASM_het_calc.rds')
#Generate bed files for genome browser
outdir='../downstream/output/dMML/'
MML_ASM=MML_allele_ASM_calc[[1]]
for (sp in unique(MML_ASM_het$Sample)){ASM_bed_gen_sp(MML_ASM_het,gr_allele_CpG,sp,'../downstream/output/dMML/bedallele/')}
#####Find the hist gram of MML
MML_allele_ASM=MML_allele[which(MML_allele$pval<=pval_cutoff)]
hist(MML_allele_ASM$Value,xlab='Mean methylation level',main='Mean methylation value at each allele at ASM')

##############################Het Cpg Analysis: Density effect #######################################################
############################## NME analysis #######################################################
#CG count enrichment
#Using all NME calculation, define ASM vs non-ASM, het CG vs non-het CG 
NME_all=NME_allele_calc[[1]] #Granges for each allele
NME_all=NME_all[!is.na(NME_all$pval)]
NME_all$CG_type=NA
NME_all$CG_type[NME_all$CGcount_diff!=0] ='Imbalanced CG'
NME_all$CG_type[NME_all$CGcount_diff==0] ='Balanced CG'
NME_all$ASM=NA
NME_all$ASM[NME_all$pval<=pval_cutoff] ='Yes'
NME_all$ASM[NME_all$pval>pval_cutoff] ='No'
NME_ASM=NME_all[NME_all$ASM=='Yes']
NME_het=NME_all[NME_all$CG_type=='Imbalanced CG']
NME_all$HetCpG=NME_all$CGcount_diff!=0
NME_ASM_het=NME_het[which(NME_het$ASM=='Yes')]
subject_old=unique(NME_all$Sample)[1:43]
subject_new=unique(NME_all$Sample)[44:49]
NME_all=NME_all[NME_all$Sample%in%c(subject_old,'merged - H1','merged - GM12878')]
#Subset for SNP-containing ranges
NME_hap=SNP_conmtaining_hap(NME_all,variant_HetCpG)
NME_SNP=NME_hap[[1]]
NME_non_SNP=NME_hap[[2]]
#Enrichement test for SNP-overlapping region
ASM_het_enrichment(NME_SNP)
ASM_het_enrich(NME_SNP,'dNME HetCpG enrichment for SNP overlapping haplotypes')

#Enrichement test for non-SNP containing het CpG
ASM_het_enrichment(NME_non_SNP)
ASM_het_enrich(NME_non_SNP,'dNME HetCpG enrichment for non SNP overlapping haplotypes')

#Enrichement test for all het CpG
ASM_het_enrichment(NME_all)
ASM_het_enrich(NME_all,'dNME HetCpG enrichment for all')

#####################################Global observation on CpG density etc #################################
#Put something into NME all
#CpG content distribution
NME_all$CG_cont_type[NME_all$ASM=='No' & NME_all$CGcount_diff!=0]='Imbalanced CG Non ASM'
NME_all$CG_cont_type[NME_all$ASM=='No' & NME_all$CGcount_diff==0]='Balanced CG Non ASM'
NME_all$CG_cont_type[NME_all$ASM=='Yes' & NME_all$CGcount_diff!=0]='Imbalanced CG ASM'
NME_all$CG_cont_type[NME_all$ASM=='Yes' & NME_all$CGcount_diff==0]='Balanced CG ASM'
#Choose from CG content, CG density, TpA/CpG, TA/CG
plot_density(data.frame(density=NME_all$TpA_count_extend/NME_all$CpG_count_extend,CG_type=NME_all$CG_cont_type),
             ylab='TpA count/CpG count',ylim=c(0,30),title='')
plot_density(data.frame(density=NME_all$TpA_count_extend,CG_type=NME_all$CG_cont_type),
             ylab='TpA count',ylim=c(0,100),title='')
plot_density(data.frame(density=NME_all$CGcount_hg19_extend,CG_type=NME_all$CG_cont_type),
             ylab='Oberseved/Expected CpG',ylim=c(0,0.6),title='')
plot_density(data.frame(density=NME_all$CpG_count_extend/NME_all$gff_size_extend,CG_type=NME_all$CG_cont_type),
             ylab='CpG density',ylim=c(0,0.03),title='')
plot_density(data.frame(density=NME_all$AT_count_extend,CG_type=NME_all$CG_cont_type),
             ylab='AT count',ylim=c(0,500),title='')
plot_density(data.frame(density=NME_all$AT_count_extend/NME_all$CG_count_extend,CG_type=NME_all$CG_cont_type),
            ylab='AT/CG count',ylim=c(0,0.3),title='')

#It looks ASM have higher AT/CG count redo analysis
plot_density(data.frame(density=NME_all$TpA_count_extend/NME_all$CpG_count_extend,CG_type=NME_all$ASM),
             ylab='TpA count/CpG count',ylim=c(0,30),title='',xlab='ASM')
plot_density(NME_CG_df=data.frame(density=NME_all$AT_count_extend/NME_all$CG_count_extend,CG_type=NME_all$ASM),
             ylab='AT/CG count',ylim=c(0,0.3),title='',xlab='ASM')
plot_density(data.frame(density=NME_all$CGcount_hg19_extend,CG_type=NME_all$ASM),
             ylab='Oberseved/Expected CpG',ylim=c(0,0.6),title='',xlab='ASM')
#Checking those regions that have no CpG number difference but with Het CpG

cor(NME_ASM_het$CGcount_diff/NME_ASM_het$CGcount_hg19_extend,NME_ASM_het$diff)
cor(NME_ASM_het$CpGdiff/width(NME_ASM_het),NME_ASM_het$diff)
cor(NME_ASM_het$CGcount_diff,NME_ASM_het$diff)

#########################################GO analysis#################################
#reference: within 5k of SNP
names(variant_HetCpG)=NULL
GO_ref=subsetByOverlaps(genomic_features$promoter,do.call('c',variant_HetCpG),maxgap = 2000)
write(unlist(GO_ref$gene_name),'../downstream/output/dNME/GO_ref_SNP.txt')
NME_ASM_het_sub=NME_ASM_het[NME_ASM_het$Sample %in% unique(NME_ASM_het$Sample)[1:43]]
NME_ASM_het_sub_promoter=subsetByOverlaps(genomic_features$promoter,NME_ASM_het_sub,maxgap = 1000)
write(unlist(NME_ASM_het_sub_promoter$gene_name),'../downstream/output/dNME/GO_het_ASM_all.txt')
#ESC
NME_ASM_het_sub=NME_ASM_het[NME_ASM_het$Sample %in% unique(NME_ASM_het$Sample)[1:5]]
NME_ASM_het_sub_promoter=subsetByOverlaps(genome_features$promoter,NME_ASM_het_sub,maxgap = 500)
write(unlist(NME_ASM_het_sub_promoter$gene_name),'../downstream/output/dNME/GO_het_ASM_ESC.txt')
#Differentiated
NME_ASM_het_sub=NME_ASM_het[NME_ASM_het$Sample %in% unique(NME_ASM_het$Sample)[6:43]]
NME_ASM_het_sub_promoter=subsetByOverlaps(genome_features$promoter,NME_ASM_het_sub,maxgap = 500)
write(unlist(NME_ASM_het_sub_promoter$gene_name),'../downstream/output/dNME/GO_het_ASM_diff.txt')

#ASM only
NME_ASM_sub=NME_ASM[NME_ASM$Sample %in% unique(NME_ASM$Sample)[1:43]]
NME_ASM_sub_promoter=subsetByOverlaps(genome_features$promoter,NME_ASM_sub,maxgap = 5000)
write(unlist(NME_ASM_sub_promoter$gene_name),'../downstream/output/dNME/GO_ASM_all_5k.txt')
###########################Plot enrichment at each genomic feature####################
#Islands
genome_feature_plot(NME_all,genomic_features$`CpG island`,'dNME','CpG Island enrichment')
#Open seas
genome_feature_plot(NME_SNP,genomic_features$`CpG open sea`,'dNME','CpG open sea enrichment')

#repetitive element
#Without het CpG
rep=import.bed('../downstream/input/simple_repeats.bed')
OR_repeats=data.frame(subject=subjects,OR=0,CpG_type='simple_repeats')
for(subj in subjects){
  OR_repeats$OR[OR_repeats$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Subject==subj],rep,'dNME')$estimate
} 
############################## MML analysis #######################################################
#CG count enrichment
#Using all MML calculation, define ASM vs non-ASM, het CG vs non-het CG 
MML_all=MML_allele_calc[[1]] #Granges for each allele
MML_all=MML_all[!is.na(MML_all$pval)]
MML_all$CG_type=NA
MML_all$CG_type[MML_all$CGcount_diff!=0] ='Imbalanced CG'
MML_all$CG_type[MML_all$CGcount_diff==0] ='Balanced CG'
MML_all$ASM=NA
MML_all$ASM[MML_all$pval<=pval_cutoff] ='Yes'
MML_all$ASM[MML_all$pval>pval_cutoff] ='No'
MML_ASM=MML_all[MML_all$ASM=='Yes']
MML_het=MML_all[MML_all$CG_type=='Imbalanced CG']
MML_all$HetCpG=MML_all$CGcount_diff!=0
MML_ASM_het=MML_het[which(MML_het$ASM=='Yes')]
subject_old=unique(MML_all$Sample)[1:43]
subject_new=unique(MML_all$Sample)[44:49]
MML_all=MML_all[MML_all$Sample%in%c(subject_old,'merged - H1','merged - GM12878')]
#Subset for SNP-containing ranges
MML_hap=SNP_conmtaining_hap(MML_all,variant_HetCpG)
MML_SNP=MML_hap[[1]]
MML_non_SNP=MML_hap[[2]]
#Enrichement test for SNP-overlapping region: Het CpG only enriched in dNME, likely MML difference are caused by dNME
ASM_het_enrichment(MML_SNP)
ASM_het_enrich(MML_SNP,'dMML HetCpG enrichment for SNP overlapping haplotypes')

#Enrichement test for non-SNP containing het CpG
ASM_het_enrichment(MML_non_SNP)
ASM_het_enrich(MML_non_SNP,'dMML HetCpG enrichment for non SNP overlapping haplotypes')

#Enrichement test for all het CpG
ASM_het_enrichment(MML_all)
ASM_het_enrich(MML_all,'dMML HetCpG enrichment for all')

###########################Plot enrichment at each genomic feature####################
genomic_features=readRDS('../downstream/input/genomic_features.rds')
#Islands
genome_feature_plot(MML_all,genomic_features$`CpG island`,'dMML','CpG Island enrichment',ylim=c(0,10))
genome_feature_plot(NME_all,genomic_features$`CpG island`,'dNME','CpG Island enrichment',ylim=c(0,10))
#Open seas
genome_feature_plot(MML_all,genomic_features$`CpG open sea`,'dMML','CpG open sea enrichment',ylim=c(0,10))
#Promoters
genome_feature_plot(MML_all,genomic_features$promoter,'dMML','dMML Promoter enrichment',ylim=c(0,6))
genome_feature_plot(NME_all,genomic_features$promoter,'dNME','dNME Promoter enrichment',ylim=c(0,6))

#enhancer
genome_feature_plot(MML_all,genomic_features$enhancer,'dMML','enhancer enrichment dMML',ylim=c(0,2))
genome_feature_plot(NME_all,genomic_features$enhancer,'dNME','enhancer enrichment dNME',ylim=c(0,2))
############################## NME analysis #######################################################


#filtering CpG islands bed
CpG_islands=read.table('../downstream/input/CpG_islands_hg19.bed',header=TRUE)
CpG_islands=makeGRangesFromDataFrame(CpG_islands,keep.extra.columns = TRUE)
genomic_features=readRDS('../downstream/input/genomic_features.rds')
promoter_island_overlap=findOverlaps(genomic_features$TSS,genomic_features$`CpG island`)
promoter_with_CpG_islands=genomic_features$promoter[queryHits(promoter_island_overlap)]
promoter_without_CpG_islands=genomic_features$promoter[-queryHits(promoter_island_overlap)]
#Find distnace to CpG island 
gr_distance(NME_all[NME_all$pval<=pval_cutoff],genomic_features$`CpG island`,'distance to CpG island (bp)','dNME ASM distance to CpG island',ylim=c(0,0.3))
gr_distance(MML_all[MML_all$pval<=pval_cutoff],genomic_features$`CpG island`,'distance to CpG island (bp)','dMML ASM distance to CpG island',ylim=c(0,0.3))
#Promoters
gr_distance(MML_all[MML_all$pval<=pval_cutoff],promoter_with_CpG_islands,'distance to TSS','dMML ASM distance to promoter with CpG islands',ylim=c(0,0.1))
gr_distance(MML_all[MML_all$pval<=pval_cutoff],promoter_without_CpG_islands,'distance to TSS','dMML ASM distance to promoter without CpG islands',ylim=c(0,0.1))

gr_distance(NME_all[NME_all$pval<=pval_cutoff],promoter_with_CpG_islands,'distance to TSS','dNME ASM distance to promoter with CpG islands',ylim=c(0,0.22))
gr_distance(NME_all[NME_all$pval<=pval_cutoff],promoter_without_CpG_islands,'distance to TSS','dNME ASM distance to promoter without CpG islands',ylim=c(0,0.22))
#Enhancers
enhancers <- readRDS("E:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/enhancers.rds") #474,004 enhancer region
gr_distance(MML_all[MML_all$pval<=pval_cutoff],enhancers,'distance to enhancer','dMML ASM distance to enhancer',ylim=c(0,0.25))
gr_distance(NME_all[NME_all$pval<=pval_cutoff],enhancers,'distance to TSS','dNME ASM distance to enhancer',ylim=c(0,0.25))


###################Find overlap between dMML, dNME and UC regions ##########################################################################
GR=readRDS('../downstream/output/GR.all.diff.H1.GM12878.rds')
##################Put 3 stat into same Granges object########################################################################################
GR_merge=GRanges()
for(sp in unique(GR$Sample)){GR_merge=c(GR_merge,stat_merge(GR[GR$Sample==sp]))}
GR_merge_exclude_GM=GR_merge[!GR_merge$Sample%in%c('merged - GM12878','1 - GM12878','2 - GM12878')]
GR_merge_exclude_GM=GR_merge_exclude_GM[width(GR_merge_exclude_GM)>10]
#GR dNME vs dMML at ASM
GR_merge_exclude_GM_sig=GR_merge_exclude_GM[GR_merge_exclude_GM$dMML_pval<=pval_cutoff |GR_merge_exclude_GM$dNME_pval<=pval_cutoff]
GR_merge_exclude_GM_sig=GR_merge_exclude_GM[!is.na(GR_merge_exclude_GM$dMML_pval) &!is.na(GR_merge_exclude_GM$dNME_pval)]
GR_merge_df=data.frame(dMML=GR_merge_exclude_GM_sig$dMML,dNME=GR_merge_exclude_GM_sig$dNME)
#GR_merge_df_agg=aggregate(GR_merge_df$dNME,by=list(GR_merge_df$dMML),FUN=mean)
#GR_merge_df_agg$sd=aggregate(GR_merge_df$dNME,by=list(GR_merge_df$dMML),FUN=sd)$x
names(GR_merge_df_agg)=c("dMML","dNME","sd")
ggplot(GR_merge_df, aes(x=dMML, y=dNME)) + geom_smooth(method = "loess", formula = y ~x, size = 1)+
  xlim(c(0,1))+ylim(c(0,0.5))+ggtitle("dMML and dNME relationship at least one ASM event region")
  
#####Heatmap trial###########
GR_merge_df_heat=data.frame(dMML=round(GR_merge_exclude_GM_sig$dMML*2,digits=1)/2,dNME=round(GR_merge_exclude_GM_sig$dNME*2,digits=1)/2)
GR_merge_agg_heat=as.data.frame(aggregate(GR_merge_df_heat,by=list(dMML=GR_merge_df_heat$dMML,dNME=GR_merge_df_heat$dNME),FUN=length))
names(GR_merge_agg_heat)=c("dMML","dNME","count1","count_normal")
GR_merge_agg_heat$dMML_total=NA
for(i in unique(GR_merge_agg_heat$dMML)){
  GR_merge_agg_heat$dMML_total[which(GR_merge_agg_heat$dMML==i)]=sum(GR_merge_agg_heat$count1[which(GR_merge_agg_heat$dMML==i)])
}
GR_merge_agg_heat$count_nromal=log(GR_merge_agg_heat$count1)/log(GR_merge_agg_heat$dMML_total)
ggplot(GR_merge_agg_heat , aes(x = dMML, y = dNME)) +
  geom_raster(aes(fill = count_nromal), interpolate=T) +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(0,1))
###########################
GR_ASM=GR[GR$pvalue<=pval_cutoff]
olap_sample=data.frame()
for(sp in unique(GR_ASM$Sample)){olap_sample=rbind(olap_sample,olap_ASM(GR_ASM[GR_ASM$Sample==sp]))}
summary_olap=colSums(olap_sample[1:43,1:10])
#dMML_dNME      dMML_UC      dNME_UC        olap3         dMML         dNME           UC dMML_nonolap dNME_nonolap   UC_nonolap 
#494         3527          494          517         9621       322840         8146         4049       320301         2574 
#############################################SNP type enrichment in ASM###############################################
variant_HetCpG_meta=readRDS('../downstream/output/ASM_enrich_meta.rds')
variant_HetCpG_meta$variants=variants_collapase(variant_HetCpG_meta)
seleted_sample= unique(variant_HetCpG_meta$Subject)[1:9]
variant_HetCpG_meta=variant_HetCpG_meta[variant_HetCpG_meta$Subject %in% seleted_sample]
#Variant frequency
variant_freq=table(variant_HetCpG_meta$variants)
variant_freq=as.data.frame(variant_freq/sum(variant_freq))
colnames(variant_freq)=c('SNP','frequency')
#Plot SNP frequency
title='SNP frequency'
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
ggplot(variant_freq,aes(x=SNP,y=frequency,fill=SNP)) +ylim(c(0,0.5))+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+theme_bar

#calculate OR for each type of SNP
OR_df_out=data.frame()
for (stat_in in unique(variant_HetCpG_meta$Statistic)){
  OR_output=list()
  for (SNP_type in unique(variant_HetCpG_meta$variants)){OR_output[[SNP_type]]=variants_OR(variant_HetCpG_meta[variant_HetCpG_meta$Statistic==stat_in],SNP_type)}
  OR_df=do.call(rbind,lapply(OR_output,function(x) data.frame(lower_CI=x$conf.int[1],upper_CI=x$conf.int[2],OR=x$estimate)))
  OR_df$stat=stat_in
  OR_df$SNP=names(OR_output)
  OR_df_out=rbind(OR_df_out,OR_df)
}
SNP_fill=data.frame(SNP=sort(unique(OR_df_out$SNP)),color=c('red','coral','darkturquoise','orchid3','dodgerblue3','seagreen3'),stringsAsFactors = F)
OR_df_out$SNP_fill=SNP_fill$color[match(OR_df_out$SNP,SNP_fill$SNP)]
title='SNP enrichement'
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom",plot.title = element_text(hjust = 0.5))
ggplot(OR_df_out,aes(x=stat,y=OR,fill=SNP)) +ylim(c(0,2))+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+theme_bar+scale_fill_manual(values=SNP_fill$color)
#dNME SNP
title='dNME SNP enrichement'
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom",plot.title = element_text(hjust = 0.5))
ggplot(OR_df_out[OR_df_out$stat=='dNME',],aes(x=SNP,y=OR,fill=SNP)) +ylim(c(0,2))+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+theme_bar+scale_fill_manual(values=SNP_fill$color)
#Trinucleotide analysis
#get masked trinucleotide
mask_tri<-function(x){
  x_sp=strsplit(x,'')[[1]]
  x_sp[2]='X'
  return(paste(x_sp,collapse = ''))
}

variant_HetCpG_meta$mask_tri=unlist(lapply(variant_HetCpG_meta$REF_tri,mask_tri))
#Plot each variant with 6 tri nucleotide
#one trinucleotide


tri_OR_out_SNP_all=list()
for (stat_type in unique(variant_HetCpG_masked$Statistic)){
  tri_OR_out_SNP=list()
  for(vari in unique(variant_HetCpG_masked$variants)){
    tri_OR_out=list()
    for(tri_mask in unique(variant_HetCpG_masked$mask_tri)){
    tri_OR_out[[tri_mask]]=tri_nucleo_OR(variant_HetCpG_masked
                                                 [variant_HetCpG_masked$Statistic==stat_type&variant_HetCpG_masked$variants==vari],tri_mask)
    }
    tri_OR_df=do.call(rbind,lapply(tri_OR_out,function(x) data.frame(lower_CI=x$conf.int[1],upper_CI=x$conf.int[2],OR=x$estimate)))
    tri_OR_df$tri=rownames(tri_OR_df)
    tri_OR_df$SNP=vari
    tri_OR_df=tri_OR_df[order(tri_OR_df$OR,decreasing = T),]
    tri_OR_out_SNP[[vari]]=tri_OR_df
  }
  tri_OR_out_SNP_all[[stat_type]]=tri_OR_out_SNP
}
saveRDS(tri_OR_out_SNP_all,'../downstream/output/tri_OR_out_SNP_all.rds')
dNME_tri_OR_all=tri_OR_out_SNP_all$dNME
dNME_tri_plot_ls=list()
for(vari in names(dNME_tri_OR_all)){
title=vari
  
  theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="",plot.title = element_text(hjust = 0.5))
  dNME_tri_plot_ls[[vari]]=ggplot(dNME_tri_OR_all[[vari]],aes(x=reorder(tri,-OR),y=OR,fill=SNP)) +ylim(c(0,4))+
    geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
    geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                  position=position_dodge(.9))+theme_bar+xlab('')+scale_fill_manual(values=SNP_fill$color[SNP_fill$SNP==vari])
}
dNME_tri_plot_ls=dNME_tri_plot_ls[order(names(dNME_tri_plot_ls))]
grid.arrange(grobs=dNME_tri_plot_ls, nrow = 2)
###############################E14 analysis############################################
E14=read.table('../Data/E14_ASM.vcf',skip=30,header = T)
colnames(E14)=c('chr','start','filter','Ref','Alt','QUAL','Het','Info','E14')
E14$end=E14$start
E14=makeGRangesFromDataFrame(E14)

##############################Motifbreak_R analysis###################################
if (!requireNamespace("MotifDb", quietly = TRUE))
{BiocManager::install("MotifDb")}
library(MotifDb)
if (!requireNamespace("motifbreakR", quietly = TRUE))
{BiocManager::install("motifbreakR")}
library(motifbreakR)
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP142.GRCh37", quietly = TRUE))
{BiocManager::install("SNPlocs.Hsapiens.dbSNP142.GRCh37")}
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
{BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
library("BSgenome.Hsapiens.UCSC.hg19")
library(BiocParallel)
library(motifbreakR)
pval_cutoff=0.1
variant_HetCpG_meta=readRDS('../downstream/output/ASM_enrich_meta.rds')
variant_HetCpG_NME=variant_HetCpG_meta[!variant_HetCpG_meta$Subject=='GM12878' & variant_HetCpG_meta$Statistic=='dNME']
#All SNP
motif=motif_break(unique(variant_HetCpG_NME))
#Unique SNP: do it tonight
#Do it for all motifs
variant_HetCpG_NME=variant_HetCpG_meta[variant_HetCpG_meta$Statistic=='dNME' & variant_HetCpG_meta$pvalue<=pval_cutoff]
for(sp in unique(variant_HetCpG_NME$Sample)){saveRDS(variant_HetCpG_NME[variant_HetCpG_NME$Sample==sp],paste('../downstream/motif/',sp,'_dNME_SNP.rds',sep=''))}
#Create a list of files
variant_HetCpG_NME_uq=unique(variant_HetCpG_NME)
variant_HetCpG_NME_uq=variant_HetCpG_NME_uq[!variant_HetCpG_NME_uq$Subject=='GM12878']
variant_HetCpG_NME_uq=variant_HetCpG_NME_uq[variant_HetCpG_NME_uq$HetCpG]
motif=motif_break(variant_HetCpG_NME_uq)
motif_break_noH1GM <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/motif_break_noH1GM.rds")
#Looing for strong motif
motif_break_noH1GM_strong=motif_break_noH1GM[motif_break_noH1GM$effect=='strong']
#This is the location of SNP with significant pval Need to do it for each sample
motif_break<-function(gr_in){
  gr_in_gr=granges(gr_in)
  strand(gr_in_gr)='+'
  attributes(gr_in_gr)$genome.package="BSgenome.Hsapiens.UCSC.hg19"
  #Make sure ref and alt agree with BS genome #may not necessary
  ref_BS=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,gr_in_gr))
  alt_BS=gr_in$ALT
  switch_nucleotide=which(gr_in$REF!=ref_BS)
  alt_BS[switch_nucleotide]=gr_in$REF[switch_nucleotide]
  #Get the necessary information for program running
  gr_in_gr$REF=unlist(DNAStringSetList(ref_BS))
  gr_in_gr$ALT=unlist(DNAStringSetList(alt_BS))
  names(gr_in_gr)=gr_in$snpId
  HS.SELEX=subset (MotifDb, organism=='Hsapiens'&dataSource=="jolma2013")
  results <- motifbreakR(snpList = gr_in_gr, filterp = TRUE,
                         pwmList = HS.SELEX,
                         threshold = 1e-4,
                         method = "ic",
                         bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),verbose=T,
                         BPPARAM = SnowParam(workers=22))

  return(results)
}
##############################Onuchic paper###################################
if (!requireNamespace("jsonlite", quietly = TRUE))
{BiocManager::install("jsonlite")}
library(jsonlite)
Onuchic_SNP=fromJSON('../downstream/input/AllelicEpigenome-sigOnly-AllDocs.jsonld',flatten=FALSE)
#Extract data: in paper 241360 loci, here 240811 loci, there's pval missing
Onuchic_SNP_df=data.frame(REF=Onuchic_SNP$`Reference Allele`,
                          is_enc=Onuchic_SNP$`Is In Enhancer Region`,
                          GWAS=Onuchic_SNP$`Is Near GWAS Variant`,
                          chr=Onuchic_SNP$Chromosome,
                          start=Onuchic_SNP$Position,
                          ALT=Onuchic_SNP$`Alternative Allele`,
                          REF_met=Onuchic_SNP$`Combined Data Analysis`$`Ref Allele Methylated CpG Count`,
                          REF_unmet=Onuchic_SNP$`Combined Data Analysis`$`Ref Allele Unmethylated CpG Count`,
                          ALT_met=Onuchic_SNP$`Combined Data Analysis`$`Alt Allele Methylated CpG Count`,
                          ALT_unmet=Onuchic_SNP$`Combined Data Analysis`$`Alt Allele Unmethylated CpG Count`,
                          stringsAsFactors = F)
Onuchic_SNP_df$end=Onuchic_SNP_df$start
#After combining data, much more CpG in REF than in ALT, given more homozygous than heterozygous
Onuchic_SNP_df$REF_CpG=Onuchic_SNP_df$REF_met+Onuchic_SNP_df$REF_unmet
Onuchic_SNP_df$ALT_CpG=Onuchic_SNP_df$ALT_met+Onuchic_SNP_df$ALT_unmet
hist(log10(Onuchic_SNP_df$REF_CpG/Onuchic_SNP_df$ALT_CpG))
Onuchic_SNP_gr=makeGRangesFromDataFrame(Onuchic_SNP_df,keep.extra.columns = TRUE)
GR=readRDS('../downstream/output/GR.all.diff.H1.GM12878.rds')
GR=GR[!GR$Subject %in% c("H1","GM12878")]
GR_UC=GR[GR$Statistic=='UC' & GR$pvalue<=pval_cutoff]
GR_dMML=GR[GR$Statistic=='dMML' & GR$pvalue<=pval_cutoff]
GR_dNME=GR[GR$Statistic=='dNME']
#Find overlapped regions
Onuchic_SNP_gr=subsetByOverlaps(Onuchic_SNP_gr,GR)
Onuchic_SNP_gr$diff=abs(Onuchic_SNP_gr$REF_met/(Onuchic_SNP_gr$REF_met+Onuchic_SNP_gr$REF_unmet)-
                          Onuchic_SNP_gr$ALT_met/(Onuchic_SNP_gr$ALT_met+Onuchic_SNP_gr$ALT_unmet))

Onuchic_SNP_gr_dMML=subsetByOverlaps(Onuchic_SNP_gr,GR_merge_exclude_GM[GR_merge_exclude_GM$dMML_pval<=pval_cutoff])
Onuchic_SNP_gr_dMML=Onuchic_SNP_gr_dMML[Onuchic_SNP_gr_dMML$REF_CpG>100 &Onuchic_SNP_gr_dMML$ALT_CpG>100]
Onuchic_SNP_gr_dMML=Onuchic_SNP_gr_dMML[order(Onuchic_SNP_gr_dMML$diff,decreasing = T)]
GR_Onuchic_olap=findOverlaps(GR_merge_exclude_GM[GR_merge_exclude_GM$dMML_pval<=pval_cutoff],Onuchic_SNP_gr_dMML)
GR_dMML_Onuchic=GR_merge_exclude_GM[GR_merge_exclude_GM$dMML_pval<=pval_cutoff][queryHits(GR_Onuchic_olap)]
GR_dMML_Onuchic$Onuchic_diff=Onuchic_SNP_gr_dMML$diff[subjectHits(GR_Onuchic_olap)]
GR_dMML_Onuchic_top_dMML=GR_dMML_Onuchic[GR_dMML_Onuchic$Onuchic_diff>0.8 &GR_dMML_Onuchic$dMML>0.8]
unique(GR_dMML_Onuchic_top_dMML)[1:6]
#High dNME region can have little dMML


#Find overlaps 
length(unique(GR_dMML))#7699
length(unique(subsetByOverlaps(GR_UC,subsetByOverlaps(GR_dMML,GR_dNME,type='equal'),type='equal')))#2243
length(unique(subsetByOverlaps(GR_dMML,GR_dNME,type='equal')))#3254
length(unique(subsetByOverlaps(GR_dMML,GR_UC,type='equal')))#3749
length(unique(subsetByOverlaps(GR_dNME,GR_UC,type='equal')))#3063
length(unique(GR_UC))#7207
length(unique(GR_dNME))#207626
#Find overlap with Onuchic SNP
length(subsetByOverlaps(Onuchic_SNP_gr,GR_dMML))#1883
length(subsetByOverlaps(Onuchic_SNP_gr,GR_dNME))#5617
length(subsetByOverlaps(Onuchic_SNP_gr,GR_UC))#1677
length(subsetByOverlaps(Onuchic_SNP_gr,GR,maxgap = 100))
#Find dMML data
GR_dMML=GR[GR$Statistic=='dMML']
GR_dNME=GR[GR$Statistic=='dNME']
GR_dMML_SD_ASM=subsetByOverlaps(GR_dMML,Onuchic_SNP_gr)

#Perchr data summary
fp='../../../../../../allele_specific/AllelicEpigenome-AllChrs-AllDocs.jsonld/'
chr_df=list(Subject=c(),Chromosome=c(),Position=c(),`Reference Allele`=c(),`Alternative Allele`=c(),
            `Is In Enhancer Region`=c(), `Is Near GWAS Variant`=c(),`Is On Heterogenous CpG`=c(),
            `Is In Promoter Region`=c(),`Tissue Specific Analysis`=list(),
            `Transcription Factor`=c(),`1000 Genomes Allele Frequency`=c(),
            REF_met=c(), REF_unmet=c(),ALT_met=c(),ALT_unmet=c())
for (fin in dir(fp)){
  tt1=proc.time()[[1]]
  cat("Reading", fin,'\n')
  gc()
  chr_in=fromJSON(paste(fp,fin,sep=''),flatten=FALSE)
  cat("processing",fin,'\n')
  chr_in_list=chr_in[c("Subject","Chromosome","Position","Reference Allele","Alternative Allele",
           "Is In Enhancer Region", "Is Near GWAS Variant","Is On Heterogenous CpG",
           "Is In Promoter Region","Tissue Specific Analysis",
           "Transcription Factor","1000 Genomes Allele Frequency")]
  chr_in_list$REF_met=chr_in$`Combined Data Analysis`$`Ref Allele Methylated CpG Count`
  chr_in_list$REF_unmet=chr_in$`Combined Data Analysis`$`Ref Allele Unmethylated CpG Count`
  chr_in_list$ALT_met=chr_in$`Combined Data Analysis`$`Alt Allele Methylated CpG Count`
  chr_in_list$ALT_unmet=chr_in$`Combined Data Analysis`$`Alt Allele Unmethylated CpG Count`
  chr_df=mapply(c,chr_df,chr_in_list,SIMPLIFY=F)
  cat(tail(chr_df$REF_met),'\n')
  cat("Finish processing",fin,"in",proc.time()[[1]]-tt1,'\n')
}
saveRDS(chr_df,"D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/Onuchic_SNP.rds")
Onuchic_SNP <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/Onuchic_SNP.rds")
#Make granges over Onuchic SNP
Onuchic_SNP_df=data.frame(REF=Onuchic_SNP$`Reference Allele`,
                          GWAS=Onuchic_SNP$`Is Near GWAS Variant`,
                          chr=Onuchic_SNP$Chromosome,
                          start=Onuchic_SNP$Position,
                          ALT=Onuchic_SNP$`Alternative Allele`,
                          REF_met=Onuchic_SNP$REF_met,
                          REF_unmet=Onuchic_SNP$REF_unmet,
                          ALT_met=Onuchic_SNP$ALT_met,
                          ALT_unmet=Onuchic_SNP$ALT_unmet,
                          stringsAsFactors = F)
Onuchic_SNP_df=Onuchic_SNP_df[!is.na(Onuchic_SNP_df$REF_met),]
Onuchic_SNP_df$end=Onuchic_SNP_df$start
#After combining data, much more CpG in REF than in ALT, given more homozygous than heterozygous
Onuchic_SNP_df$REF_CpG=Onuchic_SNP_df$REF_met+Onuchic_SNP_df$REF_unmet
Onuchic_SNP_df$ALT_CpG=Onuchic_SNP_df$ALT_met+Onuchic_SNP_df$ALT_unmet
hist(log10(Onuchic_SNP_df$REF_CpG/Onuchic_SNP_df$ALT_CpG))
Onuchic_SNP_gr=makeGRangesFromDataFrame(Onuchic_SNP_df,keep.extra.columns = TRUE)
Onuchic_SNP_gr$diff=abs(Onuchic_SNP_gr$REF_met/Onuchic_SNP_gr$REF_CpG-Onuchic_SNP_gr$ALT_met/Onuchic_SNP_gr$ALT_CpG)
#Find overlap with dNME region
dNME_Onuchic=subsetByOverlaps(GR_merge_exclude_GM[GR_merge_exclude_GM$dNME_pval<=pval_cutoff],Onuchic_SNP_gr)
dNME_Onuchic$N=NA
dNME_Onuchic$Subject=unlist(lapply(dNME_Onuchic$Sample,function(x) strsplit(x,' - ')[[1]][2]))
for (sp in unique(dNME_Onuchic$Subject)){
  olap=findOverlaps(dNME_Onuchic[dNME_Onuchic$Subject==sp],gff_in[gff_in$Subject==sp])
  dNME_Onuchic$N[dNME_Onuchic$Subject==sp][queryHits(olap)]=gff_in$N[gff_in$Subject==sp][subjectHits(olap)]
  
  
}
dNME_Onuchic=dNME_Onuchic[dNME_Onuchic$N>=10]
dNME_Onuchic_high_dNME=dNME_Onuchic[dNME_Onuchic$dNME>=0.8]
unique(dNME_Onuchic_high_dNME[order(dNME_Onuchic_high_dNME$dMML,decreasing = F)])#2nd one
subsetByOverlaps(Onuchic_SNP_gr, unique(dNME_Onuchic_high_dNME[order(dNME_Onuchic_high_dNME,decreasing = F)])[2])
dNME_Onuchic_low_dMML=dNME_Onuchic[dNME_Onuchic$dMML<=0.1 &dNME_Onuchic$dNME>=0.6]
dNME_Onuchic_low_dMML[order(dNME_Onuchic_low_dMML$dMML)]

#True is no data, false is have data
Onuchic_SNP_gr_raw$Tissue_data=is.na(Onuchic_SNP$`Tissue Specific Analysis`) | unlist(lapply(Onuchic_SNP$`Tissue Specific Analysis`,is.null))
SNP_sig=readRDS('../downstream/output/Onuchic_SNP_gr.rds')
#Number of SD-ASM overlapping Onuchic SNP
Onuchic_SNP_gr_sig=subsetByOverlaps(Onuchic_SNP_gr_raw,SNP_sig)
variant_HetCpG_new <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/variant_HetCpG_new.rds")
variant_HetCpG_new=variant_HetCpG_new[1:7]
names(variant_HetCpG_new)=NULL
variant_HetCpG_new=do.call(c,variant_HetCpG_new)
olap_SNP=findOverlaps(Onuchic_SNP_gr_raw,variant_HetCpG_new)
subsetByOverlaps(Onuchic_SNP_gr_raw,variant_HetCpG_new)
non_olap=which(!1:length(Onuchic_SNP_gr_raw) %in% queryHits(olap_SNP))
non_olap_not_na=non_olap[!non_olap %in% which(Onuchic_SNP_gr_raw$Tissue_data)]
#Which patient are those
non_olap_patient=do.call("c",lapply(Onuchic_SNP$`Tissue Specific Analysis`[non_olap_not_na],function(x) x$Patient))
#GR analyzed
GR <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/GR.all.diff.H1.GM12878.rds")
GR=GR[!GR$Subject %in% c("H1","GM12878")]
Onuchic_SNP_gr_vcf=subsetByOverlaps(Onuchic_SNP_gr_raw,variant_HetCpG_new)
SNP_with_stat=subsetByOverlaps(Onuchic_SNP_gr_vcf,GR,maxgap = 100)#Might be some double counts, maxgap=100
#Check SNP for those we don't have stat
SNP_with_stat_olap=findOverlaps(Onuchic_SNP_gr_vcf,GR,maxgap = 100)
Onuchic_SNP_gr_vcf_noolap=Onuchic_SNP_gr_vcf[!(1:length(Onuchic_SNP_gr_vcf) %in% queryHits(SNP_with_stat_olap))]
Onuchic_SNP_gr_vcf_noolap$SNP=variants_collapase(Onuchic_SNP_gr_vcf_noolap)
SNP_with_stat$SNP=variants_collapase(SNP_with_stat)

SNP_are_sig=subsetByOverlaps(Onuchic_SNP_gr_sig,Onuchic_SNP_gr_vcf)
subsetByOverlaps(SNP_are_sig,SNP_with_stat)
sum(Onuchic_SNP_gr_sig$Tissue_data)
SNP_tissue_no_WGS=subsetByOverlaps(Onuchic_SNP_gr_sig[!Onuchic_SNP_gr_sig$Tissue_data],Onuchic_SNP_gr_raw[non_olap_not_na])
subsetByOverlaps(Onuchic_SNP_gr_vcf,GR.all.diff)
sum(1:length(Onuchic_SNP_gr_vcf) %in% queryHits(findOverlaps(Onuchic_SNP_gr_vcf,GR.all.diff)))
sum(!1:length(Onuchic_SNP_gr_vcf) %in% queryHits(findOverlaps(Onuchic_SNP_gr_vcf,GR.all.diff)))
#Finding SNPs that have Stats
Onuchic_SNP_gr_olap=findOverlaps(Onuchic_SNP_gr_raw,GR[GR$Statistic=='dMML'],maxgap = 100)
Onuchic_SNP_qt=unique(queryHits(Onuchic_SNP_gr_olap))#1163614, maxgap=200:2314666
Onuchic_SNP_gr_df=as.data.frame(Onuchic_SNP_gr_raw[Onuchic_SNP_qt])
#Getting ASM stats
library(data.table)
Onuchic_SNP_stats=lapply(Onuchic_SNP$`Tissue Specific Analysis`[Onuchic_SNP_qt],extract_stats)
Onuchic_SNP_gr_stats=lapply(1:length(Onuchic_SNP_stats), function(x) merge_stat(Onuchic_SNP_stats[[x]],Onuchic_SNP_gr_df[x,]))
Onuchic_SNP_gr_stats=rbindlist(Onuchic_SNP_gr_stats)
Onuchic_SNP_gr_stats=as.data.frame(Onuchic_SNP_gr_stats)
Onuchic_SNP_gr_stats$end=Onuchic_SNP_gr_stats$start
Onuchic_SNP_gr_stats_sig=Onuchic_SNP_gr_stats[Onuchic_SNP_gr_stats$pval<=pval_cutoff,]
dMML=GR[GR$Statistic=='dMML']
Onuchic_SNP_gr_stats=makeGRangesFromDataFrame(Onuchic_SNP_gr_stats,keep.extra.columns = T)

#Find number of overlap with each tissue-sample combination, both have data
olap_single_sample=data.frame()

for (subj in unique(dMML$Subject)){
  for (ts in unique(dMML$Tissue)){
    cat("processing", subj,ts,'\n')
    dMML_sp=dMML[dMML$Subject==subj & dMML$Tissue==ts]
    Onuchic_sp=Onuchic_SNP_gr_stats[Onuchic_SNP_gr_stats$subject == subj & Onuchic_SNP_gr_stats$sample==ts]
    if (length(dMML_sp)!=0){
      cat("Counting", subj,ts,'\n')
      olap_sig=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
      olap_sig_cpel=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
      olap_sig_Onuchic=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
      olap_nonsig=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
      olap_single_sample=rbind(olap_single_sample,data.frame(sp=paste(subj,ts,sep='-'),olap_sig=olap_sig,
                                          olap_sig_cpel=olap_sig_cpel,olap_sig_Onuchic=olap_sig_Onuchic,olap_nonsig=olap_nonsig))
    }
  }
}


saveRDS(Onuchic_SNP_gr_stats,'../downstream/output/Onuchic_SNP_df_within_GR.rds')
#reform all regions
Onuchic_SNP_gr_df=as.data.frame(Onuchic_SNP_gr_raw)
Onuchic_SNP_stats=lapply(Onuchic_SNP$`Tissue Specific Analysis`,extract_stats)
Onuchic_SNP_gr_stats=lapply(1:length(Onuchic_SNP_stats), function(x) merge_stat(Onuchic_SNP_stats[[x]],Onuchic_SNP_gr_df[x,]))
Onuchic_SNP_gr_stats=rbindlist(Onuchic_SNP_gr_stats)
Onuchic_SNP_gr_stats$end=Onuchic_SNP_gr_stats$start
Onuchic_SNP_gr_stats=makeGRangesFromDataFrame(as.data.frame(Onuchic_SNP_gr_stats),keep.extra.columns = T)
saveRDS(Onuchic_SNP_gr_stats,'../downstream/output/Onuchic_SNP_df_all.rds')
#Check if regions having stats
#Get rid of regions with unavailable dataset
Onuchic_SNP_gr_stats_ft=Onuchic_SNP_gr_stats[!(Onuchic_SNP_gr_stats$subject %in% c("112","149","150","HuFGM02") | 
                                               Onuchic_SNP_gr_stats$sample == "Penis, Foreskin, Fibroblast Primary Cells")]
dMML_sp=single_sample(GR,Onuchic_SNP_gr_stats_ft,'dMML')
dNME_sp=single_sample(GR,Onuchic_SNP_gr_stats_ft,'dNME')
single_sample<-function(CPEL_in,Onuchic_in,stat_in){
  #Matching names
  GR_no_olap=GRanges()
  CPEL_in=CPEL_in[CPEL_in$Statistic==stat_in]
  CPEL_in$Tissue[CPEL_in$Tissue=='Adipose']="Adipose Tissue"
  CPEL_in$Tissue[CPEL_in$Tissue=='Embryonic Stem Cell']="H9 Cell Line"
  CPEL_in$Tissue[CPEL_in$Tissue=='Foreskin Melanocyte']="Penis, Foreskin, Melanocyte Primary Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Foreskin Keratinocyte']="Penis, Foreskin, Keratinocyte Primary Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Embyonic Stem Cell']="HUES64 Cell Line"
  CPEL_in$Tissue[CPEL_in$Tissue=='Ectoderm']="hESC Derived CD56+ Ectoderm Cultured Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Mesoderm']="hESC Derived CD56+ Mesoderm Cultured Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Endoderm']="hESC Derived CD184+ Endoderm Cultured Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Liver']="Adult Liver"
  olap_single_sample=data.frame()
  for (subj in unique(CPEL_in$Subject)){
    for (ts in unique(CPEL_in$Tissue)){
      cat("processing", subj,ts,'\n')
      CPEL_in_sp=CPEL_in[CPEL_in$Subject==subj & CPEL_in$Tissue==ts]
      Onuchic_sp=Onuchic_in[Onuchic_in$subject == subj & Onuchic_in$sample==ts]
      if (length(CPEL_in_sp)!=0){
        cat("Counting", subj,ts,'\n')
        olap=findOverlaps(Onuchic_sp,CPEL_in_sp,maxgap = 100)
        GR_no_olap=c(GR_no_olap,Onuchic_sp[!(1:length(Onuchic_sp) %in% queryHits(olap))])
        olap_GR=length(subsetByOverlaps(Onuchic_sp,CPEL_in_sp,maxgap = 100))
        olap_sig=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
        olap_sig_cpel=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
        olap_sig_Onuchic=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
        olap_nonsig=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
        olap_single_sample=rbind(olap_single_sample,data.frame(olap=olap_GR,non_olap=length(Onuchic_sp)-olap_GR,sp=paste(subj,ts,sep='-'),olap_sig=olap_sig,
                                                               olap_sig_cpel=olap_sig_cpel,olap_sig_Onuchic=olap_sig_Onuchic,olap_nonsig=olap_nonsig))
      }
    }
  }
  return(list(olap_single_sample,GR_no_olap))
}
Onuchic_SNP_gr_stats_ft_with_stat=subsetByOverlaps(Onuchic_SNP_gr_stats_ft,GR,maxgap = 100)
sum(Onuchic_SNP_gr_stats$subject %in% c("112","149","150","HuFGM02") | Onuchic_SNP_gr_stats$sample == "Penis, Foreskin, Fibroblast Primary Cells" )
#Extract stats from Onuchic paper
extract_stats<-function(dat_in){
  dat_in=dat_in[dat_in$Experiment$`foaf:name`=="Bisulfite-Seq",]
  return(data.frame(subject=dat_in$Patient,sample=dat_in$Tissue$`rdfs:label`,diff=dat_in$`Methylation Difference`,pval=dat_in$`FDR P-value`,
                    REF_unmet=dat_in$`Ref Allele Unmethylated CpG Count`,
                    REF_met=dat_in$`Ref Allele Methylated CpG Count`,
                    ALT_unmet=dat_in$`Alt Allele Unmethylated CpG Count`,
                    ALT_met=dat_in$`Alt Allele Methylated CpG Count`,
                    stringsAsFactors = F))
}
#Extract location for each stat
merge_stat<-function(stat_df,gr_df){
  if(nrow(stat_df)>0){
  stat_df$chr=gr_df$seqnames
  stat_df$start=gr_df$start
  return(stat_df)
  }
  
}
