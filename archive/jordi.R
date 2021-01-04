###############################################################################################################
# Deps
###############################################################################################################
# Genomics
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(coMET)

# # Parallelization
# numCores <- detectCores()
# library(foreach)
# library(doParallel)

# Plots
library(ggplot2)
library(RColorBrewer)

# Source main functions
setwd("/Users/jordiabante/Documents/code/R/HASM-MetaAnalysis/")
source("mainFunctions.R")

###############################################################################################################
# Read in CpelAsm results
###############################################################################################################

# Directories
inDir <- "/Users/jordiabante/Documents/CloudStation/ASM/Data/Onuchic 2018/"
# inDir <- "/cis/home/jabante/Desktop/Differential Analysis/"
CpGdir <- "/Users/jordiabante/Documents/CloudStation/ASM/Data/CpG sites hg19/"
# CpGdir <- "/cis/home/jabante/Desktop/Differential Analysis/"
GTEx <- "/Users/jordiabante/Documents/CloudStation/ASM/Data/GTEx/"
# GTEx <- "/cis/home/jabante/Desktop/Differential Analysis/"

# Read in CpelAsm data
# cpelGRs <- resultsCpelAsm(inDir)
# saveRDS(cpelGRs,file="/Users/jordiabante/Documents/CloudStation/ASM/Data/Onuchic 2018/cpelGRs.rds")
cpelGRs <- readRDS(file="/Users/jordiabante/Documents/CloudStation/ASM/Data/Onuchic 2018/cpelGRs.rds")

###############################################################################################################
# Satacked Histogram of ASM in each sample
###############################################################################################################

# dMML plot
prop_table <- 100*prop.table(xtabs(~ASM+Sample,data=cpelGRs[cpelGRs$Statistic=="dMML"]),margin=2)
pdf(file=paste(inDir,"dMML - Proportion of haplotype ASM events per sample.pdf",sep=""),width=8,height=6)
barplot(prop_table,beside=FALSE,col = c("lightblue","mistyrose"),legend=TRUE,las=2,
        cex.names=.5,args.legend = list(x = "topright", y=1.5, bty = "n"),
        main="dMML: Proportion of haplotype ASM events per sample")
dev.off()

# dNME plot
prop_table <- 100*prop.table(xtabs(~ASM+Sample,data=cpelGRs[cpelGRs$Statistic=="dNME"]),margin=2)
pdf(file=paste(inDir,"dNME - Proportion of haplotype ASM events per sample.pdf",sep=""),width=8,height=6)
barplot(prop_table,beside=FALSE,col = c("lightblue","mistyrose"),legend=TRUE,las=2,
        cex.names=.5,args.legend = list(x = "topright", y=1.5, bty = "n"),
        main="dNME: Proportion of haplotype ASM events per sample")
dev.off()

# UC plot
prop_table <- 100*prop.table(xtabs(~ASM+Sample,data=cpelGRs[cpelGRs$Statistic=="UC"]),margin=2)
pdf(file=paste(inDir,"UC - Proportion of haplotype ASM events per sample.pdf",sep=""),width=8,height=6)
barplot(prop_table,beside=FALSE,col = c("lightblue","mistyrose"),legend=TRUE,las=2,
        cex.names=.5,args.legend = list(x = "topright", y=1.5, bty = "n"),
        main="UC: Proportion of haplotype ASM events per sample")
dev.off()

# Clean environment
rm(list = c("prop_table"))

###############################################################################################################
# Genomic Features dMML, dNME, and UC break-down
###############################################################################################################

# Get all genomic features
genomicFeatures <- getGeneralFeats(CpGdir)

# Get overlaps
GRs_feats <- GRanges()
for(feature in names(genomicFeatures)){
  # Get feature GR
  GRs_tmp <- subsetByOverlaps(cpelGRs,genomicFeatures[[feature]],type="within")
  GRs_tmp$Feature <- feature
  GRs_feats <- append(GRs_feats,GRs_tmp)
}
GRs_feats_df <- as.data.frame(GRs_feats)
GRs_feats_df <- GRs_feats_df[,6:ncol(GRs_feats_df)]
GRs_feats_df$Feature <- factor(GRs_feats_df$Feature,levels=c("genome-wide","CpG island","CpG shore",
                              "CpG shelf","CpG open sea","gene body","exon","intron","intergenic","promoter"))
rm(list = c("genomicFeatures","GRs_tmp","GRs_feats"))

# Boxplots for all haplotypes
pdf(file=paste(inDir,"Observed Statistics in All Haplotypes - Genomic Features Boxplot.pdf",sep=""),width=8,height=6)
ggplot(GRs_feats_df,aes(x=Feature,y=Value,fill=Statistic)) + 
  geom_boxplot(alpha=0.5,position=position_dodge(1),outlier.shape=NA) + ylim(c(0,1)) +
  stat_boxplot(geom ='errorbar',position=position_dodge(1)) +
  theme(legend.position="top",panel.border=element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour="black")) +
  ggtitle("Observed Statistics in All Haplotypes")
dev.off()

# Boxplots for ASM haplotypes
pdf(file=paste(inDir,"Observed Statistics in ASM Haplotypes - Genomic Features Boxplot.pdf",sep=""),width=8,height=6)
ggplot(GRs_feats_df[GRs_feats_df$ASM=="Yes",],aes(x=Feature,y=Value,fill=Statistic)) + 
  geom_boxplot(alpha=0.5,position=position_dodge(1),outlier.shape=NA) + ylim(c(0,1)) +
  stat_boxplot(geom ='errorbar',position=position_dodge(1)) +
  theme(legend.position="top",panel.border=element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour="black")) +
  ggtitle("Observed Statistics in ASM Haplotypes")
dev.off()

# Boxplots for non-ASM haplotypes
pdf(file=paste(inDir,"Observed Statistics in Non-ASM Haplotypes - Genomic Features Boxplot.pdf",sep=""),width=8,height=6)
ggplot(GRs_feats_df[GRs_feats_df$ASM=="No",],aes(x=Feature,y=Value,fill=Statistic)) + 
  geom_boxplot(alpha=0.5,position=position_dodge(1),outlier.shape=NA) + ylim(c(0,1)) +
  stat_boxplot(geom ='errorbar',position=position_dodge(1)) +
  theme(legend.position="top",panel.border=element_blank(),
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour="black")) +
  ggtitle("Observed Statistics in Non-ASM Haplotypes")
dev.off()

# Clean environment
# rm(list = c("GRs_feats_df"))

###############################################################################################################
# Enrichment/depletion in different genomic features (constant across tissues)
###############################################################################################################
cpgIslandEnrichment <- testEnrichment(GRs_feats_df,"CpG island")
promoterEnrichment <- testEnrichment(GRs_feats_df,"promoter")

###############################################################################################################
# Relationship between N in haplotypes and proportion of differentially methylated haplotypes
###############################################################################################################

# dMML plot considering only those haplotypes with data
gffGR <- gffFormat(inDir,cpelGRs,"dMML")
prop_table <- 100*prop.table(xtabs(~ASM_bin+N,data=gffGR[gffGR$Data>0]),margin=2)
pdf(file=paste(inDir,"dMML - Percentage of haplotypes with at least one ASM event VS N.pdf",sep=""),width=10,height=6)
barplot(prop_table,beside=FALSE,col=c("lightblue","mistyrose"),legend=TRUE,las=1,
        cex.names=1,args.legend=list(x="topright",y=1.5,bty="n"),
        xlab="N",main="dMML: Percentage of haplotypes with at least one ASM event")
dev.off()

# dNME plot considering only those haplotypes with data
gffGR <- gffFormat(inDir,cpelGRs,"dNME")
prop_table <- 100*prop.table(xtabs(~ASM_bin+N,data=gffGR[gffGR$Data>0]),margin=2)
pdf(file=paste(inDir,"dNME - Percentage of haplotypes with at least one ASM event VS N.pdf",sep=""),width=10,height=6)
barplot(prop_table,beside=FALSE,col=c("lightblue","mistyrose"),legend=TRUE,las=1,
        cex.names=1,args.legend = list(x = "topright", y=1.5, bty = "n"),
        xlab="N",main="dNME: Percentage of haplotypes with at least one ASM event")
dev.off()

# UC plot considering only those haplotypes with data
gffGR <- gffFormat(inDir,cpelGRs,"UC")
prop_table <- 100*prop.table(xtabs(~ASM_bin+N,data=gffGR[gffGR$Data>0]),margin=2)
pdf(file=paste(inDir,"UC - Percentage of haplotypes with at least one ASM event VS N.pdf",sep=""),width=10,height=6)
barplot(prop_table,beside=FALSE,col=c("lightblue","mistyrose"),legend=TRUE,las=1,
        cex.names=1,args.legend = list(x = "topright", y=1.5, bty = "n"),
        xlab="N",main="UC: Percentage of haplotypes with at least one ASM event")
dev.off()

# Clean environment
rm(list=c("gffGR","prop_table"))

###############################################################################################################
# Relationship between CpG density (N/size) and proportion of differentially methylated haplotypes
###############################################################################################################

# dMML
gffGR <- gffFormat(inDir,cpelGRs,"dMML")
gffGR <- cpgDensity(gffGR)
gffGR <- gffGR[gffGR$Data>0]
gffGR$CpgDensBin <- cut(gffGR$CpgDens,seq(0,0.3,0.1),right=FALSE)
prop_table <- 100*prop.table(xtabs(~ASM_bin+CpgDensBin,data=gffGR),margin=2)
pdf(file=paste(inDir,"dMML - Percentage of haplotypes with at least one ASM event VS CpG density.pdf",sep=""),width=6,height=4)
barplot(prop_table,beside=FALSE,col=c("lightblue","mistyrose"),legend=TRUE,las=1,
        cex.names=1,args.legend = list(x = "topright", y=1.5, bty = "n"),
        xlab="CpG density",main="dMML: Percentage of haplotypes with at least one ASM event")
dev.off()

# dNME
gffGR <- gffFormat(inDir,cpelGRs,"dNME")
gffGR <- cpgDensity(gffGR)
gffGR <- gffGR[gffGR$Data>0]
gffGR$CpgDensBin <- cut(gffGR$CpgDens,seq(0,0.3,0.1),right=FALSE)
prop_table <- 100*prop.table(xtabs(~ASM_bin+CpgDensBin,data=gffGR),margin=2)
pdf(file=paste(inDir,"dNME - Percentage of haplotypes with at least one ASM event VS CpG density.pdf",sep=""),width=6,height=4)
barplot(prop_table,beside=FALSE,col=c("lightblue","mistyrose"),legend=TRUE,las=1,
        cex.names=1,args.legend = list(x = "topright", y=1.5, bty = "n"),
        xlab="CpG density",main="dNME: Percentage of haplotypes with at least one ASM event")
dev.off()

# UC
gffGR <- gffFormat(inDir,cpelGRs,"UC")
gffGR <- cpgDensity(gffGR)
gffGR <- gffGR[gffGR$Data>0]
gffGR$CpgDensBin <- cut(gffGR$CpgDens,seq(0,0.3,0.1),right=FALSE)
prop_table <- 100*prop.table(xtabs(~ASM_bin+CpgDensBin,data=gffGR),margin=2)
pdf(file=paste(inDir,"UC - Percentage of haplotypes with at least one ASM event VS CpG density.pdf",sep=""),width=6,height=4)
barplot(prop_table,beside=FALSE,col=c("lightblue","mistyrose"),legend=TRUE,las=1,
        cex.names=1,args.legend = list(x = "topright", y=1.5, bty = "n"),
        xlab="CpG density",main="UC: Percentage of haplotypes with at least one ASM event")
dev.off()

# Clean environment
rm(list=c("gffGR","prop_table"))

###############################################################################################################
# Relationship between number of SNPs in haplotypes and ASM
###############################################################################################################


###############################################################################################################
# Enrichment of ASM events in imprinted regions
###############################################################################################################

# dMML
dmmlImprinting <- imprintingEnrichment(GTEx,CpGdir,cpelGRs,testStat="dMML")
dmmlImprintingTest <- dmmlImprinting["test"]
dmmlImprintingGR <- dmmlImprinting["GR"]

# dNME
dnmeImprinting <- imprintingEnrichment(GTEx,CpGdir,cpelGRs,testStat="dNME")
dnmeImprintingTest <- dnmeImprinting["test"]
dnmeImprintingGR <- dnmeImprinting["GR"]

# UC
ucImprinting <- imprintingEnrichment(GTEx,CpGdir,cpelGRs,testStat="UC")
ucImprintingTest <- ucImprinting["test"]
ucImprintingGR <- ucImprinting["GR"]

###############################################################################################################
# Plot some imprinting control regions
###############################################################################################################

# Get CpG islands
genomicFeatures <- getGeneralFeats(CpGdir)
cpgIslands <- genomicFeatures[["CpG island"]]
cpgSites <- getCpgSitesH19()

# Plot H19
plotGR(asmInProms[asmInProms$geneSymbol=="H19"],chr="chr11",lim=c(1995176,2001470))

###############################################################################################################
# Co-ocurrence within subject across tissues
###############################################################################################################

# Do test on co-occurrence on dMML for each subject (donor or cell line)
perms <- 1000
uniqueSubjectLabels <- unique(subject_labels)
cooccurrence_table <- data.frame(Subject=uniqueSubjectLabels,Statistic_dMML=NA,pValue_dMML=NA,
                                 Statistic_dNME=NA,pValue_dNME=NA,Statistic_UC=NA,pValue_UC=NA)
for(i in 1:length(uniqueSubjectLabels)) {
  
  # Do test for subject 
  subject <- uniqueSubjectLabels[i]
  GRs_subject <- GRs[(GRs$Subject==subject)&(GRs$Statistic=="dMML")]
  
  # Skip if only one tissue
  if(length(unique(GRs_subject$Tissue))==1){
    next
  }
  
  # Prepare table
  GR_haps <- unique(GRs_subject[,c()])
  GR_haps$Haplotype <- paste("Hap-",seq(1:length(GR_haps)),sep="")
  olaps <- findOverlaps(GRs_subject,GR_haps,type="equal",ignore.strand=TRUE,select="all")
  GRs_subject$Haplotype <- NA
  GRs_subject[queryHits(olaps)]$Haplotype <- GR_haps[subjectHits(olaps)]$Haplotype
  subject_table <- xtabs(~ASM+Haplotype,data=GRs_subject)
  subject_table <- as.data.frame.matrix(t(subject_table))
  subject_table$Data <- subject_table$No+subject_table$Yes
  subject_table$No <- NULL

  # Compute statistic for permutation test
  obsStat <- testStat(subject_table)
 
  # Generate null stats permutation test
  nullStats <- c()
  nullStats <- foreach (j=1:perms,.combine=c) %dopar% {
    testStat(permuteTable(subject_table))
  }

  # Compute p-value permutation test
  pVal <- sum(nullStats>obsStat)/length(nullStats)
  
  # Store result permutation test
  cooccurrence_table[i,2] <- obsStat
  cooccurrence_table[i,3] <- pVal
}

# Do test on co-occurrence on dNME for each subject (donor or cell line)
for(i in 1:length(uniqueSubjectLabels)) {
  # Do test for subject 
  subject <- uniqueSubjectLabels[i]
  GRs_subject <- GRs[(GRs$Subject==subject)&(GRs$Statistic=="dNME")]
  
  # Skip if only one tissue
  if(length(unique(GRs_subject$Tissue))==1){
    next
  }
  
  # Prepare table
  GR_haps <- unique(GRs_subject[,c()])
  GR_haps$Haplotype <- paste("Hap-",seq(1:length(GR_haps)),sep="")
  olaps <- findOverlaps(GRs_subject,GR_haps,type="equal",ignore.strand=TRUE,select="all")
  GRs_subject$Haplotype <- NA
  GRs_subject[queryHits(olaps)]$Haplotype <- GR_haps[subjectHits(olaps)]$Haplotype
  subject_table <- xtabs(~ASM+Haplotype,data=GRs_subject)
  subject_table <- as.data.frame.matrix(t(subject_table))
  subject_table$Data <- subject_table$No+subject_table$Yes
  subject_table$No <- NULL
  
  # Compute statistic for permutation test
  obsStat <- testStat(subject_table)
  
  # Generate null stats permutation test
  nullStats <- c()
  nullStats <- foreach (j=1:perms,.combine=c) %dopar% {
    testStat(permuteTable(subject_table))
  }
  
  # Compute p-value permutation test
  pVal <- sum(nullStats>obsStat)/length(nullStats)
  
  # Store result permutation test
  cooccurrence_table[i,4] <- obsStat
  cooccurrence_table[i,5] <- pVal
}

# Do test on co-occurrence on UC for each subject (donor or cell line)
for(i in 1:length(uniqueSubjectLabels)) {
  # Do test for subject 
  subject <- uniqueSubjectLabels[i]
  GRs_subject <- GRs[(GRs$Subject==subject)&(GRs$Statistic=="UC")]
  
  # Skip if only one tissue
  if(length(unique(GRs_subject$Tissue))==1){
    next
  }
  
  # Prepare table
  GR_haps <- unique(GRs_subject[,c()])
  GR_haps$Haplotype <- paste("Hap-",seq(1:length(GR_haps)),sep="")
  olaps <- findOverlaps(GRs_subject,GR_haps,type="equal",ignore.strand=TRUE,select="all")
  GRs_subject$Haplotype <- NA
  GRs_subject[queryHits(olaps)]$Haplotype <- GR_haps[subjectHits(olaps)]$Haplotype
  subject_table <- xtabs(~ASM+Haplotype,data=GRs_subject)
  subject_table <- as.data.frame.matrix(t(subject_table))
  subject_table$Data <- subject_table$No+subject_table$Yes
  subject_table$No <- NULL
  
  # Compute statistic for permutation test
  obsStat <- testStat(subject_table)
  
  # Generate null stats permutation test
  nullStats <- c()
  nullStats <- foreach (j=1:perms,.combine=c) %dopar% {
    testStat(permuteTable(subject_table))
  }
  
  # Compute p-value permutation test
  pVal <- sum(nullStats>obsStat)/length(nullStats)
  
  # Store result permutation test
  cooccurrence_table[i,6] <- obsStat
  cooccurrence_table[i,7] <- pVal
}

###############################################################################################################
# Co-ocurrence across subjects
###############################################################################################################
