###############################################################################################################
# Deps
###############################################################################################################
# Genomics
setwd("E:/Dropbox/JHU/Projects/Allele-spcific/code/downstream")
library(ggplot2)
library(topGO)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(reshape2)
library(VIM)
library(FactoMineR)
library(missMDA)
library(ggbiplot)
# Parallelization
numCores <- detectCores()
library(foreach)
library(doParallel)
library(VariantAnnotation)

# Source
source("../../../git_rep/HASM-MetaAnalysis/mainFunctions.R")
#Test message: to delete
###############################################################################################################
# Read in
###############################################################################################################

# Directories
#inDir <- "/Users/jordiabante/Documents/CloudStation/ASM/Data/Onuchic 2018/"
inDir <- "../downstream/input/Julia_bedfile/"
# CpGdir <- "/Users/jordiabante/Documents/CloudStation/ASM/Data/CpC sites h19/"
#CpGdir <- "/cis/home/jabante/Desktop/Differential Analysis/"

GRs <- resultsCpelAsm (inDir)
GRs <- readRDS("E:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/GRs.rds")
#Generating PCA plot

shared_regions_STL<-function(cpelGRs,stats){
  # Get ranges
  rangesStl001 <- unique(cpelGRs[(cpelGRs$Subject=="STL001"),c()])
  rangesStl002 <- unique(cpelGRs[(cpelGRs$Subject=="STL002"),c()])
  rangesStl003 <- unique(cpelGRs[(cpelGRs$Subject=="STL003"),c()])
  rangesInt <- subsetByOverlaps(rangesStl001,rangesStl002,type="any")
  rangesInt <- subsetByOverlaps(rangesInt,rangesStl003,type="any")
  rangesInt$name <- paste("sharedHap-",seq(1:length(rangesInt)),sep="")
# Get data
  olaps <- findOverlaps(cpelGRs,rangesInt)
  cpelGRs$HapName <- NA
  cpelGRs$HapName[queryHits(olaps)] <- rangesInt$name[subjectHits(olaps)]
  pcaData <- cpelGRs[!is.na(cpelGRs$HapName)]
  olaps <- findOverlaps(cpelGRs,rangesInt)
  cpelGRs$HapName <- NA
  cpelGRs$HapName[queryHits(olaps)] <- rangesInt$name[subjectHits(olaps)]
  pcaData <- cpelGRs[!is.na(cpelGRs$HapName)]
  pcaDataDf <- as.data.frame(pcaData[,c("HapName","Statistic","Value","ASM","Sample")])[,c(6,7,8,9,10)]
  # Get matrix
  pcaDataDfstats <- pcaDataDf[pcaDataDf$Statistic==stats,][,c(1,3,4,5)]
  sigHapNames <- unique(pcaDataDfstats$HapName[pcaDataDfstats$ASM=="Yes"])
  pcaDataDfstats <- pcaDataDfstats[pcaDataDfstats$HapName %in% sigHapNames,]
  print(pcaDataDfstats)
  pca_decast=dcast(data=pcaDataDfstats,formula=Sample~HapName,fun.aggregate=mean,value.var="Value")
  rownames(pca_decast)=pca_decast$Sample
  return(pca_decast[,-1])
}

pca_temp=shared_regions_STL(GRs,'dMML')
#Generate sample name
sample=do.call('c',lapply(rownames(pca_temp),function(x) strsplit(x,'- ')[[1]][2]))
#Filter samples
STL_PCA=pca_temp[sample %in% c('STL001','STL002','STL003'),]
pca_na_select=apply(STL_PCA,2,function(x) sum(!is.na(x))>=30)
pca_decast_select=STL_PCA[pca_na_select]
#Perfrom PCA analysis
res.comp <- imputePCA(pca_decast_select, ncp = 2)
res.pca <- prcomp(res.comp$completeObs,center=TRUE,scale. = TRUE)
tt='dMML PCA'
PCA_plot=ggbiplot(res.pca,var.axes = FALSE,groups = sample[sample %in% c('STL001','STL002','STL003')],scale=1)+theme(legend.position="top")+
  ggtitle(tt) +
  theme(panel.border=element_blank(),panel.background=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"))

#for a given sample, calcualte enrichment of the enhancer region
CpG_sites=getCpgSitesH19()
dMML=GRs[GRs$Statistic=='dMML']
dMML_sig=dMML[dMML$ASM=="Yes"]
load("C:/Users/vince/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/enhancers_intersect.RData")
enhancer_gr_all=do.call('c',lapply(rownames(max_states),rownames2Granges))
#Coverting rownames to granges
rownames2Granges<-function(rn){
  stp=strsplit(rn,':')[[1]]
  chr=stp[1]
  start_end=strsplit(stp[2],'-')[[1]]
  start=as.numeric(start_end[1])
  end=as.numeric(start_end[2])
  return(GRanges(seqnames=chr,ranges=IRanges(start,end)))
}

enhancer_OR<-function(gr){
  CpG_enhancer_op=findOverlaps(CpG_sites,enhancer_gr_all)
  CpG_ASM_op=findOverlaps(CpG_sites,gr[gr$ASM=='Yes'])
  
  CpG_enhancer=CpG_sites[queryHits(CpG_enhancer_op)]
  CpG_not_enhancer=CpG_sites[-queryHits(CpG_enhancer_op)]
  CpG_ASM=CpG_sites[queryHits(CpG_ASM_op)]
  CpG_not_ASM=CpG_sites[-queryHits(CpG_ASM_op)]
  
  n_ASM_enhancer=length(subsetByOverlaps(CpG_enhancer,CpG_ASM))
  n_ASM_not_enhancer=length(subsetByOverlaps(CpG_not_enhancer,CpG_ASM))
  n_not_ASM_enhancer=length(subsetByOverlaps(CpG_enhancer,CpG_not_ASM))
  n_not_ASM_not_enhancer=length(subsetByOverlaps(CpG_not_enhancer,CpG_not_ASM))
  
  cont_table=matrix(c(n_ASM_enhancer,n_ASM_not_enhancer,n_not_ASM_enhancer,n_not_ASM_not_enhancer),ncol=2)
  # overflow_prevent=10^6
  # OR=(cont_table[1,1]*(cont_table[2,2]/overflow_prevent))/(cont_table[1,2]*(cont_table[2,1]/overflow_prevent))
  
  ft=fisher.test(cont_table)
  return(c(ft$estimate,ft$p.value))
}
stat_OR_calc<-function(GR_in){
  all_sample=unique(GR_in$Sample)
  out_sum = data.frame(sample=all_sample)
  out_sum$OR=NULL
  out_sum$pval=NULL
  rownames(out_sum)=all_sample
  for(sp in out_sum$sample){
    out_sum[sp,c(2,3)]=enhancer_OR(GR_in[GR_in$Sample==sp])
  }
  colnames(out_sum)=c('sample','OR','p-val')
  return(out_sum)
}
dMML_enrich=stat_OR_calc(dMML)
dNME_enrich=stat_OR_calc(dNME)
UC_enrich=stat_OR_calc(UC)


dMML=GRs[GRs$Statistic=='dMML']

gene_ranking<-function(gene_gr,GR_stat){
  GR_stat_all=subsetByOverlaps(gene_gr,GR_stat)
  GR_stat_all_gene=GR_stat_all$gene_name
  GR_stat_ASM=GR_stat[GR_stat$ASM=='Yes']
  GR_stat_ASM=GR_stat_ASM[order(GR_stat_ASM$Value,decreasing=TRUE)]
  GR_stat_ASM$rank=1:length(GR_stat_ASM)
  GR_stat_ASM_op=findOverlaps(gene_gr,GR_stat_ASM)
  #for each gene sum their GR_stat value
  gene_df=data.frame(gene=NA,value_sum=NA)
  for(i in unique(queryHits(GR_stat_ASM_op))){
    gf=gene_gr[i]
    val_sum=sum(GR_stat_ASM$Value[subjectHits(GR_stat_ASM_op)[queryHits(GR_stat_ASM_op)==i]])
    gene_df=rbind(gene_df,data.frame(gene=gf$gene_name,value_sum=val_sum))
  }
  gene_df=gene_df[order(gene_df$value_sum,decreasing = TRUE),]
  return(list(gene_df,GR_stat_all_gene))
}
GO_anno<-function(myInterestingGenes,geneNames,topNodes=200){
  #Make factor
  geneList=factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)=geneNames
  #find annotation orgs annot=annFUN.org,mapping="org.Hs.eg.db"
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, mapping="org.Hs.eg.db",ID="symbol")
  resultF <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultF.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  resultF.weight <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  list(GenTable(GOdata,classicF=resultF,elimF=resultF.elim,weightF=resultF.weight,orderBy="weightF",topNodes=topNodes),GOdata)
}
dMML_rank=gene_ranking(genome_features[[11]],dMML) #genome_features is at output
dMML_GO=GO_anno(dMML_rank[[1]]$gene[dMML_rank[[1]]$value_sum>1],dMML_rank[[2]])
dMML_GO_table=dMML_GO[[1]]
dMML_GO_table[dMML_GO_table$Expected>1,]
GO_genes<-function(goID,Go_dat,sig=TRUE){
  gene_term=genesInTerm(Go_dat,whichGO=goID)
  gene_sig=sigGenes(Go_dat)
  gene_sig_GO=gene_sig[unlist(gene_sig)%in%unlist(gene_term)]
  if (sig){gene_out=gene_sig_GO}else{gene_out=gene_term}
  return(gene_anno(gene_out))
  #get gene description

}
gene_anno<-function(gene_in){
  ensembl=biomaRt::useMart("ensembl",host="http://grch37.ensembl.org") #hg19
  ensembl = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
  gene_des=biomaRt::getBM(attributes=c('hgnc_symbol',"description"),filters='hgnc_symbol',values=gene_in,mart=ensembl)
  gene_des$description=gsub("\\[Source:HGNC Symbol;Acc:.*","",gene_des$description)
  return(gene_des)
}
#get variant frequency
read_vcf<-function(vcf_file){
  vcf=readVcf(vcf_file)
  return(rowRanges(vcf))
}
readin_all_vcf<-function(vcfs){
  varsDiff=do.call('c',lapply(vcfs,read_vcf))
  variants <- paste(as.character(varsDiff$REF),as.character(unlist(varsDiff$ALT)),sep="-")
  variants[variants %in% c("A-C","C-A")] <- "A-C"
  variants[variants %in% c("A-G","G-A")] <- "A-G"
  variants[variants %in% c("A-T","T-A")] <- "A-T"
  variants[variants %in% c("C-G","G-C")] <- "C-G"
  variants[variants %in% c("C-T","T-C")] <- "C-T"
  variants[variants %in% c("G-T","T-G")] <- "G-T"
  varsDiff$merge=variants
  return(varsDiff)
}
het_CpG<-function(var_gr){
  #Check het CpG separately
  add_var=c('C-T','A-C')
  minus_var=c('A-G','G-T')
  GC_var='C-G'
  var_gr$Het=NA
  #Check CpG
  var_gr$Het[var_gr$merge %in% add_var]=het_CpG_check(var_gr[var_gr$merge %in% add_var],'add')
  var_gr$Het[var_gr$merge %in% minus_var]=het_CpG_check(var_gr[var_gr$merge %in% minus_var],'minus')
  var_gr$Het[var_gr$merge %in% GC_var]=het_CpG_check(var_gr[var_gr$merge %in% GC_var],'GC')
  var_gr$Het[var_gr$merge =='A-T']=FALSE
  return(var_gr)
}
het_CpG_check<-function(var_gr,var_type){
  if (var_type=='add'){#check XG type 
    nucleotide=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,1)))
    pos=start(matchPattern('G',nucleotide))
    out=rep(FALSE,length(nucleotide))
    out[pos]=TRUE
    return(out) #If XG, heterogzyous nulceotide
  }else if(var_type=='minus'){#Check CX type
    nucleotide=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,-1)))
    pos=start(matchPattern('C',nucleotide))
    out=rep(FALSE,length(nucleotide))
    out[pos]=TRUE
    return(out)
  }else if(var_type=='GC'){#check C-G nucleotide separately
    nucleotide_before=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,-1)))
    nucleotide_after=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,1)))
    print(length(nucleotide_before))
    print(length(nucleotide_after))
    pos=unique(c(start(matchPattern('C',nucleotide_before)),start(matchPattern('G',nucleotide_after))))
    out=rep(FALSE,length(var_gr))
    out[pos]=TRUE
    return(out)
    
  }else{return('error type')}
  
}

getTriNucleotide<-function(x,var_gr){
  n_next=Hsapiens[[x]][start(var_gr[seqnames(var_gr)==x])+1]
  n_before=Hsapiens[[x]][start(var_gr[seqnames(var_gr)==x])-1]
  return(do.call('c',lapply(1:length(n_before),
                            function(x) paste(n_before[x],n_next[x],sep='X'))))
}

TriNucleotide_freq<-function(varsDiff,variant,GRin){
 varsDiff=varsDiff[varsDiff$merge==variant]
 varsDiff_gr=subsetByOverlaps(varsDiff,GRin,type='within')
  chrs=paste('chr',1:22,sep='')
  vars=do.call('c',lapply(chrs,getTriNucleotide,var_gr=varsDiff_gr))
  table(vars)
}
bar_plot_all<-function(bar_df){
  ggplot(data=bar_df,aes(x=vars,y=Freq,label=Freq))+
    geom_bar(stat='identity',fill='steelblue')
}
plot_het_CpG<-function(df_in,het_CpG){
  #Here frequency is average number of trinucleotide that discrupts CpG, i.e. numberof variants/number of context
  CG_het=mean(df_in$Freq[df_in$vars %in% het_CpG])
  CG_hom=mean(df_in$Freq[!(df_in$vars %in% het_CpG)])
  df_bar=data.frame(vars=c('Het CpG','No Het CpG'),Freq=c(CG_het,CG_hom))
  print(df_bar)
  bar_plot_all(df_bar)
}
calc_het_CpG<-function(df_in,het_CpG){
  #Here frequency is average number of trinucleotide that discrupts CpG, i.e. numberof variants/number of context
  CG_het=sum(df_in$Freq[df_in$vars %in% het_CpG])
  CG_hom=sum(df_in$Freq[!(df_in$vars %in% het_CpG)])
  df_bar=data.frame(vars=c('Het CpG','No Het CpG'),Freq=c(CG_het,CG_hom))
  return(df_bar)
}
DIR='../downstream/vcfFiles/'
vcfs=dir(DIR)
vcfs=gsub('.tbi','',vcfs)
vcfs=unique(vcfs)
vcf_in=readin_all_vcf(paste(DIR,vcfs,sep=''))
genome <- BSgenome.Hsapiens.UCSC.hg19
GRs=readRDS('../downstream/output/GRs_no_brain.rds')
dNME=GRs[GRs$Statistic=='dNME']

dNME_ASM_CT=as.data.frame(sort(TriNucleotide_freq(vcf_in,'C-T',dNME[dNME$ASM=='Yes']),decreasing=TRUE))
dNME_ASM_AG=as.data.frame(sort(TriNucleotide_freq(vcf_in,'A-G',dNME[dNME$ASM=='Yes']),decreasing=TRUE))
dNME_all_CT=as.data.frame(sort(TriNucleotide_freq(vcf_in,'C-T',dNME),decreasing=TRUE))
dNME_all_AG=as.data.frame(sort(TriNucleotide_freq(vcf_in,'A-G',dNME),decreasing=TRUE))
#Plot of all
bar_plot_all(dNME_ASM_AG)+ylab("Number of variant") +ggtitle('A-G variant at dNME ASM')
bar_plot_all(dNME_ASM_CT)+ylab("Number of variant") +ggtitle('C-T variant at dNME ASM')
bar_plot_all(dNME_all_AG)+ylab("Number of variant") +ggtitle('All A-G variant')
bar_plot_all(dNME_all_CT)+ylab("Number of variant") +ggtitle('All C-T variant')
#Plot of Het CpG or not
#Het CpG for CT AXG, CXG, TXG,GXG
het_CpG_CT=c('AXG','CXG','TXG','GXG')
plot_het_CpG(dNME_ASM_CT,het_CpG_CT)+ 
  xlab("Variant type") + ylab("Average number of variant") +ggtitle('C-T variant at dNME ASM')

plot_het_CpG(dNME_all_CT,het_CpG_CT)+ 
  xlab("Variant type") + ylab("Average number of variant") +ggtitle('all C-T variant at dNME analyzed')
#Het CpG for GA CXA, CXT, CXC,CXG
het_CpG_AG=c('CXA','CXT','CXC','CXG')
plot_het_CpG(dNME_ASM_AG,het_CpG_AG)+ 
  xlab("Variant type") + ylab("Average number of variant") +ggtitle('A-G variant at dNME ASM')

plot_het_CpG(dNME_all_AG,het_CpG_AG)+ 
  xlab("Variant type") + ylab("Average number of variant") +ggtitle('all G-A variant at dNME analyzed')
ASM_enrich<-function(ASM_var,het_CpG,dNME_all){
  het_calc_ASM=calc_het_CpG(ASM_var,het_CpG)
  het_calc_notASM=data.frame(vars=c('Het CpG','No Het CpG'),Freq=calc_het_CpG(dNME_all,het_CpG)$Freq-het_calc_ASM$Freq)
  ASM_enrich=cbind(het_calc_ASM$Freq,het_calc_notASM$Freq)
  colnames(ASM_enrich)=c('ASM','Not ASM')
  rownames(ASM_enrich)=c('Het','Not Het')
  print(fisher.test(ASM_enrich))
  return((ASM_enrich))
  
}
AG_enrich=ASM_enrich(dNME_ASM_AG,het_CpG_AG,dNME_all_AG)
CT_enrich=ASM_enrich(dNME_ASM_CT,het_CpG_CT,dNME_all_CT)
#CT AG enrichment
CT_het_calc=calc_het_CpG(dNME_all_CT,het_CpG_CT)
AG_het_calc=calc_het_CpG(dNME_all_AG,het_CpG_AG)
AG_CT_enrich=cbind(CT_het_calc$Freq,AG_het_calc$Freq)
colnames(AG_CT_enrich)=c('CT','AG')
rownames(AG_CT_enrich)=c('Het','Not Het')
fisher.test(AG_CT_enrich)

#Genome wide test
genome_wide_stat<-function(GR_in,vcf_het){
  het_size=length(unique(vcf_het[vcf_het$Het]$merge))
  all_var_length=length(unique(vcf_het))
  #Get numbers
  het_stat=length(subsetByOverlaps(vcf_het[vcf_het$Het],GR_in[GR_in$ASM=='Yes']))
  not_het_stat=length(subsetByOverlaps(vcf_het[-which(vcf_het$Het)],GR_in[GR_in$ASM=='Yes']))
  het_not_stat=length(subsetByOverlaps(vcf_het[which(vcf_het$Het)],GR_in[GR_in$ASM=='No']))
  not_het_not_stat=length(subsetByOverlaps(vcf_het[-which(vcf_het$Het)],GR_in[GR_in$ASM=='No']))
  cont_table=matrix(c(het_stat,not_het_stat,het_not_stat,not_het_not_stat),nrow=2)
  rownames(cont_table)=c('Het CpG','Not Het CpG')
  colnames(cont_table)=c('ASM','Not ASM')
  print(fisher.test(cont_table))
  print(cont_table)
  bar_plot_df=data.frame(vars=c('Het CpG','Not Het CpG'),Freq=c(het_stat,not_het_stat))
  bar_plot_all(bar_plot_df)
}
genome_wide_stat(dNME,vcf_het)+
  xlab("Variant type") + ylab("Average number of variant") +ggtitle('All variant at dNME ASM')
genome_wide_stat(GRs[GRs$Statistic=='dMML'],vcf_het)+
  xlab("Variant type") + ylab("Average number of variant") +ggtitle('All variant at dMML ASM')

genome_wide_stat(GRs[GRs$Statistic=='UC'],vcf_het)+
  xlab("Variant type") + ylab("Average number of variant") +ggtitle('All variant at UC ASM')

GR_dNME_het=subsetByOverlaps(dNME,vcf_het[vcf_het$Het])
dNM_het_rank=gene_ranking(genome_features[[11]],GR_dNME_het)
dNME_het_GO=GO_anno(dNM_het_rank[[1]]$gene,dNM_het_rank[[2]])
dNME_het_GO=dNME_het_GO[[1]]
dNME_het_GO[dNME_het_GO$Expected>=10,c(1,2,9)]

dNME_rank=gene_ranking(genome_features[[11]],dNME)
dNME_GO=GO_anno(dNME_rank[[1]]$gene,dNME_rank[[2]])
dNME_GO=dNME_GO[[1]]
dNME_GO[dNME_GO$Expected>=10,c(1,2,9)]


dMML_rank=gene_ranking(genome_features[[11]],GRs[GRs$Statistic=='dMML'])
dMML_GO=GO_anno(dMML_rank[[1]]$gene,dMML_rank[[2]])
dMML_GO=dMML_GO[[1]]
dMML_GO[dMML_GO$Expected>=10,c(1,2,9)]


UC_rank=gene_ranking(genome_features[[11]],GRs[GRs$Statistic=='UC'])
UC_GO=GO_anno(UC_rank[[1]]$gene,UC_rank[[2]])
UC_GO=UC_GO[[1]]
UC_GO[UC_GO$Expected>=10,c(1,2,9)]


#Venn diagram of all samples overlapping:

#GR from all Granges:
region_olap<-function(gr,subject){
  gr_in=granges(gr)
  gr_in$Subject=elementMetadata(gr)[,subject]
  gr_in_unique=unique(gr_in)
  gr_in_unique$CpG=paste('CpG:',1:length(gr_in_unique),sep='')
  olap=findOverlaps(gr_in,gr_in_unique)
  gr_in$CpG[queryHits(olap)]=gr_in_unique$CpG[subjectHits(olap)]
  gr_in$value=1
  gr_in_meta=elementMetadata(gr_in)
  gr_in_table=dcast(data=as.data.frame(gr_in_meta),formula=CpG ~ Subject,fun.aggregate =mean,value.var = "value",fill=0)
  gr_in_table$CpG=NULL
  gr_in_table$olap_number=rowSums(gr_in_table)
  sum_gr_in=table(gr_in_table$olap_number)
  return(sum_gr_in)
}
gff_in_sum=region_olap(gff_in,'Subject')
gff_in_sum_df=data.frame(overlap=0:7,rn= as.numeric(gff_in_sum))
ggplot(data=gff_in_sum_df,aes(x=overlap,y=rn))+geom_bar(stat='identity',fill="steelblue")+
  xlab('Number of overlaps with other sample')+ylab('Number of regions')+
  geom_text(aes(label=rn), vjust=-0.3, color="black",
            position = position_dodge(0.9), size=3.5)

GR_sum=region_olap(GR,'Sample')
GR_sum_df=data.frame(overlap=0:41,rn=as.numeric(GR_sum))

ggplot(data=GR_sum_df,aes(x=overlap,y=rn))+geom_bar(stat='identity',fill="steelblue")+
  xlab('Number of overlaps with other sample')+ylab('Number of regions')+
  geom_text(aes(label=rn), vjust=-0.3, color="black",
            position = position_dodge(0.9), size=3.5)

GR_H9_MML=GR[GR$Subject=='H9' & GR$ASM=='Yes' & GR$Statistic=='dMML']
GR_HUES64_MML=GR[GR$Subject=='HUES64' & GR$ASM=='Yes' & GR$Statistic=='dMML']

H9_MML_anno=subsetByOverlaps(genomic_features[['promoter']],GR_H9_MML)

HUES64_MML_anno=subsetByOverlaps(genomic_features[['promoter']],GR_HUES64_MML)
HUES64_MML_olap=findOverlaps(HUES64_MML_anno,GR_HUES64_MML)
HUES64_MML_anno$dMML[queryHits(HUES64_MML_olap)]=GR_HUES64_MML$Value[subjectHits(HUES64_MML_olap)]
HUES64_MML_anno=HUES64_MML_anno[order(HUES64_MML_anno$dMML,decreasing=TRUE)]
#Look at het CpG annotation
dat_check_MML_HUES64=dat_check_MML[dat_check_MML$sub=='HUES64']
dat_check_MML_olap_HUES64=findOverlaps(dat_check_MML_HUES64,HUES64_MML_anno,type='within')
HUES64_MML_anno_het=HUES64_MML_anno[subjectHits(dat_check_MML_olap_HUES64)]
dat_check_MML_olap_HUES64=findOverlaps(dat_check_MML_HUES64,HUES64_MML_anno_het,type='within')
HUES64_MML_anno_het$dMML_het[subjectHits(dat_check_MML_olap_HUES64)]=dat_check_MML_HUES64$diff[queryHits(dat_check_MML_olap_HUES64)]
HUES64_MML_anno_het=HUES64_MML_anno_het[order(abs(HUES64_MML_anno_het$dMML_het),decreasing=T)]


subsetByOverlaps(GR_HUES64_MML,HUES64_MML_anno[HUES64_MML_anno$gene_name=='SNURF'])
subsetByOverlaps(GR_HUES64_MML,HUES64_MML_anno[HUES64_MML_anno$gene_name=='PIWIL1'])
#Imprinted region + het CpG
subsetByOverlaps(GR_HUES64_MML,HUES64_MML_anno[HUES64_MML_anno$gene_name=='GNAS'])
subsetByOverlaps(GR_HUES64_MML,HUES64_MML_anno[HUES64_MML_anno$gene_name=='NDN'])
#Only Het CpG
subsetByOverlaps(GR_HUES64_MML,HUES64_MML_anno[HUES64_MML_anno$gene_name=='SOX30'])

#NME analysis

GR_H9_NME=GR[GR$Subject=='H9' & GR$ASM=='Yes' & GR$Statistic=='dNME']
GR_HUES64_NME=GR[GR$Subject=='HUES64'&GR$Tissue=='Embyonic Stem Cell' & GR$ASM=='Yes' & GR$Statistic=='dNME']

H9_NME_anno=subsetByOverlaps(genomic_features[['promoter']],GR_H9_NME)
HUES64_NME_anno=subsetByOverlaps(genomic_features[['promoter']],GR_HUES64_NME)


HUES64_NME_olap=findOverlaps(HUES64_NME_anno,GR_HUES64_NME)
HUES64_NME_anno$dNME[queryHits(HUES64_NME_olap)]=GR_HUES64_NME$Value[subjectHits(HUES64_NME_olap)]
HUES64_NME_anno=HUES64_NME_anno[order(HUES64_NME_anno$dNME,decreasing=TRUE)]
#Look at het CpG annotation
dat_check_NME_HUES64=dat_check_NME[dat_check_NME$sub=='HUES64']
dat_check_NME_olap_HUES64=findOverlaps(dat_check_NME_HUES64,HUES64_NME_anno,type='within')
HUES64_NME_anno_het=HUES64_NME_anno[subjectHits(dat_check_NME_olap_HUES64)]
dat_check_NME_olap_HUES64=findOverlaps(dat_check_NME_HUES64,HUES64_NME_anno_het,type='within')
HUES64_NME_anno_het$dNME_het[subjectHits(dat_check_NME_olap_HUES64)]=dat_check_NME_HUES64$diff[queryHits(dat_check_NME_olap_HUES64)]
HUES64_NME_anno_het=HUES64_NME_anno_het[order(abs(HUES64_NME_anno_het$dNME_het),decreasing=T)]

#H9

H9_NME_olap=findOverlaps(H9_NME_anno,GR_H9_NME)
H9_NME_anno$dNME[queryHits(H9_NME_olap)]=GR_H9_NME$Value[subjectHits(H9_NME_olap)]
H9_NME_anno=H9_NME_anno[order(H9_NME_anno$dNME,decreasing=TRUE)]
#Look at het CpG annotation
dat_check_NME_H9=dat_check_NME[dat_check_NME$sub=='H9']
dat_check_NME_olap_H9=findOverlaps(dat_check_NME_H9,H9_NME_anno,type='within')
H9_NME_anno_het=H9_NME_anno[subjectHits(dat_check_NME_olap_H9)]
dat_check_NME_olap_H9=findOverlaps(dat_check_NME_H9,H9_NME_anno_het,type='within')
H9_NME_anno_het$dNME_het[subjectHits(dat_check_NME_olap_H9)]=dat_check_NME_H9$diff[queryHits(dat_check_NME_olap_H9)]
H9_NME_anno_het=H9_NME_anno_het[order(abs(H9_NME_anno_het$dNME_het),decreasing=T)]


#NME ASM

subsetByOverlaps(GR_HUES64_NME,HUES64_NME_anno[which(HUES64_NME_anno$gene_name=='BST1')])
subsetByOverlaps(GR_HUES64_NME,HUES64_NME_anno[which(HUES64_NME_anno$gene_name=='TFCP2L1')])


subsetByOverlaps(GR_H9_NME,H9_NME_anno[H9_NME_anno$gene_name=='ACVR1B'])

# gene_sum=lapply(unique(queryHits(dMML_ASM_op)),gene_df_gen)
# gene_df_gen<-function(i){
#   gf=genome_features[[11]][i]
#   val_sum=sum(dMML$Value[subjectHits(dMML_ASM_op)[queryHits(dMML_ASM_op)==i]])
#   return(data.frame(gene=gf$gene_name,value_sum=val_sum))
#   
# }