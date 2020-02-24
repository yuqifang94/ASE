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
source("mainFunctions_sub.R")
############################## Motif analysis: Finding if there's other motif affecting dNME/dMML ##################################################
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003","STL011","GM12878","H1")
variant_motif=lapply(subjects,function(x) extractmotif('../downstream/vcfFiles/',x))
#add dMML and dNME to the variant, we need to do it for each sample
gr_allele=readRDS('../downstream/output/GR.all.allele.H1.GM12878.rds')
GR=readRDS('../downstream/output/GR.diff.allele.H1.GM12878.rds')
#variant_motif=readRDS('../downstream/output/variant_motif.rds')
gr_sample=unique(GR$Sample)

variant_values=lapply(gr_sample,function(x) extract_diff_values(x,GR[GR$Sample ==x],variant_motif))
names(variant_values)=NULL
variant_values=do.call('c',variant_values)
saveRDS(variant_values,'../downstream/output/variant_values_all.rds')
names(variant_values)=gr_sample
#Reshape variant motif to allele, keep extra column like dNME,dMML and pval
variant_values_allele=lapply(variant_values,reshape_sample_variant)

#Add MML and NME information
variant_values_allele=lapply(gr_sample, function(x) extract_allele_value(variant_values_allele[[x]],gr_allele[gr_allele$Sample==x]))
names(variant_values_allele)=gr_sample
saveRDS(variant_values_allele,'../downstream/output/variant_values_allele_all.rds')
variant_values_allele=do.call("c",variant_values_allele)
saveRDS(variant_values_allele,'../downstream/output/motif/variant_values_allele_merge_all.rds')
#We need some enrichment test
#Enrichment test: ASM vs nonASM, nucleo vs non nucleo
#Start with 7 nucleotide
#dNME
dNME_7_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_7_X)){
  fisher_motif=motif_enrichment(variant_values_allele,pval_cutoff=0.1,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_7_X')
  dNME_7_motif=rbind(dNME_7_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_7_motif=dNME_7_motif[order(dNME_7_motif$OR,decreasing = T),]
saveRDS(dNME_7_motif,'../downstream/output/dNME_7_motif.rds')
#dMML
dMML_7_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_7_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_7_X')
  dMML_7_motif=rbind(dMML_7_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_7_motif=dMML_7_motif[order(dMML_7_motif$OR,decreasing = T),]
saveRDS(dMML_7_motif,'../downstream/output/motif/dMML_7_motif.rds')

#Start with 5 nucleotide
#dNME
dNME_5_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_5_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_5_X')
  dNME_5_motif=rbind(dNME_5_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_5_motif=dNME_5_motif[order(dNME_5_motif$OR,decreasing = T),]
saveRDS(dNME_5_motif,'../downstream/output/motif/dNME_5_motif.rds')
#dMML
dMML_5_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_5_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_5_X')
  dMML_5_motif=rbind(dMML_5_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_5_motif=dMML_5_motif[order(dMML_5_motif$OR,decreasing = T),]
saveRDS(dMML_5_motif,'../downstream/output/motif/dMML_5_motif.rds')

#Start with 3 nucleotide
#dNME
dNME_3_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_3_X)){
  fisher_motif=motif_enrichment(variant_values_allele,pval_cutoff=0.1,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_3_X')
  dNME_3_motif=rbind(dNME_3_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_3_motif=dNME_3_motif[order(dNME_3_motif$OR,decreasing = T),]
saveRDS(dNME_3_motif,'../downstream/output/motif/dNME_3_motif.rds')
#dMML
dMML_3_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_3_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_3_X')
  dMML_3_motif=rbind(dMML_3_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_3_motif=dMML_3_motif[order(dMML_3_motif$OR,decreasing = T),]
saveRDS(dMML_3_motif,'../downstream/output/motif/dMML_3_motif.rds')


#Start with 9 nucleotide
#dNME
dNME_9_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_9_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_9_X')
  dNME_9_motif=rbind(dNME_9_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_9_motif=dNME_9_motif[order(dNME_9_motif$OR,decreasing = T),]
saveRDS(dNME_9_motif,'../downstream/output/motif/dNME_9_motif.rds')
#dMML
dMML_9_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_9_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_9_X')
  dMML_9_motif=rbind(dMML_9_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_9_motif=dMML_9_motif[order(dMML_9_motif$OR,decreasing = T),]
saveRDS(dMML_9_motif,'../downstream/output/motif/dMML_9_motif.rds')


#Start with 11 nucleotide
#dNME
dNME_11_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_11_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_11_X')
  dNME_11_motif=rbind(dNME_11_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_11_motif=dNME_11_motif[order(dNME_11_motif$OR,decreasing = T),]
saveRDS(dNME_11_motif,'../downstream/output/motif/dNME_11_motif.rds')
#dMML
dMML_11_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_11_X)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_11_X')
  dMML_11_motif=rbind(dMML_11_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_11_motif=dMML_11_motif[order(dMML_11_motif$OR,decreasing = T),]
saveRDS(dMML_11_motif,'../downstream/output/motif/dMML_11_motif.rds')

###################Without X###########################
#Start with 7 nucleotide
#dNME
dNME_7_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_7)){
  fisher_motif=motif_enrichment(variant_values_allele,pval_cutoff=0.1,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_7')
  dNME_7_motif=rbind(dNME_7_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_7_motif=dNME_7_motif[order(dNME_7_motif$OR,decreasing = T),]
saveRDS(dNME_7_motif,'../downstream/output/motif/dNME_7_motif_noX.rds')
#dMML
dMML_7_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_7)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_7')
  dMML_7_motif=rbind(dMML_7_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_7_motif=dMML_7_motif[order(dMML_7_motif$OR,decreasing = T),]
saveRDS(dMML_7_motif,'../downstream/output/motif/dMML_7_motif_noX.rds')

#Start with 5 nucleotide
#dNME
dNME_5_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_5)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_5')
  dNME_5_motif=rbind(dNME_5_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_5_motif=dNME_5_motif[order(dNME_5_motif$OR,decreasing = T),]
saveRDS(dNME_5_motif,'../downstream/output/motif/dNME_5_motif_noX.rds')
#dMML
dMML_5_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_5)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_5')
  dMML_5_motif=rbind(dMML_5_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_5_motif=dMML_5_motif[order(dMML_5_motif$OR,decreasing = T),]
saveRDS(dMML_5_motif,'../downstream/output/motif/dMML_5_motif_noX.rds')

#Start with 3 nucleotide
#dNME
dNME_3_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_3)){
  fisher_motif=motif_enrichment(variant_values_allele,pval_cutoff=0.1,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_3')
  dNME_3_motif=rbind(dNME_3_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_3_motif=dNME_3_motif[order(dNME_3_motif$OR,decreasing = T),]
saveRDS(dNME_3_motif,'../downstream/output/motif/dNME_3_motif_noX.rds')
#dMML
dMML_3_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_3)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_3')
  dMML_3_motif=rbind(dMML_3_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_3_motif=dMML_3_motif[order(dMML_3_motif$OR,decreasing = T),]
saveRDS(dMML_3_motif,'../downstream/output/motif/dMML_3_motif_noX.rds')


#Start with 9 nucleotide
#dNME
dNME_9_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_9)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_9')
  dNME_9_motif=rbind(dNME_9_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_9_motif=dNME_9_motif[order(dNME_9_motif$OR,decreasing = T),]
saveRDS(dNME_9_motif,'../downstream/output/motif/dNME_9_motif_noX.rds')
#dMML
dMML_9_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_9)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_9')
  dMML_9_motif=rbind(dMML_9_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_9_motif=dMML_9_motif[order(dMML_9_motif$OR,decreasing = T),]
saveRDS(dMML_9_motif,'../downstream/output/motif/dMML_9_motif_noX.rds')


#Start with 11 nucleotide
#dNME
dNME_11_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_11)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dNME_pval',motif=mf,motif_type = 'nucleo_11')
  dNME_11_motif=rbind(dNME_11_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dNME_11_motif=dNME_11_motif[order(dNME_11_motif$OR,decreasing = T),]
saveRDS(dNME_11_motif,'../downstream/output/motif/dNME_11_motif_noX.rds')
#dMML
dMML_11_motif=data.frame(motif=NULL,OR=NULL,pval=NULL)
for (mf in unique(variant_values_allele$nucleo_11)){
  fisher_motif=motif_enrichment(variant_values_allele,p_stat='dMML_pval',motif=mf,motif_type = 'nucleo_11')
  dMML_11_motif=rbind(dMML_11_motif,data.frame(motif=mf,OR=fisher_motif$estimate,pval=fisher_motif$p.value))
}
dMML_11_motif=dMML_11_motif[order(dMML_11_motif$OR,decreasing = T),]
saveRDS(dMML_11_motif,'../downstream/output/motif/dMML_11_motif_noX.rds')

