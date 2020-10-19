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
if (!requireNamespace("BiocParallel", quietly = TRUE))
{BiocManager::install("BiocParallel")}
library(BiocParallel)
motif_break<-function(gr_in,motif_list){
  gr_in_gr=granges(gr_in)
  strand(gr_in_gr)='*'
  attributes(gr_in_gr)$genome.package="BSgenome.Hsapiens.UCSC.hg19"
  #Make sure ref and alt agree with BS genome #may not necessary
  ref_BS=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,gr_in_gr))
  alt_BS=gr_in$ALT
  switch_nucleotide=which(gr_in$REF!=ref_BS)
  print(switch_nucleotide)
  alt_BS[switch_nucleotide]=gr_in$REF[switch_nucleotide]
  #Get the necessary information for program running
  gr_in_gr$REF=unlist(DNAStringSetList(ref_BS))
  gr_in_gr$ALT=unlist(DNAStringSetList(alt_BS))
  names(gr_in_gr)=gr_in$snpId
  results <- motifbreakR(snpList = gr_in_gr, filterp = TRUE,
                         pwmList = motif_list,
                         threshold = 1e-4,
                         method = "log",
                         bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),verbose=T,
                         BPPARAM = MulticoreParam(workers=22,progressbar=TRUE))
  
  return(results)
}

#Using all SNP in vcf to do motif
HS.SELEX=subset (MotifDb, organism=='Hsapiens'&dataSource=="jolma2013")
variant_in=readRDS('../downstream/input/variant_HetCpG_new.rds')
names(variant_in)=NULL
variant_in=do.call('c',variant_in)
elementMetadata(variant_in)=elementMetadata(variant_in)[1:4]
enchancer_DNase_gr=readRDS("../downstream/input/enchancer_DNase.rds")
prom_DNase_gr=readRDS("../downstream/input/prom_DNase.rds")
DNase_all=c(prom_DNase_gr,enchancer_DNase_gr)
variant_in_dnase=subsetByOverlaps(variant_in,DNase_all)
rm(variant_in)
rm(DNase_all)
motif_gene <- readRDS("../downstream/output/motif_break_all_unique_SNP.rds")
olap_done=findOverlaps(variant_in_dnase,motif_gene)
motif_out_dnase=motif_break(unique(variant_in_dnase[-queryHits(olap_done)]),HS.SELEX)
saveRDS(motif_out_dnase,'../downstream/output/motif_out_new_dnase.rds')
olap_done=findOverlaps(variant_in,motif_gene)
motif_out=motif_break(unique(variant_in[-queryHits(olap_done)]),HS.SELEX)
saveRDS(motif_out,'../downstream/output/motif_out_new.rds')
variant1=readRDS('../downstream/output/ASM_enrich_meta2.rds')
variant2=readRDS('../downstream/output/ASM_enrich_meta_new.rds')
olap=findOverlaps(variant1,variant2)
variant2=unique(variant2[-subjectHits(olap)])
elementMetadata(variant2)=elementMetadata(variant2)[1:4]
variant2=subsetByOverlaps(variant2,DNase_all)
motif_out=motif_break(variant2,HS.SELEX)
saveRDS(motif_out_dnase,'../downstream/output/motif_out_new_v2_diff.rds')

genomic_features <- readRDS("../downstream/input/genomic_features2020.rds")
HS.SELEX=subset (MotifDb, organism=='Hsapiens'&dataSource=="jolma2013")
CTCF= query(HS.SELEX,"CTCF")
gene_range_body_promoter=c(genomic_features$promoter[genomic_features$promoter$gene_name%in%HS.SELEX@elementMetadata$geneSymbol],
                           genomic_features$`gene body`[genomic_features$`gene body`$gene_name%in%HS.SELEX@elementMetadata$geneSymbol])
gene_range_body_promoter=c(genomic_features$promoter, genomic_features$`gene body`)
pval_cutoff=0.1
variant_HetCpG_meta=readRDS('../downstream/output/ASM_enrich_meta.rds')
variant_HetCpG_NME=variant_HetCpG_meta[!variant_HetCpG_meta$Subject=='GM12878' & variant_HetCpG_meta$Statistic=='dNME']
variant_HetCpG_NME=variant_HetCpG_meta[!variant_HetCpG_meta$Subject=='GM12878' & variant_HetCpG_meta$Statistic=='dMML']
GR_merge=readRDS("../downstream/output/gr_merge.H1.GM12878.rds")
GR_merge_exclude_GM=GR_merge[!GR_merge$Sample%in%c('merged - GM12878','1 - GM12878','2 - GM12878')]
GR_merge_exclude_GM=GR_merge_exclude_GM[GR_merge_exclude_GM$N>1]
variant_sub_dNME=GRanges()
for(sp in unique(variant_HetCpG_NME$Sample)){
  variant_sub_dNME=c(variant_sub_dNME,
                     subsetByOverlaps(variant_HetCpG_NME[variant_HetCpG_NME$Sample==sp],GR_merge_exclude_GM[GR_merge_exclude_GM$Sample==sp]))
}
#All SNP
rm(variant_HetCpG_meta)
rm(GR_merge)
rm(GR_merge_exclude_GM)
motif=motif_break(unique(variant_sub_dNME,HS.SELEX))#Maybe we can downscale the SNPs
motif_sig=motif_break(unique(variant_sub_dNME[variant_sub_dNME$pvalue<=pval_cutoff]))#Significant SNPs
SNP_sub=unique(subsetByOverlaps(variant_HetCpG_meta,gene_range_body_promoter,maxgap = 5000))
motif_sub=motif_break(SNP_sub,HS.SELEX)
#Downscale SNPs using SNPs near those gene
#Find SNPs that have at least one dNME event in all samples
variant_HetCpG_NME_ASM=unique(variant_HetCpG_NME[variant_HetCpG_NME$pvalue<=pval_cutoff])
#Remove those SNP
olap=findOverlaps(variant_HetCpG_NME,variant_HetCpG_NME_ASM)
variant_HetCpG_NME_no_ASM
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
