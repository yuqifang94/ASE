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
  results <- motifbreakR(snpList = gr_in_gr, filterp = FALSE,
                         pwmList = HS.SELEX,
                         threshold = 1e-4,
                         method = "ic",
                         bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),verbose=T,
                         BPPARAM = (workers=22))
  
  return(results)
}

motif_num<-function(number_in,pval_cutoff=0.05){
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
  samples_input=c("42_embryonic_stem_cell_single - H9","stem_27_undifferentiated_paired - HUES64",
  "ectoderm_paired - HUES64","endoerm_27_paired - HUES64","mesoderm_23_paired - HUES64","foreskin_keratinocyte_paired - skin03",
  "foreskin_melanocyte_paired - skin03","Adipose_single - STL001","Small_Intestine_single - STL001",
  "Bladder_single - STL001","Gastric_single - STL001","Left_Ventricle_single - STL001",
  "Lung_single - STL001","Psoas_Muscle_single - STL001","Right_Ventricle_single - STL001",
  "Sigmoid_Colon_single - STL001","Spleen_single - STL001","Thymus_single - STL001",
  "Adipose_single - STL002","Adrenal_Gland_single - STL002","Aorta_single - STL002",
  "Esophagus_single - STL002","Gastric_single - STL002","Lung_single - STL002",
  "Ovary_single - STL002","Pancreas_single - STL002","Psoas_Muscle_single - STL002",
  "Small_Intestine_single - STL002","Spleen_single - STL002","Adipose_Tissue_single - STL003",
  "Adrenal_Gland_single - STL003","Aorta_single - STL003","Esophagus_single - STL003","Gastric_single - STL003",
  "Left_Ventricle_single - STL003","Pancreas_single - STL003","Psoas_Muscle_single - STL003","Right_Atrium_single - STL003",
  "Right_Ventricle_single - STL003","Sigmoid_Colon_single - STL003","Small_Intestine_single - STL003","Spleen_single - STL003",
  "Liver_single - STL011","1 - GM12878","2 - GM12878","merged - GM12878","rep1 - H1","rep2 - H1","merged - H1")
  tt1=proc.time()[[1]]
  sp=samples_input[as.numeric(number_in)]
  cat(paste("Processing",sp,'\n'))
  SNP_list=readRDS('../downstream/output/ASM_enrich_meta.rds')
  SNP_list_NME=SNP_list[SNP_list$Statistic=='dNME']
  SNP_list_NME_sig=SNP_list_NME[SNP_list_NME$pvalue<=pval_cutoff & SNP_list_NME$Sample==sp]
  SNP_list_NME_non_sig=SNP_list_NME[SNP_list_NME$pvalue>pval_cutoff& SNP_list_NME$Sample==sp]
  #SNP_list=readRDS(paste('../downstream/motif/',sp,'_dNME_SNP.rds',sep=''))
  cat('Processing significant motif \n')
  motif_sp_sig=motif_break(SNP_list_NME_sig)
  saveRDS(motif_sp_sig,paste('../downstream/motif/',gsub(' - ','_',sp),'_dNME_motif_sig.rds',sep=''))
  cat('Processing nonsignificant motif \n')
  motif_sp_non_sig=motif_break(SNP_list_NME_sig)
  saveRDS(motif_sp_non_sig,paste('../downstream/motif/',gsub(' - ','_',sp),'_dNME_motif_nonsig.rds',sep=''))
  cat("Process finish in",proc.time()[[1]],'\n')
}
args <- commandArgs(trailingOnly = TRUE)
motif_num(args)
