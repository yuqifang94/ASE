source('mainFunctions_sub.R')
RNA_seq_process<-function(dir_in,fds,condition_rep=3){
  files=paste(dir_in,fds,"/t_data.ctab",sep="")
  tmp=fread(files[1])
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
GR_merge=readRDS(GR_merge_file)
###Reading in RNA-seq data
#"../downstream/input/RNA_seq/"
RNA_dir='../downstream/input/human_analysis/MML_expression/RNA_seq/'
fds_H1=c("H1_rep1_genome1","H1_rep2_genome1","H1_rep3_genome1","H1_rep1_genome2","H1_rep2_genome2","H1_rep3_genome2")
res_RNA_H1=RNA_seq_process(RNA_dir,fds_H1)
H1_df=RNA_df(GR_merge[GR_merge$Subject=='H1'],res_RNA_H1[!is.na(res_RNA_H1$log2FoldChange),])

fds_HUES64_ectoderm=c("HUES64_ectoderm_genome1","HUES64_ectoderm_genome2")
res_RNA_HUES64_ectoderm=RNA_seq_process(RNA_dir,fds_HUES64_ectoderm,condition_rep = 1)
HUES64_ectoderm_df=RNA_df(GR_merge[GR_merge$Sample=='ectoderm_paired - HUES64'],res_RNA_HUES64_ectoderm[!is.na(res_RNA_HUES64_ectoderm$log2FoldChange),])

fds_HUES64_endoerm=c("HUES64_endoderm_genome1","HUES64_endoderm_genome2")
res_RNA_HUES64_endoerm=RNA_seq_process(RNA_dir,fds_HUES64_endoerm,condition_rep = 1)
HUES64_endoerm_df=RNA_df(GR_merge[GR_merge$Sample=='endoerm_27_paired - HUES64'],res_RNA_HUES64_endoerm[!is.na(res_RNA_HUES64_endoerm$log2FoldChange),])

fds_HUES64_stem=c("HUES64_stem_genome1","HUES64_stem_genome2")
res_RNA_HUES64_stem=RNA_seq_process(RNA_dir,fds_HUES64_stem,condition_rep = 1)
HUES64_stem=RNA_df(GR_merge[GR_merge$Sample=='stem_27_undifferentiated_paired - HUES64'],res_RNA_HUES64_stem[!is.na(res_RNA_HUES64_stem$log2FoldChange),])

fds_HUES64_mesoderm=c("HUES64_mesoderm_genome1","HUES64_mesoderm_genome2")
res_RNA_HUES64_mesoderm=RNA_seq_process(RNA_dir,fds_HUES64_mesoderm,condition_rep = 1)
HUES64_mesoderm=RNA_df(GR_merge[GR_merge$Sample=='mesoderm_23_paired - HUES64'],res_RNA_HUES64_mesoderm[!is.na(res_RNA_HUES64_mesoderm$log2FoldChange),])
###At dMML hap, large correlation with RNA-seq, data too few not significant
RNA_seq_df=rbind(H1_df,HUES64_ectoderm_df,HUES64_endoerm_df,HUES64_stem,HUES64_mesoderm)
cor.test(RNA_seq_df$dMML_relative[RNA_seq_df$dMML_pval<=pval_cutoff],RNA_seq_df$ASE_log2FC[RNA_seq_df$dMML_pval<=pval_cutoff])
ggplot(RNA_seq_df[RNA_seq_df$dMML_pval<=pval_cutoff,],aes(x=dMML_relative,y=ASE_log2FC))+geom_smooth(method='lm')+geom_point()+
  xlab('relative dMML (genome2-genome1)')+ylab('log2Fc (genome2/genome1)')

cor.test(RNA_seq_df$dNME[RNA_seq_df$dNME_pval<=pval_cutoff],abs(RNA_seq_df$ASE_log2FC[RNA_seq_df$dNME_pval<=pval_cutoff]))
ggplot(RNA_seq_df[RNA_seq_df$dNME_pval<=pval_cutoff,],aes(x=dNME_relative,y=ASE_log2FC))+geom_smooth(method='lm')+geom_point()+
  xlab('relative dNME (genome2-genome1)')+ylab('log2Fc (genome2/genome1)')