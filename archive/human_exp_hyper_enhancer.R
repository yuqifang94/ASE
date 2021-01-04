# Enhancer with nearest gene ----------------------------------------------
ah = AnnotationHub()
ENCODE_name=ENCODE_to_sample(unique(GR_merge$Sample))
chromHMM_dMML_all_ls=list()
ah_gr=GRanges()
chromHMM=list()
#Do it for all available data
suppressMessages({for (sp in ENCODE_name$sample[!is.na(ENCODE_name$ENCODE)]){
  ah_num=names(query(ah, c("chromhmmSegmentations", ENCODE_name$ENCODE[ENCODE_name$sample==sp])))
  chromHMM[[sp]]=ah[[ah_num]]
  
  # }
}
})
#Look for "7_Enh"  "6_EnhG"    "12_EnhBiv"
chromHMM=lapply(chromHMM,function(x) x[x$abbr %in% c("7_Enh","6_EnhG")])
JASPAR_motif=readRDS('../downstream/output/motif_JASPAR_hg19.rds')
hyper_var_all_enhancer_NME=GRanges()
hyper_var_all_enhancer_MML=GRanges()
for (sp in names(chromHMM)){
  hyper_var_all_enhancer_NME=c(hyper_var_all_enhancer_NME,subsetByOverlaps(hyper_var_all$NME_hypervar_calc[hyper_var_all$NME_hypervar_calc$Sample==sp],chromHMM[[sp]]))
  hyper_var_all_enhancer_MML=c(hyper_var_all_enhancer_MML,subsetByOverlaps(hyper_var_all$MML_meanvar_calc[hyper_var_all$MML_meanvar_calc$Sample==sp],chromHMM[[sp]]))
  cat("Percent enhancer covered for",sp,":",
      length(subsetByOverlaps(chromHMM[[sp]],hyper_var_all$MML_meanvar_calc[hyper_var_all$MML_meanvar_calc$Sample==sp]))/length(chromHMM[[sp]]),'\n')
}
hyper_var_all_enhancer_MML=hyper_var_all_enhancer_MML[abs(hyper_var_all_enhancer_MML$dist)>=3000]
hyper_var_all_enhancer_NME=hyper_var_all_enhancer_NME[abs(hyper_var_all_enhancer_NME$dist)>=3000]
MML_region_select=data.table(MML=hyper_var_all_enhancer_MML$score,mean_exp=hyper_var_all_enhancer_MML$exp_stat)
MML_region_select$region=paste0(seqnames(hyper_var_all_enhancer_MML),':',start(hyper_var_all_enhancer_MML),'-',end(hyper_var_all_enhancer_MML))
region_enough_sample=table(MML_region_select$region)

MML_region_select=MML_region_select[region %in% names(region_enough_sample[region_enough_sample>=20]),
                                    list(cor=cor(MML,mean_exp),pval=cor.test(MML,mean_exp)$p.value),by=list(region)]
MML_region_select$qval=p.adjust(MML_region_select$pval,method="BH")
NME_hypervar_calc=data.table(NME=hyper_var_all_enhancer_NME$score,hypervarquant=hyper_var_all_enhancer_NME$hypervarquant001,
                             Sample=hyper_var_all_enhancer_NME$Sample,gene_exp=hyper_var_all_enhancer_NME$exp_stat)
NME_hypervar_calc$region=paste0(seqnames(hyper_var_all_enhancer_NME),':',start(hyper_var_all_enhancer_NME),'-',end(hyper_var_all_enhancer_NME))

NME_hypervar_calc_sel=NME_hypervar_calc[region %in% MML_region_select[ qval<=0.1]$region]

NME_hypervar_calc=NME_hypervar_calc[order(NME_hypervar_calc$cor,decreasing = F),]
NME_hypervar_calc$Sample=factor(NME_hypervar_calc$Sample,levels = unique(NME_hypervar_calc$Sample))
NME_hypervar_calc[NME_hypervar_calc$region %in% NME_hypervar_calc_agg[qvalue<=0.1&cor>0]$region]

#CTCF?
motif_cor=lapply(names(JASPAR_motif),function(x){
  NME_motif=subsetByOverlaps(hyper_var_all_enhancer_MML,JASPAR_motif[[x]])
  NME_motif_dt=data.table(NME=NME_motif$score,exp_hypervar=NME_motif$exp_stat,N=NME_motif$N,
                          sample=NME_motif$Sample,gene=NME_motif$gene)
  
  NME_motif_dt_sp=NME_motif_dt[,list(cor=cor(NME,exp_hypervar)),by=sample]
  NME_motif_dt_sp$motif=x
  NME_motif_dt_gene=NME_motif_dt[,list(cor=cor(NME,exp_hypervar)),by=gene]
  NME_motif_dt_gene$motif=x
  return(list(cor.test(NME_motif_dt$NME,NME_motif_dt$exp_hypervar),NME_motif_dt_sp,NME_motif_dt_gene))
})

#Max correlation=1.2
names(motif_cor)=names(JASPAR_motif)
motif_correlation=lapply(motif_cor,function(x) x[[1]]$estimate)
motif_df=data.table(motif=names(motif_correlation),correlation=unlist(motif_correlation))
motif_df$pvalue=unlist(lapply(motif_cor,function(x) x[[1]]$p.value))
hist(unlist(motif_correlation),xlab="NME hypervariability enrichment in motif")
motif_cor_motif=unlist(lapply(motif_cor,function(x) x[[1]]$motif[x[[1]]$cor>=0.5]))
motif_cor_gene=unlist(lapply(motif_cor,function(x) x[[2]]$gene[x[[2]]$cor>=0.9]))
motif_cor_gene[!is.na(motif_cor_gene)]


# FANTOM enhancer ---------------------------------------------------------
#An atlas of active enhancers across human cell types and tissues
#
FANTOM=import.bed('../downstream/input/enhancer_tss_associations_hg19.bed')
FANTOM$gene=unlist(lapply(strsplit(FANTOM$name,';'),function(x) x[[3]]))
FANTOM$enhancer=unlist(lapply(strsplit(FANTOM$name,';'),function(x) x[[1]]))
FANTOM$cor=as.numeric(sub('R:','',unlist(lapply(strsplit(FANTOM$name,';'),function(x) x[[4]]))))
FANTOM$FDR=as.numeric(sub('FDR:','',unlist(lapply(strsplit(FANTOM$name,';'),function(x) x[[5]]))))
FANTOM_enhancer=convert_GR(FANTOM$enhancer)
mcols(FANTOM_enhancer)=mcols(FANTOM)[,c("gene","cor","FDR")]

FANTOM=fread('../downstream/input/hg19_enhancer_promoter_correlations_distances_organ.txt')
FANTOM$promoter=sub('\\..','-',FANTOM$promoter)
FANTOM_promoter_gene=convert_GR(FANTOM$promoter)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- GenomicFeatures::genes(txdb)
TSS<-promoters(genes,upstream=0,downstream=0)
TSS$gene_name=AnnotationDbi::select(Homo.sapiens,key=as.character(TSS$gene_id),
                                    keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
FANTOM_promoter_gene=dist_calc(FANTOM_promoter_gene,TSS)
FANTOM$gene=FANTOM_promoter_gene$gene
FANTOM$dist_promoter=FANTOM_promoter_gene$dist
FANTOM_enhancer=convert_GR(FANTOM$enhancer)
mcols(FANTOM_enhancer)=FANTOM[,c("gene","correlation","FDR")]
NME_in=readRDS("../downstream/input/allele_agnostic_hg19_cov10_3kb_FANTOM_NME.rds")
MML_in=readRDS("../downstream/input/allele_agnostic_hg19_cov10_3kb_FANTOM_MML.rds")
NME_olap=findOverlaps(NME_in,FANTOM_enhancer,maxgap = 0)
NME_in=NME_in[queryHits(NME_olap)]
NME_in$gene=FANTOM_enhancer$gene[subjectHits(NME_olap)]
MML_olap=findOverlaps(MML_in,FANTOM_enhancer,maxgap = 0)
MML_in=MML_in[queryHits(MML_olap)]
MML_in$gene=FANTOM_enhancer$gene[subjectHits(MML_olap)]
NME_calc=GRanges()
MML_calc=GRanges()


#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){
    
    sp_hyper_var=read_hypervar(hyper_var_file)
    #scRNA_result=rbind(scRNA_result,sp_hyper_var)
    sp_hyper_var$log2mean=log2(sp_hyper_var$mean)
    #Add hypervaribility inforamtion
    #Add hypervaribility inforamtion
    NME_in_sp=NME_in[NME_in$Sample==sp]
    mcols(NME_in_sp)=cbind(mcols(NME_in_sp),sp_hyper_var[match(NME_in_sp$gene,sp_hyper_var$gene_name)])
    NME_calc=c(NME_calc,NME_in_sp)
    MML_in_sp=MML_in[MML_in$Sample==sp]
    mcols(MML_in_sp)=cbind(mcols(MML_in_sp),sp_hyper_var[match(MML_in_sp$gene,sp_hyper_var$gene_name)])
    MML_calc=c(MML_calc,MML_in_sp)
    
  }else{cat("file not exist for:",sp,'\n')}
}

NME_calc=NME_calc[!is.na(NME_calc$gene_name)]
NME_calc_dt=as.data.table(mcols(NME_calc))
NME_calc_dt=NME_calc_dt[,list(mean_NME=mean(NME)),by=list(Sample,gene,mean,hypervar_logvar)]

NME_calc_dt[,list(cor=cor.test(mean_NME,hypervar_logvar,method="pearson")$estimate,
                  pvalue=cor.test(mean_NME,hypervar_logvar,method="pearson")$p.value),by=list(Sample)]

#CTCF?
motif_cor=lapply(names(JASPAR_motif),function(x){
  NME_motif=subsetByOverlaps(NME_calc,JASPAR_motif[[x]])
  NME_motif_dt=data.table(NME=NME_motif$NME,exp_hypervar=NME_motif$hypervar_logvar,N=NME_motif$N,
                          sample=NME_motif$Sample,gene=NME_motif$gene)
  
  NME_motif_dt_sp=NME_motif_dt[,list(cor=cor(NME,exp_hypervar)),by=sample]
  NME_motif_dt_sp$motif=x
  NME_motif_dt_gene=NME_motif_dt[,list(cor=cor(NME,exp_hypervar)),by=gene]
  NME_motif_dt_gene$motif=x
  return(list(cor.test(NME_motif_dt$NME,NME_motif_dt$exp_hypervar),NME_motif_dt_sp,NME_motif_dt_gene))
})
names(motif_cor)=names(JASPAR_motif)
motif_correlation=lapply(motif_cor,function(x) x[[1]]$estimate)
motif_df=data.table(motif=names(motif_correlation),correlation=unlist(motif_correlation))
motif_df$pvalue=unlist(lapply(motif_cor,function(x) x[[1]]$p.value))

NME_hypervar_calc_dt=as.data.table(mcols(liver))
NME_hypervar_calc_dt=NME_hypervar_calc_dt[!is.na(gene_name)]
NME_hypervar_calc_dt=NME_hypervar_calc_dt[,list(mean_NME=mean(NME),total_N=sum(N)),by=list(Sample,gene,hypervar_logvar)]
NME_hypervar_calc_dt[,list(cor=cor.test(mean_NME,hypervar_logvar,method="pearson")$estimate,pvalue=cor.test(mean_NME,hypervar_logvar,method="pearson")$p.value),by=list(Sample)]
liver_enhancer=import.bed('../downstream/input/liver_enhancer.bed')
liver=subsetByOverlaps(NME_calc[NME_calc$Sample%in%c("Liver_single - STL011")],liver_enhancer)

saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc,
             MML_meanvar_calc=MML_meanvar_calc,
             NME_meanvar_calc=NME_meanvar_calc),'../downstream/output/allele_agnostic_var_FANTOM.rds')


