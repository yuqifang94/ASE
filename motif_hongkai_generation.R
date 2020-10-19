source('mainFunctions_sub.R')
GR_merge=readRDS(GR_merge_file)
genomic_features=readRDS(genomic_features_file)
GR_merge$fn=paste(GR_merge$Subject,GR_merge$tissue,sep='_')
#Imprinted_Genes <- as.data.frame(read_excel("../downstream/input/Imprinted Genes.xlsx"))
#Imprinted_Genes=Imprinted_Genes[Imprinted_Genes$`Expressed Allele`%in% c('Maternal','Paternal'),]
in_dir='../downstream/data/allele_agnostic/'
out_dir='../downstream/output/motif_hongkai/'
#sp='stem_27_undifferentiated_paired - HUES64'
log_output=data.frame()
for (sp in unique(GR_merge$Sample)){

  GR_merge_sp=GR_merge[GR_merge$Sample==sp]
  fn=unique(GR_merge$fn[GR_merge$Sample==sp])
  #allele-agnositc model
  fn_nme=paste(in_dir,fn,'_phased_allele_agnostic_nme.bedGraph',sep='')
  if( file.exists(fn_nme)){
  ag_nme=import.bedGraph(fn_nme)
  
  elementMetadata(ag_nme)=elementMetadata(ag_nme)[,'score']
  colnames(elementMetadata(ag_nme))='NME'
  ag_nme=ag_nme[!is.infinite(ag_nme$NME)]
  if(all(seqlevels(ag_nme)==gsub('chr','',seqlevels(ag_nme)))){seqlevels(ag_nme)=paste('chr',seqlevels(ag_nme),sep='')}
  ag_nme=dist_calc(ag_nme,genomic_features$TSS)
  ag_nme=ag_nme[!ag_nme$gene %in% Imprinted_Genes$Gene]
  export_bedGraph(ag_nme[order(ag_nme$NME,decreasing = T)][1:10000],paste(out_dir,fn,'-top10K_allele_agnostic_NME.bedGraph',sep=''))
  export_bedGraph(ag_nme[order(ag_nme$NME,decreasing = F)][1:10000],paste(out_dir,fn,'-bottom10K_allele_agnostic_NME.bedGraph',sep=''))
  ag_nme_TSS=ag_nme[abs(ag_nme$dist)<=2000]
  if(length(ag_nme_TSS)<10000){
    top_length=round(length(ag_nme_TSS)/3,digits=-3)
    TSS_10K=FALSE
  }else{
      top_length=10000
      TSS_10K=TRUE}
  export_bedGraph(ag_nme_TSS[order(ag_nme_TSS$NME,decreasing = T)][1:top_length],paste(out_dir,fn,'-top10K_allele_agnostic_NME_2k_TSS.bedGraph',sep=''))
  export_bedGraph(ag_nme_TSS[order(ag_nme_TSS$NME,decreasing = F)][1:top_length],paste(out_dir,fn,'-bottom10K_allele_agnostic_NME_2k_TSS.bedGraph',sep=''))

  #dNME model
  GR_merge_sp$score=GR_merge_sp$dNME
  cat(sp,':',sum(GR_merge_sp$dNME_pval<=pval_cutoff),'\n')
  dNME_regions=GR_merge_sp[GR_merge_sp$dNME_pval<=pval_cutoff]
  non_dNME_regions=GR_merge_sp[GR_merge_sp$dNME_pval>=0.8&GR_merge_sp$dNME<=min(dNME_regions$dNME)/2]
  non_dNME_regions_idx=sample(1:length(non_dNME_regions),length(dNME_regions),replace = FALSE)
  export.bedGraph(dNME_regions,paste(out_dir,fn,'-dNME.bedGraph',sep=''))
  export.bedGraph(non_dNME_regions[non_dNME_regions_idx],paste(out_dir,fn,'-non_dNME.bedGraph',sep=''))
  log_output=rbind(log_output,data.frame(diff_length=length(GR_merge_sp),
                                         dNME_length=sum(GR_merge_sp$dNME_pval<=pval_cutoff),
                                         TSS_10K=TSS_10K,
                                         ag_length=length(ag_nme),
                                         min_ag_top10K=min(ag_nme$NME[order(ag_nme$NME,decreasing = T)][1:10000]),
                                         max_ag_bot10K=max(ag_nme$NME[order(ag_nme$NME,decreasing = F)][1:10000]),
                                         min_diff_top=min(GR_merge_sp$dNME[GR_merge_sp$dNME_pval<=pval_cutoff]),
                                         max_diff_bot=max(non_dNME_regions$dNME[non_dNME_regions_idx]),
                                         sample=sp))
  }
}
#Think about how to use the jordi's new method on my analysis: calculating the differential statistics on given genes/cluster of genes etc
export_bedGraph<-function(gr_in,fn){
  gr_in=as.data.frame(gr_in)
  colnames(gr_in)=c('chr','start','end','width','strand','NME','dist_to_TSS','gene_at_TSS','dist_round')
  gr_in=gr_in[,c('chr','start','end','NME','dist_to_TSS','gene_at_TSS')]
  gr_in$start=gr_in$start-1
  write.table(gr_in,file=fn,quote=F,sep='\t',row.names = F,col.names = F)
}
