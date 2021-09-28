rm(list=ls())
source('mainFunctions_sub.R')
# reading in the Ken's file -----------------------------------------------
Ken_motif_in='../downstream/input/mouse_analysis/motif_enrichment/'
motif_enrichment_out=data.table()
for(fn in dir(Ken_motif_in,pattern='NME_only.csv')){
  Ken_motif=fread(paste0(Ken_motif_in,fn))
  Ken_motif$tissue=gsub('_.*','',fn)
  Ken_motif=Ken_motif[human_high_NME==TRUE]
  Ken_motif$TF=gsub('.*_','',Ken_motif$motif)
  motif_enrichment_out=rbind(motif_enrichment_out,Ken_motif[,list(motif,TF,tissue)])
  
  
}

motif_human_mouse=motif_enrichment_out[,list(nTS=length(tissue),TS=paste(tissue,collapse = ';')),by=TF]
motif_human_mouse_sub=motif_human_mouse[grepl('heart|forebrain|limb',TS)]
tissue_sel=c('forebrain','heart','limb')
key_terms=c(
  #"((craniofacial[Title/Abstract]) OR (face[Title/Abstract]) OR (head[Title/Abstract]))",
  "((forebrain[Title/Abstract]) OR (brain[Title/Abstract]))",
  "(heart[Title/Abstract])",
  #"(hindbrain[Title/Abstract])",
  "(limb[Title/Abstract])"
  #"(liver[Title/Abstract])",
  #"(midbrain[Title/Abstract])"
  
  
  
)
names(key_terms)=tissue_sel
PMC_out=data.table()
for(ts in tissue_sel){
  cat('Processing:',ts, 'with key:',key_terms[ts],'\n')
  tt1=proc.time()[[3]]
  selected_genes=motif_human_mouse_sub[grepl(ts,TS)]$TF
  gene_pub_med=do.call(rbind,lapply(selected_genes,pubmed_rec,keys=key_terms[ts]))
  if(!is.null(gene_pub_med)){
    #maximum countof 300
    PMID_unique=unique(gene_pub_med$PMID)
    
    if(length(PMID_unique)>300){
    #Counting number of citations
    breaks=ceiling(length(PMID_unique)/300)
    PMID_cut=cut(1:length(PMID_unique),breaks,label=FALSE)
    Sys.sleep(300)
    pmc_ref_count=c()
    for(idx in 1:breaks){

      pmc_ref_count=c(pmc_ref_count,do.call(c,lapply(entrez_summary(db="pubmed", id=PMID_unique[PMID_cut==idx]),function(x) x$pmcrefcount)))
      Sys.sleep(120)
    }
    }else( pmc_ref_count=do.call(c,lapply(entrez_summary(db="pubmed", id=PMID_unique),function(x) x$pmcrefcount)))

    gene_pub_med$pmc_count=pmc_ref_count[gene_pub_med$PMID]
    
    gene_pub_med=cbind(gene_pub_med,motif_human_mouse_sub[match(gene_pub_med$gene,TF)])
    gene_pub_med$tissue_PMID=ts
    PMC_out=rbind(PMC_out,gene_pub_med)
  }

  
  
}
saveRDS(PMC_out,'../downstream/output/mouse_analysis/motif_analysis/PMC_out_enrichment_human_mouse.rds')
