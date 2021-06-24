source('mainFunctions_sub.R')
#create pubMed searching queue based on tissue
#This needs to be run in Windows environment
enhancer_regions=readRDS(enhancer_region_fn)
#table for key items
#In format of  "Dlx3[Title/Abstract] AND craniofacial[Title/Abstract]AND (Human OR Mouse)"
#"forebrain" "heart"     "limb"      "midbrain"  "hindbrain" "liver"     "EFP"

tissue_sel=names(enhancer_regions)
key_terms=c(
  "((craniofacial[Title/Abstract]) OR (face[Title/Abstract]) OR (head[Title/Abstract]))",
  "(forebrain[Title/Abstract])",
  "(heart[Title/Abstract])",
  "(hindbrain[Title/Abstract])",
  "(limb[Title/Abstract])",
  "(liver[Title/Abstract])",
  "(midbrain[Title/Abstract])"
 
  
  
)
names(key_terms)=tissue_sel
enhancer_regions_motif_dNME_all=list()
enhancer_regions_motif_dMML_all=list()
for(ts in tissue_sel){
  cat('Processing:',ts, 'with key:',key_terms[ts],'\n')
  tt1=proc.time()[[3]]
  GO_in=fread(paste0(gene_example_dir,unique(ts),'.csv'))
  enhancer_regions[[ts]]$GO_ID=GO_in[match(enhancer_regions[[ts]]$region,region)]$GO_ID
  enhancer_regions_motif=enhancer_regions[[ts]][(!is.na(dNME_motif)|!is.na(dMML_motif))&(!is.na(GO_ID))]
  enhancer_regions_motif=enhancer_regions_motif[order(dNME_max_pair,decreasing=T)]
  selected_genes=unique(enhancer_regions_motif$gene)
  gene_pub_med=do.call(rbind,lapply(selected_genes,pubmed_rec,keys=key_terms[ts]))
  if(!is.null(gene_pub_med)){
    #maximum countof 300
    PMID_unique=unique(gene_pub_med$PMID)

    gene_pub_med_collapse=gene_pub_med[,list(PMID=paste(PMID,collapse=";"),title=paste(title,collapse=' ; ')),by=gene]
    #Only genes with Pubmed literature
    enhancer_regions_motif_dNME=enhancer_regions_motif[!is.na(dNME_motif)]
    enhancer_regions_motif_dMML=enhancer_regions_motif[!is.na(dMML_motif)]
    enhancer_regions_motif_dNME=cbind(enhancer_regions_motif_dNME,gene_pub_med_collapse[match(enhancer_regions_motif_dNME$gene,gene)])
    enhancer_regions_motif_dMML=cbind(enhancer_regions_motif_dMML,gene_pub_med_collapse[match(enhancer_regions_motif_dMML$gene,gene)])
  }
  else{
  enhancer_regions_motif_dNME=enhancer_regions_motif[!is.na(dNME_motif)]
  enhancer_regions_motif_dMML=enhancer_regions_motif[!is.na(dMML_motif)]
  enhancer_regions_motif_dMML$PMID=NA
  enhancer_regions_motif_dNME$PMID=NA
  enhancer_regions_motif_dMML$title=NA
  enhancer_regions_motif_dNME$title=NA
  }

    cat("Writing result\n")
    write.csv(enhancer_regions_motif_dNME[,list(region,tissue,dNME_motif,cluster,gene,region_type,
                                                dNME_max_pair ,dMML_max_pair,UC_max_time,
                                               UC_max_pair,dNME_cor,
                                               dMML_cor,PMID,title,GO_ID)],
              paste0(region_motif_dir,ts,'_dNME_Pubmed_annotated.csv'),
              row.names = F,quote=F)
    write.csv(enhancer_regions_motif_dMML[,list(region,tissue,dNME_motif,cluster,gene,region_type,
                                                dNME_max_pair ,dMML_max_pair,UC_max_time,
                                                UC_max_pair,dNME_cor,
                                                dMML_cor,PMID,title,GO_ID)],
              paste0(region_motif_dir,ts,'_dMML_Pubmed_annotated.csv'),
              row.names = F,quote=F)
    enhancer_regions_motif_dNME_all[[ts]]=enhancer_regions_motif_dNME
    enhancer_regions_motif_dMML_all[[ts]]=enhancer_regions_motif_dMML
  
 
}
#Identical chack passed in clean run
saveRDS(enhancer_regions_motif_dNME_all,paste0(mouse_motif_dir,'enhancer_regions_motif_dNME_all2.rds'))
saveRDS(enhancer_regions_motif_dMML_all,paste0(mouse_motif_dir,'enhancer_regions_motif_dMML_all2.rds'))