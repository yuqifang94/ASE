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
saveRDS(enhancer_regions_motif_dNME_all,enhancer_motif_all_dNME_fn)
saveRDS(enhancer_regions_motif_dMML_all,enhancer_motif_all_dMML_fn)
enhancer_regions_motif_dNME_all=readRDS(enhancer_motif_all_dNME_fn)
enhancer_regions_motif_dMML_all=readRDS(enhancer_motif_all_dMML_fn)
#creating merged csv file #Supplmentary data
enhancer_regions_motif_dNME_all_dNME=do.call(rbind,enhancer_regions_motif_dNME_all)
enhancer_regions_motif_dNME_all_dNME=enhancer_regions_motif_dNME_all_dNME[region_type=='NME only',list(tissue,dNME_motif,cluster,gene,dNME_max_pair,dNME_cor,dMML_max_pair,dMML_cor,UC_max_time,PMID)]
colnames(enhancer_regions_motif_dNME_all_dNME)=c('Tissue','Motif','Cluster','Gene','dNME','dNME-UC correlation','dMML','dMML-UC correlation','Stage','PMID')
enhancer_regions_motif_dMML_all_dMML=do.call(rbind,enhancer_regions_motif_dMML_all)
enhancer_regions_motif_dMML_all_dMML=enhancer_regions_motif_dMML_all_dMML[region_type=='MML only',list(tissue,dMML_motif,cluster,gene,dNME_max_pair,dNME_cor,dMML_max_pair,dMML_cor,UC_max_time,PMID)]
colnames(enhancer_regions_motif_dMML_all_dMML)=c('Tissue','Motif','Cluster','Gene','dNME','dNME-UC correlation','dMML','dMML-UC correlation','Stage','PMID')
write.csv(enhancer_regions_motif_dNME_all_dNME,paste0(mouse_motif_dir,'merged_motif_dNME.csv',row.names = F)
write.csv(enhancer_regions_motif_dMML_all_dMML,paste0(mouse_motif_dir,'merged_motif_dMML.csv',row.names = F)
# Reformat motif from Ken -------------------------------------------------

dMML_motifs=data.table()
for(fn in dir(Ken_motif_folder,pattern='dMML')){
  ts=gsub('_.*','',fn)
  csv_in=fread(paste0(Ken_motif_folder,fn))
  csv_in$tissue=ts
  csv_in$motif=gsub('.*_','',csv_in$motif)
  dMML_motifs=rbind(dMML_motifs,csv_in[,list(tissue,motif,human_high_NME,residual,log_OR_dMML,normalized_log_OR_dMML,FDR_dMML)])
}
write.csv(dMML_motifs,paste0(mouse_motif_dir,'Ken_dMML_all.csv')
dNME_motifs=data.table()
for(fn in dir(Ken_motif_folder,pattern='dNME')){
  ts=gsub('_.*','',fn)
  csv_in=fread(paste0(dir_Ken,fn))
  csv_in$tissue=ts
  csv_in$motif=gsub('.*_','',csv_in$motif)
  dNME_motifs=rbind(dNME_motifs,csv_in[,list(tissue,motif,human_high_NME,residual,log_OR_dNME,normalized_log_OR_dNME,FDR_dNME)])
}
write.csv(dNME_motifs,paste0(mouse_motif_dir,'Ken_dNME_all.csv')