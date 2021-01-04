#GO annotation
GO_run<-function(gl,back,motif){
  #GO object generation
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                                  annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")})
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  #Filter out the GO terms
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  #maybe remove Expected filter, ask Jason if it's a good idea add Expected filter, look for standard way, Expected is more strict, and redndent, keep one 
  sigres <- sigres[(as.numeric(sigres$FDR) <= 0.1)&(sigres$Expected>=5),] 
  #Calculate FC based
  fc <- (sigres[,"Significant"]/sum(GOdata@allScores[GOdata@feasible]==1))/(sigres[,"Annotated"]/sum(GOdata@feasible))
  sigres <- data.frame(sigres,FC=fc)
  #Order terms
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
  sigres=as.data.table(sigres)
  #Find significant genes with rank, GenTable function
  siggene_forID=lapply(sigres$GO.ID,function(x,GOdata){
    gene=sigGenes(GOdata)[sigGenes(GOdata)%in%unlist(genesInTerm(GOdata, x))]
    gl_dt=data.table(rank=1:length(gl),gene=gl)
    mt=match(gl_dt$gene,gene)
    mt=mt[!is.na(mt)]
    highest_rank=min(gl_dt$rank[gl_dt$gene %in% gene])
    
    motif_in=motif[motif%in%gene]
    if(length(motif_in)==0){motif_in=NA}
    return(list(paste(gene[mt],collapse =";"),highest_rank,paste(motif_in,collapse = ";")))
    
    
  },GOdata=GOdata)
  siggene=unlist(lapply(siggene_forID,function(x) x[[1]]))
  max_rank=unlist(lapply(siggene_forID,function(x) x[[2]]))
  motif_in=unlist(lapply(siggene_forID,function(x) x[[3]]))
  if(nrow(sigres)>0){ 
    sigres$genes=siggene
    sigres$higest_ranks=max_rank
    sigres$motif=motif_in
  }
  return(sigres[sigres$FC>=1.5])
}

#CSV files are from Jason

csv_files=dir('../downstream/input/mm10_cluster/',pattern="csv")
tissue=unique(sub(".csv*","",csv_files))
#jsd_enhancer=readRDS('../downstream/input/chromHMM_enhancer.rds')
for (ts in tissue){
  fn=paste0(ts,'.csv')
  #read in csv file for given tissue
  csv_in_ts=fread(paste0('../downstream/input/mm10_cluster/',fn))
  csv_in_ts=csv_in_ts[order(dNME_maxJSD_rank,decreasing = F)]
  # Getting enhancer based on choice
  if(enc_type=="chromHMM_enhancer"){csv_in_ts=csv_in_ts[csv_in_ts$chromHMM_enhancer]}else
    if(enc_type=="bin_enhancer"){
      
      csv_in_gr=convert_GR(csv_in_ts$region)
      mcols(csv_in_gr)=csv_in_ts
      olap=findOverlaps(csv_in_gr,enhancer)
      csv_in_gr=csv_in_gr[queryHits(olap)]
      csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
      csv_in_gr$distance=NA
      csv_in_ts=as.data.table(mcols(csv_in_gr))
    }else
      if(enc_type=="TSS"){csv_in_ts=csv_in_ts[abs(distance)<=500]}
  
  if(nrow(csv_in_ts)>1){
    #GO annotation for each cluster
    csv_out=fastDoCall("rbind",lapply(1:10,function(clu){
      sp=paste0(ts,'-',clu)
      csv_in_ts_clu=csv_in_ts[cluster==clu]
      #rank csv based on dNME_maxJSD
      csv_in_ts_clu=csv_in_ts_clu[order(dNME_maxJSD_rank,decreasing=F)]
      
      if(nrow(csv_in_ts_clu)>1){
        
        #GO annotation for each cluster
        GO_out_cluster=GO_run(csv_in_ts_clu$gene,unique(csv_in_ts$gene),motif=motif_all)
        if(any(GO_out_cluster$motif!="NA")){cat(ts,clu,"Have motif in GO\n")}
        write.csv(GO_out_cluster,
                  file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GO/',sp,'_cluster_GO.csv'),row.names = F,quote = T)
        csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
        return(csv_in_ts_clu)
        
      }
      
    }))
    write.csv(csv_out,file=paste0('../downstream/output/mm10_result/',enc_type,'/all_gene_list/',ts,'_all.csv'))
  }
  
}
