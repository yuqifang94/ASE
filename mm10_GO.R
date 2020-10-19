rm(list=ls())
library(data.table)
suppressMessages(library(GenomicRanges))
suppressMessages(library(topGO))
library(RColorBrewer)
library(ggplot2)
enc_type="chromHMM_enhancer"
corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}
#GO annotation
GO_run<-function(gl,back,motif){
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")})
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  sigres <- sigres[as.numeric(sigres$FDR) <= 0.1&sigres$Expected>=5,]
  fc <- (sigres[,"Significant"]/sum(GOdata@allScores[GOdata@feasible]==1))/(sigres[,"Annotated"]/sum(GOdata@feasible))
  sigres <- data.frame(sigres,FC=fc)
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
  sigres=as.data.table(sigres)
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

# mml=readRDS("../downstream/output/mml_matrix_DNase.rds")
# nme=readRDS("../downstream/output/nme_matrix_DNase.rds")
# uc=readRDS('../downstream/output/uc_matrix_DNase.rds')

#get chromHMM data
if(enc_type=="chromHMM_enhancer"){enhancer_file="../downstream/output/chromHMM_enhancer.rds"}else 
  if(enc_type=="FeDMR"){enhancer_file="../downstream/output/FeDMR.rds"}
enhancer=readRDS(enhancer_file)

#Initialize data
# GO_out_cluster=list()
# GO_out_dNME_only=list()
# csv_in_ts_out=list()
# csv_in_ts_out_clu_gene=list()
cutoff_out=data.table()
csv_files=dir('../downstream/input/mm10_cluster/',pattern="csv")
tissue=unique(sub("_.*","",csv_files))
jsd_enhancer=readRDS('../downstream/input/chromHMM_enhancer.rds')
for (ts in tissue){

  #read in csv file for given tissue
  csv_in_ts=data.table()

  for(clu_in in 1:10){
    fn=paste0(ts,'_',clu_in,'.csv')
    csv_in=fread(paste0('../downstream/input/mm10_cluster/',fn))
  
    sn=sub('.csv','',fn)
    csv_in=csv_in[order(dNME_maxJSD_rank,decreasing = F)]
    setkeyv(csv_in,enc_type)
    
    csv_in=csv_in[.(TRUE)]
    if(nrow(csv_in)>1){
    csv_in$cluster=clu_in
    # csv_in$correlation_nme_mml=
    #   corfunc(nme[csv_in$region,sub('-.*','',colnames(nme))==ts],
    #           mml[csv_in$region,sub('-.*','',colnames(mml))==ts])
    # write.csv(csv_in,file=paste0('../downstream/output/mm10_result/',enc_type,'/gene_list_cluster/',ts,'_',clu_in,'.csv'))
     }
        csv_in_ts=rbind(csv_in_ts,csv_in)
        print(nrow(csv_in))
  }
  #reorder cluster section
  #csv_in_ts$cluster=jsd_enhancer[[ts]][csv_in_ts$region]
  cat(head( csv_in_ts$cluster),'\n')
  #csv file enhancer check
  csv_in_gr=GRanges(seqnames=sub(':.*','',csv_in_ts$region),
                                     IRanges(start=as.numeric(sub('-.*','',sub('.*:','',csv_in_ts$region))),
                                             end=as.numeric(sub('.*-','',csv_in_ts$region))))
  cat("enhancer check:",
      length(subsetByOverlaps(csv_in_gr,enhancer))==length(csv_in_gr),'\n')
  
  
  #Calculate nme mml correlation
  write.csv(csv_in_ts,file=paste0('../downstream/output/mm10_result/',enc_type,'/all_gene_list/',ts,'_all.csv'))
  if(nrow(csv_in_ts)>1){
  #GO annotation for each cluster
  cutoff_out=rbind(cutoff_out,do.call(rbind,lapply(1:10,function(clu){
    sp=paste0(ts,'-',clu)
    csv_in_ts_clu=csv_in_ts[cluster==clu]
    csv_in_ts_clu=csv_in_ts_clu[order(dNME_maxJSD_rank,decreasing=F)]
    csv_in_motif=fread(paste0('../downstream/input/mouse_motif_enrichment_enhancer/',ts,'/','motif_',ts,'_cluster_',clu,'_enhancer.csv')) 
    
    if(nrow(csv_in_ts_clu)>1){
    
    #GO annotation for each cluster
    motif_all=sub('\\(.*','',sub('.*_','',csv_in_motif$motif[csv_in_motif$FDR<=0.05]))
    motif_all=unlist(strsplit(motif_all,'::'))
    if (length(motif_all)==0){cat('No motif for',ts,clu,'\n')}
    GO_out_cluster=GO_run(csv_in_ts_clu$gene,unique(csv_in_ts$gene),motif=motif_all)
    if(any(GO_out_cluster$motif!="NA")){cat(ts,clu,"Have motif in GO\n")}
    write.csv(GO_out_cluster,
              file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GO/',sp,'_cluster_GO.csv'),row.names = F,quote = T)
    #get genes for each clusters with high dNME max and dMML max
    #cat('GO for high dNME\n')
    #NME and MML correlation cutoff
    # cor_cutoff=quantile(csv_in_ts_clu$correlation_nme_mml,prob=0.1,na.rm=T)
    # #NME and MML ratio cutoff
    # ratio_cutoff=quantile(csv_in_ts_clu$maxJSD_rankratio,prob=0.10,na.rm=T)
    # csv_in_ts_clu=csv_in_ts_clu[order(maxJSD_rankratio,decreasing=F)]
    # #select high dNME genes
    # csv_in_ts_clu_ft=csv_in_ts_clu[maxJSD_rankratio<=ratio_cutoff]
    # gene_dNME=unique(csv_in_ts_clu_ft$gene)
    # if(length(gene_dNME)>0){
    #  
    #   GO_out_dNME_only=GO_run(gene_dNME,unique(csv_in_ts$gene))
    #   cat("writing:",sp,'\n')
    #   write.csv(GO_out_dNME_only,
    #             file=paste0('../downstream/output/mm10_result/',enc_type,'/dNME_GO/',sp,'_dNME_dMML_high_ratio.csv'),row.names = F,quote = T)
    #   return( data.table(cluster=clu,tissue=ts,cutoff_dMML_dNME_ratio=ratio_cutoff,
    #                      min_dNME=min(csv_in_ts_clu_ft$dNME_maxpair),
    #                      min_dMML=min(csv_in_ts_clu_ft$dMML_maxpair)))
    #   }
    #csv_in_ts_out_clu_gene[[sp]]=csv_in_ts_clu_ft
   
    }
    
  })))
  }
  
}
#write.csv(cutoff_out,paste0('../downstream/output/mm10_result/',enc_type,'/dNME_GO/ratio_cutoff.csv'),row.names = F,quote = T)
