#GO annotation
source('mainFunctions_sub.R')
#Only save the output, not saving the csv file etc
GO_run_tissue_perm<-function(ts,dir_in,enc_type,dist_cutoff=NA,permute=F){
  GO_out_all=list()
  csv_files=dir(paste0('../downstream/input/',dir_in),pattern="csv")
  print(csv_files)
  cat("Processing:",ts,'\n')
   fn=paste0(ts,'.csv')
    #read in csv file for given tissue
    csv_in_ts=fread(paste0('../downstream/input/',dir_in,'/',fn))
    
    #Note some times Jason use dNME_maxJSD_rank
    csv_in_ts=csv_in_ts[order(dNME_maxUC_rank,decreasing = F)]
    # Getting enhancer
    if(enc_type=="chromHMM_enhancer"){csv_in_ts=csv_in_ts[csv_in_ts$chromHMM_enhancer]}else
      if(enc_type=="non_chromHMM_enhancer"){csv_in_ts=csv_in_ts[!csv_in_ts$chromHMM_enhancer]}else 
        if(enc_type=="promoter"){csv_in_ts=csv_in_ts}else 
          if(enc_type=="all_regions"){csv_in_ts=csv_in_ts}else{
            if(enc_type=="bin_enhancer"){
              enhancer=readRDS("../downstream/output/bin_enhancer.rds")
              csv_in_gr=convert_GR(csv_in_ts$region)
              mcols(csv_in_gr)=csv_in_ts
              olap=findOverlaps(csv_in_gr,enhancer)
              csv_in_gr=csv_in_gr[queryHits(olap)]
              csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
              csv_in_gr$distance=NA
              csv_in_ts=as.data.table(mcols(csv_in_gr))
            }
            
          }
    if(!is.na(dist_cutoff)){csv_in_ts=csv_in_ts[abs(distance)<=dist_cutoff]}
    if(permute){csv_in_ts$cluster=sample(csv_in_ts$cluster,length(csv_in_ts$cluster))}
    #GO annotation
    if(nrow(csv_in_ts)>1){
      #GO annotation for each cluster
      csv_out=mclapply(1:10,function(clu){
        sp=paste0(ts,'-',clu)
        csv_in_ts_clu=csv_in_ts[cluster==clu]
        csv_in_ts_clu=csv_in_ts_clu[order(dNME_maxUC_rank,decreasing=F)]
        #Add NME and mml cor
        # csv_in_ts_clu$nme_cor=nme_cor[[ts]][match(csv_in_ts_clu$region,names(nme_cor[[ts]]))]
        # csv_in_ts_clu$mml_cor=mml_cor[[ts]][match(csv_in_ts_clu$region,names(mml_cor[[ts]]))]
        
        if(nrow(csv_in_ts_clu)>1){
          
          #GO annotation for chromHMM 
          GO_out_cluster=GO_run(csv_in_ts_clu$gene,unique(csv_in_ts$gene),cluster=clu)
          csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
          
          
          write.csv(GO_out_cluster,row.names = F,quote = T,
                    file=paste0(GO_out,sp,'_cluster_GO.csv'))
          #return(list(csv_in_ts_clu=csv_in_ts_clu,GO_out_cluster=GO_out_cluster))
          return(list(GO_out_cluster_all=GO_out_cluster,csv_in_ts_clu=csv_in_ts_clu))
        }
        
      })
    }
    
  return(csv_out)
}
#calling code
#tissue, enhancer type,cutoff=3 permutation array number
args = commandArgs(trailingOnly=TRUE)
dir_in='zero_cutoff_mm10'
tissue=args[1]
enc_type=args[2]
cutoff=args[3]
array_id=args[4]
