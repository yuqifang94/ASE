GO_run_tissue_perm<-function(ts,dir_in,enc_type,dist_cutoff=NA,permute=F,bg=NULL,extend=0,GO_type="TopGO",ranking_stat=NA){
  #ranking_stat = "dNME_maxUC" or "dMML_maxUC"
  print(ts)
  print(enc_type)
  print(dist_cutoff)
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
        if(enc_type=="all_regions"){csv_in_ts=csv_in_ts}else
          if(enc_type=="bin_enhancer"){
            enhancer=readRDS("../downstream/output/bin_enhancer.rds")
            enhancer=resize(enhancer,width = width(enhancer)+extend*2,fix='center')
            csv_in_gr=convert_GR(csv_in_ts$region)
            mcols(csv_in_gr)=csv_in_ts
            olap=findOverlaps(csv_in_gr,enhancer)
            csv_in_gr=csv_in_gr[queryHits(olap)]
            csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
            csv_in_gr$distance=NA
            csv_in_ts=as.data.table(mcols(csv_in_gr))
          }else
            if(enc_type=="FANTOM5"){
              #heart 47831 regions
              enhancer=import.bed('../downstream/input/F5.mm10.enhancers.bed.gz')
              csv_in_gr=convert_GR(csv_in_ts$region)
              olap=findOverlaps(csv_in_gr,enhancer)
              print("Filtering")
              csv_in_ts=csv_in_ts[queryHits(olap)]
            }else
              if(enc_type=="chromHMM_enhancer_union"){
                csv_in_ts=csv_in_ts[csv_in_ts$chromHMM_enhancer_union]
                
              }else
                if(enc_type=="chromHMM_enhancer_tissuespecific"){
                  csv_in_ts=csv_in_ts[csv_in_ts$chromHMM_enhancer_tissuespecific]
                  
                }else
                  if(enc_type=="bin_enhancer_promoter"){
                    enhancer=readRDS("../downstream/output/bin_enhancer.rds")
                    enhancer=resize(enhancer,width = width(enhancer)+extend*2,fix='center')
                    csv_in_gr=convert_GR(csv_in_ts$region)
                    mcols(csv_in_gr)=csv_in_ts
                    olap=findOverlaps(csv_in_gr,enhancer)
                    csv_in_gr=csv_in_gr[queryHits(olap)]
                    csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
                    csv_in_gr$distance=NA
                    csv_in_ts=as.data.table(mcols(csv_in_gr))
                    csv_in_ts=rbind(csv_in_ts,csv_in_ts[distance<=dist_cutoff])
                  }
  if(GO_type=="GSEA"){
    library(qusage)
    library(msigdbr)
    library(fgsea)
    all_gene_sets = msigdbr(species = "Mus musculus")
    h_gene_sets = msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:BP")
    h_gene_sets_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
  }
  if(dist_cutoff>0&enc_type=="promoter"){csv_in_ts=csv_in_ts[abs(distance)<=dist_cutoff]}
  if(permute){csv_in_ts$cluster=sample(csv_in_ts$cluster,length(csv_in_ts$cluster))}
  print(csv_in_ts)
  #GO annotation
  if(nrow(csv_in_ts)>1){
    if(is.null(bg)){bg=unique(csv_in_ts$gene)}
    #GO annotation for each cluster
    print(nrow(csv_in_ts))
    print(csv_in_ts)
    csv_out=lapply(1:max(csv_in_ts$cluster),function(clu){
      sp=paste0(ts,'-',clu)
      csv_in_ts_clu=csv_in_ts[cluster==clu]
      #csv_in_ts_clu=csv_in_ts_clu[order(dNME_maxUC_rank,decreasing=F)]
      #Add NME and mml cor
      # csv_in_ts_clu$nme_cor=nme_cor[[ts]][match(csv_in_ts_clu$region,names(nme_cor[[ts]]))]
      # csv_in_ts_clu$mml_cor=mml_cor[[ts]][match(csv_in_ts_clu$region,names(mml_cor[[ts]]))]
      
      if(nrow(csv_in_ts_clu)>1){
        #Sys.sleep(clu*60)
        #GO annotation for chromHMM
        cat('start processing cluster:',clu,'\n')
        if(GO_type=="TopGO"){
        tt1=proc.time()[[3]]
        cat('length of background gene:',length(bg),'\n')
        GO_out_cluster=GO_run(unique(csv_in_ts_clu$gene),bg,cluster=clu)
        csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
        
        cat('Finish processing cluster:',clu,'in:',proc.time()[[3]]-tt1,'\n')
        # write.csv(GO_out_cluster,row.names = F,quote = T,
        #           file=paste0(GO_out,sp,'_cluster_GO.csv'))
        #return(list(csv_in_ts_clu=csv_in_ts_clu,GO_out_cluster=GO_out_cluster))
        return(list(GO_out_cluster_all=GO_out_cluster,csv_in_ts_clu=csv_in_ts_clu))
        }else 
          if(GO_type=="GSEA"){
            if( ranking_stat == "dNME_maxUC"){
                csv_in_ts_clu_max=csv_in_ts_clu[,list(stat_max=max(dNME_maxUC)),by=list(gene)]}else
              if(ranking_stat=="dMML_maxUC"){
                csv_in_ts_clu_max=csv_in_ts_clu[,list(stat_max=max(dMML_maxUC)),by=list(gene)]
                
              }
            ranks=csv_in_ts_clu_max$stat_max
            names(ranks)=csv_in_ts_clu_max$gene
            gsea_out<- fgsea(h_gene_sets_list, ranks, minSize=10, maxSize = 500, nperm=100000)
            
            gsea_out$FDR=gsea_out$padj
            gsea_out$FC=gsea_out$NES
            gsea_out$`GO.ID`=gsea_out$pathway
            gsea_out$classicFisher=gsea_out$pval
            gsea_out$cluster=clu
            return(list(GO_out_cluster_all=gsea_out))
            
          }
      }
      
    })
  }
  
  return(csv_out)
}

plot_GO_heatmap<-function(tissue,GO_anno,GO_in,enc_type,FC_cutoff=1.5,FDR_cutoff=0.05,clu_in=NULL){
  
  
  
  GO_in=fastDoCall('rbind',lapply(GO_in,function(x) {
    x=x[[GO_anno]]
    x$sig_num=sum(x$FC>=FC_cutoff&x$FDR<=FDR_cutoff)
    return(x)
  }))
  total_clu=max(GO_in$cluster)
  if(is.null(clu_in)){clu_in=1:total_clu}
  if(nrow(GO_in)>0&!is.null(GO_in)){
    
    GO_in=GO_in[,.(GO.ID,Term,classicFisher,FDR,FC,cluster,sig_num)]
    
    GO_in_top=do.call(c,lapply(clu_in,function(x) {
      
      return(GO_in[cluster==x][FC>=FC_cutoff&FDR<=FDR_cutoff][order(-FC,decreasing=F)][1:5]$GO.ID)
    }))
    GO_in=GO_in[GO.ID%in%GO_in_top&cluster%in%clu_in]
    GO_in$log10FDR=-log10(GO_in$FDR)
    GO_in_main= dcast_matrix(GO_in,"FC")
    GO_in_main=GO_in_main[order(max.col(GO_in_main),decreasing = F),]
    GO_in_FDR= dcast_matrix(GO_in,"FDR")
    #GO_in_FDR_log10= dcast_matrix(GO_in,"log10FDR")
    GO_in_FDR[GO_in_FDR<=0.1]="*"
    GO_in_FDR[GO_in_FDR>0.1]=""
    GO_in_FDR=GO_in_FDR[rownames(GO_in_main),]
    col_label=unique(GO_in[,.(cluster,sig_num)])
    col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    #c2 <- brewer.pal(total_clu,'Set3')
    c2=sample(color,length(clu_in))
    names(c2) <- clu_in
    breaksList = seq(-1, 1, by = 0.01)
    colann= data.frame(cluster=as.character(clu_in))
    #pdf(paste0('../downstream/output/graphs/Figure6/GO_', tissue,'_',GO_anno,'_FC.pdf'),width=25,height=14)
    pheatmap(scalematrix(GO_in_main),cluster_rows =F,cluster_cols = F,
             show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR,border_color = NA,
             color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(200),
             filename=paste0('../downstream/output/GO_', tissue,'_',GO_anno,'_FC_',enc_type,'.pdf'),
             cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
             fontsize=30,legend = F)#,labels_col=col_label$sig_num)
    #dev.off()
  }else{print(GO_in)}
  #,breaks=breaksList,color=col(length(breaksList))
}

plot_GO_heatmap_variable<-function(fn_in,FC_cutoff=1.5,FDR_cutoff=0.1,clu_in=NULL){
  cat("Processing:",fn_in,'\n')
  GO_in=readRDS(fn_in)
  GO_anno="GO_out_cluster_all"
  
  GO_in=fastDoCall('rbind',lapply(GO_in,function(x) {
    if(!is.null(x[[GO_anno]])){
      if(nrow(x[[GO_anno]])>0){
      x=x[[GO_anno]]
      x=x[Annotated>=10]
      x$FDR=p.adjust(x$classicFisher,method='fdr')
      x$sig_num=sum(x$FC>=FC_cutoff&x$FDR<=FDR_cutoff)
      return(x)
      }
    }
  }))
  if(grepl("GSEA",fn_in)){GO_in$Term=GO_in$GO.ID}
  if(!is.null(GO_in)){
    if(nrow(GO_in)>0){
      #print(GO_in)
      total_clu=max(GO_in$cluster)
      if(is.null(clu_in)){clu_in=1:total_clu}
   
      GO_in=GO_in[,.(GO.ID,Term,classicFisher,FDR,FC,cluster,sig_num)]
      
      GO_in_top=do.call(c,lapply(clu_in,function(x) {
        
        return(GO_in[cluster==x][FC>=FC_cutoff&FDR<=FDR_cutoff][order(FDR,-FC,decreasing=F)][1:10]$GO.ID)
      }))
      GO_in=GO_in[GO.ID%in%GO_in_top&cluster%in%clu_in]
      if(length(GO_in_top)>1&nrow(GO_in)>0){
       
        GO_in$log10FDR=-log10(GO_in$FDR)
        GO_in_main= dcast_matrix(GO_in,"FC")
        GO_in_FDR= dcast_matrix(GO_in,"FDR")
        if(nrow(GO_in_main)>1){
        GO_in_main=GO_in_main[order(max.col(GO_in_main),decreasing = F),]
        
        
        #GO_in_FDR_log10= dcast_matrix(GO_in,"log10FDR")
        GO_in_FDR[GO_in_FDR<=FDR_cutoff]="*"
        GO_in_FDR[GO_in_FDR>FDR_cutoff]=" "
        GO_in_FDR=GO_in_FDR[rownames(GO_in_main),]
        }
        col_label=unique(GO_in[,.(cluster,sig_num)])
        col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
        color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
        #c2 <- brewer.pal(total_clu,'Set3')
        c2=sample(color,length(clu_in))
        names(c2) <- clu_in
        breaksList = seq(-1, 1, by = 0.01)
        colann= data.frame(cluster=as.character(clu_in))
        GO_in_main[is.na(GO_in_main)]=1
        GO_in_FDR[is.na(GO_in_FDR)]=" "
        #pdf(paste0('../downstream/output/graphs/Figure6/GO_', tissue,'_',GO_anno,'_FC.pdf'),width=25,height=14)
        pheatmap(scalematrix(GO_in_main),cluster_rows =F,cluster_cols = F,
                 show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR,border_color = NA,
                 color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(200),
                 filename=gsub('.rds','.pdf',fn_in),
                 cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
                 fontsize=30,legend = F)#,labels_col=col_label$sig_num)
      }
      #dev.off()
    }else{print(GO_in)}
  }
  #,breaks=breaksList,color=col(length(breaksList))
}



GO_all_variabile<-function(dir_target,tissue_all=c("heart","limb","forebrain"),GO_type="TopGO",
                           ranking_stat=NA,bg_enhancer=NA,bg_promoter=NA){
  cutoff_UC=dir(paste0('../downstream/input/',dir_target))
  for(tissue in tissue_all){
    for(cutoff_uc_in in cutoff_UC){
      dir_in=paste0(dir_target,cutoff_uc_in,'/full')
      #Setting variables
     
     
      #Loop over 2 enhancer setting and two promoter settings
      for(enc_type in c("bin_enhancer","promoter")){
        if(enc_type =="bin_enhancer"){
          dist_cutoff=-1
          bg=bg_enhancer
        }else
          if(enc_type == "promoter"){
            dist_cutoff=2000
            bg=bg_promoter
          }
        # for(bg_setting in c("all_bg","ts_bg")){
        #   if(bg_setting=="all_bg"){
        #     dir_all_ts=paste0('../downstream/input/',dir_in,'/')
        #     csv_all=data.table()
        #     for (fn in dir(dir_all_ts)){
        #       csv_all=rbind(csv_all,fread(paste0(dir_all_ts,fn)))
        #     }
        #     if(enc_type=="bin_enhancer"){
        #     bin_enhancer_in=readRDS('../downstream/output/bin_enhancer.rds')
        #     bg=subsetByOverlaps(bin_enhancer_in,convert_GR(csv_all$region))$`Target Gene`
        #     }else
        #       if(enc_type=="promoter"){
        #         bg=unique(csv_all[distance<=dist_cutoff]$gene)
        #         
        #         
        #       }
        #   }else
        #     if(bg_setting=="ts_bg"){
        #       bg=NULL
        #     }
        bg_setting="enhancer_bg"
          cat("Processing:",tissue,"in",dir_in,"with",enc_type,',',bg_setting,'\n')
          GO_out_all=GO_run_tissue_perm(tissue,dir_in,enc_type,dist_cutoff=dist_cutoff,
                                        permute=F,bg=bg,extend=0,GO_type=GO_type,ranking_stat=ranking_stat)
          fn_out=paste0('../downstream/output/',dir_target,tissue,'_',bg_setting,'_',cutoff_uc_in,'_',enc_type,'_',GO_type,'_',ranking_stat,'.rds')
          cat("saving:",fn_out,'\n')
          saveRDS(GO_out_all,fn_out)
          
        #}
        
      }
      
    }
    
  }
}

#Run GO analysis for each region type
GO_run_tissue<-function(ts,dir_in,region_type_sel=NA,bg=NULL){
  #ranking_stat = "dNME_maxUC" or "dMML_maxUC"
  GO_out_all=list()
  csv_files=dir(dir_in,pattern="csv")
  cat("Processing:",ts,'\n')
  fn=paste0(ts,'.csv')
  #read in csv file for given tissue
  csv_in_ts=fread(paste0(dir_in,fn))
  #Note some times Jason use dNME_maxJSD_rank
  csv_in_ts=csv_in_ts[order(dNME_maxUC_rank,decreasing = F)]
  # Getting enhancer
  enhancer=readRDS("../downstream/output/bin_enhancer.rds")
  csv_in_gr=convert_GR(csv_in_ts$region)
  mcols(csv_in_gr)=csv_in_ts
  olap=findOverlaps(csv_in_gr,enhancer)
  csv_in_gr=csv_in_gr[queryHits(olap)]
  csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
  csv_in_gr$distance=NA
  csv_in_ts=as.data.table(mcols(csv_in_gr))
  if(!is.na(region_type_sel)){
    csv_in_ts=csv_in_ts[region_type==region_type_sel]
    
  }
  #GO annotation
  if(nrow(csv_in_ts)>1){
    if(is.null(bg)){bg=unique(csv_in_ts$gene)}
    #GO annotation for each cluster
    print(nrow(csv_in_ts))
    print(csv_in_ts)
    csv_out=lapply(1:max(csv_in_ts$cluster),function(clu){
      sp=paste0(ts,'-',clu)
      csv_in_ts_clu=csv_in_ts[cluster==clu]
      if(nrow(csv_in_ts_clu)>1){
        cat('start processing cluster:',clu,'\n')
          tt1=proc.time()[[3]]
          cat('length of background gene:',length(bg),'\n')
          GO_out_cluster=GO_run(unique(csv_in_ts_clu$gene),bg,cluster=clu)
          csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
          
          cat('Finish processing cluster:',clu,'in:',proc.time()[[3]]-tt1,'\n')
          return(list(GO_out_cluster_all=GO_out_cluster,csv_in_ts_clu=csv_in_ts_clu))
   
      
      
    }
      })
  }
  
  return(csv_out)
}
