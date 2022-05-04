source('mainFunctions_sub.R')
dir_out='../downstream/output/mouse_analysis/cell_population/'
GO_run_tissue_remove_gene<-function(ts,dir_in,enc_type,geneToRemove,region_type_sel=NA,bg=NULL,
                         active_enc=F,enc_cor=NA){
  #ranking_stat = "dNME_maxUC" or "dMML_maxUC"
  GO_out_all=list()
  cat("Processing:",ts,'\n')
  fn=paste0(ts,'.csv')
  #read in csv file for given tissue
  csv_in_ts=fread(paste0(dir_in,fn))

  #Note some times Jason use dNME_maxJSD_rank
  csv_in_ts=csv_in_ts[order(dNME_max_UC_pair_adj,decreasing = T)]

  # Getting enhancer
  print(enc_type)
  if(enc_type=="enhancer"&(!active_enc)){
    enhancer=readRDS(bin_enhancer_rds)
    csv_in_gr=convert_GR(csv_in_ts$regions)
   
    mcols(csv_in_gr)=csv_in_ts
    olap=findOverlaps(csv_in_gr,enhancer)
    csv_in_gr=csv_in_gr[queryHits(olap)]
    csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
    csv_in_gr$distance=NA
    csv_in_ts=as.data.table(mcols(csv_in_gr))

  }else 
    if(enc_type=="enhancer"&active_enc){
      cat("Analyzing active enhancer\n")
      enc_cor_ts=enc_cor[tissue==ts&cor_FDR<=0.1]
      olap=findOverlaps(convert_GR(csv_in_ts$regions,direction='GR'),
                      convert_GR(enc_cor_ts$region,direction='GR'))
      csv_in_ts=csv_in_ts[queryHits(olap)]
      csv_in_ts$gene=enc_cor_ts[subjectHits(olap)]$target_gene
  }else
    if(enc_type=="promoter"){
      csv_in_ts=csv_in_ts[abs(distance)<=2000]
      
      
    }

  if(region_type_sel!="all"){
    csv_in_ts=csv_in_ts[region_type==region_type_sel]
    print(csv_in_ts)
  }
  #GO annotation
  if(nrow(csv_in_ts)>1){
    #GO annotation for each cluster
    print(nrow(csv_in_ts))
  
    csv_out=lapply(1:max(csv_in_ts$cluster),function(clu){
      sp=paste0(ts,'-',clu)
      csv_in_ts_clu=csv_in_ts[cluster==clu]
      if(nrow(csv_in_ts_clu)>1){
        cat('start processing cluster:',clu,'\n')
        tt1=proc.time()[[3]]
        cat('length of background gene:',length(bg),'\n')
        geneIn=unique(csv_in_ts_clu$gene)
        geneIn=geneIn[!geneIn%in%geneToRemove]
        GO_out_cluster=GO_run(geneIn,bg,cluster=clu)
        csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
        
        cat('Finish processing cluster:',clu,'in:',proc.time()[[3]]-tt1,'\n')
        return(list(GO_out_cluster_all=GO_out_cluster,csv_in_ts_clu=csv_in_ts_clu))
        
        
        
      }
    })
  }
  print(csv_out[[1]]$GO_out_cluster_all)
  return(csv_out)
}

UC_merge=readRDS(UC_merge_file)#Define all analyzed regions, were using UC_merge_max_loc_cluster01.rds,4626
cutoff_fn='01'
#Runnning
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
#prepare enhancer background gene list
uc_gr=lapply(UC_merge,function(x) rownames(x))
uc_gr=Reduce(intersect,uc_gr)
uc_gr=convert_GR(uc_gr)
enhancer=readRDS(bin_enhancer_rds)#21441
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)
heart_markerGene=c("Myh6","Ttn","Myh7","Tnnt2","Actc1","Tnnc1","Smpx","Tnni3","Nexn","Sh3bgr",#CM
                    "Cd93","Gpr116","Pecam1","Klhl4","Rasgrp3","Emcn","Cav2","Cdh5","Ctla2b",#Endothelial
                    "Col1a2","Col3a1","Col1a1","Pde1a","Col5a1","Sox9","Thbs1","Dcn","Postn","Tcf21")#Fibroblast
bg_enhancer = bg_enhancer[!tolower(bg_enhancer)%in%tolower(heart_markerGene)]#23/30 are in enhancer list
heat_GO_gene_remove=GO_run_tissue_remove_gene("heart",dir_out_cluster01,geneToRemove=heart_markerGene,enc_type="enhancer",region_type_sel="all",bg=bg_enhancer)
heat_GO_gene_remove=lapply(heat_GO_gene_remove,function(x){
    return(list(GO_out_cluster_all=x$GO_out_cluster_all,
                csv_in_ts_clu=cbind(x$csv_in_ts_clu,as.data.table(UC_merge[[ts]][x$csv_in_ts_clu$region,!grepl('max',colnames(UC_merge[[ts]]))]))))
    
    })
GO_gene_remove=list()
GO_gene_remove[["heart"]]=heat_GO_gene_remove
select_top_GO_out=select_top_GO(GO_gene_remove,"heart",ptcount=0,FDR_cutoff=0.2,FC_cutoff=1.5)
  plot_GO_heatmap_all("heart",GO_gene_remove,region_type="enhancer",enc_type="enhancer_removeMarker",ptcount=0,
                        FDR_cutoff=0.2,
                      dir_plot="../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run_01/remove_heart_marker_gene/")
#Check the marker gene UC and NME distribution
NME_in=readRDS(NME_matrix_file)
#Extract enhancer of marker gene
enhancer_marker=enhancer[tolower(enhancer$`Target Gene`) %in% tolower(heart_markerGene)]
#Extract UC>0.1 enhancer
tissue='heart'
UC_ts=UC_merge[[tissue]]#4876367
# UC_ts=UC_ts[rowSums(UC_ts[,grepl('UC-',colnames(UC_ts))]>0.1)>0,]#994001
# #Extract UC>0.1 in all other tissues
# highUC_regions=lapply(UC_merge[which(names(UC_merge)!=tissue)],function(x) names(which(rowSums(x[,grepl('UC-',colnames(x))]>0.1)>0)))
# UC_ts=UC_ts[!rownames(UC_ts)%in%unlist(highUC_regions),]#51534
#Asign clusters
cluster_result=readRDS(cluster_01_region_out_fn)
UC_ts=UC_ts[rownames(UC_ts)%in%cluster_result[[ts]]$regions,!grepl("max",colnames(UC_ts))]
UC_ts$cluster=as.numeric(cluster_result[[ts]][match(rownames(UC_ts),regions)]$cluster)
UC_ts=UC_ts[order(UC_ts$cluster,decreasing=F),]
UC_ts_marker_olap=findOverlaps(convert_GR(rownames(UC_ts),direction="GR"),enhancer_marker)
UC_ts_nonMarker=UC_ts[-queryHits(UC_ts_marker_olap),]#49891
UC_ts_marker=UC_ts[queryHits(UC_ts_marker_olap),]#25
UC_ts_marker$marker=enhancer_marker$`Target Gene`[subjectHits(UC_ts_marker_olap)]
plot_heatmap_cluster_tissue(UC_ts_marker,ts,paste0(dir_out,ts,"_markerGene_dNME.png"),"dNME",marker=T)
plot_heatmap_cluster_tissue(UC_ts_nonMarker,ts,paste0(dir_out,ts,"_nonMarkerGene_dNME.png"),"dNME")
plot_heatmap_cluster_tissue(UC_ts_marker,ts,paste0(dir_out,ts,"_markerGene_dMML.png"),"dMML",marker=T)
plot_heatmap_cluster_tissue(UC_ts_nonMarker,ts,paste0(dir_out,ts,"_nonMarkerGene_dMML.png"),"dMML")
plot_heatmap_cluster_tissue(UC_ts_marker,ts,paste0(dir_out,ts,"_markerGene_UC.png"),"UC",marker=T)
plot_heatmap_cluster_tissue(UC_ts_nonMarker,ts,paste0(dir_out,ts,"_nonMarkerGene_UC.png"),"UC")
NME_in=readRDS(NME_matrix_file)
plotClusterMarker_NMEMML(NME_in,"NME",UC_ts_marker,ts)
MML_in=readRDS(MML_matrix_file)
plotClusterMarker_NMEMML(MML_in,"MML",UC_ts_marker,ts)
plotClusterMarker_NMEMML<-function(datIn,statIn,UC_ts_marker,ts){
  datIn_ts=as.data.frame(convert_GR(datIn,direction="matrix"))
  datIn_ts=datIn_ts[cluster_result[[ts]]$regions,grepl(ts,colnames(datIn_ts))]
  colnames(datIn_ts)=paste0(statIn,"-",colnames(datIn_ts))
  datIn_ts$cluster=as.numeric(cluster_result[[ts]][match(rownames(datIn_ts),regions)]$cluster)
  datIn_ts_marker=datIn_ts[rownames(UC_ts_marker),]#25
  datIn_ts_nonMmarker=datIn_ts[rownames(UC_ts_nonMarker),]#49891
  datIn_ts_marker$marker=UC_ts_marker$marker#25
  plot_heatmap_cluster_tissue(datIn_ts_marker,ts,paste0(dir_out,ts,"_markerGene_",statIn,".png"),statIn,difference=F,marker=T)
  plot_heatmap_cluster_tissue(datIn_ts_nonMmarker,ts,paste0(dir_out,ts,"_nonMarkerGene_",statIn,".png"),statIn,difference=F)
}


#Plotting UC, dNME, dMML heatmap 
plot_heatmap_cluster_tissue<-function(d,tissue_all,figure_name,stat_in,figure_width=1000,figure_height=1000,res=200,difference=T,scale=T,marker=F){
    cat('Plotting heatmap\n')
    d=d[order(d$cluster,decreasing=F),]
    if(marker){markerGene=d$marker}
    cl=as.character(d$cluster)
    d=d[,grepl(stat_in,colnames(d))]
    library(gplots)
    timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))      
    colnames(d)=gsub(paste0(stat_in,'-|-all|',paste(tissue_all,'-',sep='',collapse = '|')),'',colnames(d))
    
    d <- d[rowSums(d) > 0,]
    if(difference){
    d <- d[,colnames(d) %in% timeorder]
    d <- d[,order(match(colnames(d),timeorder))] 
    }
    #d <- scalematrdx(d)
    d <- d[complete.cases(d),]
    mat <-d
    rm(d)
    na_ma=-which(rowSums(is.na(mat))>0)
    if(length(na_ma)>0){
      mat= mat[na_ma,]
      rowann <- data.frame(cluster=cl[na_ma],
                          #dMMLJSDcor=dmml[[n]][rownames(mat)],
                          #dNMEJSDcor=dnme[[n]][rownames(mat)],
                          stringsAsFactors = F)
      rowann=rowann[na_ma,]
    }else{
    
      rowann <- data.frame(cluster=cl,
                          #dMMLJSDcor=dmml[[n]][rownames(mat)],
                          #dNMEJSDcor=dnme[[n]][rownames(mat)],
                          stringsAsFactors = F)
    }
    
    
      #Note for non-tissue specific UC 01, add tissue name to mat
    rownames(rowann)=rownames(mat)
    #Refine plotting parameters
   
    colann <- data.frame(time=sub('.*:','',colnames(mat)),stringsAsFactors = F)

    rownames(colann) <- colnames(mat)
   
    c2 <- brewer.pal(10,'Set3')
    names(c2) <- 1:10
    c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
    names(c4) <- sort(unique(colann[,1]))
    anno_color=list(cluster=c2,time=c4)
    if(marker){
      rowann$marker=markerGene
      coul <- brewer.pal(8, "Set2")
      marker_col = colorRampPalette(coul)(length(unique(markerGene)))
      names(marker_col)=unique(markerGene)
      anno_color$marker=marker_col

    }
    #remove row with all NA 

    #sub_sp=sort(sample(1:nrow(mat),round(nrow(mat))))
    if(scale){
    mat_sc=scalematrix(mat)
    }else{mat_sc=mat}
    png(figure_name,
        width=figure_width,
        height=figure_height,
        res=res,
        type='cairo')
    if(marker){
    pheatmap(mat_sc,cluster_rows = F,
            annotation_row = rowann,
            cluster_cols = F,
            annotation_col = colann,show_colnames = F,show_rownames = F,
            #gaps_row = row_gap[-1],
            display_numbers=round(mat,digits=2),
            #filename=figure_name,
            annotation_colors =anno_color
                                      #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10))
                                      
           
    )}else{

        pheatmap(mat_sc,cluster_rows = F,
            annotation_row = rowann,
            cluster_cols = F,
            annotation_col = colann,show_colnames = F,show_rownames = F,
            #gaps_row = row_gap[-1],
            
            #filename=figure_name,
            annotation_colors =anno_color
                                      #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10))
                                      
           
    )
    }
    dev.off()
  
}
