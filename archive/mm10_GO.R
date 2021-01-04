rm(list=ls())
ptcount=0
library(data.table)
suppressMessages(library(GenomicRanges))
suppressMessages(library(topGO))
library(RColorBrewer)
library(ggplot2)
#library(rGREAT)
library(Gmisc)
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
GO_run<-function(gl,back){
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")})
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  #sigres <- sigres[as.numeric(sigres$FDR) <= 0.1,]
  fc <- ((sigres[,"Significant"]+ptcount)/(sum(GOdata@allScores[GOdata@feasible]==1)+ptcount))/((sigres[,"Annotated"]+ptcount)/(sum(GOdata@feasible)+ptcount))
  sigres <- data.frame(sigres,FC=fc)
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
  sigres=as.data.table(sigres)
  siggene_forID=lapply(sigres$GO.ID,function(x,GOdata){
    gene=sigGenes(GOdata)[sigGenes(GOdata)%in%unlist(genesInTerm(GOdata, x))]
    gl_dt=data.table(rank=1:length(gl),gene=gl)
    mt=match(gl_dt$gene,gene)
    mt=mt[!is.na(mt)]
    #highest_rank=min(gl_dt$rank[gl_dt$gene %in% gene])
    highest_rank=NA
    #motif_in=paste(motif[motif%in%gene],collapse = ";")
    #if(length(motif_in)==0){motif_in=NA}
    return(list(paste(gene[mt],collapse =";"),highest_rank))

    
  },GOdata=GOdata)
  siggene=unlist(lapply(siggene_forID,function(x) x[[1]]))
  max_rank=unlist(lapply(siggene_forID,function(x) x[[2]]))
   #motif_in=unlist(lapply(siggene_forID,function(x) x[[3]]))
  if(nrow(sigres)>0){ 
    sigres$genes=siggene
    sigres$higest_ranks=max_rank
   # sigres$motif=motif_in
   }
  return(sigres)
}

# mml=readRDS("../downstream/output/mml_matrix_DNase.rds")
# nme=readRDS("../downstream/output/nme_matrix_DNase.rds")
# uc=readRDS('../downstream/output/uc_matrix_DNase.rds')
# enhancer=fread('../downstream/input/s8C_enhancer_bin.csv',skip=1)
# enhancer=makeGRangesFromDataFrame(enhancer,seqnames.field = "chrom",keep.extra.columns = T)
# saveRDS(enhancer,"../downstream/output/bin_enhancer.rds")
#get chromHMM data
if(enc_type=="chromHMM_enhancer"){enhancer=readRDS("../downstream/output/chromHMM_enhancer.rds")}else
  if(enc_type=="FeDMR"){enhancer=readRDS("../downstream/output/FeDMR.rds")}else
    if(enc_type=="bin_enhancer"){enhancer=readRDS("../downstream/output/bin_enhancer.rds")}
nme_cor <- readRDS('../downstream/input/dnmecor.rds')
mml_cor <- readRDS('../downstream/input/dmmlcor.rds')
GO_out_all=list()

dir_in='mm10_cluster_all'
#Initialize data
# GO_out_cluster=list()
# GO_out_dNME_only=list()
# csv_in_ts_out=list()
# csv_in_ts_out_clu_gene=list()
#cutoff_out=data.table()
csv_files=dir(paste0('../downstream/input/',dir_in),pattern="csv")
#tissue=unique(sub(".csv*","",csv_files))
tissue=c("forebrain","heart","limb")
#jsd_enhancer=readRDS('../downstream/input/chromHMM_enhancer.rds')
for (ts in tissue){
  cat("Processing:",ts,'\n')
  fn=paste0(ts,'.csv')
  #read in csv file for given tissue
  csv_in_ts=fread(paste0('../downstream/input/',dir_in,'/',fn))
  csv_in_ts=csv_in_ts[order(dNME_maxJSD_rank,decreasing = F)]
# Getting enhancer
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
      if(enc_type=="TSS"){csv_in_ts=csv_in_ts[abs(distance)<=500]}else
  if(enc_type=="non_chromHMM_enhancer"){csv_in_ts=csv_in_ts[!csv_in_ts$chromHMM_enhancer]}else 
  if(enc_type=="all_regions"){csv_in_ts=csv_in_ts}
#   for(clu_in in 1:10){
#    
#     csv_in=fread(paste0('../downstream/input/mm10_cluster/',fn))
#   
#     sn=sub('.csv','',fn)
#     csv_in=csv_in[order(dNME_maxJSD_rank,decreasing = F)]
#     setkeyv(csv_in,enc_type)
#     
#     csv_in=csv_in[.(TRUE)]
#     if(nrow(csv_in)>1){
#     csv_in$cluster=clu_in
#     # csv_in$correlation_nme_mml=
#     #   corfunc(nme[csv_in$region,sub('-.*','',colnames(nme))==ts],
#     #           mml[csv_in$region,sub('-.*','',colnames(mml))==ts])
#     # write.csv(csv_in,file=paste0('../downstream/output/mm10_result/',enc_type,'/gene_list_cluster/',ts,'_',clu_in,'.csv'))
#      }
#         csv_in_ts=rbind(csv_in_ts,csv_in)
#         print(nrow(csv_in))
#   }
  #reorder cluster section
  #csv_in_ts$cluster=jsd_enhancer[[ts]][csv_in_ts$region]
 # cat(head( csv_in_ts$cluster),'\n')
  #csv file enhancer check
  # csv_in_gr=convert_GR(csv_in_ts$region)
  # 
  # cat("enhancer check:",
  #     length(subsetByOverlaps(csv_in_gr,enhancer))==length(csv_in_gr),'\n')
  
  
  #Calculate nme mml correlation
 
  if(nrow(csv_in_ts)>1){
  #GO annotation for each cluster
  csv_out=lapply(1:10,function(clu){
    sp=paste0(ts,'-',clu)
    csv_in_ts_clu=csv_in_ts[cluster==clu]
    csv_in_ts_clu=csv_in_ts_clu[order(dNME_maxJSD_rank,decreasing=F)]
    #Add NME and mml cor
    csv_in_ts_clu$nme_cor=nme_cor[[ts]][match(csv_in_ts_clu$region,names(nme_cor[[ts]]))]
    csv_in_ts_clu$mml_cor=mml_cor[[ts]][match(csv_in_ts_clu$region,names(mml_cor[[ts]]))]
    #csv_in_motif=fread(paste0('../downstream/input/mouse_motif_enrichment_enhancer/',ts,'/','motif_',ts,'_cluster_',clu,'_enhancer.csv')) 
    
    if(nrow(csv_in_ts_clu)>1){
    
    #GO annotation for each cluster
    # motif_all=sub('\\(.*','',sub('.*_','',csv_in_motif$motif[csv_in_motif$FDR<=0.05]))
    # motif_all=unlist(strsplit(motif_all,'::'))
    #if (length(motif_all)==0){cat('No motif for',ts,clu,'\n')}
    GO_out_cluster=GO_run(csv_in_ts_clu$gene,unique(csv_in_ts$gene))
    
   # if(any(GO_out_cluster$motif!="NA")){cat(ts,clu,"Have motif in GO\n")}
    write.csv(GO_out_cluster[FC>=1.5],
              file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GO/',dir_in,'/',sp,'_cluster_GO.csv'),row.names = F,quote = T)
    GO_out_cluster_NME=GO_run(csv_in_ts_clu[nme_cor>=0.7]$gene,unique(csv_in_ts$gene))
    GO_out_cluster_NME_only=GO_run(csv_in_ts_clu[nme_cor>=0.7&mml_cor<0.7]$gene,unique(csv_in_ts$gene))
    GO_out_cluster_MML=GO_run(csv_in_ts_clu[mml_cor>=0.7]$gene,unique(csv_in_ts$gene))
    GO_out_cluster_MML_only=GO_run(csv_in_ts_clu[mml_cor>=0.7&nme_cor<0.7]$gene,unique(csv_in_ts$gene))
    GO_out_cluster_non_MML=GO_run(csv_in_ts_clu[mml_cor<0.7]$gene,unique(csv_in_ts$gene))
    GO_out_cluster_NME_MML=GO_run(csv_in_ts_clu[nme_cor>=0.7&mml_cor>=0.7]$gene,unique(csv_in_ts$gene))
    GO_out_cluster_non_NME_non_MML=GO_run(csv_in_ts_clu[nme_cor<0.7&mml_cor<0.7]$gene,unique(csv_in_ts$gene))
    # write.csv(GO_out_cluster_NME,
    #           file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GO_high_NME_cor/',sp,'_cluster_GO.csv'),row.names = F,quote = T)
    csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
    
    GO_out_cluster$cluster=clu
    GO_out_cluster_NME$cluster=clu
    GO_out_cluster_NME_only$cluster=clu
    GO_out_cluster_MML$cluster=clu
    GO_out_cluster_MML_only$cluster=clu
    GO_out_cluster_non_MML$cluster=clu
    #return(list(csv_in_ts_clu=csv_in_ts_clu,GO_out_cluster=GO_out_cluster))
    return(list(GO_out_cluster_all=GO_out_cluster,GO_out_cluster_NME=GO_out_cluster_NME,GO_out_cluster_non_MML=GO_out_cluster_non_MML,csv_in_ts_clu=csv_in_ts_clu,
                GO_out_cluster_NME_only=GO_out_cluster_NME_only,GO_out_cluster_MML=GO_out_cluster_MML,GO_out_cluster_MML_only=GO_out_cluster_MML_only,
                GO_out_cluster_non_NME_non_MML=GO_out_cluster_non_NME_non_MML,GO_out_cluster_NME_MML=GO_out_cluster_NME_MML))
# GREAT analysis ----------------------------------------------------------
    # bed=convert_GR(csv_in_ts_clu$region)
    # bg=convert_GR(csv_in_ts$region)
    # job=submitGreatJob(bed,bg=bg,rule="oneClosest",adv_oneDistance= 1000.0,species="mm10")#default adv_oneDistance = 1000.0
    # tb=getEnrichmentTables(job, download_by = "tsv")
    # saveRDS(job,paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GREAT/job/',sp,'_GREAT_job.rds'))
    # write.csv(tb$`GO Biological Process`[tb$`GO Biological Process`$HyperP<=0.05&tb$`GO Biological Process`$RegionFoldEnrich>=1.5,],
    #           file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GREAT/table/',sp,'_GREAT_table.csv'),row.names = F,quote = T)
    # write.csv(tb$`GO Biological Process`[tb$`GO Biological Process`$HyperFdrQ<=0.1&tb$`GO Biological Process`$RegionFoldEnrich>=1.5,],
    #           file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GREAT/table/',sp,'_GREAT_table_FDR.csv'),row.names = F,quote = T)
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
    
  })
  GO_out_all[[ts]]=csv_out
  #write.csv(csv_in_ts,file=paste0('../downstream/output/mm10_result/',enc_type,'/all_gene_list/',ts,'_all.csv'))
  write.csv(fastDoCall('rbind',lapply(csv_out,function(x) x$csv_in_ts_clu))[order(dNME_maxJSD,decreasing=T)],
            file=paste0('../downstream/output/mm10_result/',enc_type,'/enhancer_gene_list/',dir_in,'/',ts,'_all.csv'))
  #write.csv(fastDoCall('rbind',lapply(csv_out,function(x) x[[2]])),file=paste0('../downstream/output/mm10_result/',enc_type,'/all_gene_list/',ts,'_all_GO_unft.csv'))
  }
 
}
# #write.csv(cutoff_out,paste0('../downstream/output/mm10_result/',enc_type,'/dNME_GO/ratio_cutoff.csv'),row.names = F,quote = T)
saveRDS(GO_out_all,'../downstream/output/GO_out_all_chromHMM_all_stat.rds')
GO_out=readRDS(paste0('../downstream/output/GO_out_all_chromHMM_all_stat.rds'))
#Plot GO terms with heatmap

dcast_matrix<-function(dt_in,value_in,colnames_order=colnames_order){
  dt_in=dcast.data.table(dt_in,Term~cluster,value.var  = value_in)
  dt_in_mt=as.matrix(dt_in[,-1])
  rownames(dt_in_mt)=dt_in$Term
  dt_in_mt=dt_in_mt[order(max.col(dt_in_mt,)),]
  dt_in_mt=dt_in_mt[,colnames_order]
  return(dt_in_mt)
}
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
#Using raw values to define a set of GO terms here, mainly tissue specific

library(pheatmap)
plot_GO_heatmap<-function(selected_tissue,GO_anno,GO_out){
   GO_tissue=list()
   for(tissue in selected_tissue){
    GO_in=GO_out[[tissue]]
    GO_in=fastDoCall('rbind',lapply(GO_in,function(x) x[[GO_anno]]))
    
    GO_in=GO_in[,.(GO.ID,Term,classicFisher,FDR,FC,cluster)]
    GO_in_top=do.call(c,lapply(1:10,function(x) {
  
      return(GO_in[cluster==x][FC>=1.5&FDR<=0.1][order(FDR,-FC,decreasing=F)][1:5]$GO.ID)
    }))
    GO_in=GO_in[GO.ID%in%GO_in_top]
    #GO_in=GO_in[,list(Term,FC,FDR,cluster)]
    #if(tissue=="heart"){GO_in=GO_in[Term %in% unique(GO_in$Term)[-c(2,3,8,12,13,14,15,16,18,20,21,25)]]}
    GO_in$log10FDR=-log10(GO_in$FDR)
    GO_in_main= dcast_matrix(GO_in,"FC")
    GO_in_main=GO_in_main[order(max.col(GO_in_main),decreasing = F),]
    GO_in_FDR= dcast_matrix(GO_in,"FDR")
    #GO_in_FDR_log10= dcast_matrix(GO_in,"log10FDR")
    GO_in_FDR[GO_in_FDR<=0.1]="*"
    GO_in_FDR[GO_in_FDR>0.1]=""
    GO_in_FDR=GO_in_FDR[rownames(GO_in_main),]
    col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
    c2 <- brewer.pal(10,'Set3')
    names(c2) <- 1:10
    breaksList = seq(-1, 1, by = 0.01)
    colann= data.frame(cluster=as.character(1:10))
    #pdf(paste0('../downstream/output/graphs/Figure6/GO_', tissue,'_',GO_anno,'_FC.pdf'),width=25,height=14)
   pheatmap(scalematrix(GO_in_main),cluster_rows =F,cluster_cols = F,
                   show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR,border_color = NA,
                   color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),
                   filename=paste0('../downstream/output/graphs/Figure6/all_regions/',GO_anno,'/GO_', tissue,'_',GO_anno,'_FC_chromHMM.pdf'),
                  cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
                  fontsize=30,legend = F)
    #dev.off()
    GO_tissue[[tissue]]=GO_in_main
    #,breaks=breaksList,color=col(length(breaksList))
  }

}
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME_only",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_MML",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_MML_only",GO_out)
GO_merge_terms=list()
for(tissue in selected_tissue){
  GO_in=GO_out[[tissue]]
  GO_in=fastDoCall('rbind',lapply(GO_in,function(x) x[[GO_anno]]))
  GO_in=GO_in[,.(GO.ID,Term,classicFisher,FDR,FC,cluster)]
  GO_in_top=do.call(c,lapply(1:10,function(x) {

    return(GO_in[cluster==x][FC>=1.5&FDR<=0.1][1:5]$Term)
  }))
  GO_merge_terms[[tissue]]=GO_in_top
}
#Forebrain manual selection
GO_merge_terms[["forebrain"]]=GO_merge_terms[["forebrain"]][c(-4,-5,-9,-10,-34,-46,-47)]
#Heart manual selection
GO_merge_terms[["heart"]]=GO_merge_terms[["heart"]][c(-6,-7,-9,-10,-14,-15,-40,-44,-45)]
#limb manual selection
GO_merge_terms[["limb"]]=GO_merge_terms[["limb"]][c(-16,-17,-18,-20,-22,-23,-25,-40,-50)]

# #Forebrain manual selection
# GO_merge_terms[["forebrain"]]=GO_merge_terms[["forebrain"]][c(-4,-5,-9,-10,-34,-46,-47)]
# #Heart manual selection
# GO_merge_terms[["heart"]]=GO_merge_terms[["heart"]][c(-6,-7,-9,-10,-14,-15,-40,-44,-45)]
# #limb manual selection
# GO_merge_terms[["limb"]]=GO_merge_terms[["limb"]][c(-16,-17,-18,-20,-22,-23,-25,-40,-50)]
saveRDS(GO_merge_terms,'../downstream/output/GO_terms_raw.rds')

GO_qc<-function(pc,scale){
  
  GO_out=readRDS(paste0('../downstream/output/GO_out_all_pcount',pc,'.rds'))
  #Get result for each tissue selected terms
  GO_tissue_selected_tissue=list()
  for(ts in selected_tissue){
    GO_tissue_selected=data.table()
    for(tissue in selected_tissue){
      GO_in=GO_out[[tissue]]
      GO_in=fastDoCall('rbind',lapply(GO_in,function(x) x[[GO_anno]]))
      GO_in=GO_in[,.(Term,classicFisher,FDR,FC,cluster)]
      GO_in$cluster=paste0(tissue,'-',GO_in$cluster)
      GO_tissue_selected=rbind(GO_tissue_selected,GO_in[Term %in% GO_merge_terms[[ts]]])
    }
    colnames_order=c()
    selected_tissue_ts=c(ts,selected_tissue[selected_tissue!=ts])
    for(tissue in selected_tissue_ts){colnames_order=c(colnames_order,paste0(tissue,'-',1:10))}
    FC_tissue_selected_mt=dcast_matrix(GO_tissue_selected,"FC")
    FDR_tissue_selected_mt=dcast_matrix(GO_tissue_selected,"FDR")
    FDR_tissue_selected_mt[FDR_tissue_selected_mt<=0.1]="*"
    FDR_tissue_selected_mt[FDR_tissue_selected_mt>0.1]=""
    FC_tissue_selected_mt=FC_tissue_selected_mt[,colnames_order]
    FDR_tissue_selected_mt= FDR_tissue_selected_mt[,colnames_order]
    FC_tissue_selected_mt_order=FC_tissue_selected_mt
    #FC_tissue_selected_mt_order[FDR_tissue_selected_mt==""]=1
    FC_tissue_selected_mt_order= FC_tissue_selected_mt_order[,gsub('-.*','',colnames(FC_tissue_selected_mt))==ts]
    
    
    FC_tissue_selected_mt=FC_tissue_selected_mt[order(max.col(FC_tissue_selected_mt_order)),]
    #scale by tissue
    FC_tissue_selected_mt_scale=NULL
    if(scale=='scale_tissue'){
      for(tissue_scale in selected_tissue_ts){
        FC_tissue_selected_mt_scale=cbind(FC_tissue_selected_mt_scale,
                                          scalematrix(FC_tissue_selected_mt[,gsub('-.*','',colnames(FC_tissue_selected_mt))==tissue_scale]))
        
        
      }
    }else if(scale=='scaled'){FC_tissue_selected_mt_scale=scalematrix(FC_tissue_selected_mt)}else if(scale=='unscaled'){
      FC_tissue_selected_mt_scale=FC_tissue_selected_mt
      
    }else if(scale=="scaled_FDR"){
      FC_tissue_selected_mt_scale=FC_tissue_selected_mt
     FC_tissue_selected_mt_scale[FDR_tissue_selected_mt==""]=1
     FC_tissue_selected_mt_scale=scalematrix(FC_tissue_selected_mt_scale)
    }
    #Plotting heatmap
    colann <- data.frame(tissue=sub('-.*','',colnames(FC_tissue_selected_mt_scale)),cluster=sub('.*-','',colnames(FC_tissue_selected_mt_scale)),stringsAsFactors = F)
    rownames(colann) <- colnames(FC_tissue_selected_mt_scale)
    col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
    breaksList = seq(-3, 3, by = 0.01)
    tiff(paste0('../downstream/output/GO_qc/',ts,'_GO_all_',scale,'_pc',pc,'.tiff'),width=2200,height=1000)
    print(pheatmap(t(FC_tissue_selected_mt_scale),cluster_rows =F,cluster_cols = F,display_numbers=t(FDR_tissue_selected_mt),border_color = NA,fontsize =24,
                   show_colnames = T,show_rownames = T,gaps_row  = cumsum(rle(colann[,1])$lengths),breaks=breaksList,color=col(length(breaksList))))
    dev.off()
    
    #,breaks=breaksList,color=col(length(breaksList))
  }

}
GO_qc(0,'unscaled')
GO_qc(0,'scaled')
GO_qc(0,'scale_tissue')
GO_qc(0,'scaled_FDR')
GO_qc(10,'unscaled')
GO_qc(10,'scaled')
GO_qc(10,'scale_tissue')
GO_qc(100,'unscaled')
GO_qc(100,'scaled')
GO_qc(100,'scale_tissue')




#breaksList = seq(-3, 3, by = 0.01)

pdf('GO_nopc.pdf',width=20,height=7)
print(pheatmap(scalematrix(FC_tissue_selected_mt),cluster_rows =F,cluster_cols = F,display_numbers=FDR_tissue_selected_mt,
               show_colnames = T,show_rownames = T,gaps_col = cumsum(rle(colann[,1])$lengths)))
dev.off()
ggplot(GO_tissue$heart,aes(x=cluster,y=Term,fill=-log10(FDR)))+geom_tile()+ scale_fill_distiller(palette = "RdBu")
