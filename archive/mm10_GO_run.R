rm(list=ls())
source('mainFunctions_sub.R')




# Plotting GO with heatmap ------------------------------------------------


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
library(pheatmap)
plot_GO_heatmap<-function(selected_tissue,GO_anno,GO_out){
  GO_tissue=list()
  for(tissue in selected_tissue){
    GO_in=GO_out[[tissue]]
    GO_in=fastDoCall('rbind',lapply(GO_in,function(x) {
      x=x[[GO_anno]]
      x$sig_num=sum(x$FDR<=0.1)
      return(x)
    }))

    
    GO_in=GO_in[,.(GO.ID,Term,classicFisher,FDR,FC,cluster,sig_num)]
    GO_in_top=do.call(c,lapply(1:10,function(x) {
      
      return(GO_in[cluster==x][FC>=1.5&FDR<=0.1][order(FDR,-FC,decreasing=F)][1:5]$GO.ID)
    }))
    GO_in=GO_in[GO.ID%in%GO_in_top]
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
    c2 <- brewer.pal(10,'Set3')
    names(c2) <- 1:10
    breaksList = seq(-1, 1, by = 0.01)
    colann= data.frame(cluster=as.character(1:10))
    #pdf(paste0('../downstream/output/graphs/Figure6/GO_', tissue,'_',GO_anno,'_FC.pdf'),width=25,height=14)
    pheatmap(scalematrix(GO_in_main),cluster_rows =F,cluster_cols = F,
             show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR,border_color = NA,
             color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),
             filename=paste0('../downstream/output/graphs/Figure6/all_regions_chromHMM/',GO_anno,'/GO_', tissue,'_',GO_anno,'_FC_chromHMM.pdf'),
             cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
             fontsize=30,legend = F,labels_col=col_label$sig_num)
    #dev.off()
    GO_tissue[[tissue]]=GO_in_main
    #,breaks=breaksList,color=col(length(breaksList))
  }
  
}
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME_only",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_MML",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_MML_only",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME_MML",GO_out)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_non_NME_non_MML",GO_out)
