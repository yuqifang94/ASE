#Search to negative

#Search to positive
# der=100
# pos_idx=which(den_der$diff_round==0)
# der_max=max(den_der$der_round)
# proportion_cutoff=0.05
# #Search to negative
# while(der>proportion_cutoff&pos_idx<(length(den_der$diff_round)-smooth_window)){
#   der=mean(abs(den_der$der[pos_idx:(pos_idx+smooth_window)]/der_max))
#   pos_idx=pos_idx+1
#   
# }
# 
# 
# 
# tissue_in[cor_mean_round==cutoff]$classification=mod1$classification
# plot(mod1, what = "classification",xlab=paste0(ts,":mean_correlation:",cutoff),main=ts)


#mod1 <- Mclust(tissue_in[cor_mean_round==cutoff]$cor_diff,G=3,modelNames="V")
#plot(fpc::dbscan(tissue_in_sig[cor_mean_round==cutoff]$cor_diff, eps = 0.15, MinPts = 10),df, main = "DBSCAN", frame = FALSE)
# tissue_in$classification=-1

#tissue_in_sig=tissue_in[dNME_FDR<=FDR_cutoff|dMML_FDR<FDR_cutoff]
#tissue_in_non=tissue_in[dNME_FDR>FDR_cutoff&dMML_FDR>FDR_cutoff]
# minimums <- function(x) which(x - data.table::shift(x, 1) < 0  & x - data.table::shift(x, 1, type='lead') < 0)
pdf('../downstream/output/cor_cutoff_all_001_heart_all_V.pdf')


classification_str=c("dMML only","Both","dNME only")
tissue_in$type_out="NA"
tissue_in[classification!=-1]$type_out=classification_str[tissue_in[classification!=-1]$classification]
tissue_in[dNME_FDR>FDR_cutoff&dMML_FDR>FDR_cutoff]$type_out="None"


png('../downstream/output/cor_heart_cluster_MA_all_V.png')
ggplot(tissue_in,aes(x=cor_mean ,y=cor_diff    ,color=type_out,fill=type_out))+geom_point(alpha=0.1)
dev.off()


# tissue_in_non$type_out="None"
# tissue_in_sig$type_out="NA"
# tissue_in_sig[classification!=-1]$type_out=classification_str[tissue_in_sig[classification!=-1]$classification]
# tissue_out=rbind(tissue_in_non,tissue_in_sig)
png('../downstream/output/cor_heart_cluster.png')
ggplot(tissue_out,aes(x=dMML_cor,y=dNME_cor,color=type_out,fill=type_out))+geom_point(alpha=0.1)
dev.off()

png('../downstream/output/cor_heart_cluster_MA.png')
ggplot(tissue_out,aes(x=cor_mean ,y=cor_diff    ,color=type_out,fill=type_out))+geom_point(alpha=0.1)
dev.off()
#For non cluster do it manually



dev.off()



density_smooth=matrixSmooth(fhat.exp$z_trans,passes=1)
#Brute force scanning for valleys
#Define valleys: sudden change to neighboring elements
nr=nrow(density_smooth)
nc=ncol(density_smooth)

gtp=2#proportion larger
valley_out=matrix(data=0,nrow=nr,ncol=nc)
valley_out_loc=list()
for(i in 1:nr){
  for (j in 1:nc){
    #Points of interest
    pt_in=density_smooth[i,j]
    #Looking for surrounding points
    i_sur=c(i-1,i,i+1)
    j_sur=c(j-1,j,j+1)
    i_sur=i_sur[i_sur>0&i_sur<=nr]
    j_sur=j_sur[j_sur>0&j_sur<=nc]
    #include itself is ok
    pt_sur=density_smooth[i_sur,j_sur]
    diff_sur=pt_sur-pt_in
    if(any(diff_sur>=0.3)){
      #potential valley
      valley_out[i,j]=1
      valley_out_loc=c(valley_out_loc,list(c(i,j)))
    }
  }
  
  
}



jpeg('../downstream/output/dmml_dnme_cor/dmml_dnme_cor_density.jpg')
valley_out_loc_x=fhat.exp$x[unlist(lapply(valley_out_loc,function(x) x[1]))]
valley_out_loc_y=fhat.exp$y[unlist(lapply(valley_out_loc,function(x) x[2]))]
filled.contour(x = fhat.exp$x, y = fhat.exp$y, z = density_smooth,
               plot.axes={points(valley_out_loc_x,valley_out_loc_y)}) 
dev.off()
#Using mle on all data below FDR=0.25
set.seed(12345)
FDR_cor_cutoff=0.25
contour_data_out=list()
estimage_out=list()
x_seq <- seq(-1, 1, by=0.001)
y_seq <- seq(-1, 1, by=0.001)
for(cor_mv_norm in cor_dt){
  tissue=unique(cor_mv_norm$tissue)
  cor_mv_norm_trunc=cor_mv_norm[(dMML_FDR>FDR_cor_cutoff)&(dNME_FDR>FDR_cor_cutoff)]
  #cor_mv_norm_trunc=cor_mv_norm[(dMML_cor<=0)&(dNME_cor<=0)]
  #Use range instead of max due to decimal problem
  upper=c(range(cor_mv_norm_trunc$dMML_cor)[2],range(cor_mv_norm_trunc$dNME_cor)[2])
  #There's some digit problem?
  cor_mv_norm_trunc[dMML_cor<(-1)]$dMML_cor=-1
  cor_mv_norm_trunc[dNME_cor<(-1)]$dNME_cor=-1
  cor_mv_norm_trunc[dMML_cor>(1)]$dMML_cor=1
  cor_mv_norm_trunc[dNME_cor>(1)]$dNME_cor=1
  lower=c(-1,-1)
  trunc_matrix=as.matrix(cor_mv_norm_trunc[,list(dMML_cor,dNME_cor)])
  mle.fit1 <- mle.tmvnorm(trunc_matrix, lower=lower, upper=upper,method = "L-BFGS-B",lower.bounds=c(-1,-1),upper.bounds = c(1,1))
  mle.fit1_sum=summary(mle.fit1)
  
  mu= mle.fit1_sum@coef[c('mu_1','mu_2'),"Estimate"]
  sigma= matrix(mle.fit1_sum@coef[c('sigma_1.1','sigma_1.2',"sigma_1.2","sigma_2.2"),"Estimate"],nrow=2,byrow=T)
  
  contour_data=as.data.table(expand.grid(x_seq,y_seq))
  colnames(contour_data)=c("x_seq","y_seq")
  contour_data$MD=mahalanobis(as.matrix(contour_data[,.(x_seq,y_seq)]),center=mu,cov=sigma)
  contour_data$density=NA
  # for (i in 1:nrow(contour_data)){
  # 
  #     #prepare for contour
  #   contour_data$density[i] = dmvnorm(contour_data[i,.(x_seq,y_seq)], mean=mu,sigma=sigma)
  # 
  # }
  # jpeg(paste0('../downstream/output/dmml_dnme_cor/truncated_nrom_small_0_',tissue,'.jpg'))
  # print(ggplot(cor_mv_norm[(dMML_FDR>FDR_cor_cutoff)&(dNME_FDR>FDR_cor_cutoff)],aes(x=dMML_cor,y=dNME_cor))+ geom_point(alpha=0.01)+
  #   geom_vline(xintercept = upper[1],size=1)+geom_hline(yintercept = upper[2],size=1)+
  #   geom_contour(mapping=aes(x=x_seq,y=y_seq,z=MD),data=contour_data,breaks=seq(1,5,1),color='red')+ 
  #   geom_label_contour(mapping=aes(x=x_seq,y=y_seq,z=MD),data=contour_data,breaks=seq(1,5,1))+
  #   geom_point(mapping=aes(x=mu1,y=mu2),data=data.frame(mu1=mu[1],mu2=mu[2]),size=2,color='blue'))
  #  dev.off()
  contour_data_out[[tissue]]=contour_data
  estimage_out[[tissue]]=list(mu=mu,sigma=sigma)
  gc()
  #geom_contour_fill()+ scale_fill_gradient2(high = "blue",low="white")+
}
#Use distance cutoff: 4
MD_cutoff=4
FDR_cutoff=0.25
cor_mv_norm_MD_out=list()
col_theme=c("red","blue","green","purple")
sig_levels=c("dMML_only","Both","dNME_only","None")
for(cor_mv_norm_MD in cor_dt){
  tissue=unique(cor_mv_norm_MD$tissue)
  estimage_tissue=estimage_out[[tissue]]
  cor_mv_norm_MD$MD=mahalanobis(as.matrix(cor_mv_norm_MD[,.(dMML_cor,dNME_cor)]),center=estimage_tissue$mu,cov=estimage_tissue$sigma)
  cor_mv_norm_MD$cor_type="NA"
  cor_mv_norm_MD[MD<=MD_cutoff&(dMML_FDR<=FDR_cutoff|dNME_FDR<=FDR_cutoff)]$cor_type="Both"
  cor_mv_norm_MD[(dMML_FDR>FDR_cutoff&dNME_FDR>FDR_cutoff)]$cor_type="None"
  cor_mv_norm_MD[MD>MD_cutoff&dNME_FDR<=FDR_cutoff]$cor_type="dNME_only"
  cor_mv_norm_MD[MD>MD_cutoff&dMML_FDR<=FDR_cutoff]$cor_type="dMML_only"
  dNME_cutoff=min(cor_mv_norm_MD[dNME_FDR<=FDR_cutoff]$dNME_cor)
  dMML_cutoff=min(cor_mv_norm_MD[dMML_FDR<=FDR_cutoff]$dMML_cor)
  
  cor_mv_norm_MD$cor_type=factor(cor_mv_norm_MD$cor_type,levels=sig_levels)
  
  
  png(paste0('../downstream/output/dmml_dnme_cor/cor_type_',tissue,'_025_4_md.png'),width=7,height=7,units = "in",res=96)
  print(
    ggplot(cor_mv_norm_MD,aes(x=dMML_cor,y=dNME_cor))+
      geom_point(aes(color=cor_type),alpha=0.1)+
      theme(legend.position = "bottom")+ 
      scale_color_manual(breaks= sig_levels, values = col_theme)+
      geom_hline(yintercept = dNME_cutoff,size=1,color="black")+geom_vline(xintercept = dMML_cutoff,size=1,color="black")+
      geom_contour(mapping=aes(x=x_seq,y=y_seq,z=MD),data=contour_data_out[[tissue]],breaks=MD_cutoff,color='black')
  )
  
  dev.off()
  cor_mv_norm_MD_out[[tissue]]=cor_mv_norm_MD
}

cor_mv_norm_MD_out=fastDoCall('rbind',cor_mv_norm_MD_out)
cor_mv_norm_MD_out_plot=cor_mv_norm_MD_out[,.N,by=c("tissue","cor_type")][,list(prop=N/sum(N),cor_type=cor_type),by=tissue]
pdf('../downstream/output/dmml_dnme_cor/bar_plot_025_4_trunc.pdf')
ggplot(cor_mv_norm_MD_out_plot,aes(x=tissue,y=prop))+geom_bar(stat="identity",aes(group=cor_type,fill=cor_type))+theme(legend.position = 'bottom')+
  scale_fill_manual(breaks= sig_levels, values = col_theme)
dev.off()
library(fpc)
for(i in 1:length(cor_dt)){
  df=cor_dt[[i]][,.(dMML_cor,dNME_cor)]
  tissue=unique(cor_dt[[i]]$tissue)
  for(eps in seq(0.001,0.01,0.001)){
    png(paste0('../downstream/output/dmml_dnme_cor/DBscan_',eps,'_',tissue,'.png'),width=7,height=7,units = "in",res=96)
    print(plot(fpc::dbscan(df, eps = eps, MinPts = 10),df, main = "DBSCAN", frame = FALSE))
    dev.off()
    
  }
}
library(mclust)
cor_dt_mc_out_all=list()
set.seed(123)
for(i in 1:length(cor_dt)){
  cor_dt_ts=cor_dt[[i]]
  cor_dt_mt=as.matrix(cor_dt_ts[,.(dMML_cor,dNME_cor)])
  tissue=unique(cor_dt_ts$tissue)
  cor_dt_mc_out=list()
  for(g in 3:12){
    cor_dt_mc_out[[as.character(g)]]=Mclust(cor_dt_mt,G=g)
    png(paste0('../downstream/output/dmml_dnme_cor/mcluster_',tissue,'_',g,'.png'),width=7,height=7,units = "in",res=96)
    plot( cor_dt_mc_out[[as.character(g)]], what = 'classification')
    dev.off()
    
  }
  cor_dt_mc_out_all[[tissue]]=cor_dt_mc_out
  
  
  
  
  plot_GO_heatmap<-function(selected_tissue,GO_out,region_type,enc_type,GO_anno="GO_out_cluster_all",
                            FDR_cutoff=0.05,FC_cutoff=1.5,ptcount=0){
    GO_tissue=list()
    for(tissue in selected_tissue){
      cat('Plotting:',tissue,'\n')
      GO_in=GO_out[[tissue]]
      
      if(!is.null(GO_in)){
        GO_in=fastDoCall('rbind',lapply(GO_in,function(x) {
          if(!is.null(x[[GO_anno]])){
            x=x[[GO_anno]]
            #((sigres[,"Significant"]+ptcount)/(sum(GOdata@allScores[GOdata@feasible]==1)+ptcount))/((sigres[,"Annotated"]+ptcount)/(sum(GOdata@feasible)+ptcount))
            x$FC_raw=x$FC
            x$FC=((x$Significant+ptcount)/(x$feasible_allscore+ptcount))/((x$Annotated+ptcount)/(x$feasible+ptcount))
            x$sig_num=sum(x$FDR<=FDR_cutoff&x$FC>=FC_cutoff)
            return(x)
          }
        }))
        
        
        GO_in=GO_in[,.(GO.ID,Term,classicFisher,FDR,FC,FC_raw,cluster,sig_num,p_cond)]
        GO_in_top=do.call(rbind,lapply(1:10,function(x) {
          
          return(GO_in[cluster==x][FC>=FC_cutoff][order(FDR,p_cond,-FC,decreasing=F)][1:5])
        }))
        GO_in=GO_in[GO.ID%in%GO_in_top$GO.ID]
        
        if(nrow(GO_in)>0){
          
          GO_in_main= dcast_matrix(GO_in,"FC",order=T)
          GO_in_FDR= dcast_matrix(GO_in,"FDR")
          # # GO_in$log10FDR=-log10(GO_in$FDR)
          # # GO_in_FDR_log10= dcast_matrix(GO_in,"log10FDR")
          
          GO_in_FDR[GO_in_FDR<=FDR_cutoff]="*"
          GO_in_FDR[GO_in_FDR>FDR_cutoff]=""
          if(nrow(GO_in_FDR)>1){GO_in_FDR=GO_in_FDR[rownames(GO_in_main),]}
          col_label=unique(GO_in[,.(cluster,sig_num)])
          col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
          c2 <- brewer.pal(10,'Set3')
          names(c2) <- 1:10
          #breaksList = seq(-1, 1, by = 0.01)
          colann= data.frame(cluster=as.character(1:10))
          #Using all cluster use GO_in_main[GO_in_top$Term,]
          #if(nrow(GO_in_main)>1){GO_in_main=scalematrix(GO_in_main)}
          #pdf(paste0('../downstream/output/graphs/Figure6/GO_', tissue,'_',GO_anno,'_FC.pdf'),width=25,height=14)
          # #Manually change the color theme
          # pheatmap(GO_in_main[GO_in_top$Term,],cluster_rows =F,cluster_cols = F,
          #          show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR[GO_in_top$Term,],border_color = NA,
          #          color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),number_color = "white",
          #          filename=paste0('../downstream/output/graphs/GO_figures/',tissue,'_',GO_anno,'_FC_pval_raw_all_cluster_manbreak_',region_type,'_',enc_type,'.pdf'),
          #          cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
          #          fontsize=30,legend = T,labels_col=col_label$cluster,gaps_row = seq(5,50,5),
          #          breaks=quantile_breaks(GO_in_main[GO_in_top$Term,]))
          pheatmap(scalematrix(GO_in_main),cluster_rows =F,cluster_cols = F,
                   show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR,border_color = NA,
                   color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),number_color = "white",
                   filename=paste0('../downstream/output/graphs/GO_figures/single_sample/',tissue,'_FC_pval_raw_all_cluster_',region_type,'_',enc_type,'.pdf'),
                   cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
                   fontsize=30,legend = T,labels_col=col_label$cluster)
          #dev.off()
          GO_tissue[[tissue]]=GO_in_main
          #,breaks=breaksList,color=col(length(breaksList))
        }
      }
    }
  }
  plot_GO_heatmap_perm<-function(selected_tissue,GO_anno,GO_out,FC_cutoff=1.5){
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
        
        return(GO_in[cluster==x][FC>=FC_cutoff&FDR<=0.1][order(FDR,-FC,decreasing=F)][1:5]$GO.ID)
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
               filename=paste0('../downstream/output/graphs/Figure6/all_regions_chromHMM_perm/',GO_anno,'/GO_', tissue,'_',GO_anno,'_FC_chromHMM.pdf'),
               cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
               fontsize=30,legend = F,labels_col=col_label$sig_num)
      #dev.off()
      GO_tissue[[tissue]]=GO_in_main
      #,breaks=breaksList,color=col(length(breaksList))
    }
    
  }