cor_dt_all=lapply(names(dNME_cor),cor_dt_preprocessing,dMML_cor=dMML_cor,
                  dNME_cor=dNME_cor,dmml_perm=dmml_perm,dnme_perm=dnme_perm,
                  filtered=FALSE)
names(cor_dt_all)=names(dNME_cor)
saveRDS(cor_dt_all,'../downstream/output/correlation/correlation_dt_unfiltered.rds')
cor_dt_all=readRDS('../downstream/output/correlation/correlation_dt_unfiltered.rds')
tissue_out_all=lapply(names(cor_dt_all),correlation_processing,cor_dt=cor_dt_all,filtered=F)
names(tissue_out_all)=names(cor_dt_all)
saveRDS(tissue_out_all,'../downstream/output/correlation/tissue_out_unfiltered.rds')   

correlation_processing<-function(ts,cor_dt,filtered=F,density_plot=T,FDR_cutoff=0.2,quant=0.25){
  cat('Processing:',ts,'\n')
  tissue_in=cor_dt[[ts]]
  #generate MA plot
  tissue_in$cor_diff=tissue_in$dNME_cor -tissue_in$dMML_cor
  tissue_in$cor_mean=(tissue_in$dNME_cor +tissue_in$dMML_cor)/2
  
  #There might be some floating point issue showing mean ==1 but stored as >1
  tissue_in[cor_mean>1]$cor_mean=1
  
  if(density_plot==TRUE){
    cat("Generating density plot\n")
    tt1=proc.time()[[3]]
    tissue_in_mt=as.matrix(tissue_in[,list(dMML_cor,dNME_cor)])
    print(head(tissue_in_mt))
    if(filtered==T){
      
      row_select=1:nrow(tissue_in_mt)
    }else{
      
      row_select=sample(1:nrow(tissue_in_mt),round(1/5*nrow(tissue_in_mt)),replace = F)
      
    }
    bw_opt_dmml_cor=ucv(tissue_in_mt[row_select,'dMML_cor'])
    bw_opt_dnme_cor=ucv(tissue_in_mt[row_select,'dNME_cor'])
    
    # saveRDS(bw_opt_dmml_cor,paste0('../downstream/output/correlation/bw_opt_dmml_cor_',ts,'_',filtered,'.rds'))
    # saveRDS(bw_opt_dnme_cor,paste0('../downstream/output/correlation/bw_opt_dnme_cor_',ts,'_',filtered,'.rds'))
    fhat.exp.dMML.dNME=kde2d(tissue_in_mt[,'dMML_cor'], tissue_in_mt[,'dNME_cor'],
                             h=c(bw_opt_dmml_cor,bw_opt_dnme_cor),
                             n=100,
                             lims=c(c(-1,1),c(-1,1)))
    fhat.exp.dMML.dNME$z_trans=log( fhat.exp.dMML.dNME$z+0.1)
    
    pdf(paste0('../downstream/output/correlation/',ts,'_ggplot_raw_density_',filtered,'.pdf'))
    filled.contour(x = fhat.exp.dMML.dNME$x, y = fhat.exp.dMML.dNME$y, z = fhat.exp.dMML.dNME$z_trans,
                   xlab='dMML cor',ylab='dNME cor')
    dev.off()
    
    
    tissue_in_mt=as.matrix(tissue_in[,list(cor_mean,cor_diff)])
    print(head(tissue_in_mt))
    bw_opt_mean=ucv(tissue_in_mt[row_select,'cor_mean'])
    bw_opt_diff=ucv(tissue_in_mt[row_select,'cor_diff'])
    # saveRDS(bw_opt_mean,paste0('../downstream/output/correlation/bw_opt_mean_',ts,'_',filtered,'.rds'))
    # saveRDS(bw_opt_diff,paste0('../downstream/output/correlation/bw_opt_diff_',ts,'_',filtered,'.rds'))
    fhat.exp.MA=kde2d(tissue_in_mt[,'cor_mean'], tissue_in_mt[,'cor_diff'],h=c(bw_opt_mean,bw_opt_diff),
                      lims=c(c(-1,1),c(-2,2)),n=100)
    
    fhat.exp.MA$z_trans=log(fhat.exp.MA$z+0.1)
    pdf(paste0('../downstream/output/correlation/',ts,'_raw_MA_',filtered,'.pdf'))
    filled.contour(x = fhat.exp.MA$x, y = fhat.exp.MA$y, z = fhat.exp.MA$z_trans,xlab='mean correlation',ylab='correlation difference')
    dev.off()
    cat("Finish generating density plot in:",proc.time()[[3]]-tt1,'\n')
  }
  
  cat("Finding cutoffs\n")
  tissue_in$cor_mean_round=round(tissue_in$cor_mean*2,digits=1)/2
  cor_cutoffs=data.table()
  pdf(paste0('../downstream/output/correlation/cor_cutoff_all_05_',ts,'_all_',filtered,'.pdf'))
  for(cutoff in unique(tissue_in$cor_mean_round)){
    if(sum(tissue_in$cor_mean_round==cutoff)>=10){
      #Find saddle point for density difference
      cor_diff_den=density(tissue_in[cor_mean_round==cutoff]$cor_diff,bw="SJ")
      den_der=der_calc(cor_diff_den$x,cor_diff_den$y)
      den_der$diff=den_der$x
      den_der$density=den_der$y
      
      
      #Function to calculate derivatives
      
      neg_cutoff=der_flat_finder(den_der$der,den_der$diff,den_der$density,direction=-1,quant=0.25)
      pos_cutoff=der_flat_finder(den_der$der,den_der$diff,den_der$density,direction=1,quant=0.25)
      plot(den_der$diff,den_der$der,xlab="dNME-dMML cor",ylab="1st derivitative",main=cutoff)
      abline(v=pos_cutoff$x_out)
      abline(v=neg_cutoff$x_out)
      plot(cor_diff_den$x,cor_diff_den$y,xlab="dNME-dMML cor",ylab="density",main=cutoff)
      abline(v=pos_cutoff$x_out)
      abline(v=neg_cutoff$x_out)
      
      
      #Assigning the catogries based on cutoff
      diff_cutoff_dMML_only=neg_cutoff$x_out
      diff_cutoff_dNME_only=pos_cutoff$x_out
      cor_cutoffs=rbind(cor_cutoffs,data.table(cor_mean_round=cutoff,
                                               diff_cutoff_dMML_only=diff_cutoff_dMML_only,
                                               diff_cutoff_dNME_only=diff_cutoff_dNME_only))
      
    }
    
  }
  dev.off()
  #Smoothing cutoffs using weighted loess
  cor_cutoffs=cor_cutoffs[order(cor_mean_round)]
  
  weight_loess= ecdf(tissue_in$cor_mean)(cor_cutoffs$cor_mean_round)
  dMML_fit=loess(diff_cutoff_dMML_only  ~ cor_mean_round,data=cor_cutoffs,span=0.75,
                 weights=weight_loess)
  dNME_fit=loess(diff_cutoff_dNME_only ~ cor_mean_round,data=cor_cutoffs,span=0.75,
                 weights=weight_loess)
  tissue_in$diff_cutoff_dMML_only=predict(dMML_fit,data.table(cor_mean_round=tissue_in$cor_mean))
  tissue_in$diff_cutoff_dNME_only=predict(dNME_fit, data.table(cor_mean_round=tissue_in$cor_mean))
  #smoothing plot
  png(paste0('../downstream/output/correlation/cor_',ts,'_cluster_all_cor_dMML_sm_quant_weighted_',filtered,'.png'))
  plot(cor_cutoffs$cor_mean_round,cor_cutoffs$diff_cutoff_dMML_only,xlab="mean correlation",ylab="correlation differnce cutoff",xlim=c(-1,1))
  lines(cor_cutoffs$cor_mean_round,predict(dMML_fit,data=cor_cutoffs))
  dev.off()
  png(paste0('../downstream/output/correlation/cor_',ts,'_cluster_all_cor_dNME_sm_quant_weighted_',filtered,'.png'))
  plot(cor_cutoffs$cor_mean_round,cor_cutoffs$diff_cutoff_dNME_only,xlab="mean correlation",ylab="correlation differnce cutoff")
  lines(cor_cutoffs$cor_mean_round,predict(dNME_fit,data=cor_cutoffs))
  dev.off()
  #Assigning regions
  rm(diff_cutoff_dNME_only)
  rm(diff_cutoff_dMML_only)
  tissue_in$region_type="NA"
  #Both
  tissue_in[cor_diff<diff_cutoff_dNME_only&
              cor_diff>diff_cutoff_dMML_only&
              (dNME_FDR<=FDR_cutoff|dMML_FDR<=FDR_cutoff)]$region_type="Both"
  #One significant the other one <0 are also quantified as only
  #dMML_only
  tissue_in[(dNME_FDR<=FDR_cutoff|dMML_FDR<=FDR_cutoff)&
              (cor_diff<=diff_cutoff_dMML_only|
                 (dMML_cor>0&dNME_cor<0))]$region_type="dMML_only"
  #dNME_only
  tissue_in[(dNME_FDR<=FDR_cutoff|dMML_FDR<=FDR_cutoff)&
              (cor_diff>=diff_cutoff_dNME_only|
                 (dNME_cor>0&dMML_cor<0))]$region_type="dNME_only"
  
  
  #Neither
  tissue_in[dNME_FDR>FDR_cutoff&dMML_FDR>FDR_cutoff]$region_type="Neither"
  
  #dev.off()
  #Smooth cutoffs
  if(filtered==TRUE){alpha_plot=0.1}else{alpha_plot=0.01}
  png(paste0('../downstream/output/correlation/cor_',ts,'_cluster_all_cor_smoothed_quant_weighted_',filtered,'.png'))
  print(ggplot(tissue_in,aes(x=dMML_cor,y=dNME_cor,color=region_type,fill=region_type))+geom_point(alpha=alpha_plot)+
          guides(colour = guide_legend(override.aes = list(alpha = 1))))
  dev.off()
  #MA plot
  png(paste0('../downstream/output/correlation/cor_',ts,'_cluster_all_cor_MA_smoothed_quant_',filtered,'.png'))
  print(ggplot(tissue_in,aes(x=cor_mean ,y=cor_diff ,color=region_type,fill=region_type))+geom_point(alpha=alpha_plot)+
          guides(colour = guide_legend(override.aes = list(alpha = 1))))
  dev.off()
  
  
  if(density_plot==TRUE){
    cat("Generating density plot after finish\n")
    cutoff_dt=dMML_dNME_cutoff_dt(tissue_in,FDR_cutoff)
    if(diff(range(cutoff_dt[region_type=='dNME only']$dMML))>=0.05&diff(range(cutoff_dt[region_type=='dMML only']$dMML))>=0.05){
      cutoff_dt_dNME_only_sm=loess(dNME~dMML,data=cutoff_dt[region_type=='dNME only'])
      cutoff_dt_dMML_only_sm=loess(dNME~dMML,data=cutoff_dt[region_type=='dMML only'])
      pdf(paste0('../downstream/output/correlation/',ts,'_catogry_raw_density_',filtered,'.pdf'))
      filled.contour(x = fhat.exp.dMML.dNME$x, y = fhat.exp.dMML.dNME$y, z = fhat.exp.dMML.dNME$z_trans,xlab='dMML cor',ylab='dNME cor',
                     plot.axes = { 
                       axis(1); 
                       axis(2);
                       lines(cutoff_dt[region_type=='dMML FDR cutoff']$dMML, cutoff_dt[region_type=='dMML FDR cutoff']$dNME, col = "black",lwd=2);
                       lines(cutoff_dt[region_type=='dNME FDR cutoff']$dMML, cutoff_dt[region_type=='dNME FDR cutoff']$dNME, col = "black",lwd=2);
                       lines(cutoff_dt[region_type=='dNME only']$dMML,
                             predict(cutoff_dt_dNME_only_sm,data.table(dMML=cutoff_dt[region_type=="dNME only"]$dMML)),
                             col = "black",lwd=2);
                       lines(cutoff_dt[region_type=='dMML only']$dMML,
                             predict(cutoff_dt_dMML_only_sm,data.table(dMML=cutoff_dt[region_type=="dMML only"]$dMML)),
                             col = "black",lwd=2);
                       points(cutoff_dt[region_type=='dNME only']$dMML,
                              cutoff_dt[region_type=="dNME only"]$dNME,
                              col = "black",lwd=2);
                       points(cutoff_dt[region_type=='dMML only']$dMML,
                              cutoff_dt[region_type=="dMML only"]$dNME,
                              col = "black",lwd=2);
                     }
      )
      dev.off()
      
      cutoff_dt$mean=(cutoff_dt$dMML+cutoff_dt$dNME)/2
      cutoff_dt$diff=(cutoff_dt$dNME-cutoff_dt$dMML)
      cutoff_dt=cutoff_dt[order(mean)]
      pdf(paste0('../downstream/output/correlation/',ts,'_catogry_raw_MA_',filtered,'.pdf'))
      filled.contour(x = fhat.exp.MA$x, y = fhat.exp.MA$y, z = fhat.exp.MA$z_trans,xlab='mean correlation',ylab='correlation difference',
                     plot.axes = { 
                       axis(1); 
                       axis(2);
                       lines(cutoff_dt[region_type=='dMML FDR cutoff']$mean, cutoff_dt[region_type=='dMML FDR cutoff']$diff, col = "black",lwd=2);
                       lines(cutoff_dt[region_type=='dNME FDR cutoff']$mean, cutoff_dt[region_type=='dNME FDR cutoff']$diff, col = "black",lwd=2);
                       lines(cutoff_dt[region_type=='dNME only']$mean, 
                             predict(dNME_fit,data.table(cor_mean_round=cutoff_dt[region_type=="dNME only"]$mean)), col = "black",lwd=2);
                       lines(cutoff_dt[region_type=='dMML only']$mean, 
                             predict(dMML_fit,data.table(cor_mean_round=cutoff_dt[region_type=="dMML only"]$mean)), col = "black",lwd=2);
                     }
      )
      dev.off()
      cat("Finish generating density plot in:",proc.time()[[3]]-tt1,'\n')
    }else{cat("high FDR threshold for density plot")}
  }
  return(tissue_in)
  
}

ggplot(tissue_in[,list(dMML_cor,dNME_cor)],aes(x=dMML_cor,y=dNME_cor))+
  stat_density_2d( geom = "raster",  aes(fill = log(after_stat(density)+0.1)), contour = FALSE,n=200)+
  scale_fill_viridis_c()+geom_abline(slope=1,intercept=0,color='red')
