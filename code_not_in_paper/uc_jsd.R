source('mainFunctions_sub.R')
jsd=readRDS('../downstream/output/jsd_matrix_DNase.rds')
uc=readRDS('../downstream/output/uc_matrix_DNase.rds')
cor_output=list()
for(ts in names(uc)){
  cat('processing:',ts,'\n')
  tt1=proc.time()[[3]]
  jsd_ts=jsd[[ts]]
  uc_ts=uc[[ts]]
  region_inter=intersect(rownames(jsd_ts),rownames(uc_ts))
  cut_num=2000
  split_data=cut(1:length(region_inter),cut_num,label=FALSE)
  tt2=proc.time()[[3]]
  cor_output[[ts]]=do.call(c,lapply(1:cut_num,function(x){
    if(round(x/cut_num,digits=1)==(x/cut_num)){
      cat("Finished:",x/cut_num*100,'%','in ',proc.time()[[3]]-tt2,'\n')
      }
    return(diag(cor(t(jsd_ts[region_inter[split_data==x],colnames(uc_ts)]),t(uc_ts[region_inter[split_data==x],]),method="spearman")))
}))
  
  #pvalue: rcorr 
  cat("Finish processing:",ts,'in ',proc.time()[[3]]-tt1,'\n')
  
}
