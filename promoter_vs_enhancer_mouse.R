source('mainFunctions_sub.R')
# Promoter vs enhancer by tissue ------------------------------------------
#max dMML vs dNME in promoter vs enhancer
UC_merge=readRDS('../downstream/output/UC_merge_max_loc.rds')
folder_in="../downstream/input/ts_cluster_0_1/"
folder_out=paste0(folder_in,'region_assigned/')

enhancer=readRDS('../downstream/output/bin_enhancer.rds')
csv_out=data.table()
for(fn in dir(folder_out,pattern="csv")){
  csv_in=fread(paste0(folder_out,fn))
  tissue=gsub('.csv','',fn)
  if(tissue !="NT"){
    csv_in=cbind(csv_in,UC_merge[[tissue]][csv_in$region,c("dMML_max_pair","dNME_max_pair","dMML_max_time","dNME_max_time","UC_max_pair","UC_max_time")])
    olap=findOverlaps(convert_GR(csv_in$region),enhancer)
    csv_in$enhancer=F
    csv_in[queryHits(olap)]$enhancer=TRUE
    csv_in$enhancer_gene="NA"
    csv_in[queryHits(olap)]$enhancer_gene=enhancer$`Target Gene`[subjectHits(olap)]
    csv_in$tissue=tissue
    csv_out=rbind(csv_out,csv_in)
  }
}
#csv_out[enhancer==TRUE &abs(distance)<=2000]
csv_out$states="NA"
csv_out[enhancer==TRUE]$states="enhancers"
csv_out[abs(distance)<=2000]$states="promoters"
#boxplot compare max dNME and dMML
#maximum value

ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dMML_max_pair,fill=states))+geom_boxplot()+ylab("dMML")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dNME_max_pair,fill=states))+geom_boxplot()+ylab("dNME")

#at max UC
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dMML_maxUC ,fill=states))+geom_boxplot()+ylab("dMML")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dNME_maxUC ,fill=states))+geom_boxplot()+ylab("dNME")

#Max UC

#at max UC
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=UC_max_pair ,fill=states))+geom_boxplot()+ylab("UC")

#max dMML/dNME in non-adjacent vs adjacent 

csv_out$dMML_max_time=gsub('E','',gsub('\\.5','',csv_out$dMML_max_time))
csv_out$dMML_max_time_diff=abs(as.numeric(gsub('-.*','',csv_out$dMML_max_time))-as.numeric(gsub('.*-','',csv_out$dMML_max_time)))
csv_out$dNME_max_time=gsub('E','',gsub('\\.5','',csv_out$dNME_max_time))
csv_out$dNME_max_time_diff=abs(as.numeric(gsub('-.*','',csv_out$dNME_max_time))-as.numeric(gsub('.*-','',csv_out$dNME_max_time)))
csv_out$UC_max_time=gsub('E','',gsub('\\.5','',csv_out$UC_max_time))
csv_out$UC_max_time_diff=abs(as.numeric(gsub('-.*','',csv_out$UC_max_time))-as.numeric(gsub('.*-','',csv_out$UC_max_time)))
#boxplot time of maximum value
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dMML_max_time_diff,fill=states))+geom_boxplot()+ylab("dMML time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=dNME_max_time_diff,fill=states))+geom_boxplot()+ylab("dNME time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=tissue,y=UC_max_time_diff,fill=states))+geom_boxplot()+ylab("UC time change")

ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=dMML_max_time_diff,color=states))+geom_density()+ylab("dMML time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=dNME_max_time_diff,color=states))+geom_density()+ylab("dNME time change")
ggplot(csv_out[states%in%c("enhancers","promoters")],aes(x=UC_max_time_diff,color=states))+geom_density()+ylab("dNME time change")