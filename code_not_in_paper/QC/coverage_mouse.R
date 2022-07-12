source('mainFunctions_sub.R')
read_bedGraph<-function(fnIn){
  cat("Processing: ",fnIn,"\n")
  bedIn=fread(fnIn)
  bedIn=bedIn[,c(1:3,7)]
  colnames(bedIn)=c("chr","start","end","coverage")
  bedIn$start=bedIn$start+1
  bedIn$Sample=gsub('-5','.5',gsub('day1','E1',gsub('_','-',gsub('mm10_|_all','',gsub('\\..*','',gsub('.*/','',fnIn))))))
  return(bedIn)
}
avgCpG_cov_tissue<-function(CpGCov,regionIn,grIn=F){
  sp=unique(CpGCov$Sample)
  cat(format(Sys.time(), "%a %b %d %X %Y"),": Start processing coverage for ",sp ,"\n")
  if(!grIn){CpGCov=makeGRangesFromDataFrame(CpGCov,keep.extra.columns=T)}
   #CpGCov = CpGCov[1:100000]#26752510
  olap=as.data.table(findOverlaps(regionIn,CpGCov,minoverlap=2))
  olap$CpGCov=CpGCov[olap$subjectHits]$coverage
  olap = olap[,list(meanCov=mean(CpGCov),N=length(CpGCov)),by=list(queryHits)]
  regionIn=regionIn[olap$queryHits]
  regionIn$coverage=olap$meanCov
  regionIn$NCpG=olap$N
  regionIn$Sample=sp
  cat(format(Sys.time(), "%a %b %d %X %Y"),": Finish processing coverage for ",sp,"\n")
  return(convert_GR(regionIn,direction="DT"))
}
#Loading coverage per CpG
dir_Cov='../downstream/data/coverage_mouse/'
coverage_mouse=mclapply(paste0(dir_Cov,dir(dir_Cov,pattern="day1")),read_bedGraph,mc.cores=20)
saveRDS(coverage_mouse,'../downstream/output/mouse_analysis/QC/mouse_coverage_CpG.rds')
NME_in=readRDS(NME_matrix_file)
GR=granges(NME_in)
rm(NME_in)
#Average coverage per region
coverage_mouse_avgCpG=do.call(rbind,mclapply(coverage_mouse,avgCpG_cov_tissue,regionIn=GR,grIn=F,mc.cores=20))
fwrite(coverage_mouse_avgCpG[,list(coverage=mean(coverage)),by=list(Sample)],'../downstream/output/mouse_analysis/QC/coverage_mouse.csv')

#Coverage check
hist(coverage_mouse_avgCpG$cov_mean,main='coverage of all sample',xlab='mean coverage')
hist(coverage_mouse_avgCpG$cov_sd,main='coverage of all sample',xlab='sd coverage')