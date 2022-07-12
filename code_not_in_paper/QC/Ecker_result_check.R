source('mainFunctions_sub.R')
#Checking result with Ecker's group
cor_rep<-function(gr_in,stat){
  cor_out=list()
  gr_in$tissue_stage=paste(gr_in$tissue,gr_in$stage,sep='-')
  for(ts in unique(gr_in$tissue_stage)){
    gr_rep1=gr_in[gr_in$tissue_stage==ts&gr_in$bioreplicate==1]
    gr_rep2=gr_in[gr_in$tissue_stage==ts&gr_in$bioreplicate==2]
    olap=findOverlaps(gr_rep1,gr_rep2,type='equal')
    gr_rep1=gr_rep1[queryHits(olap)]
    gr_rep2=gr_rep2[subjectHits(olap)]
    cor_out[[ts]]=cor.test(mcols(gr_rep1)[,stat],mcols(gr_rep2)[,stat])}
  return(cor_out)
}
MML_rename<-function(fn){
  fn_split=strsplit(gsub('.*/|.bed.gz','',fn),'_')[[1]]
  tissue=fn_split[2]
  stage=fn_split[1]
  rep=fn_split[3]
  if(stage=='P0'){stage='day0'}else{
    stage=paste(gsub('E','day',stage),'_5',sep='')
    
  }
  return(paste(tissue,stage,rep,sep='-'))

}
MML_Ecker_read<-function(fn){
  sample_in=MML_rename(fn)
  bed_in=fread(fn)
  bed_in=bed_in[,c(1,2,3,11)]
  colnames(bed_in)=c('chr','start','end','MML')
  bed_in=bed_in[,list(region=paste0(chr,':',start,'-',end),MML=MML/100)]
  colnames(bed_in)[2]=sample_in
  
  return(bed_in[grepl(paste(paste0('chr',c(1:20,"X"),':'),collapse='|'),region)])
}
MML_fn_check<-function(sample_in,MML){
  Ecker_MML_dir='../downstream/input/mouse_analysis/MML_bed_Ecker/'
  tt1=proc.time()[[3]]
  cat('start processing:',sample_in,'\n')
  sample_sp=unlist(strsplit(sample_in,'\\.'))
  fn=dir(Ecker_MML_dir,pattern=paste0(sample_sp[2],"_",sample_sp[1]))
  if(length(fn)>0){
  #Read in files
  bed_in=lapply(paste0(Ecker_MML_dir,fn),MML_Ecker_read)
  bed_in=lapply(bed_in,function(x) {
      colnames(x)[2]="MML"
      return(x)
  })
   bed_in=fastDoCall('rbind',bed_in)
  olap=findOverlaps(convert_GR(bed_in$region,direction="GR")  ,convert_GR(MML$region,direction="GR"))
  olap_df=data.table(region=MML$region[subjectHits(olap)],bed_in_M=bed_in$MML[queryHits(olap)],MML=MML[[sample_in]][subjectHits(olap)])
  olap_df_agg=olap_df[,list(MML_Ecker=mean(bed_in_M,na.rm=T)),by=list(region,MML)]
  olap_df_agg$Sample=sample_in
  cat('finish processing:',sample_in,'in:',proc.time()[[3]]-tt1,'\n')
  return(olap_df_agg)
  }
  
}
MML_regionMean<-function(dtIn,olapIn){
  dtIn=dtIn[subjectHits(olapIn)]
  dtIn$qt=queryHits(olapIn)
  sp=colnames(dtIn)[2]
  dtIn=dtIn[,list(meanMML=mean(get(sp))),by=list(qt)]
  colnames(dtIn)[2]=sp
  return(dtIn)
}

#check MML for each tissue between Ecker and ours
MML_in=readRDS(MML_matrix_file)
MML_in=convert_GR(MML_in,direction="DT")
Ecker_MML_dir='../downstream/input/mouse_analysis/MML_bed_Ecker/'
MML_in_sp=lapply(colnames(MML_in),MML_fn_check,MML=MML_in)
saveRDS(MML_in_sp,'../downstream/output/mouse_analysis/QC/MML_merged_Ecker.rds')
correlation = lapply(MML_in_sp[!unlist(lapply(MML_in_sp,is.null))],function(x) {
  out=cor.test(x$MML,x$MML_Ecker)
  return(list(cor_out=out,dt_out=data.table(Sample=x$Sample[1],correlation=out$estimate,p_value=out$p.value)))
})
saveRDS(correlation,'../downstream/output/mouse_analysis/QC/cor_out_MML_Ecker.rds')
correlation=readRDS('../downstream/output/mouse_analysis/QC/cor_out_MML_Ecker.rds')
correlation_df=do.call(rbind,lapply(correlation,function(x) x[[2]]))
pdf('../downstream/output/mouse_analysis/QC/correlation_Ecker.pdf')
ggplot(correlation_df,aes(x=Sample,y=correlation ))+geom_point()+ylab('correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
hist(correlation_df$correlation,main='correlation between estimated MML and paper MML',xlab='cor')
dev.off()

pdf('../downstream/output/mouse_analysis/QC/correlation_Ecker_0_1.pdf')
ggplot(correlation_df,aes(x=Sample,y=correlation ))+geom_point()+ylab('correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylim(c(0.9,1))
hist(correlation_df$correlation,main='correlation between estimated MML and paper MML',xlab='cor',xlim=c(0.9,1))
dev.off()


#check DMV
library(readxl)

read_supp<-function(file_in){
  sheets=readxl::excel_sheets(file_in)
  df_out <- lapply(sheets, function(X) {sheets_in=readxl::read_excel(file_in, sheet = X,skip=2)
  colnames(sheets_in)=c('chr','start','end')
  sheets_in=makeGRangesFromDataFrame(sheets_in)
  sheets_in$sample=X
  return(sheets_in)})
  df_out=do.call('c',df_out)
  df_out$sample[df_out$sample=="craniofacial"]='EFP'
  df_out$sample[df_out$sample=="neural tube"]='NT'
  df_out$sample[df_out$sample=="lung"]='Lung'
  return(df_out)
}


# MML_in_DMV=lapply(unique(DMV$sample),function(x){ 
#   MML_in_DMV_sub=subsetByOverlaps(MML_in,DMV[DMV$sample ==x],type='within')
#   MML_in_DMV_sub=data.frame(tissue=MML_in_DMV_sub$tissue,MML=MML_in_DMV_sub$MML,hypo_tissue=x)
#   return(ggplot(MML_in_DMV_sub,aes(x=tissue,y=MML))+geom_boxplot(outlier.shape = NA)+ggtitle(paste('DMV in:',x))+
#            xlab('tissue')+ylab('MML'))
#   
# })
check_agg<-function(gr_in,MML_matrix){
  #gr_in=reduce(gr_in)
  #MML_in_gr=subsetByOverlaps(MML_matrix,gr_in,type='within')
  MML_in_gr_olap=findOverlaps(MML_matrix,gr_in,type='within')
  MML_in_dt=cbind(data.table(idx=convert_GR(gr_in,direction="DT")$region[subjectHits(MML_in_gr_olap)]),as.data.table(mcols(MML_matrix)[queryHits(MML_in_gr_olap),]))
  #Select region with at least 50% region covered
  MML_in_dt_ft =data.table(idx=convert_GR(gr_in,direction="DT")$region[subjectHits(MML_in_gr_olap)],
                width_region=width(gr_in)[subjectHits(MML_in_gr_olap)],width_CPEL=width(MML_matrix)[queryHits(MML_in_gr_olap)])
  MML_in_dt_ft=MML_in_dt_ft[,list(covered_percent=sum(width_CPEL)/mean(width_region)),by=list(idx)]
  MML_in_dt_agg=MML_in_dt[idx %in% MML_in_dt_ft[covered_percent>=0.5]$idx,lapply(.SD,mean,na.rm=T),by=idx]
  return(list(MML_in_dt_agg,MML_in_gr_olap,gr_in))
}
col_fun <- colorRampPalette(
  c(
    ##rev(brewer.pal(8,"RdYlBu"))[1:4],
    rev(brewer.pal(8,"RdYlBu"))[1:2],
    #"white",
    rev(brewer.pal(8,"RdYlBu"))[5:8]
  )
)
MML_check_heatmap<-function(MML_in_agg,fn,cluster_cols=FALSE,annotation_col=NA){
  MML_in_agg=MML_in_agg[,-1]
  coln=rownames(MML_in_agg)
  MML_in_agg_mt=t(as.matrix(MML_in_agg))
  colnames(MML_in_agg_mt)=coln
  MML_in_agg_rn=strsplit(rownames(MML_in_agg_mt),'\\.')
  MML_in_agg_rn_df=data.frame(tissue=unlist(lapply(MML_in_agg_rn,function(x) x[1])),
                              stage=unlist(lapply(MML_in_agg_rn,function(x) x[2])),stringsAsFactors = F)
  rownames(MML_in_agg_rn_df)=rownames(MML_in_agg_mt)
  tissue_col= rainbow(length(unique(MML_in_agg_rn_df$tissue)))
  names(tissue_col)=unique(MML_in_agg_rn_df$tissue)
  stage_col= heat.colors(length(unique(MML_in_agg_rn_df$stage)))
  names(stage_col)=unique(MML_in_agg_rn_df$stage)
  biorep_col= gray.colors(2)
  names(biorep_col)=c(1,2)
  ann_colors = list(
    tissue = tissue_col,
    stage =stage_col,
    biorep =biorep_col
  )
  
  pheatmap(MML_in_agg_mt, labels_col = "", labels_row = "",col=col_fun(1024),cluster_rows =FALSE,annotation_colors = ann_colors,
           Rowv=NA,colv=NA,annotation_row =MML_in_agg_rn_df, na.rm = T,scale="none",cluster_cols =cluster_cols,filename=fn,annotation_col=annotation_col
  )
  
}
MML_in_matrix=readRDS(MML_matrix_file)
DMV=read_supp('../downstream/input/mouse_analysis/QC/DMV.xls')#Table S5
MML_in_DMV_agg=check_agg(DMV,MML_in_matrix)
MML_in_DMV_agg[[1]] = MML_in_DMV_agg[[1]][,!grepl('P0',colnames(MML_in_DMV_agg[[1]])),with=F]
MML_check_heatmap(MML_in_DMV_agg[[1]],fn="../downstream/output/mouse_analysis/QC/Ecker_DMV.pdf")
# MML_in_DMV_agg[[1]][which(rowMeans(MML_in_DMV_agg[[1]][,-1])>0.8),]
# #hist(subjectHits(MML_in_DMV_agg[[2]]),xlab='Number of regioins',main='Number of regions')
# qt_tb=table(subjectHits(MML_in_DMV_agg[[2]]))
# hist((qt_tb*250)/width(MML_in_DMV_agg[[3]][as.numeric(names(qt_tb))]),xlab='%',main='Percent DMV covered')
# hist(width(MML_in_DMV_agg[[3]][as.numeric(names(qt_tb))]),xlab='bp',main='DMV size')
# dev.off()

#check large hypo-CG DMR
hypo_DMR=read_supp('../downstream/input/mouse_analysis/QC/large_hypo_met.xls')#table S3
hypo_DMR_agg=check_agg(hypo_DMR,MML_in_matrix)
hypo_DMR_dt=convert_GR(hypo_DMR,direction="DT")
hypo_DMR_dt=hypo_DMR_dt[order(sample)]
hypo_DMR_dt=hypo_DMR_dt[sample%in% gsub('\\..*','',colnames(hypo_DMR_agg[[1]]))]
hypo_DMR_agg[[1]]=hypo_DMR_agg[[1]][match(hypo_DMR_dt$region,idx)]
hypo_DMR_agg[[1]]=hypo_DMR_agg[[1]][,!grepl("liver|P0",colnames(hypo_DMR_agg[[1]])),with=F]
rownames(hypo_DMR_agg[[1]])=1:nrow(hypo_DMR_agg[[1]])
col_ann=data.frame(tissue=hypo_DMR_dt$sample)
rownames(col_ann)=rownames(hypo_DMR_agg[[1]])
MML_check_heatmap(hypo_DMR_agg[[1]],fn="../downstream/output/mouse_analysis/QC/large_hypo_DMR_heatmap.pdf",annotation_col=col_ann).

#hist(subjectHits(MML_in_DMV_agg[[2]]),xlab='Number of regioins',main='Number of regions')
qt_tb=table(subjectHits(hypo_DMR_agg[[2]]))
pdf("../downstream/output/mouse_analysis/QC/large_hypo_DMR_coverage.pdf")
  percent_covered=(qt_tb*250)/width(hypo_DMR_agg[[3]][as.numeric(names(qt_tb))])*100
  percent_covered[percent_covered>100]=100
  hist(percent_covered,xlab='%',main='Percent hypo DMR covered')
  hist(width(hypo_DMR_agg[[3]][as.numeric(names(qt_tb))]),xlab='bp',main='hypo DMR size')
dev.off()
#Plot global methylation
mean_MML=colMeans(as.matrix(mcols(MML_in_matrix)),na.rm=T)
metadata=strsplit(gsub('-all','',names(mean_MML)),'\\-')
mean_MML_df=data.table(MML=mean_MML,tissue=unlist(lapply(metadata,function(x) x[1])),
                       stage=unlist(lapply(metadata,function(x) x[2])),
                       stringsAsFactors = F)
mean_MML_df=mean_MML_df[stage!="P0"]
mean_MML_df$stage=factor(mean_MML_df$stage,levels=c("E10.5","E11.5","E12.5","E13.5","E14.5","E15.5","E16.5"))
pdf("../downstream/output/mouse_analysis/QC/MML_stage.pdf")
ggplot(mean_MML_df,aes(x=stage,y=MML,color=tissue,group=tissue))+geom_line(size=1)+ylim(c(0.5,0.85))
dev.off()
# ref=as.data.frame(readxl::read_xlsx('../downstream/output/global_dat.xlsx',skip=2))
# ref$tissue=tolower(ref$tissue)   
# ref$tissue[ref$tissue=='craniofacial']='EFP'
# ref$stage[ref$stage=='P0']='day0'
# ref$tissue[ref$tissue=='neural tube']='NT'
# ref$tissue[ref$tissue=='lung']='Lung'
# ref$stage=gsub('\\.','_',ref$stage)
# ref$stage=gsub('E','day',ref$stage)
# ref=ref[,c(2,3,15,16)]
# ref_mt=melt(ref,id=c('tissue','stage'))
# colnames(ref_mt)[c(3,4)]=c('rep','MML')
# ref_mt$rep=as.numeric(gsub('mCG level r','',ref_mt$rep))
# ref_mt$sample=paste(ref_mt$tissue,ref_mt$stage,ref_mt$rep,sep='.')
# MML_compare=data.frame(sample=rownames(mean_MML_df),MML_DNase=mean_MML_df$MML)
# MML_compare$MML_paper=ref_mt$MML[match(MML_compare$sample,ref_mt$sample)]

#Check coverage
refDat=as.data.table(readxl::read_xlsx('../downstream/input/mouse_analysis/QC/Ecker_coverage.xlsx',skip=2))#sup 1 from Ecker
mouse_coverage=fread('../downstream/output/mouse_analysis/QC/coverage_mouse.csv')
refDat$mergedCov=refDat$`coverage r1` + refDat$`coverage r2`
mouse_coverage$Sample=tolower(mouse_coverage$Sample)
refDat$Sample=tolower(paste0(refDat$tissue,'-',refDat$stage))
refDat$Sample=gsub('neural tube','nt',refDat$Sample)
refDat$Sample=gsub('craniofacial','efp',refDat$Sample)
refDatCov=refDat[match(mouse_coverage$Sample,Sample),list(mergedCov,Sample)]
mouse_coverage$Ecker_coverage=refDatCov$mergedCov
#Loading biological replicates
#MML
dir_comp='../downstream/data/mouse_with_rep/'
MML_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*mml"),
                               read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
MML_in$bioreplicate=unlist(mclapply(strsplit(MML_in$Sample,'_'),function(x) gsub('merged','',x[2]),mc.cores=10))
MML_in$MML=MML_in$score
MML_in$score=NULL
saveRDS(MML_in,'../downstream/output/mouse_analysis/QC/MML_with_replicates.rds')
#NME
NME_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*nme"),read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
NME_in$bioreplicate=unlist(mclapply(strsplit(NME_in$Sample,'_'),function(x) gsub('merged','',x[2]),mc.cores=10))
NME_in$NME=NME_in$score
NME_in$score=NULL
saveRDS(NME_in,'../downstream/output/mouse_analysis/QC/NME_with_replicates.rds')



replication_cor<-function(datIn,statIn){
  datIn=convert_GR(datIn,direction="DT")
  datIn$Sample=gsub("_merged1|merged2_","",datIn$Sample)
  datIn=datIn[N>=2]
  return(do.call(rbind,lapply(unique(datIn$Sample),function(x){
      datCor=dcast.data.table(datIn[Sample==x],region~bioreplicate,value.var=statIn)
      cor=cor.test(datCor$`1`,datCor$`2`)
      return(data.table(Sample=x,cor=cor$estimate,Pval=cor$p.value,lowerCI=cor$conf.int[1],upperCI=cor$conf.int[2]))

      })
    )
  )
}
MML_in=readRDS('../downstream/output/mouse_analysis/QC/MML_with_replicates.rds')
NME_in=readRDS('../downstream/output/mouse_analysis/QC/NME_with_replicates.rds')

MML_in_cor=replication_cor(MML_in,statIn="MML")
NME_in_cor=replication_cor(NME_in,statIn="NME")
saveRDS(MML_in_cor,'../downstream/output/mouse_analysis/QC/MML_in_cor.rds')
saveRDS(NME_in_cor,'../downstream/output/mouse_analysis/QC/NME_in_cor.rds')
fwrite(MML_in_cor,'../downstream/output/mouse_analysis/QC/MML_in_cor.csv')
fwrite(NME_in_cor,'../downstream/output/mouse_analysis/QC/NME_in_cor.csv')

