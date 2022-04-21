
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
MML_Ecker_read<-function(fn){
  fn_split=strsplit(gsub('.*/|.bed.gz','',fn),'_')[[1]]
  tissue=fn_split[2]
  stage=fn_split[1]
  rep=fn_split[3]
  if(stage=='P0'){stage='day0'}else{
    stage=paste(gsub('E','day',stage),'_5',sep='')
    
  }
  sample_in=paste(tissue,stage,rep,sep='-')
  bed_in=fread(fn)
  bed_in=bed_in[,c(1,2,3,11)]
  colnames(bed_in)=c('chr','start','end','MML')
  bed_in=bed_in[,list(region=paste0(chr,':',start,'-',end),MML=MML/100)]
  bed_in$Sample=sample_in
  return(bed_in[grepl(paste0('chr',c(1:20,"X"),':'),region)])
}
MML_fn_check<-function(fn,MML){
  tt1=proc.time()[[3]]

  cat('start processing:',sample_in,'\n')
  MMLIn=MML[,sample_in]
  MML=convert_GR(ronames(MML),direction="GR")  
  olap=findOverlaps(bed_in,MML)
  olap_df=data.frame(MML_loc=subjectHits(olap),bed_in_M=bed_in$M[queryHits(olap)],MML_M=MMLIn[subjectHits(olap)])
  olap_df_agg=aggregate(olap_df[,c(2,3)],by=list(olap_df$MML_loc),mean)
  cat('finish processing:',sample_in,'in:',proc.time()[[3]]-tt1,'\n')
  return(list(sample_in,cor.test(olap_df_agg$bed_in_M,olap_df_agg$MML_M)))
  
}
#checking methylation value if they're same from bed file
MML_in=readRDS(MML_matrix_file)
#check MML for each tissue between Ecker and ours
Ecker_MML_dir='../downstream/input/mouse_analysis/MML_bed_Ecker/'
MML_Ecker=mclapply(paste0(Ecker_MML_dir,dir(Ecker_MML_dir,pattern='bed.gz')),MML_Ecker_read,mc.cores=20)
saveRDS("../downstream/output/mouse_analysis/QC/Ecker_MML_raw.rds")
MML_Ecker=do.call(rbind,MML_Ecker)
MML_Ecker = dcast.data.table(MML_Ecker,region~Sample,fun.aggregate=mean)
cor_out=mclapply(dir('MML_bed',pattern='bed.gz'),MML_fn_check,MML=MML_in,mc.cores=20)
saveRDS('../downstream/output/mouse_analysis/QC/cor_out_pearson_DNase_dedup.rds')
correlation=readRDS('../downstream/output/mouse_analysis/QC/cor_out_pearson_DNase_dedup.rds')
correlation_df=do.call(rbind,lapply(correlation,function(x){
  data.frame(sample=x[[1]],correlation=x[[2]]$estimate,p_value=x[[2]]$p.value)
})
)
ggplot(correlation_df,aes(x=sample,y=correlation ))+geom_point()+ylab('correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
hist(correlation_df$correlation,main='correlation between estimated MML and paper MML',xlab='cor')

#Coverage check
hist(coverage$cov_mean,main='coverage of all sample',xlab='mean coverage')
hist(coverage$cov_sd,main='coverage of all sample',xlab='sd coverage')

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

DMV=read_supp('../downstream/input/mouse_analysis/QC/DMV.xls')#Table S5
# MML_in_DMV=lapply(unique(DMV$sample),function(x){ 
#   MML_in_DMV_sub=subsetByOverlaps(MML_in,DMV[DMV$sample ==x],type='within')
#   MML_in_DMV_sub=data.frame(tissue=MML_in_DMV_sub$tissue,MML=MML_in_DMV_sub$MML,hypo_tissue=x)
#   return(ggplot(MML_in_DMV_sub,aes(x=tissue,y=MML))+geom_boxplot(outlier.shape = NA)+ggtitle(paste('DMV in:',x))+
#            xlab('tissue')+ylab('MML'))
#   
# })
check_agg<-function(gr_in,MML_matrix){
  gr_in=reduce(gr_in)
  MML_in_gr=subsetByOverlaps(MML_matrix,gr_in,type='within')
  MML_in_gr_olap=findOverlaps(MML_in_gr,gr_in,type='within')
  MML_in_gr=cbind(data.frame(idx=subjectHits(MML_in_gr_olap)),mcols(MML_in_gr))
  MML_in_gr_agg=aggregate(MML_in_gr[,-1],by=list(MML_in_gr$idx),function(x) mean(x,na.rm=T))
  
  return(list(MML_in_gr_agg,MML_in_gr_olap,gr_in))
}
col_fun <- colorRampPalette(
  c(
    ##rev(brewer.pal(8,"RdYlBu"))[1:4],
    rev(brewer.pal(8,"RdYlBu"))[1:2],
    #"white",
    rev(brewer.pal(8,"RdYlBu"))[5:8]
  )
)
MML_check_heatmap<-function(MML_in_agg){
  MML_in_agg=MML_in_agg[,-1]
  MML_in_agg_mt=t(as.matrix(MML_in_agg))
  MML_in_agg_rn=strsplit(rownames(MML_in_agg_mt),'\\.')
  MML_in_agg_rn_df=data.frame(tissue=unlist(lapply(MML_in_agg_rn,function(x) x[1])),
                              stage=unlist(lapply(MML_in_agg_rn,function(x) x[2])),
                              biorep=unlist(lapply(MML_in_agg_rn,function(x) x[3])),stringsAsFactors = F)
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
           Rowv=NA,colv=NA,annotation_row =MML_in_agg_rn_df, na.rm = T,scale="none",cluster_cols =FALSE
  )
  
}
MML_in_DMV_agg=check_agg(DMV,MML_in_matrix)

MML_in_DMV_agg=readRDS('../downstream/output/MML_in_DMV_agg_DNase.rds')
pdf('../downstream/output/DMV_heatmap.pdf',width=15)
MML_check_heatmap(MML_in_DMV_agg[[1]])
#hist(subjectHits(MML_in_DMV_agg[[2]]),xlab='Number of regioins',main='Number of regions')
qt_tb=table(subjectHits(MML_in_DMV_agg[[2]]))
hist((qt_tb*250)/width(MML_in_DMV_agg[[3]][as.numeric(names(qt_tb))]),xlab='%',main='Percent DMV covered')
hist(width(MML_in_DMV_agg[[3]][as.numeric(names(qt_tb))]),xlab='bp',main='DMV size')
dev.off()

#check large hypo-CG DMR
hypo_DMR=read_supp('../downstream/input/mouse_analysis/QC/large_hypo_met.xls')#table S3
hypo_DMR_agg=check_agg(hypo_DMR,MML_in_matrix)
hypo_DMR_agg=readRDS('../downstream/output/hypo_DMR_agg.rds')
pdf('../downstream/output/large_hypo_DMR_heatmap.pdf',width=15)
MML_check_heatmap(hypo_DMR_agg[[1]])
#hist(subjectHits(MML_in_DMV_agg[[2]]),xlab='Number of regioins',main='Number of regions')
qt_tb=table(subjectHits(hypo_DMR_agg[[2]]))
hist((qt_tb*250)/width(hypo_DMR_agg[[3]][as.numeric(names(qt_tb))]),xlab='%',main='Percent hypo DMR covered')
hist(width(hypo_DMR_agg[[3]][as.numeric(names(qt_tb))]),xlab='bp',main='hypo DMR size')
dev.off()
#Plot global methylation
mean_MML=readRDS('../downstream/output/tissue_mean_MML.rds')
metadata=strsplit(names(mean_MML),'\\.')
mean_MML_df=data.frame(MML=mean_MML,tissue=unlist(lapply(metadata,function(x) x[1])),
                       stage=unlist(lapply(metadata,function(x) x[2])),
                       replicates=unlist(lapply(metadata,function(x) x[3])),stringsAsFactors = F)
mean_MML_df$stage=factor(mean_MML_df$stage,levels=c("day10_5","day11_5","day12_5","day13_5","day14_5","day15_5","day16_5","day0"))
ggplot(mean_MML_df,aes(x=stage,y=MML,color=tissue,group=tissue))+geom_line(size=1)+ylim(c(0.5,0.8))
ref=as.data.frame(readxl::read_xlsx('../downstream/output/global_dat.xlsx',skip=2))
ref$tissue=tolower(ref$tissue)   
ref$tissue[ref$tissue=='craniofacial']='EFP'
ref$stage[ref$stage=='P0']='day0'
ref$tissue[ref$tissue=='neural tube']='NT'
ref$tissue[ref$tissue=='lung']='Lung'
ref$stage=gsub('\\.','_',ref$stage)
ref$stage=gsub('E','day',ref$stage)
ref=ref[,c(2,3,15,16)]
ref_mt=melt(ref,id=c('tissue','stage'))
colnames(ref_mt)[c(3,4)]=c('rep','MML')
ref_mt$rep=as.numeric(gsub('mCG level r','',ref_mt$rep))
ref_mt$sample=paste(ref_mt$tissue,ref_mt$stage,ref_mt$rep,sep='.')
MML_compare=data.frame(sample=rownames(mean_MML_df),MML_DNase=mean_MML_df$MML)
MML_compare$MML_paper=ref_mt$MML[match(MML_compare$sample,ref_mt$sample)]