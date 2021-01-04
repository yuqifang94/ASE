


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

#checking methylation value if they're same from bed file
MML_in=readRDS('../../../../allele_agnostic_DNase_dedup/MML_agnostic_mouse_DNAase_all_DNase_dedup.rds')
#check MML for each tissue
cor_out_pearson_DNase=mclapply(dir('MML_bed',pattern='bed.gz'),MML_fn_check,MML=MML_in,mc.cores=20)
MML_fn_check<-function(fn,MML){
  tt1=proc.time()[[3]]
  fn_split=strsplit(gsub('.bed.gz','',fn),'_')[[1]]
  tissue=fn_split[2]
  stage=fn_split[1]
  rep=fn_split[3]
  if(stage=='P0'){stage='day0'}else{
    stage=paste(gsub('E','day',stage),'_5',sep='')
    
  }
  sample_in=paste(tissue,stage,rep,sep='-')
  cat('start processing:',sample_in,'\n')
  MML=MML[MML$Sample==sample_in]
  bed_in=read.table(paste('MML_bed/',fn,sep=''))
  bed_in=bed_in[,c(1,2,3,11)]
  colnames(bed_in)=c('chr','start','end','M')
  bed_in=makeGRangesFromDataFrame(bed_in,keep.extra.columns = T)
  bed_in$M=bed_in$M/100
  olap=findOverlaps(bed_in,MML)
  olap_df=data.frame(MML_loc=subjectHits(olap),bed_in_M=bed_in$M[queryHits(olap)],MML_M=MML$MML[subjectHits(olap)])
  olap_df_agg=aggregate(olap_df[,c(2,3)],by=list(olap_df$MML_loc),mean)
  cat('finish processing:',sample_in,'in:',proc.time()[[3]]-tt1,'\n')
  return(list(sample_in,cor.test(olap_df_agg$bed_in_M,olap_df_agg$MML_M)))
  
}

coverage=readRDS('../downstream/output/coverage_mean_sd.rds')
ggplot(coverage,aes(x=sample,y=cov_mean))+geom_point()+ylab('mean_cov')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(coverage,aes(x=sample,y=cov_sd))+geom_point()+ylab('sd(coverage)')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

correlation=readRDS('../downstream/output/cor_out_pearson_DNase_dedup.rds')
correlation_df=do.call(rbind,lapply(correlation,function(x){
  data.frame(sample=x[[1]],correlation=x[[2]]$estimate,p_value=x[[2]]$p.value)
})
)
ggplot(correlation_df,aes(x=sample,y=correlation ))+geom_point()+ylab('correlation')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
hist(correlation_df$correlation,main='correlation between estimated MML and paper MML',xlab='cor')
hist(coverage$cov_mean,main='coverage of all sample',xlab='mean coverage')
hist(coverage$cov_sd,main='coverage of all sample',xlab='sd coverage')

#read in deduplication report in bam/dedup_report
dup_out=data.frame()
for (fn in dir(pattern='txt')){
  report_in=read.table(fn,skip=7)
  sample=gsub('merged','',strsplit(fn,'\\.')[[1]][1])
  dup_out=rbind(dup_out,data.frame(sample=sample,dup=report_in[,10]*100))
  
}
dup_out=readRDS('../downstream/output/dup_percent.rds')
ggplot(dup_out,aes(x=sample,y=dup))+geom_point()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  xlab('Sample')+ylab("duplication percent")+geom_text(aes(label=sample),hjust=0.5, vjust=0,angle=90)+ylim(c(0,75))

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

DMV=read_supp('DMV.xls')
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
hypo_DMR=read_supp('large_hypo_met.xls')
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


#Check overlap between filtered regions

region_filter_sp=function(x){
  UC_comp=unique(unlist(lapply(strsplit(colnames(mcols(x)),'-'),
                               function(x) paste(x[1],x[2],x[3],sep='-'))))
  
  region_idx=unique(unlist(mclapply(UC_comp,region_filter_comp,matrix_in=mcols(x),mc.cores=20)))
  region_idx=sort(region_idx)
  return(x[region_idx])
}
region_filter_comp<-function(comp_in,matrix_in){
  matrix_in_comp=matrix_in[,c(paste(comp_in,'1',sep='-'),paste(comp_in,'2',sep='-'))]
  return(which(apply(matrix_in_comp,1,function(x) sum(x>=0.1)==2)))
  
  
}
overlap_tissue_ts<-function(x,UC_in){
  #for each element,count overlap between others
  olap_tissue=lapply(UC_in,function(UC_ls) {
    tissue=strsplit(colnames(mcols(UC_ls)),'-')[[1]][1]
    olap_len=length(findOverlaps(x,UC_ls,type='equal'))
    tissue_olap=data.table(tissue1=tissue,olap=olap_len,olap_percent=olap_len/length(x),
                           tissue2=strsplit(colnames(mcols(x)),'-')[[1]][1])
    
    return(tissue_olap)
  })
  return(do.call(rbind,olap_tissue))
  
}
UC_in_filter=lapply(UC_in,region_filter_sp)
olap_UC_tissue=mclapply(UC_in_filter,overlap_tissue_ts,UC_in=UC_in_filter,mc.cores=6)
olap_UC_tissue=do.call(rbind,olap_UC_tissue)
olap_UC_tissue_ds=dcast(olap_UC_tissue,formula=tissue1~tissue2,value.var="olap_percent")


# PCA_plot ----------------------------------------------------------------

UC_in_Epiblast=readRDS('../downstream/output/UC_agnostic_mouse_all_matrix_dedup_N2_all_merged_Epiblast.rds')
UC_in_Epiblast= UC_in_Epiblast[which(apply(mcols(UC_in_Epiblast),1,function(x) !anyNA(x)))]




UC_in_Epiblast_large_JSD=UC_in_Epiblast[which(apply(mcols(UC_in_Epiblast),1,function(x) any(x>=0.1)))]
quantile(unlist(mcols(UC_in_Epiblast)),prob=0.95)

mcols(UC_in_Epiblast_large_JSD)=mcols(UC_in_Epiblast_large_JSD[,which(unlist(lapply(colnames(mcols(UC_in_Epiblast_large_JSD)),
                                                                                    function(x) strsplit(x,'-')[[1]][1]!='liver')))])
UC_in_Epiblast_PCA=PCA_df_prep(UC_in_Epiblast_large_JSD)
pdf('../downstream/output/UC_BL6DBA_PCA_Epiblast.pdf')
autoplot(UC_in_Epiblast_PCA[[1]],data=UC_in_Epiblast_PCA[[2]],colour='sample')+geom_text(label=UC_in_Epiblast_PCA[[2]]$sample)+
  theme(legend.position = 'none')
dev.off()


# Heatmap  ---------------------------------------------------
#Dealing with the UC matrix
#dNME between epiblast and NME
differential_calc<-function(ref_in,UC_gr,NME_matrix){
  ref_score=import.bedGraph(ref_in)
  start(ref_score)=start(ref_score)-1
  olap=findOverlaps(ref_score,UC_gr,type='equal')
  ref_score=ref_score[queryHits(olap)]
  olap=findOverlaps(NME_matrix,UC_gr,type='equal')
  NME_matrix=NME_matrix[queryHits(olap)]
  diff_matrix=apply(mcols(NME_matrix),2,function(x) abs(x-ref_score$score))
  diff_gr=granges(NME_matrix)
  mcols(diff_gr)=diff_matrix
  return(diff_gr)
  
}
NME_matrix=readRDS('../downstream/output/NME_matrix_mouse_all_dedup_N2.rds')
NME_matrix=NME_matrix[which(apply(mcols(NME_matrix),1,function(x) !anyNA(x)))]
#dNME/dMML calculation
dNME_gr=differential_calc('BL6DBA_Epiblast_all_allele_agnostic_nme.bedGraph',UC_in_Epiblast,NME_matrix)
saveRDS(dNME_gr,'../downstream/output/dNME_Epiblast_N2.rds')
MML_matrix=readRDS('../downstream/output/MML_matrix_mouse_all_dedup_N2.rds')
MML_matrix=MML_matrix[which(apply(mcols(MML_matrix),1,function(x) !anyNA(x)))]
dMML_gr=differential_calc('BL6DBA_Epiblast_all_allele_agnostic_mml.bedGraph',UC_in_Epiblast,MML_matrix)
dNME_gr=differential_calc('BL6DBA_Epiblast_all_allele_agnostic_mml.bedGraph',UC_in_Epiblast,MML_matrix)
saveRDS(dMML_gr,'../downstream/output/dMML_Epiblast_N2.rds')

heatmap_tissue_dNME_dMML<-function(tissue,UC_in,dNME_in,dMML_in,kmeans_k=10,
                                   region_select=NULL,sd_cutoff=0.05,JSD_cutoff=0.1,plot_scale=F,
                                   dist_scale=T,pic_header_diff='../downstream/output/dNME_dMML_UC_scale_',
                                   pic_header_all='../downstream/output/all_UC_scale_'){
  #filter regions with high row sd
  cat('Processing:',tissue,'\n')
  tt1=proc.time()[[3]]
  if(!is.null(region_select)){
    UC_in=subsetByOverlaps(UC_in,region_select,type='equal')
  }
  #Find all tissue in matrix
  tissue_all=unlist(lapply(strsplit(colnames(mcols(UC_in)),'-'),function(x) x[1]))
  
  #getting tissue for UC, dNME, dMML
  UC_in_tissue=granges(UC_in)
  mcols(UC_in_tissue)=mcols(UC_in)[,unlist(lapply(strsplit(colnames(mcols(UC_in)),'-'),function(x) x[1]))==tissue]
  mcols(dNME_in)=mcols(dNME_in)[,unlist(lapply(strsplit(colnames(mcols(dNME_in)),'-'),function(x) x[1]))==tissue]
  mcols(dMML_in)=mcols(dMML_in)[,unlist(lapply(strsplit(colnames(mcols(dMML_in)),'-'),function(x) x[1]))==tissue]
  #Filter regions for UC_in or UC_in_tissue
  UC_in=UC_in[which(rowSds(as.matrix(mcols(UC_in_tissue)))>=sd_cutoff&apply(as.matrix(mcols(UC_in_tissue)),1,function(x) any(x>=JSD_cutoff))),]
  UC_in_tissue=subsetByOverlaps(UC_in_tissue,UC_in,type='equal')
  sample_number=ncol(mcols(UC_in_tissue))#sample number for subsetting calculation matrix
  if(plot_scale){
    mcols(UC_in_tissue)= t(scale(t(as.matrix(mcols(UC_in_tissue)))))
    mcols(dNME_in)=t(scale(t(as.matrix(mcols(dNME_in)))))
    mcols(dMML_in)=t(scale(t(as.matrix(mcols(dMML_in)))))
    mcols(UC_in)=fastDoCall("cbind",lapply(unique(tissue_all),function(x) {
      out=t(scale(t(as.matrix(mcols(UC_in))[,which(unlist(lapply(strsplit(colnames(as.matrix(mcols(UC_in))),'-'),function(y) y[1]))==x)])))
      attributes(out)[c('scaled:center','scaled:scale')] <- NULL
      out=apply(out,2,function(y) {y[which(is.na(y))]=0
      return(y)})
      return(out)}))
    cat('finishing scaling in:',proc.time()[[3]]-tt1,'\n')
  }
  #subset overlapped regions
  dNME_olap=findOverlaps(dNME_in,UC_in_tissue,type='equal')
  UC_in_tissue=UC_in_tissue[subjectHits(dNME_olap)]
  colnames(mcols(dNME_in))=paste(colnames(mcols(dNME_in)),'-dNME',sep='')
  mcols(UC_in_tissue)=cbind(mcols(UC_in_tissue),mcols(dNME_in)[queryHits(dNME_olap),])
  dMML_olap=findOverlaps( dMML_in,UC_in_tissue,type='equal')
  colnames(mcols(dMML_in))=paste(colnames(mcols(dMML_in)),'-dMML',sep='')
  UC_in_tissue=UC_in_tissue[subjectHits(dMML_olap)]
  mcols(UC_in_tissue)=cbind(mcols(UC_in_tissue),mcols(dMML_in)[queryHits(dMML_olap),])
  UC_in=subsetByOverlaps(UC_in,UC_in_tissue,type='equal')
  #Fine calculation matrix
  UC_in_mt_all=as.matrix(mcols(UC_in))
  UC_in_mt_diff=as.matrix(mcols(UC_in_tissue))
  #subset dNME and dMML regions
  
  #matrix for calculating kmeans
  UC_in_calc=UC_in_mt_diff[,1:sample_number]
  if(dist_scale&!plot_scale){UC_in_calc=t(scale(t(UC_in_calc)))}
  kmeans_out=kmeans(UC_in_calc,center=kmeans_k,iter.max = 10000,nstart=20,algorithm ="MacQueen")
  
  #Define row gap and column gap for both heatmap
  row_gap_dNME_dMML=as.vector(cumsum(table(sort(kmeans_out$cluster,decreasing=F))))
  row_gap_all_tissue=as.vector(cumsum(table(sort(kmeans_out$cluster,decreasing=F))))
  col_break_dNME_dMML=c(sample_number,sample_number*2)
  col_break_tissue=cumsum(table(tissue_all)[unique(tissue_all)])
  
  #Plotting with dNME and dMML
  cat('Plotting start in:',proc.time()[[3]]-tt1,'\n')
  png(paste(pic_header_diff,tissue,'.png',sep=''),width=700,height=700)
  heatmap_out_dNME_dMML=pheatmap(UC_in_mt_diff[order(kmeans_out$cluster,decreasing=F),],cluster_cols = FALSE,scale='none',col=col_fun(1024),
                                 cluster_rows = FALSE,gaps_col=col_break_dNME_dMML,main=tissue,rasterized =T,
                                 gaps_row = row_gap_dNME_dMML)
  print(heatmap_out_dNME_dMML)
  dev.off()
  #Plotting with all tissues
  png(paste(pic_header_all,tissue,'.png',sep=''),width=700,height=700)
  heatmap_out_all_tissue=pheatmap(UC_in_mt_all[order(kmeans_out$cluster,decreasing=F),],cluster_cols = FALSE,scale='none',col=col_fun(1024),
                                  cluster_rows = FALSE,gaps_col=col_break_tissue,main=tissue,rasterized =T,
                                  gaps_row = row_gap_all_tissue)
  print(heatmap_out_all_tissue)
  dev.off()
  
  
  cat('Finish plotting start in:',proc.time()[[3]]-tt1,'\n')
  return(list(granges(UC_in),kmeans_out,heatmap_out_dNME_dMML,heatmap_out_all_tissue,UC_in_tissue))
}

dNME_gr=readRDS('../downstream/output/dNME_Epiblast_N2.rds')
mcols(dNME_gr)=mcols(dNME_gr)[,which(unlist(lapply(strsplit(colnames(mcols(dNME_gr)),'-'),function(x) x[2]!='day0')))]
dMML_gr=readRDS('../downstream/output/dMML_Epiblast_N2.rds')
mcols(dMML_gr)=mcols(dMML_gr)[,which(unlist(lapply(strsplit(colnames(mcols(dMML_gr)),'-'),function(x) x[2]!='day0')))]
DNase=readRDS('../downstream/output/mm10_DNase.rds')
UC_in_Epiblast=readRDS('../downstream/output/UC_agnostic_mouse_all_matrix_dedup_N2_all_merged_Epiblast.rds')
mcols(UC_in_Epiblast)=mcols(UC_in_Epiblast)[,which(unlist(lapply(strsplit(colnames(mcols(UC_in_Epiblast)),'-'),function(x) x[2]!='day0')))]
tissue_all=list('EFP','liver','Lung','midbrain','hindbrain','forebrain','kidney','heart','stomach','limb','NT','intestine')
dNME_dMML_all=mclapply(tissue_all,heatmap_tissue_dNME_dMML,UC_in=UC_in_Epiblast,dNME_in=dNME_gr,dMML_in=dMML_gr,
                       region_select=DNase[DNase$region_type=='DNase'],sd_cutoff=0.05,JSD_cutoff = 0.15,plot_scale=T,dist_scale=T,mc.cores=4)

dNME_dMML_all=mclapply(tissue_all,heatmap_tissue_dNME_dMML,UC_in=UC_in_Epiblast,dNME_in=dNME_gr,dMML_in=dMML_gr,
                       region_select=DNase[DNase$region_type=='DNase'],sd_cutoff=0.05,JSD_cutoff = 0.15,plot_scale=F,dist_scale=T,
                       pic_header_diff='../downstream/output/dNME_dMML_UC_no_scale_',
                       pic_header_all='../downstream/output/all_UC_no_scale_',mc.cores=4)

names(dNME_dMML_all)=tissue_all
saveRDS(dNME_dMML_all,'../downstream/output/mm10/dNME_dMML_all_scaled_no_day0.rds')
for(ts in tissue_all){saveRDS(dNME_dMML_all[[ts]],paste('../downstream/output/dNME_dMML_all_scaled_no_day0',ts,'.rds'))}
#Forebrain example
dNME_dMML_in=readRDS('../downstream/output/dNME_dMML_all_scaled_no_day0.rds')
for(tissue in names(dNME_dMML_in)){
  png(paste('../downstream/output/mm10/mm10_',tissue,'_scaled_tissue.png'),width=700,height=700)
  print(dNME_dMML_in[[tissue]][[4]])
  dev.off()
  png(paste('../downstream/output/mm10/mm10_',tissue,'_scaled_dNME_dMML.png'),width=700,height=700)
  print(dNME_dMML_in[[tissue]][[3]])
  dev.off()
  
}

kmeans_forebrain=dNME_dMML_in$forebrain
forbrain_all_regions=kmeans_forebrain[[1]]
forbrain_all_clusters=kmeans_forebrain[[2]]$cluster



gene_cluster(DNase,forbrain_all_regions,forbrain_all_clusters,4)
gene_cluster(DNase,forbrain_all_regions,forbrain_all_clusters,7)
gene_cluster(DNase,forbrain_all_regions,forbrain_all_clusters,8)
gene_cluster(DNase,forbrain_all_regions,forbrain_all_clusters,9)
gene_cluster(DNase,forbrain_all_regions,forbrain_all_clusters,10)
gene_cluster<-function(annote_in,region_in,cluster_in,cluster,output_head='../downstream/output/forebrain_C'){
  region_in=region_in[cluster_in==cluster]
  genes_anote=subsetByOverlaps(annote_in, region_in,type='equal')
  print(genes_anote)
  gene_sig=unique(genes_anote$gene[abs(genes_anote$dist)<=3000])
  write(gene_sig,paste(output_head,cluster,'_JSD.txt',sep=''))
}
gene_all=unique(DNase$gene[abs(DNase$dist)<=3000&DNase$region_type=='DNase'])
write(gene_all,'../downstream/output/forebrain_JSD_all.txt')

#Sample specific analysis
#Filtering the regions with all samples available
UC_in_Epiblast=readRDS('../downstream/output/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_Epiblast.rds')
UC_in_Epiblast=lapply(UC_in_Epiblast,function(x) x[which(apply(mcols(x),1,function(y) !anyNA(y)))])
#UC_90_quantile=lapply(UC_in_Epiblast,function(x) quantile(unlist(mcols(x)),prob=0.9))
UC_in_Epiblast_large_JSD=lapply(UC_in_Epiblast,function(x) x[which(apply(mcols(x),1,function(y) any(y>=0.25)))])
lapply(UC_in_Epiblast_large_JSD,length)
#For a single tissue
UC_tissue=UC_in_Epiblast_large_JSD[[1]]
UC_tissue_mt=as.matrix(mcols(UC_tissue))
UC_tissue_mt=UC_tissue_mt[order(rowSds(UC_tissue_mt),decreasing=T),][1:30000,] 
#Plot heatmap to show clusters
pheatmap(UC_tissue_mt,cluster_cols = FALSE,scale='row',cutree_rows=10,col=col_fun(1024),treeheight_row = FALSE)

#motif enrichment analysis

motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_OR=do.call(rbind,lapply(c('forebrain','midbrain','hindbrain','liver','heart','limb','EFP','NT'),function(y) {
  do.call(rbind,lapply(1:10,function(x) {
    motif_tissue_enrich(y,x,motif_dir[motif_dir$prob<0.5,])
  }))
  
  
}))
motif_OR$cluster=as.factor(motif_OR$cluster)
motif_OR$sig=''
motif_OR$qval=p.adjust(motif_OR$pval,method='BH')
motif_OR$sig[motif_OR$qval<=0.05]='*'
ggplot(motif_OR,aes(x=tissue,y=OR,fill=cluster,group=cluster,label=sig))+geom_bar(stat="identity",position=position_dodge())+
  geom_text(position = position_dodge(width = 0.9))

write.csv(motif_OR,'../downstream/output/motif_prefer_less_ent_enrich_mm10.csv')

sig_TF=unique(unlist(strsplit(motif_OR$overlap_TF[motif_OR$qval<=0.05&motif_OR$OR>1],',')))

motif_tissue_enrich<-function(tissue,cluster,motif_dir){
  motif_tissue=read.csv(paste('../downstream/input/mm10_motif_development/',tissue,'/motif_',tissue,'_cluster_',cluster,'.csv',sep=''))
  motif_tissue$motif_short=gsub('.*\\_','',motif_tissue$motif)
  OR=motif_enriched_mm10(motif_dir,motif_tissue,cutoff=0.1)
  return(cbind(data.frame(tissue=tissue,cluster=cluster),OR))
}
motif_enriched_mm10<-function(motif_dNME,motif_tissue,cutoff=0.1,pseduo=10){
  dNME_motif=motif_dNME[motif_dNME$qval_binom<=cutoff,]
  non_dNME_motif=motif_dNME[motif_dNME$qval_binom>cutoff,]
  tissue_motif=motif_tissue[motif_tissue$FDR<=cutoff,]
  nontissue_motif=motif_tissue[motif_tissue$FDR>cutoff,]
  OR=c(sum(dNME_motif$TF %in% tissue_motif$motif_short)+pseduo,
       sum(dNME_motif$TF %in% nontissue_motif$motif_short)+pseduo,
       sum(non_dNME_motif$TF %in% tissue_motif$motif_short)+pseduo,
       sum(non_dNME_motif$TF %in% nontissue_motif$motif_short)+pseduo)
  OR_matrix=matrix(OR,nrow=2)
  OR=fisher.test(OR_matrix)
  out=data.frame(overlap=OR_matrix[1,1], OR=OR$estimate,pval=OR$p.value,low_CI=OR$conf.int[1],up_CI=OR$conf.int[2])
  out$overlap_TF=paste(dNME_motif$TF[dNME_motif$TF %in% tissue_motif$motif_short],collapse = ',')
  return(out)
}

#Looking for regions driven by dNME or dMML
UC=mcols(UC_in_tissue)[,unlist(lapply(strsplit(colnames(mcols(UC_in_tissue)),'-'),function(x) x[length(x)]=='all'))]
colnames(UC)=unlist(lapply(strsplit(colnames(UC),'-'),function(x) x[[2]]))
dNME=mcols(UC_in_tissue)[,unlist(lapply(strsplit(colnames(mcols(UC_in_tissue)),'-'),function(x) x[length(x)]=='dNME'))]
colnames(dNME)=unlist(lapply(strsplit(colnames(dNME),'-'),function(x) x[[2]]))
dNME=dNME[,match(colnames(UC),colnames(dNME))]

#read in NME, MML file etc
source('mainFunctions_sub.R')
subjects="BL6DBA"
tissue="Epiblast_merged_paired"
in_dir="../downstream/output/"
setGenomeLengths <- function(GR,chrsOfInterest=paste0("chr",c(1:19))){
  # Get genome info
  GR=chr_check(GR)
  mm10<-getBSgenome("mm10")
  genome.seqinfo <- seqinfo(mm10)
  genome.seqinfo <- genome.seqinfo[chrsOfInterest]
  GR <- GR[seqnames(GR) %in% chrsOfInterest]
  genome(GR) <- genome(genome.seqinfo)
  seqlevels(GR) <- seqlevels(genome.seqinfo)
  seqlengths(GR) <- seqlengths(genome.seqinfo)
  
  return(GR)
}
#Only need ref and alt NME which in this case is the NME1 and NME2
mm10_Epiblast_allele=read.alleleGR(subjects,tissue,in_dir,chrsOfInterest=paste0("chr",1:19))
mm10_Epiblast_diff=read.diffGR(subjects,tissue,in_dir,chrsOfInterest=paste0("chr",1:19))
#Merge
mm10_Epiblast_allele$Genome_char=c("ref","alt")[as.numeric(mm10_Epiblast_allele$Genome)]
mm10_Epiblast_allele$stat_gene=paste0(mm10_Epiblast_allele$Genome_char,mm10_Epiblast_allele$Statistic)
mm10_Epiblast_diff$stat_gene=paste0(mm10_Epiblast_diff$Statistic,"_pval")
mm10_Epiblast_diff$abs_diff=mm10_Epiblast_diff$Value
mm10_Epiblast_diff$Value=mm10_Epiblast_diff$pvalue
mm10_Epiblast_all=c(mm10_Epiblast_allele,mm10_Epiblast_diff)
mm10_Epiblast_gr=granges(unique(mm10_Epiblast_all))
olap=findOverlaps(mm10_Epiblast_all,mm10_Epiblast_gr)
mm10_Epiblast_dt=as.data.table(mcols(mm10_Epiblast_all))
mm10_Epiblast_dt$idx=subjectHits(olap)
mm10_Epiblast_dc=dcast.data.table(mm10_Epiblast_dt,idx~stat_gene,value.var="Value")
mm10_Epiblast_gr=mm10_Epiblast_gr[mm10_Epiblast_dc$idx]
mcols(mm10_Epiblast_gr)=mm10_Epiblast_dc
mm10_Epiblast_gr=mm10_Epiblast_gr[!is.na(mm10_Epiblast_gr$dNME_pval)]
variant_in=readRDS("../downstream/output/SNP_DBA.rds")
olap=findOverlaps(variant_in,mm10_Epiblast_gr)
variant_in=variant_in[queryHits(olap)]
mcols(variant_in)=mcols(mm10_Epiblast_gr[subjectHits(olap)])
motif_gene=readRDS('../downstream/motif/JASPAR_mm10_default_out_savedpwd.rds')
direction_enriched_sample<-function(tf,variant_gene,motif_gene_subj,pval_cutoff,nperm=0,stablen=50){
  motif_gene_subj=motif_gene_subj[motif_gene_subj$geneSymbol==tf]
  variant_gene=variant_gene[variant_gene$dNME_pval<=pval_cutoff]
  olap=findOverlaps(variant_gene,motif_gene_subj)
  variant_gene=variant_gene[queryHits(olap)]
  variant_gene$alleleDiff=motif_gene_subj$alleleDiff[subjectHits(olap)]
  variant_gene=variant_gene[!is.na(variant_gene$alleleDiff)]
  
  #alleleDiff is calculated use ref - alt, prefer low ent ones
  variant_gene$NME_diff=variant_gene$altNME-variant_gene$refNME
  variant_gene$MML_diff=variant_gene$altMML-variant_gene$refMML
  same_dir=sum(sign(variant_gene$alleleDiff)== sign(variant_gene$NME_diff),na.rm = TRUE)+stablen
  opposite_dir=sum(sign(variant_gene$alleleDiff)!= sign(variant_gene$NME_diff),na.rm = TRUE)+stablen
  #same_dir=sum(sign(variant_gene$alleleDiff)== sign(variant_gene$MML_diff),na.rm = TRUE)
  # opposite_dir=sum(sign(variant_gene$alleleDiff)!= sign(variant_gene$MML_diff),na.rm = TRUE)
  total_data=same_dir+opposite_dir
  
  variant_gene_df=data.frame(alleleDiff=sign(variant_gene$alleleDiff),NME_diff=sign(variant_gene$NME_diff))
  len_x=nrow(variant_gene_df)
  if(nperm>0){
    same_dir_perm=replicate(nperm,
                            sum(sample(variant_gene_df$alleleDiff,len_x,replace = F)==sample(variant_gene_df$NME_diff,len_x,replace = F)))
    same_dir_perm_prob=same_dir_perm/total_data
  }else{same_dir_perm_prob=-1}
  if(same_dir >0 &opposite_dir>0){
    
    binom=binom.test(same_dir,(same_dir+opposite_dir),0.5)
    #print(binom)
    # binom=summary(lm(abs(variant_gene$alleleDiff)~abs(variant_gene$dMML)))
    #binom$p.value=pf(binom$fstatistic[1],df1 = binom$fstatistic[2],df2 = binom$fstatistic[3],lower.tail = F)
    # return(data.frame(TF=unique(motif_gene_subj$geneSymbol),total_data=same_dir+opposite_dir,same_dir=same_dir,opposite_dir=opposite_dir,
    #                   binom.pval=binom$p.value,prob=binom$estimate[[1]],NSNP=length(variant_gene),stringsAsFactors = F))
    prob_binom=binom$estimate[[1]]
    if(nperm>0){binom.pval=sum(abs(same_dir_perm_prob-0.5)>=abs(prob_binom-0.5))/nperm}else{binom.pval=NA}
    # if(prob_binom>0.5){
    #   binom.pval=sum(same_dir_perm_prob>=(same_dir/total_data)|same_dir_perm_prob<=(1-(same_dir/total_data)))/nperm
    #   
    # }else if(prob_binom<0.5){binom.pval=sum(same_dir_perm_prob<=(same_dir/total_data)|same_dir_perm_prob>=(1-(same_dir/total_data)))/nperm}
    #cat(tf,':',binom.pval,'\n')
    return(data.frame(TF=unique(motif_gene_subj$geneSymbol),total_data=total_data,same_dir=same_dir,opposite_dir=opposite_dir,
                      binom.pval_perm=binom.pval,binom.pval=binom$p.value,prob=prob_binom,NSNP=length(variant_gene),stringsAsFactors = F))
  }
}
motif_dir=direction_calc_enriched_subj(motif_gene,variant_in,unique(motif_gene$geneSymbol),pval_cutoff=0.1)
motif_dir_hg19=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_dir_mm10=readRDS('../downstream/output/motif_dir_mm10_stable50.rds')
motif_dir_df=motif_dir_hg19[,c("TF","prob","qval_binom")]
motif_dir_df$prob_hg19=motif_dir_df$prob
motif_dir_df$prob_mm10=motif_dir_mm10$prob[match(motif_dir_hg19$TF,motif_dir_mm10$TF)]
sum(motif_dir_df$qval_binom<=0.1&motif_dir_df$prob_hg19>0.5&motif_dir_df$prob_mm10>0.5)
sum(motif_dir_df$qval_binom<=0.1&motif_dir_df$prob_hg19>0.5&motif_dir_df$prob_mm10<0.5)
plot(motif_dir_df$prob_hg19,motif_dir_df$prob_mm10,xlab='Human',ylab='Mouse')
abline(h=0.5)
plot(motif_dir_df$prob_hg19[motif_dir_df$qval_binom<=0.1],motif_dir_df$prob_mm10[motif_dir_df$qval_binom<=0.1],xlab='Human',ylab='Mouse')
abline(h=0.5)
