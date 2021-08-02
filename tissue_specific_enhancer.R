source('mainFunctions_sub.R')
library(biomaRt)
RNA_mouse_dir='../downstream/data/mouse_RNA_tsv/'
RNA_mouse_out_rds='../downstream/output/mouse_analysis/tissue_specific_enhancer/RNA_out.rds'
RNA_mouse_out_sub_rds='../downstream/output/mouse_analysis/tissue_specific_enhancer/RNA_out_sub.rds'
H3K27AC_mouse_out_rds='../downstream/output/mouse_analysis/tissue_specific_enhancer/H3K27AC_RPKM.rds'
analyzed_tissue=c("forebrain","heart","hindbrain","liver","midbrain","EFP","limb")

#prepare ensembl conversion
# ensembl <- useEnsembl(biomart = "ensembl",version=97)
# searchDatasets(mart = ensembl, pattern = "mmusculus")#mmusculus_gene_ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",version=97)#check version before conversion
searchFilters(mart = ensembl, pattern = ".*name")#external_gene_name
#att=listAttributes(ensembl)
attributes_BM=c('ensembl_gene_id','ensembl_gene_id_version','external_gene_name')

#Reading in RNA data

name_conversion=data.table(RNA_tissue=c('embryonic facial prominence','neural tube'),paper_tissue=c("EFP",'NT'))
RNA_out=data.table()
for(fn in dir(RNA_mouse_dir,pattern='.tsv')){
    cat("Processing:",fn,'\n')
    tt1=proc.time()[[3]]
  RNA_in=fread(paste0(RNA_mouse_dir,fn))
  RNA_in$gene_id_no_version=gsub('\\..*','',RNA_in$gene_id)

  #Getting stage and tissue
  fn_sp=unlist(strsplit(gsub('\\.tsv','',fn),'_'))
  #remove E0 to P0
  if(fn_sp[1]=="E0"){fn_sp[1]="P0"}
  if(fn_sp[2] %in% name_conversion$RNA_tissue){fn_sp[2]=name_conversion$paper_tissue[which(name_conversion$RNA_tissue==fn_sp[2])]}
  RNA_in$tissue=fn_sp[2]
  RNA_in$stage=fn_sp[1]
  RNA_in$replicate=fn_sp[3]
  RNA_out=rbind(RNA_out,RNA_in)
  cat("Finishing processing in:",proc.time()[[3]]-tt1,'\n')
}
#Convert enesembl id to gene name
   #use useCache = FALSE 
#   gene_id_conv_nonid=getBM(mart=ensembl,filters='ensembl_gene_id',attributes=attributes_BM,values=unique(RNA_out$gene_id_no_version),useCache = FALSE)
#   cat('Percent name have gene name:',sum(RNA_out$gene_id_no_version%in%gene_id_conv_nonid$ensembl_gene_id)/nrow(RNA_out),'\n')#0.0.6775687 
   
  #This is more accurate
   gene_id_conv=getBM(mart=ensembl,filters='ensembl_gene_id_version',attributes=attributes_BM,values=unique(RNA_out$gene_id),useCache = FALSE)
   sum(RNA_out$gene_id%in%gene_id_conv$ensembl_gene_id_version)/nrow(RNA_out)#0.6739048
  #Filter RNA ensembl id if they have name
    RNA_out=RNA_out[RNA_out$gene_id%in%gene_id_conv$ensembl_gene_id_version]
  #Debug
#   RNA_out=RNA_out[sample(1:nrow(RNA_out),replace=F)]
  
#   head(gene_id_conv[match(RNA_out$gene_id,gene_id_conv$ensembl_gene_id_version),'ensembl_gene_id'])
#   head(gene_id_conv$ensembl_gene_id)
#    head(RNA_out$gene_id_no_version)
  RNA_out$gene_name=gene_id_conv[match(RNA_out$gene_id,gene_id_conv$ensembl_gene_id_version),'external_gene_name']
  saveRDS(RNA_out,RNA_mouse_out_rds)
  RNA_out=readRDS(RNA_mouse_out_rds)
  
  RNA_out_sub=RNA_out[(tissue%in%analyzed_tissue)&(stage!="P0"),list(tissue,stage,replicate,gene_name,FPKM)]
  #Replace 0 with smallest value
  smallest_value=min(RNA_out_sub$FPKM[RNA_out_sub$FPKM!=0])
  RNA_out_sub$FPKM[RNA_out_sub$FPKM==0]=smallest_value
  RNA_out_sub$log2FPKM=log2(RNA_out_sub$FPKM)
  Bin_enhancer=readRDS(bin_enhancer_rds)
  RNA_out_sub=RNA_out_sub[gene_name%in%Bin_enhancer$`Target Gene`]#before 507650,after491648
  saveRDS(RNA_out_sub,RNA_mouse_out_sub_rds)
#Get coverage
#in data collection folder
#sbatch coverage_calc.sh ../../downstream/input/mouse_analysis/enhancer_selection/bin_enhancer.bed ../../downstream/data/mouse_ChIP/bam_files/E16_5_midbrain_H3K27ac_2.bam ../../downstream/data/mouse_ChIP/cov_files/E16_5_midbrain_H3K27ac_2.cov
#for fn in ../../downstream/data/mouse_ChIP/bam_files/*bam; do echo sbatch coverage_calc.sh ../../downstream/input/mouse_analysis/enhancer_selection/bin_enhancer.bed $fn $fn.cov; done
#cp  ../../downstream/data/mouse_ChIP/bam_files/*.bam.cov ../../downstream/data/mouse_ChIP/cov_files/
#rm  ../../downstream/data/mouse_ChIP/bam_files/*.bam.cov
chip_cov_dir='../downstream/data/mouse_ChIP/cov_files/'
chip_bam_dir='../downstream/data/mouse_ChIP/bam_files/'
H3K27ac_output=GRanges()
for (fn in dir(chip_cov_dir,pattern='.cov')){
  cat("Processing: ", fn, '\n')
  tt1=proc.time()[[3]]
  coverage_sp=fread(paste0(chip_cov_dir,fn))
  colnames(coverage_sp)=c('chr','start','end','name','V5','V6','coverage')
  coverage_sp=makeGRangesFromDataFrame(coverage_sp[,c(1,2,3,7)],keep.extra.columns=T)
  fn_bam=gsub('.cov','',fn)
  chip_read_info=system(paste0('samtools flagstat -@24 ',chip_bam_dir,fn_bam),intern=T)
  total_reads=as.numeric(gsub(' \\+.*','',chip_read_info[1]))
  coverage_sp$peakLength=width(coverage_sp)
  coverage_sp$RPKM=(coverage_sp$coverage)/(coverage_sp$peakLength/1000*total_reads/10^6)
  fn_sp=unlist(strsplit(gsub('_5','.5',gsub('.bam.cov','',fn)),'_'))
  coverage_sp$tissue=fn_sp[2]
  coverage_sp$stage=fn_sp[1]
  coverage_sp$marker=fn_sp[3]
  coverage_sp$replicates=fn_sp[4]
  coverage_sp$total_reads=total_reads
  H3K27ac_output=c(H3K27ac_output,coverage_sp)
  cat("Finishing processing: ",fn,'in', proc.time()[[3]]-tt1,'\n')
}
#Replace 0 by "(zeros were replaced by the smallest detectable value larger than zero)"
smallest_value=min(H3K27ac_output$RPKM[H3K27ac_output$RPKM!=0])
H3K27ac_output$RPKM[H3K27ac_output$RPKM==0]=smallest_value
H3K27ac_output$log2RPKM=log2(H3K27ac_output$RPKM)
pdf(paste0(figure_path,'H3K27ac_RPKM_dist.pdf'))
hist(H3K27ac_output$RPKM)
dev.off()
pdf(paste0(figure_path,'H3K27ac_RPKM_dist_log2.pdf'))
hist(H3K27ac_output$log2RPKM)
dev.off()
saveRDS(H3K27ac_output,H3K27AC_mouse_out_rds)

#Merge RNA data and ChIP data
RNA_out_sub=readRDS(RNA_mouse_out_sub_rds)
H3K27ac_output=readRDS(H3K27AC_mouse_out_rds)
Bin_enhancer=readRDS(bin_enhancer_rds)
#Need plus 1 from the start due to bed conversion
start(H3K27ac_output)=start(H3K27ac_output)+1
H3K27ac_output_dt=convert_GR(H3K27ac_output,direction='DT')
Bin_enhancer_dt=convert_GR(Bin_enhancer,direction="DT")
#get gene
H3K27ac_output_dt$target_gene=Bin_enhancer_dt[match(H3K27ac_output_dt$region,region)]$`Target.Gene`
H3K27ac_output_dt$sample=paste0(H3K27ac_output_dt$tissue,'-',H3K27ac_output_dt$stage,'-',H3K27ac_output_dt$replicate)
H3K27ac_output_dt$sample_gene=paste0(H3K27ac_output_dt$sample,'-',H3K27ac_output_dt$target_gene)
RNA_out_sub$sample_gene=paste0(RNA_out_sub$tissue,'-',RNA_out_sub$stage,'-',RNA_out_sub$replicate,'-',RNA_out_sub$gene_name)
H3K27ac_output_dt$FPKM=RNA_out_sub[match(H3K27ac_output_dt$sample_gene,sample_gene)]$FPKM
H3K27ac_output_dt$log2FPKM=RNA_out_sub[match(H3K27ac_output_dt$sample_gene,sample_gene)]$log2FPKM
#Check zero standard deviation region: they are zero expression regions in that tissue
#92 sample in total


H3K27ac_output_dt_fn='../downstream/output/mouse_analysis/tissue_specific_enhancer/H3K27ac_output_dt.rds'
saveRDS(H3K27ac_output_dt,H3K27ac_output_dt_fn)
pdf(paste0(figure_path,'H3K27ac_FPKM_dist.pdf'))
hist(H3K27ac_output_dt$FPKM)
dev.off()
pdf(paste0(figure_path,'H3K27ac_FPKM_dist_log2.pdf'))
hist(H3K27ac_output_dt$log2FPKM)
dev.off()
#Find overlap between histone marker signal vs tissue-specific regions

ts_aid_dt_fn='../downstream/output/mouse_analysis/tissue_specific_enhancer/ts_aid.rds'
ts_aid=readRDS(ts_aid_dt_fn)
H3K27ac_output_dt=readRDS(H3K27ac_output_dt_fn)
#Convert ts_aid to data.table
ts_aid_dt=do.call(rbind,lapply(names(ts_aid), function(x) data.table(region=ts_aid[[x]],tissue=x)))
H3K27ac_output_dt_region_uq=unique(H3K27ac_output_dt$region)
olap_enhancer=findOverlaps(convert_GR(H3K27ac_output_dt_region_uq,direction='GR'),convert_GR(ts_aid_dt$region,direction='GR'))
olap_enhancer=as.data.table(olap_enhancer)
olap_enhancer$tissue=ts_aid_dt$tissue[olap_enhancer$subjectHits]
ts_aid_enhancer=olap_enhancer[,list(tissue=paste(unique(tissue),collapse=','),n_tissue=length(unique(tissue))),by=list(queryHits)]
pdf(paste0(figure_path,'sample_overlap_enhancer.pdf'))
ts_aid_enhancer_tb=as.data.table(table(ts_aid_enhancer$n_tissue))
colnames(ts_aid_enhancer_tb)=c("ntissue_olap","nenhancer")
ts_aid_enhancer_tb$prop_enhancer=ts_aid_enhancer_tb$nenhancer/sum(ts_aid_enhancer_tb$nenhancer)
print(ggplot(ts_aid_enhancer_tb,aes(x=ntissue_olap,y=prop_enhancer))+geom_bar(stat='identity')+
   xlab('Number of tissue')+ylab('Proportion of enhancers')+ylim(c(0,1)))
dev.off()
#Used as background?
pdf(paste0(figure_path,'permute_tissue_enhancer.pdf'))
set.seed(123)
olap_enhancer_rand=olap_enhancer
olap_enhancer_rand$tissue=ts_aid_dt$tissue[sample(olap_enhancer_rand$subjectHits)]
ts_aid_enhancer_rand=olap_enhancer_rand[,list(tissue=paste(unique(tissue),collapse=','),n_tissue=length(unique(tissue))),by=list(queryHits)]
ts_aid_enhancer_tb_rand=as.data.table(table(ts_aid_enhancer_rand$n_tissue))
colnames(ts_aid_enhancer_tb_rand)=c("ntissue_olap","nenhancer")
ts_aid_enhancer_tb_rand$prop_enhancer=ts_aid_enhancer_tb_rand$nenhancer/sum(ts_aid_enhancer_tb_rand$nenhancer)
print(ggplot(ts_aid_enhancer_tb_rand,aes(x=ntissue_olap,y=prop_enhancer))+geom_bar(stat='identity')+xlab('Number of tissue')+
    ylab('Proportion of enhancers')+ylim(c(0,1)))
dev.off()

#Plotting raw expression and raw Histone marker signal
#Boxplot for each tissue plot ts aid region for that tissue and others
#Using tissue only enhancer
tissue_only_enhancer=H3K27ac_output_dt_region_uq[unique(ts_aid_enhancer[n_tissue==1]$queryHits)]
tissue_only_enhancer_fn='../downstream/output/mouse_analysis/tissue_specific_enhancer/tissue_only_enhancer.rds'
saveRDS(tissue_only_enhancer,tissue_only_enhancer_fn)
H3K27ac_output_dt=H3K27ac_output_dt[region %in% tissue_only_enhancer]#1944972 ->1274752 ->795800, 65%
H3K27ac_output_dt=readRDS(H3K27ac_output_dt_fn)
pdf(paste0(figure_path,'expression_histone_tissue_enhancer_boxplot_all.pdf'))
all_tissue=names(ts_aid)

#Also prepare the heatmap
library(ggpubr)
for(ts in all_tissue){

  olap=findOverlaps(convert_GR(H3K27ac_output_dt$region,dir="GR"),convert_GR(ts_aid[[ts]]))
  my_comparisons=list()
  i=1
  for(ts_others in all_tissue[all_tissue != ts]){
    my_comparisons[[i]]=c(ts,ts_others)
    i=i+1
  }
  H3K27ac_output_dt_ts=H3K27ac_output_dt[queryHits(olap)]
  H3K27ac_output_dt_ts$tissue=factor(H3K27ac_output_dt_ts$tissue,levels=c(ts,all_tissue[all_tissue!=ts]))
  plot_out_expression=ggplot(H3K27ac_output_dt_ts,aes(x=tissue,y=log2FPKM))+geom_boxplot(outlier.shape = NA)+
          ggtitle(paste0(ts,": expression"))+ylab('log2FPKM')+
          stat_compare_means(comparisons = my_comparisons,method="wilcox.test",method.args=list(alternative="greater")) 
          # Add comparisons p-value
          #stat_compare_means(label.y = 50)
          # Add global p-value

  print(plot_out_expression)
    plot_out_histone=ggplot(H3K27ac_output_dt_ts,aes(x=tissue,y=log2RPKM))+geom_boxplot(outlier.shape = NA)+
          ggtitle(paste0(ts,": histone"))+ylab('log2RPKM')+
          stat_compare_means(comparisons = my_comparisons,method="wilcox.test",method.args=list(alternative="greater")) 
          # Add pairwise comparisons p-value
          #stat_compare_means(label.y = 50)
          # Add global p-value

   print(plot_out_histone)

}

dev.off()
#Processing correlation
H3K27ac_output_dt_cor=H3K27ac_output_dt[,list(cor=cor(log2RPKM,log2FPKM,method='spearman'),std_log2RPKM=sd(log2RPKM),
                      std_log2FPKM=sd(log2FPKM)),by=list(region,tissue)]
H3K27ac_output_dt_cor=H3K27ac_output_dt_cor[std_log2RPKM!=0&std_log2FPKM!=0]
H3K27ac_output_dt_cor_dc=dcast.data.table(H3K27ac_output_dt_cor,region~tissue,value.var='cor')
H3K27ac_output_dt_cor_dc_mt=as.matrix(H3K27ac_output_dt_cor_dc[,-1])
rownames(H3K27ac_output_dt_cor_dc_mt)=H3K27ac_output_dt_cor_dc$region
#Check if tissue-specific UC also has highest value
H3K27ac_output_dt_cor_dc_rank=t(apply(H3K27ac_output_dt_cor_dc[,-1],1,function(x) rank(-x)))
rownames(H3K27ac_output_dt_cor_dc_rank)=H3K27ac_output_dt_cor_dc$region
ts_aid=readRDS(ts_aid_dt_fn)
 UC_in=readRDS(UC_in_matrix_cluster_file)
 analyzed_region=lapply(UC_in,rownames)
##Add non_tissue_specific region as control
extract_rank_dt<-function(region_in,enhancer_rank,enhancer_correlation,ts,tissue_name){
    olap=findOverlaps(convert_GR(region_in,direction="GR"),
                  convert_GR(rownames(enhancer_rank),direction="GR"))
    region_in=region_in[queryHits(olap)]
    region_in_dt=data.table(region=region_in,
                            rank_correlation=enhancer_rank[subjectHits(olap),ts],
                            enhancer_region=rownames(enhancer_rank[subjectHits(olap),]),
                            tissue=tissue_name
                            )
    region_in_dt$correlation=enhancer_correlation[region_in_dt$enhancer_region,ts]
    return(region_in_dt)
    }
ts_aid_out=data.table()
for(ts in names(ts_aid)){
   cat("Processing:",ts,'\n')
   tt1=proc.time()[[3]]
    ts_aid_out=rbind(ts_aid_out,
                      extract_rank_dt(ts_aid[[ts]],
                                      H3K27ac_output_dt_cor_dc_rank,
                                      H3K27ac_output_dt_cor_dc_mt,
                                      ts,
                                      ts),
                      extract_rank_dt(setdiff(analyzed_region[[ts]],ts_aid[[ts]]),
                                      H3K27ac_output_dt_cor_dc_rank,
                                      H3K27ac_output_dt_cor_dc_mt,
                                      ts,
                                      paste0(ts,'-control')
                                      )
    )
    cat("Finish processing:",ts, 'in:',proc.time()[[3]]-tt1,'\n')
   
}
ts_aid_out_fn='../downstream/output/mouse_analysis/tissue_specific_enhancer/tissue_specific_correlation_uc01.rds'
saveRDS(ts_aid_out,ts_aid_out_fn)
#Plotting
#Barplot for each tissue plot ts aid region for that tissue and others
ts_aid_out=readRDS(ts_aid_out_fn)

pdf(paste0(figure_path,'correlation_tissue_rank_all.pdf'))
   plot_dt=data.table(correlation=ts_aid_out$rank_correlation,cor_type='tissue_specific_UC',tissue=ts_aid_out$tissue)
    plot_dt$correlation=factor(round(plot_dt$correlation),levels=as.character(1:7))#Round tiles
    plot_dt$tissue=factor(plot_dt$tissue,levels=unique(plot_dt$tissue))
    plot_dt=plot_dt[,.N,by=list(tissue,correlation)]
    plot_dt[,prop_N:=N/sum(N),by=list(tissue)]
      print(ggplot(plot_dt,aes(x=tissue,y=prop_N,group=correlation,fill=correlation))+geom_bar(stat="identity",position="stack")+
      ggtitle("")+xlab("")+ylim(c(0,1))+ylab("")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

#Assigning cluster to each region
H3K27ac_output_dt=readRDS(H3K27ac_output_dt_fn)

tissue_out_filtered=readRDS(tissue_out_filtered_fn)