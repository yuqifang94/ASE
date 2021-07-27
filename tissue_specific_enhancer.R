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
  RNA_out_sub=RNA_out[(tissue%in%analyzed_tissue)&(stage!="P0"),list(tissue,stage,replicate,gene_name,FPKM)]
  #Replace 0 with smallest value
  smallest_value=min(RNA_out_sub$FPKM[RNA_out_sub$FPKM!=0])
  RNA_out_sub$FPKM[RNA_out_sub$FPKM==0]=smallest_value
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
pdf(paste0(figure_path,'H3K27ac_RPKM_dist.pdf'))
hist(H3K27ac_output$RPKM)
dev.off()
pdf(paste0(figure_path,'H3K27ac_RPKM_dist_log2.pdf'))
hist(log2(H3K27ac_output$RPKM))
dev.off()
H3K27ac_output$log2RPKM=log2(H3K27ac_output$RPKM)
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
H3K27ac_output_dt$log2FPKM=log2(H3K27ac_output_dt$FPKM)
#Check zero standard deviation region
H3K27ac_output_dt_cor=H3K27ac_output_dt[,list(cor=cor(log2RPKM,log2FPKM,method='spearman')),by=list(region,tissue)]
H3K27ac_output_dt_cor_fn='../downstream/output/mouse_analysis/H3K27ac_output_dt_cor.rds'
saveRDS(H3K27ac_output_dt_cor,H3K27ac_output_dt_cor_fn)
H3K27ac_output_dt_cor_dc=dcast.data.table(H3K27ac_output_dt_cor,region~tissue,value.var='cor')
#Check if tissue-specific UC also has highest value
H3K27ac_output_dt_cor_dc_rank=t(apply(H3K27ac_output_dt_cor_dc[,-1],1,function(x) rank(-x)))
rownames(H3K27ac_output_dt_cor_dc_rank)=H3K27ac_output_dt_cor_dc$region
ts_aid=readRDS(ts_aid_dt_fn)
for(ts in names(ts_aid)){
    ts_aid_ts=ts_aid[[ts]]
    olap=findOverlaps(convert_GR(ts_aid_ts,direction="GR"),
    convert_GR(rownames(H3K27ac_output_dt_cor_dc_rank),direction="GR"))
    ts_aid_ts=ts_aid_ts[queryHits(olap)]
    ts_aid_ts_dt=data.table(region=ts_aid_ts,rank=H3K27ac_output_dt_cor_dc_rank[subjectHits(olap),ts])
}