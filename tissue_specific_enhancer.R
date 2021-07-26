source('mainFunctions_sub.R')
library(biomaRt)
RNA_mouse_dir='../downstream/data/mouse_RNA_tsv/'
RNA_mouse_out_rds='../downstream/output/mouse_analysis/tissue_specific_enhancer/RNA_out.rds'
enhancer_Bin=readRDS(bin_enhancer_rds)
#prepare ensembl conversion
# ensembl <- useEnsembl(biomart = "ensembl",version=97)
# searchDatasets(mart = ensembl, pattern = "mmusculus")#mmusculus_gene_ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",version=97)#check version before conversion
searchFilters(mart = ensembl, pattern = ".*name")#external_gene_name
#att=listAttributes(ensembl)
attributes_BM=c('ensembl_gene_id','ensembl_gene_id_version','external_gene_name')
gene_mm10_conv=getBM(mart=ensembl,filters=ensembl_gene_id_version,attributes=attributes_BM,values=)
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
  gene_id_conv_nonid=getBM(mart=ensembl,filters='ensembl_gene_id',attributes=attributes_BM,values=unique(RNA_out$gene_id_no_version),useCache = FALSE)
  cat('Percent name have gene name:',sum(RNA_out$gene_id_no_version%in%gene_id_conv_nonid$ensembl_gene_id)/nrow(RNA_out),'\n')#0.0.6775687 
   
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