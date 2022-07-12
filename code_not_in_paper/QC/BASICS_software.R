#Actually run each file
# NME_agnostic=readRDS(NME_agnostic_file)
# human_hyperVar_file=gsub('.*/','',unlist(strsplit(unique(NME_agnostic$hyper_var_fn),';')))
# saveRDS(human_hyperVar_file,"../downstream/output/human_analysis/")
source('mainFunctions_sub.R')
library(BASiCS)
library(data.table)
BASiC_run_human<-function(fn,HCL_dir,output_dir,N=20000,Burn=1000,changeName=T){
  sp=unlist(strsplit(fn,";"))
  BatchInfo=NULL
  if(changeName){
    fn_out=unique(gsub("_.*","",sp))
  }else{
    fn_out=gsub(".rds","",sp)
  }
  expr=list()
  if(length(sp)>1){
    i=1
    gene=NULL
    for(sp_in in sp){
      expr[[sp_in]]=readRDS(paste0(HCL_dir,sp_in))
      BatchInfo=c(BatchInfo,rep(i,ncol(expr[[sp_in]])))
      MAV=readRDS(paste0("../downstream/input/human_analysis/NME_expression_var/scRNA/",sp_in))
      gene=unique(c(gene,rownames(MAV)))
    }
    gene_expr=Reduce(intersect,(lapply(expr,rownames)))
    gene=intersect(gene,gene_expr)
    expr=do.call(cbind,lapply(expr,function(x) x[gene,]))
    
  }else{
    expr=readRDS(paste0(HCL_dir,sp))

  }
  
  BASiC_run(expr,output_dir,N=N,Burn=Burn,fn=fn_out,BatchInfo=BatchInfo)
  return(data.table(fileIn=fn,fileOut=fn_out))
}
#Spike in concentration for mouse: 2ul
BASiC_run<-function(expr,fn,output_dir,N=20000,Burn=10000,use_spikeIn=T,BatchInfo=NULL){
  if(!is.null(BatchInfo)){
    sce <- SingleCellExperiment(list(counts=expr),
      colData=DataFrame(BatchInfo=BatchInfo),
      metadata=list(study=fn))
  }else
  {
    sce <- SingleCellExperiment(list(counts=expr),
    metadata=list(study=fn))

  }
  Chain <- BASiCS_MCMC(
    Data = sce,
    N = N, Thin = 20, Burn = Burn,
    PrintProgress = TRUE, Regression =TRUE,WithSpikes = FALSE
  )
  
  chain_summary=Summary(Chain)
  res=data.table(
    mean=displaySummaryBASiCS(chain_summary, Param = "mu")[,"median"],
    var=displaySummaryBASiCS(chain_summary, Param = "delta")[,"median"],
    hypervar_var=exp(displaySummaryBASiCS(chain_summary, Param = "epsilon")[,"median"]),
    hypervar_logvar=displaySummaryBASiCS(chain_summary, Param = "epsilon")[,"median"]
  )
  genes=names(displaySummaryBASiCS(chain_summary, Param = "epsilon")[,"median"])
  res=as.matrix(res)
  rownames(res)=genes
  saveRDS(res,paste0(output_dir,fn,".rds"))
}
# NME_in=readRDS(NME_agnostic_file)
# NME_in_dt=data.table(Sample=NME_in$Sample,scRNA_file=NME_in$hyper_var_fn)
# NME_in_dt=unique(NME_in_dt)
# saveRDS(NME_in_dt[!is.na(scRNA_file)],"../downstream/output/human_analysis/NME_MAV/human_hyperVar_file.rds")
human_hyperVar_file=readRDS("../downstream/output/human_analysis/NME_MAV/human_hyperVar_file.rds")
human_hyperVar_file$scRNA_file=gsub("../downstream/input/human_analysis/NME_expression_var/scRNA/","",human_hyperVar_file$scRNA_file)
HCL_dir='../downstream/data/HCL_scRNA/count/'
output_dir="../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/"

#file_out=mclapply(unique(human_hyperVar_file$scRNA_file),BASiC_run_human,HCL_dir=HCL_dir,output_dir=output_dir,mc.cores=10)
file_out=mclapply(unique(human_hyperVar_file$scRNA_file),BASiC_run_human,HCL_dir=HCL_dir,N=100,Burn=40,output_dir=paste0(output_dir,"small_N/"),mc.cores=10)
saveRDS(file_out,"../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/file_out.rds")
#Human separate calculation small N
mclapply(unlist(strsplit(unique(human_hyperVar_file$scRNA_file),";")),BASiC_run_human,HCL_dir=HCL_dir,N=100,Burn=40,
            output_dir=paste0("../downstream/input/human_analysis/HCL_scRNA_BASiC/small_N/"),changeName=F,mc.cores=10)

mclapply(c("AdultEsophagus_1.rds","AdultLung_1.rds"),BASiC_run_human,HCL_dir=HCL_dir,
            output_dir=paste0("../downstream/input/human_analysis/HCL_scRNA_BASiC/"),changeName=F,mc.cores=10)
#Processing mouse data
mouse_scRNA=readRDS('../downstream/input/mouse_analysis/mouse_10x/10x.rds')
time_points=unique(gsub(":.*","",colnames(mouse_scRNA)))
output_dir="../downstream/input/mouse_analysis/scRNA_BASiC/"
mclapply(time_points,function(x) BASiC_run(mouse_scRNA[,grepl(x,colnames(mouse_scRNA))],output_dir=output_dir,fn=x,use_spikeIn=F),mc.cores=length(time_points))
mclapply(time_points,function(x) BASiC_run(mouse_scRNA[,grepl(x,colnames(mouse_scRNA))],output_dir=paste0(output_dir,"small_N/"),N=100,Burn=40,fn=x,use_spikeIn=F),mc.cores=length(time_points))
#,N=100,Burn=40

# add hypervar to allele-agnostic data ------------------------------------
