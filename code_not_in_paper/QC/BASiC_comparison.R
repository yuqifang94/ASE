source('mainFunctions_sub.R')
NME_in=readRDS(NME_agnostic_file)
MML_in=readRDS(MML_agnostic_file)
genomic_features=readRDS(genomic_features_file)
#This is without combine replicates
hyper_var_dir="../downstream/input/human_analysis/HCL_scRNA_BASiC/small_N/"
scRNA_dir="../downstream/input/human_analysis/NME_expression_var/scRNA/"
GR_calc=data.frame()
scRNA_result=data.frame()
MML_hypervar_calc=GRanges()
NME_hypervar_calc=GRanges()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  hyper_var_file=gsub(scRNA_dir,hyper_var_dir,hyper_var_file)
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){
    
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    
    
  }else{cat("file not exist for:",sp,'\n')}
}
# NME_hypervar_calc_no_NA=NME_hypervar_calc[!is.na(NME_hypervar_calc$hypervar_logvar)]
# NME_hypervar_calc_no_NA_dt=convert_GR(NME_hypervar_calc_no_NA,"DT")
# NME_hypervar_calc_no_NA_dt_tss=NME_hypervar_calc_no_NA_dt[abs(dist)<=250]
#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),
             paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_small_N_single.rds'))

#This is without combine replicates
hyper_var_dir="../downstream/input/human_analysis/HCL_scRNA_BASiC/"
scRNA_dir="../downstream/input/human_analysis/NME_expression_var/scRNA/"
GR_calc=data.frame()
scRNA_result=data.frame()
MML_hypervar_calc=GRanges()
NME_hypervar_calc=GRanges()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  hyper_var_file=gsub(scRNA_dir,hyper_var_dir,hyper_var_file)
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){
    
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    
    
  }else{cat("file not exist for:",sp,'\n')}
}
# NME_hypervar_calc_no_NA=NME_hypervar_calc[!is.na(NME_hypervar_calc$hypervar_logvar)]
# NME_hypervar_calc_no_NA_dt=convert_GR(NME_hypervar_calc_no_NA,"DT")
# NME_hypervar_calc_no_NA_dt_tss=NME_hypervar_calc_no_NA_dt[abs(dist)<=250]
#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),
             paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_large_N_single.rds'))

#With combining replicates
#This is without combine replicates
hyper_var_dir="../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/small_N/"
scRNA_dir="../downstream/input/human_analysis/NME_expression_var/scRNA/"
file_out=readRDS("../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/file_out.rds")
file_out=do.call(rbind,file_out)
GR_calc=data.frame()
scRNA_result=data.frame()
MML_hypervar_calc=GRanges()
NME_hypervar_calc=GRanges()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unique(NME_in$hyper_var_fn[NME_in$Sample==sp])
  cat('Processing',sp,'\n')
  if(!is.na(hyper_var_file)){
     hyper_var_file=file_out[fileIn==gsub(scRNA_dir,"",hyper_var_file)]$fileOut
     hyper_var_file=paste0(hyper_var_dir,hyper_var_file,'.rds')
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    
    
  }else{cat("file not exist for:",sp,'\n')}
}
# NME_hypervar_calc_no_NA=NME_hypervar_calc[!is.na(NME_hypervar_calc$hypervar_logvar)]
# NME_hypervar_calc_no_NA_dt=convert_GR(NME_hypervar_calc_no_NA,"DT")
# NME_hypervar_calc_no_NA_dt_tss=NME_hypervar_calc_no_NA_dt[abs(dist)<=250]
#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),
             paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_merge_replicates_smallN.rds'))

#Large N
hyper_var_dir="../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/"
scRNA_dir="../downstream/input/human_analysis/NME_expression_var/scRNA/"
file_out=readRDS("../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/file_out.rds")
file_out=do.call(rbind,file_out)
GR_calc=data.frame()
scRNA_result=data.frame()
MML_hypervar_calc=GRanges()
NME_hypervar_calc=GRanges()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unique(NME_in$hyper_var_fn[NME_in$Sample==sp])
  cat('Processing',sp,'\n')
  if(!is.na(hyper_var_file)){
     hyper_var_file=file_out[fileIn==gsub(scRNA_dir,"",hyper_var_file)]$fileOut
     hyper_var_file=paste0(hyper_var_dir,hyper_var_file,'.rds')
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    
    
  }else{cat("file not exist for:",sp,'\n')}
}
# NME_hypervar_calc_no_NA=NME_hypervar_calc[!is.na(NME_hypervar_calc$hypervar_logvar)]
# NME_hypervar_calc_no_NA_dt=convert_GR(NME_hypervar_calc_no_NA,"DT")
# NME_hypervar_calc_no_NA_dt_tss=NME_hypervar_calc_no_NA_dt[abs(dist)<=250]
#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),
             paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_merge_replicates_largeN.rds'))

#Mouse
#Read in mouse NME and scRNA
#Need to rename the files in the folder
NME_in=readRDS(NME_matrix_file)
dir_scRNA_mouse="../downstream/input/mouse_analysis/scRNA_BASiC/small_N/"
#From JASON

mcols(NME_in)=mcols(NME_in)[,grepl('limb',colnames(mcols(NME_in)))]

gtf <- fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
genes <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
genes$gene_name <- gn
NME_in=dist_calc(NME_in,genes)
#Percent gene covered?
length(unique(NME_in[abs(NME_in$dist)<=3000]$gene))/length(genes[seqnames(genes)!="chrM"])#96%
NME_in_dt=convert_GR(NME_in,dir='DT')
NME_in_dt=melt.data.table(NME_in_dt,id.var=c('dist','gene','region'),value.name='NME',variable.name='stage')

NME_in_dt$hyper_var=-100
NME_in_dt$var=-100
NME_in_dt$mean=-100
for(st in unique(NME_in_dt$stage)){
  tt1=proc.time()[[3]]
  if(file.exists(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))){
    scRNA_in=readRDS(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))
    scRNA_in=scRNA_in[rownames(scRNA_in)%in% unique(NME_in_dt[(stage==st)]$gene),]
    NME_in_dt=NME_in_dt[gene%in%rownames(scRNA_in)]
    if(nrow(scRNA_in)>0){
      #Add hypervar to TSS 

      NME_in_dt[(stage==st)]$hyper_var=scRNA_in[NME_in_dt[(stage==st)]$gene,"hypervar_logvar"]
      NME_in_dt[(stage==st)]$var=scRNA_in[NME_in_dt[(stage==st)]$gene,"var"]
      NME_in_dt[(stage==st)]$mean=scRNA_in[NME_in_dt[(stage==st)]$gene,"mean"]
      
    }
  }else{cat("File not exist for ",st,'\n')}
  cat('Finish processing ',sub('E','',st),'in: ',proc.time()[[3]]-tt1,'\n')
  
}
 NME_in_dt[!is.na(hyper_var)&hyper_var!=-100&!is.na(NME)]
saveRDS(NME_in_dt, "../downstream/output/mouse_analysis/NME_in_limb_ENOCD3_BASiC_smallN.rds")

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
#Single samples
hyper_var_all=readRDS( paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_small_N_single.rds'))
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir=paste0(NME_MAV_human_out_dir,"BASiC_small_N_single/")) #0.1341079
hyper_var_all=readRDS( paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_large_N_single.rds'))
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir=paste0(NME_MAV_human_out_dir,"BASiC_large_N_single/")) #0.1888981


hyper_var_all=readRDS( paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_merge_replicates_smallN.rds'))
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir=paste0(NME_MAV_human_out_dir,"BASiC_small_N_batch/"))#0.176429

hyper_var_all=readRDS( paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC_merge_replicates_largeN.rds'))
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir=paste0(NME_MAV_human_out_dir,"BASiC_large_N_batch/"))#0.213196

hyper_var_all=readRDS("../downstream/output/mouse_analysis/NME_in_limb_ENOCD3_BASiC_smallN.rds")
hyper_var_all=hyper_var_all[!is.na(hyper_var)&(hyper_var!=-100)&!is.na(NME)]
hyper_var_all$Sample=hyper_var_all$stage
hyper_var_all$score=hyper_var_all$NME
dist_plot_run(hyper_var_all,theme_glob,ylab="NME",stat_in="hyper_var",dir=paste0(NME_MAV_human_out_dir,"BASiC_small_N_mouse/"))# 0.113832




#Correlation between small N, large N and MAV
library(gridExtra)
plot_dot<-function(MAV_in,BASiC_in,title,xlab="MAV",ylab="BASiC"){
    gene=intersect(rownames(MAV_in),rownames(BASiC_in))
    plotDt=data.table(MAV=MAV_in[gene,"hypervar_logvar"],BASiC=BASiC_in[gene,"hypervar_logvar"])
    plotDt=plotDt[!is.na(MAV)&!is.na(BASiC)]
    corOut=cor(plotDt$MAV,plotDt$BASiC)
    return(list(
        plotOut=ggplot(plotDt,aes(x=MAV,y=BASiC))+xlab(xlab)+ylab(ylab)+geom_bin2d(bins=100)+geom_smooth()+ggtitle(paste0(title,"\nCorrelation:",corOut)),
        corOut=corOut
      )
    )
}
BASiC_MAV_comp<-function(fn){
  fnBatch=gsub("_.*","",fn)
  batch_fn_largeN_fn=paste0("../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/",fnBatch,".rds")
  if(file.exists(batch_fn_largeN_fn)){
    MAV_fn=readRDS(paste0("../downstream/input/human_analysis/NME_expression_var/scRNA/",fn,".rds"))
    single_fn_smallN=readRDS(paste0("../downstream/input/human_analysis/HCL_scRNA_BASiC/small_N/",fn,".rds"))
    single_fn_largeN=readRDS(paste0("../downstream/input/human_analysis/HCL_scRNA_BASiC/",fn,".rds"))
    
    batch_fn_smallN=readRDS(paste0("../downstream/input/human_analysis/HCL_scRNA_BASiC_batchInfo/small_N/",fnBatch,".rds"))
    batch_fn_largeN=readRDS(batch_fn_largeN_fn)
    
    MAV_smallN=plot_dot(MAV_fn,single_fn_smallN,"MAV vs BASiC small N")
    MAV_largeN=plot_dot(MAV_fn,single_fn_largeN,"MAV vs BASiC large N")
    MAV_smallN_batch=plot_dot(MAV_fn,batch_fn_smallN,"MAV vs BASiC small N batch")
    MAV_largeN_batch=plot_dot(MAV_fn,batch_fn_largeN,"MAV vs BASiC large N batch")
    largeN_smallN=plot_dot(single_fn_largeN,single_fn_smallN,xlab="BASiC large N",ylab="BASiC small N", "BASiC not merge replicates:\nlarge N vs BASiC small N ")
    largeN_smallN_batch=plot_dot(batch_fn_largeN,batch_fn_smallN, xlab="BASiC large N",ylab="BASiC small N","BASiC merge replicates:\nlarge N vs BASiC small N ")
    pdf(paste0(NME_MAV_human_out_dir,"BASiC_QC/MAV_smallN_largeN_",fn,".pdf"))
    print(grid.arrange(MAV_smallN$plotOut,MAV_largeN$plotOut,
                       MAV_smallN_batch$plotOut,MAV_largeN_batch$plotOut,
                       largeN_smallN$plotOut,largeN_smallN_batch$plotOut,nrow=3))
    dev.off()
    return(data.table(
      fileName=fn,
      corMAV_smallN=MAV_smallN$corOut,
      corMAV_largeN=MAV_largeN$corOut,
      corMAV_smallN_batch=MAV_smallN_batch$corOut,
      corMAV_largeN_batch=MAV_largeN_batch$corOut,
      corlargeN_smallN=largeN_smallN$corOut,
      corlargeN_smallN_batch=largeN_smallN_batch$corOut

    ))
  }
}

correlation_all=do.call(rbind,
                        lapply(dir("../downstream/input/human_analysis/HCL_scRNA_BASiC/",pattern=".rds"),function(x) BASiC_MAV_comp(gsub(".rds","",x)))
)
saveRDS(correlation_all,paste0(NME_MAV_human_out_dir,"BASiC_QC/correlation_all.rds"))

colMeans(correlation_all[,-1])
  #  corMAV_smallN          corMAV_largeN    corMAV_smallN_batch 
  #            0.5978065              0.6852300              0.5679130 
  #  corMAV_largeN_batch       corlargeN_smallN corlargeN_smallN_batch 
  #            0.6262823              0.6782292              0.7984331 