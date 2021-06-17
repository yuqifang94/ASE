#GO annotation
source('mainFunctions_sub.R')
#Only save the output, not saving the csv file etc
source('GO_run_tissue_perm.R')

#calling code
#tissue, enhancer type,cutoff=3 permutation array number
args = commandArgs(trailingOnly=TRUE)
dir_in='zero_cutoff_mm10'
tissue=args[1]
enc_type=args[2]
cutoff=as.numeric(args[3])
array_id=args[4]

# chromHMM enhancer -------------------------------------------------------
tt1=proc.time()[[3]]
GO_out_all=GO_run_tissue_perm(tissue,dir_in,enc_type,dist_cutoff=cutoff,permute = T)
saveRDS(GO_out_all,paste0('../downstream/output/',tissue,"_GO_",enc_type,'_',cutoff,'_permuted_',array_id,'.rds'))
print(paste0("finish in:",proc.time()[[3]]-tt1))
dir_in="cutoff_test/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_01_cutoff_heart_pc.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_01_cutoff_forebrain_pc.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_01_cutoff_limb_pc.rds')

dir_in="cutoff_test/uc_0/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_0_cutoff_forebrain_pc.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_0_cutoff_limb_pc.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_0_cutoff_pc')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_0_cutoff_pc')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_0_cutoff_pc')

dir_in="cutoff_test/uc_0_01/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_001_cutoff_heart_pc.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_001_cutoff_forebrain_pc.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_001_cutoff_limb_pc.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_001_cutoff_pc')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_001_cutoff_pc')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_001_cutoff_pc')

dir_in='../downstream/input/cluster_all_tissue_pc/uc_0_1/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
}
bg=unique(csv_all$gene)#heart 18096, all 19044

dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_all.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_all.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_all.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(6,23,33,37,42,48,50))
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(3,7,19,25,30,43,44,47,49))
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_01_cutoff_pc_50clu_FC_sel2',
                clu_in=c(3,9,15,19,20,27,31,34,35,38,41,46,48))

GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_tissuespecific",dist_cutoff=-1,permute=F)
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_tissuespecific",dist_cutoff=-1,permute=F)
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_tissuespecific",dist_cutoff=-1,permute=F)
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_zero_cutoff')
#All gene as background
dir_in='../downstream/input/cutoff_test/uc_0_1/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bg=unique(csv_all$gene)#heart 18096, all 19044
dir_in="cutoff_test/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_01_cutoff_heart_pc_bg_all.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_01_cutoff_forebrain_pc_bg_all.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_01_cutoff_limb_pc_bg_all.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_01_cutoff_pc_all_bg')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_01_cutoff_pc_all_bg')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_01_cutoff_pc_all_bg')
#ts genes

dir_in='../downstream/input/cutoff_01_pc_ts/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bg=unique(csv_all$gene)#15603
dir_in="cutoff_01_pc_ts"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_01_cutoff_heart_pc_bg_all_ts.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_01_cutoff_forebrain_pc_bg_all_ts.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_01_cutoff_limb_pc_bg_all_ts.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_01_cutoff_pc_all_bg_ts_FC')

dir_in='../downstream/input/cutoff_01_pc_ts/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bg=unique(csv_all$gene)#15603
dir_in="cutoff_01_pc_ts"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_01_cutoff_heart_pc_bg_all_ts.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_01_cutoff_forebrain_pc_bg_all_ts.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_01_cutoff_limb_pc_bg_all_ts.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_01_cutoff_pc_all_bg_ts')

#ts bg

dir_in="cutoff_01_pc_ts"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=NULL)#4867 bg 
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_01_cutoff_heart_pc_ts.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=NULL)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_01_cutoff_forebrain_pc_ts.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=NULL)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_01_cutoff_limb_pc_ts.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_01_cutoff_pc_ts')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_01_cutoff_pc_ts')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_01_cutoff_pc_ts')

dir_in='../downstream/input/cutoff_test/uc_0/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bg=unique(csv_all$gene)#19821
dir_in="cutoff_test/uc_0/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_all.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_all.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"chromHMM_enhancer_union",dist_cutoff=-1,permute=F,bg=bg)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_all.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'chromHMM_union_0_cutoff_pc_all_bg')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'chromHMM_union_0_cutoff_pc_all_bg')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'chromHMM_union_0_cutoff_pc_all_bg')

# Bin enhancer Figure S8c -------------------------------------------------
#All tissue cluster
dir_in='../downstream/input/cluster_all_tissue_pc/uc_0_1/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
}
bin_enhancer_in=readRDS('../downstream/output/bin_enhancer.rds')
bin_enhancer_in=resize(bin_enhancer_in,width=width(bin_enhancer_in)+1000*2,fix='center')
bg=subsetByOverlaps(bin_enhancer_in,convert_GR(csv_all$region))$`Target Gene`
#heart 18096, all 19044
dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_all_bin_1k.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_all_bin_1k.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_all_bin_1k.rds')
#No result
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_01_cutoff_pc_50clu_FC_sel_bin_1k',
                clu_in=c(6,23,33,37,42,48,50))
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_01_cutoff_pc_50clu_FC_sel_bin_1k',
                clu_in=c(3,7,19,25,30,43,44,47,49))
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_01_cutoff_pc_50clu_FC_sel_bin_1k',
                clu_in=c(3,9,15,19,20,27,31,34,35,38,41,46,48))
#ts background gene

dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_ts_bin_1k.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_ts_bin_1k.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_ts_bin_1k.rds')
#ts background gene no extension

dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_ts_bin_0k.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_ts_bin_0k.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_ts_bin_0k.rds')
#Zero cutoff
dir_in='../downstream/input/cutoff_test/uc_0/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bin_enhancer_in=readRDS('../downstream/output/bin_enhancer.rds')
bin_enhancer_in=resize(bin_enhancer_in,width=width(bin_enhancer_in)+1000*2,fix='center')
bg=subsetByOverlaps(bin_enhancer_in,convert_GR(csv_all$region))$`Target Gene`
dir_in="cutoff_test/uc_0/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=500)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_all_bin_500.rds')
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_all_bin_1k.rds')#1k not significant result
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_ts_bin_1k.rds')
#Zero cutoff, ts entend 0k
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_ts_bin_0k.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_ts_bin_0k.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_ts_bin_0k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_0_cutoff_pc_ts_bg_bin_0k')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_0_cutoff_pc_ts_bg_bin_0k')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_0_cutoff_pc_ts_bg_bin_0k')

#Zero cutoff 0k
dir_in='../downstream/input/cutoff_test/uc_0/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bin_enhancer_in=readRDS('../downstream/output/bin_enhancer.rds')
bin_enhancer_in=resize(bin_enhancer_in,width=width(bin_enhancer_in)+0*2,fix='center')
bg=subsetByOverlaps(bin_enhancer_in,convert_GR(csv_all$region))$`Target Gene`
dir_in="cutoff_test/uc_0/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,extend=0)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_all_bin_0k.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_all_bin_0k.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_all_bin_0k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_0_cutoff_pc_all_bg_bin_0k')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_0_cutoff_pc_all_bg_bin_0k')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_0_cutoff_pc_all_bg_bin_0k')

#0.1 cutoff ts
dir_in='../downstream/input/cutoff_01_pc_ts/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bg=unique(csv_all$gene)#15603
dir_in="cutoff_01_pc_ts"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_01_cutoff_heart_pc_bg_all_ts.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_01_cutoff_forebrain_pc_bg_all_ts.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_01_cutoff_limb_pc_bg_all_ts.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_01_cutoff_pc_all_bg_ts')


# reading and plotting ----------------------------------------------------
#Forebrain
GO_out_all_forebrain=readRDS('../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_all_bin_0k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_0_cutoff_pc_all_bg_bin_0k')
GO_out_all_forebrain=readRDS('../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_ts_bin_1k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_1k_0_cutoff_pc_ts_bg')
GO_out_all_forebrain=readRDS('../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_ts_bin_1k_promoter.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'promoter_0_cutoff_pc_ts_bg')

plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_0k_0_1_cutoff_pc_ts_bg')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_0k_0_1_cutoff_pc_ts_bg')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_0k_0_1_cutoff_pc_ts_bg')


plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'promoter_0_1_cutoff_pc_ts_bg')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'promoter_0_1_cutoff_pc_ts_bg')
#No significance
GO_out_all_forebrain=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_all_bin_1k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_1k_50clu_0_cutoff_pc_all_bg')
#No significance
GO_out_all_forebrain=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_all_bin_1k_promoter.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_1k_promoter_50clu_01_cutoff_pc_all_bg')
GO_out_all_forebrain=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_all_promoter_2k.rds')#Typo before, this is actually tissue-specific bg
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'promoter_50clu_01_cutoff_pc_ts_bg')
#error
GO_out_all_forebrain=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_ts_bin_0k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_0k_50clu_01_cutoff_pc_ts_bg')

#Heart
GO_out_all_heart=readRDS('../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_all_bin_0k.rds')#Typo: also ts bg
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_0k_0_cutoff_pc_ts_bg')
#No significance
GO_out_all_heart=readRDS('../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_all_bin_500.rds')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_500_0_cutoff_pc_all_bg')

GO_out_all_heart=readRDS('../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_ts_bin_1k.rds')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_1k_0_cutoff_pc_ts_bg')

#No significance
GO_out_all_heart=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_all_bin_1k.rds')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_1k_50clu_01_cutoff_pc_all_bg')
#No significance
GO_out_all_heart=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_all_bin_1k_promoter.rds')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_1k_promoter_50clu_01_cutoff_pc_all_bg')

GO_out_all_heart=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_all_promoter_2k.rds')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'promoter_50clu_01_cutoff_pc_all_bg')

GO_out_all_heart=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_ts_bin_0k.rds')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_0k_50clu_01_cutoff_pc_ts_bg')

#Limb
GO_out_all_limb=readRDS('../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_all_bin_0k.rds')#ts bg typo
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_0k_0_cutoff_pc_ts_bg')
GO_out_all_limb=readRDS('../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_ts_bin_1k.rds')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_1k_0_cutoff_pc_ts_bg')

GO_out_all_limb=readRDS('../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_ts_bin_1k_promoter.rds')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_1k_promoter_0_cutoff_pc_ts_bg')
#No result
GO_out_all_limb=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_all_bin_1k.rds')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_1k_50clu_01_cutoff_pc_all_bg')
#No result
GO_out_all_limb=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_all_bin_1k_promoter.rds')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_1k_promoter_50clu_01_cutoff_pc_all_bg')

GO_out_all_limb=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_all_promoter_2k.rds')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'promoter_50clu_01_cutoff_pc_all_bg')

#error
GO_out_all_limb=readRDS('../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_ts_bin_0k.rds')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_0k_50clu_01_cutoff_pc_ts_bg')

# Promoters ---------------------------------------------------------------
#All tissue cluster
dir_in='../downstream/input/cluster_all_tissue_pc/uc_0_1/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
}
bg=unique(csv_all[distance<=2000]$gene)#15915
dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"promoter",dist_cutoff=2000,permute=F,bg=bg)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_all_promoter_2k.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"promoter",dist_cutoff=2000,permute=F,bg=bg)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_all_promoter_2k.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"promoter",dist_cutoff=2000,permute=F,bg=bg)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_all_promoter_2k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'_promoter_2k_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(6,23,33,37,42,48,50))
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'promoter_2k_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(3,7,19,25,30,43,44,47,49))
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'promoter_2k_01_cutoff_pc_50clu_FC_sel2',
                clu_in=c(3,9,15,19,20,27,31,34,35,38,41,46,48))
#Promoters ts bg
dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"promoter",dist_cutoff=2000,permute=F)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_ts_promoter_2k.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"promoter",dist_cutoff=2000,permute=F)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_ts_promoter_2k.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"promoter",dist_cutoff=2000,permute=F)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_ts_promoter_2k.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'_promoter_2k_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(6,23,33,37,42,48,50))
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'promoter_2k_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(3,7,19,25,30,43,44,47,49))
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'promoter_2k_01_cutoff_pc_50clu_FC_sel2',
                clu_in=c(3,9,15,19,20,27,31,34,35,38,41,46,48))

#Zero cutoff
dir_in='../downstream/input/cutoff_test/uc_0/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bg=unique(csv_all$gene)#19821
dir_in="cutoff_test/uc_0/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_all.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_all.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_all.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_0_cutoff_pc_all_bg')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_0_cutoff_pc_all_bg')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_0_cutoff_pc_all_bg')

#0.1 cutoff ts
dir_in='../downstream/input/cutoff_01_pc_ts/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
  
  
}
bg=unique(csv_all$gene)#15603
dir_in="cutoff_01_pc_ts"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_01_cutoff_heart_pc_bg_all_ts.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_01_cutoff_forebrain_pc_bg_all_ts.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_01_cutoff_limb_pc_bg_all_ts.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_01_cutoff_pc_all_bg_ts')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_01_cutoff_pc_all_bg_ts')


# Bin enhancer + promoter -------------------------------------------------

#All tissue cluster
dir_in='../downstream/input/cluster_all_tissue_pc/uc_0_1/full/'
csv_all=data.table()
for (fn in dir(dir_in)){
  csv_all=rbind(csv_all,fread(paste0(dir_in,fn)))
}
bin_enhancer_in=readRDS('../downstream/output/bin_enhancer.rds')
bin_enhancer_in=resize(bin_enhancer_in,width=width(bin_enhancer_in)+1000*2,fix='center')
bg=subsetByOverlaps(bin_enhancer_in,convert_GR(csv_all$region))$`Target Gene`
bg=unique(bg,csv_all[distance<=2000]$gene)#17451

dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer_promoter",dist_cutoff=2000,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_all_bin_1k_promoter.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=2000,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_all_bin_1k_promoter.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=2000,permute=F,bg=bg,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_all_bin_1k_promoter.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(6,23,33,37,42,48,50))
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_01_cutoff_pc_50clu_FC_sel',
                clu_in=c(3,7,19,25,30,43,44,47,49))
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_01_cutoff_pc_50clu_FC_sel2',
                clu_in=c(3,9,15,19,20,27,31,34,35,38,41,46,48))

dir_in="cluster_all_tissue_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer_promoter",dist_cutoff=2000,permute=F,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_50clu_01_cutoff_heart_pc_bg_ts_bin_1k_promoter.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer",dist_cutoff=2000,permute=F,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_50clu_01_cutoff_forebrain_pc_bg_ts_bin_1k_promoter.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer",dist_cutoff=2000,permute=F,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_50clu_01_cutoff_limb_pc_bg_ts_bin_1k_promoter.rds')

#0 cutoff
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer_promoter",dist_cutoff=2000,permute=F,extend=1000)
saveRDS(GO_out_all_heart,'../downstream/output/GO_out_all_0_cutoff_heart_pc_bg_ts_bin_1k_promoter.rds')
GO_out_all_forebrain=GO_run_tissue_perm('forebrain',dir_in,"bin_enhancer_promoter",dist_cutoff=2000,permute=F,extend=1000)
saveRDS(GO_out_all_forebrain,'../downstream/output/GO_out_all_0_cutoff_forebrain_pc_bg_ts_bin_1k_promoter.rds')
GO_out_all_limb=GO_run_tissue_perm('limb',dir_in,"bin_enhancer_promoter",dist_cutoff=2000,permute=F,extend=1000)
saveRDS(GO_out_all_limb,'../downstream/output/GO_out_all_0_cutoff_limb_pc_bg_ts_bin_1k_promoter.rds')
plot_GO_heatmap('forebrain',"GO_out_cluster_all",GO_out_all_forebrain,'bin_enhancer_0_cutoff_pc_all_bg_bin_1k')
plot_GO_heatmap('heart',"GO_out_cluster_all",GO_out_all_heart,'bin_enhancer_0_cutoff_pc_all_bg_bin_1k')
plot_GO_heatmap('limb',"GO_out_cluster_all",GO_out_all_limb,'bin_enhancer_0_cutoff_pc_all_bg_bin_1k')


# other test --------------------------------------------------------------


gtf <- import('../downstream/input/grcm38.gtf.gz')
gtf=gtf[which(gtf$transcript_biotype=="protein_coding"&gtf$type=="gene")]
tss <- promoters(gtf,upstream=0,downstream=1)
dir_zero='../downstream/input/zero_cutoff_union_chromHMM/'
seqlevels(tss)=paste0('chr',seqlevels(tss))
for(fn in dir(dir_zero)){
  csv_in =fread(paste0(dir_zero,fn))
  gene_check=nearest(convert_GR(csv_in$region),tss)
  csv_in$gene=tss$gene_name[gene_check]
  csv_in$dist=mcols(distanceToNearest(convert_GR(csv_in$region),tss))$distance
  write.csv(csv_in,paste0('../downstream/input/zero_cutoff_union_chromHMM_protein_coding/',fn),row.names = F)
}

genes=genes[intersect(gn,h_gene_sets$gene_symbol)]
tss <- promoters(genes,upstream=0,downstream=1)
library(qusage)
library(msigdbr)
library(fgsea)
all_gene_sets = msigdbr(species = "Mus musculus")
h_gene_sets = msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:BP")
h_gene_sets_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
csv_in=fread('../downstream/input/zero_cutoff_union_chromHMM_protein_coding/limb.csv')
csv_in_maxUC=csv_in[chromHMM_enhancer_union==TRUE,list(UC_max=max(UC_maxUC)),by=list(gene)]
ranks=csv_in_maxUC$UC_max
names(ranks)=csv_in_maxUC$gene
gsea_out_heart<- fgsea(h_gene_sets_list, ranks, minSize=10, maxSize = 500, nperm=100000)


gsea_out=vector(mode = "list", length = 10)
for(i in 1:10){
  csv_in_maxUC=csv_in[chromHMM_enhancer_union==TRUE&cluster==i,list(UC_max=mean(UC_maxUC)),by=list(gene)]
  ranks=csv_in_maxUC$UC_max
  names(ranks)=csv_in_maxUC$gene
  gsea_out_clu<- fgsea(h_gene_sets_list, ranks, minSize=10, maxSize = 500, nperm=100000)
  gsea_out_clu$cluster=i
  gsea_out[[i]]=gsea_out_clu
}
gsea_out_top=do.call('rbind',lapply(gsea_out,function(x) x[padj<=0.1][order(padj,-NES)][1:5]))

#Getting Enhancer from Bin
bin_supp=as.data.table(read_xlsx('../downstream/input/mouse_analysis/enhancer_selection//bin_supp.xlsx',sheet='S8c.Enhancer-Gene-Map-Replicatd',skip=1))
bin_supp_gr=GRanges(seqnames=bin_supp$chrom,ranges=IRanges(start=bin_supp$start,end=bin_supp$end))
mcols(bin_supp_gr)=bin_supp[,c(-1,-2,-3)]
colnames(mcols(bin_supp_gr))=gsub('\\...*','',colnames(mcols(bin_supp_gr)))
saveRDS(bin_supp_gr,'../downstream/output/bin_enhancer.rds')


#All possible combination for all tissue clustering
dir_target="cluster_all_tissue_pc/"
uc=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
uc_gr=lapply(uc,function(x) granges(x))
names(uc_gr)=NULL
uc_gr=unique(do.call(c,uc_gr))
enhancer=readRDS('../downstream/output/bin_enhancer.rds')
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)
dir_in='../downstream/output/cluster_ts_10_pc/'
gtf=fread('../downstream/input/grcm38.gtf',data.table=F)

gtf <- gtf[gtf[,3]=='gene',]
type <- sub('\".*','',sub('.*gene_type \"','',gtf[,9]))
gtf <- gtf[type=='protein_coding',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
names(gr) <- gn
tss <- promoters(gr,upstream=0,downstream=1)
bg_promoter=names(subsetByOverlaps(tss,uc_gr,maxgap = 2000))
#GO_all_variabile("cluster_all_tissue_pc/")
GO_all_variabile("cluster_ts_10_pc/",bg_enhancer=bg_enhancer)
GO_all_variabile("cluster_ts_10_pc/",tissue_all=c("heart","forebrain","limb","NT","EFP","hindbrain","liver","midbrain"),
                 bg_enhancer=bg_enhancer,bg_promoter=bg_promoter)

uc=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
uc_gr=lapply(uc,function(x) granges(x))
names(uc_gr)=NULL
uc_gr=unique(do.call(c,uc_gr))
enhancer=readRDS('../downstream/output/bin_enhancer.rds')
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)
GO_all_variabile("cluster_tissue_cutoff_pc/",bg_enhancer=bg_enhancer)
GO_all_variabile("cluster_tissue_cutoff_pc/",tissue_all=c("NT","EFP","hindbrain","liver","midbrain"),bg_enhancer=bg_enhancer)
GO_all_variabile("cluster_ts_10_pc/",GO_type="GSEA",ranking_stat="dMML_maxUC")
GO_all_variabile("cluster_ts_10_pc/",GO_type="GSEA",ranking_stat="dNME_maxUC")
#testing
dir_in="cluster_ts_10_pc/uc_0_1/full"
GO_out_all_heart=GO_run_tissue_perm('heart',dir_in,"bin_enhancer",dist_cutoff=-1,permute=F,GO_type="GSEA",ranking_stat="dNME_maxUC")
# Processing dMML and dNME ------------------------------------------------
dir_in='../downstream/output/cluster_ts_10_pc/'
for(fn in dir(dir_in,pattern='.rds')){
  
  
  plot_GO_heatmap_variable(paste0(dir_in,fn),clu_in=1:10)
}
dir_in='../downstream/output/cluster_ts_10_pc_annotated_cutoff/'
for(fn in dir(dir_in,pattern='.rds')){
  
  
  plot_GO_heatmap_variable(paste0(dir_in,fn),clu_in=1:10)
}

