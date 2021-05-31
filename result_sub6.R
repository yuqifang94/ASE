rm(list=ls())
source("mainFunctions_sub.R")
#read in JSD data
theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
# violin plot -----------------------------------------------------------
#From mm10_reformating.R
UC=readRDS('../downstream/output/mouse_analysis/CPEL_outputs/uc_matrix_DNase.rds')
dnme=readRDS('../downstream/output/mouse_analysis/CPEL_outputs/dnme_matrix_DNase.rds')
dmml=readRDS('../downstream/output/mouse_analysis/CPEL_outputs/dmml_matrix_DNase.rds')
bin_enhancer_in=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')
gtf=fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table=F)
promoter_in=gtf <- gtf[gtf[,3]=='gene',]
type <- sub('\".*','',sub('.*gene_type \"','',gtf[,9]))
gtf <- gtf[type=='protein_coding',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
names(gr) <- gn
tss <- promoters(gr,upstream=0,downstream=1)
matrix_all=list()
for(sp in names(UC)){
  dmml_in=dmml[[sp]]
  dnme_in=dnme[[sp]]
  UC_in=UC[[sp]]
  region_intersect=intersect(rownames(UC_in),intersect(rownames(dmml_in),rownames(dnme_in)))
  colnames(dmml_in)=paste0(colnames(dmml_in),'_','dMML')
  colnames(dnme_in)=paste0(colnames(dnme_in),'_','dNME')
  colnames(UC_in)=paste0(colnames(UC_in),'_','UC')
  matrix_tissue=as.data.table(cbind(UC_in[region_intersect,],dmml_in[region_intersect,],dnme_in[region_intersect,]),keep.rownames = T)
  colnames(matrix_tissue)[1]="regions"
  olap_enhancer=findOverlaps(convert_GR(matrix_tissue$regions),bin_enhancer_in)
  olap_promoter=findOverlaps(convert_GR(matrix_tissue$regions),tss,maxgap = 2000)
  matrix_tissue$promoter=FALSE
  matrix_tissue$enhancer=FALSE
  matrix_tissue$tissue=sp
  matrix_tissue$promoter[queryHits(olap_promoter)]=TRUE
  matrix_tissue$enhancer[queryHits(olap_enhancer)]=TRUE
  matrix_tissue=melt.data.table(matrix_tissue,id.vars = c("regions","promoter","enhancer","tissue"))
  matrix_tissue$stage=sub('_.*','',matrix_tissue$variable)
  matrix_tissue$stat=sub('.*_','',matrix_tissue$variable)
  #matrix_tissue=dcast.data.table(matrix_tissue,tissue+enhancer+promoter+stage+regions~stat)


  matrix_all[[sp]]=matrix_tissue
}
matrix_all=fastDoCall('rbind',matrix_all)
#matrix_all=matrix_all[matrix_all$UC!=1]
# matrix_all$UC_quant=findInterval(matrix_all$UC,quantile(matrix_all$UC,prob=c(0,0.25,0.5,0.75,1),na.rm=T))
# 
# matrix_all$UC_quant[matrix_all$UC_quant==5]=4#5th quantile is the maximum number, move to 4th
# quant_conv=c("Q1","Q2","Q3","Q4")
# matrix_all$UC_quant=quant_conv[matrix_all$UC_quant]
#Separate  based on region type
matrix_all_non_reg=matrix_all[(!promoter)&(!enhancer)]
matrix_all_non_reg$region_type="Non-regulatory"
matrix_all_enhancer=matrix_all[enhancer==TRUE]
matrix_all_enhancer$region_type="enhancer"
matrix_all_promoter=matrix_all[promoter==TRUE]
matrix_all_promoter$region_type="promoter"
matrix_all=rbind(matrix_all_non_reg,matrix_all_enhancer,matrix_all_promoter)
cat(paste0(round(table(matrix_all$region_type)/length(matrix_all$region_type)*100,digits=2),"%"),"\n")
matrix_all_max_change=matrix_all[,list(value=max(value)),by=list(stat,regions,tissue,region_type)]#add stage, it same as matrix_all
# matrix_all_agg=matrix_all[,list(dMML=median(dMML_quant),dNME=median(dNME_quant),UC=median(UC),
#                                 dMML_top25=quantile(dMML,prob=0.75),dNME_top25=quantile(dNME,prob=0.75),
#                                 dMML_bottom25=quantile(dMML,prob=0.25),dNME_bottom25=quantile(dNME,prob=0.25)),by=list(UC_quant,region_type)]
# pdf('../downstream/output/graphs/Figure5/Figure5B_dNME_dMML_enhancer_UC_quantile.pdf',width=3.5,height=3.5)
# dNME_plot=ggplot(matrix_all[region_type!="Non-regulatory"],aes(x=UC_quant,y=dNME,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
#   xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.75))
# dMML_plot=ggplot(matrix_all[region_type!="Non-regulatory"],aes(x=UC_quant,y=dMML,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
#   xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.3))
# ggarrange(dNME_plot,dMML_plot,nrow=2,ncol=1,common.legend=T)
# dev.off()

pdf('../downstream/output/graphs/Figure5/FigureS10_dNME_dMML_enhancer_boxplot_max.pdf',width=3.5,height=3.5)
ggplot(matrix_all_max_change[region_type!="Non-regulatory"&stat!="UC"],aes(x=stat,y=value,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
   xlab('')+ylab('')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.4))
dev.off()
pdf('../downstream/output/graphs/Figure5/FigureS10_promoter_enhancer_boxplot_dNME.pdf',width=3.5,height=3.5)
ggplot(matrix_all_max_change[region_type!="Non-regulatory"&stat=="dNME"],aes(x=region_type      ,y=value,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('')+ylab('')+theme(legend.position = "none",legend.title = element_blank())+ylim(c(0,0.4))
dev.off()
pdf('../downstream/output/graphs/Figure5/FigureS10_promoter_enhancer_boxplot_dMML.pdf',width=3.5,height=3.5)
ggplot(matrix_all_max_change[region_type!="Non-regulatory"&stat=="dMML"],aes(x=region_type      ,y=value,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('')+ylab('')+theme(legend.position = "none",legend.title = element_blank())+ylim(c(0,0.4))
dev.off()
pdf('../downstream/output/graphs/Figure5/FigureS10_promoter_enhancer_boxplot_UC.pdf',width=3.5,height=3.5)
ggplot(matrix_all_max_change[region_type!="Non-regulatory"&stat=="UC"],aes(x=region_type      ,y=value,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('')+ylab('')+theme(legend.position = "none",legend.title = element_blank())+ylim(c(0,0.4))
dev.off()
# Check percentage of FeDMR covered ---------------------------------------

FeDMR=readRDS('../downstream/output/FeDMR.rds')
dMML_cor=readRDS('../downstream/input/dmmlcor.rds')
dNME_cor=readRDS('../downstream/input/dNMEcor.rds')
enhancer_in=readRDS('../downstream/output/chromHMM_enhancer.rds')
#Use bin's enhancer will generate same result, use that instead if necessary
cor_OR=data.table()
for (ts in (names(dMML_cor))){
  enhancer=enhancer_in[enhancer_in$tissue==ts]
  #enhancer=enhancer_in
  dMML_region_ts=convert_GR(names(dMML_cor[[ts]]))
  dMML_region_ts$cor=dMML_cor[[ts]]
  dMML_region_ts=subsetByOverlaps(dMML_region_ts,enhancer)
  FeDMR_ts=FeDMR[rowSums(as.matrix(mcols(FeDMR)[,sub('-.*','',colnames( mcols(FeDMR))) == ts]))>=1]
  FeDMR_ts=subsetByOverlaps(FeDMR_ts,enhancer)
  dNME_region_ts=convert_GR(names(dNME_cor[[ts]]))
  dNME_region_ts$cor=dNME_cor[[ts]]
  dNME_region_ts=subsetByOverlaps(dNME_region_ts,enhancer)
  cor_OR=rbind(cor_OR,cbind(data.table(tissue=ts,statistics="dMML"),cor_dMML_dNME_enrich(dMML_region_ts,0.75,FeDMR_ts)))
  cor_OR=rbind(cor_OR,cbind(data.table(tissue=ts,statistics="dNME"),cor_dMML_dNME_enrich(dNME_region_ts,0.75,FeDMR_ts)))
}
cor_OR=cor_OR[cor_OR$tissue!="liver"]#Ecker excluded liver
cor_OR$FDR=p.adjust(cor_OR$pvalue,method='BH')
cor_OR$sig=add.significance.stars(cor_OR$FDR, cutoffs = c(0.05, 0.01, 0.001))
pdf('../downstream/output/graphs/FigureS7/OR_tissue_Ecker.pdf',width=5,height=3.5)

ggplot(cor_OR,aes(x=tissue,y=OR,fill=statistics))+geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,
                position=position_dodge(.9))+  
  geom_text(aes(y=upperCI,label=round(OR,digits = 2)),vjust=-0.5,position = position_dodge(1),size=3)+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  geom_text(aes(y=upperCI+0.1,label=sig),vjust=-0.5,position = position_dodge(1),size=2)
dev.off()
#Run for different type of analysis
dir_in='mm10_cluster_all'
nme_cor <- readRDS('../downstream/input/dnmecor.rds')
mml_cor <- readRDS('../downstream/input/dmmlcor.rds')
tissue=unique(sub(".csv*","",dir(paste0('../downstream/input/',dir_in),pattern="csv")))

#All regions cluster, using enhancer vs non enhancer
GO_out_all=GO_run_tissue(tissue,dir_in,nme_cor,mml_cor,"chromHMM_enhancer")
saveRDS(GO_out_all,'../downstream/output/GO_out_all_chromHMM_all_stat.rds')
GO_out_all=GO_run_tissue(tissue,dir_in,nme_cor,mml_cor,"non_chromHMM_enhancer")
saveRDS(GO_out_all,'../downstream/output/GO_out_all_non_chromHMM_all_stat.rds')
GO_out_all=GO_run_tissue(tissue,dir_in,nme_cor,mml_cor,"all_regions")
saveRDS(GO_out_all,'../downstream/output/GO_out_all_all_regions_all_stat.rds')


# zerocutoff --------------------------------------------------------------
# for perm_number in {1..20}; do echo ${perm_number}; echo Rscript GO_annotation.R heart FANTOM5 -1 ${perm_number}; 
# Rscript GO_annotation.R heart FANTOM5 -1 ${perm_number}; done >perm.log.FANTOM
dir_in='zero_cutoff_mm10'
dir_in="zero_cutoff_union_chromHMM"
dir_in="zero_cutoff_union_chromHMM_protein_coding"
dir_in="cutoff_test/uc_0_01/full"
tissue=c('heart','forebrain','limb')
GO_out_all=GO_run_tissue_perm(tissue,dir_in,"chromHMM_enhancer_tissuespecific")
GO_out_all=GO_run_tissue_perm(tissue,dir_in,"chromHMM_enhancer_union")
GO_out_all=GO_run_tissue_perm(tissue,dir_in,"FANTOM5")
GO_out_all=GO_run_tissue_perm(tissue,dir_in,"bin_enhancer")
GO_out_all=GO_run_tissue_perm(tissue,dir_in,"chromHMM_enhancer",dist_cutoff = 10000)
dir_in='dmml_dnme_cluster/dMML'
GO_out_all=GO_run_tissue_perm(tissue,dir_in,"chromHMM_enhancer")
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all,'chromHMM_dMML')
dir_in='dmml_dnme_cluster/dNME'
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all,'chromHMM_10k_fixed')
plot_GO_heatmap(c('forebrain','heart'),"GO_out_cluster_all",GO_out_all,'chromHMM_all')
#Permuted data
dir_in='ucrandom_all'
tissue=c('heart','forebrain','limb')
GO_out_all_perm=GO_run_tissue_perm(tissue,dir_in,"chromHMM_enhancer")
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all_perm)

# dMML cluster dNME cluster check -----------------------------------------
dMML_cluster=readRDS('../downstream/input/dmml_dnme_cluster/cluster/dMML.rds')
dNME_cluster=readRDS('../downstream/input/dmml_dnme_cluster/cluster/dNME.rds')
dir_in='../downstream/input/mm10_cluster_all/'
dir_out_dNME='../downstream/input/dmml_dnme_cluster/dNME/'
dir_out_dMML='../downstream/input/dmml_dnme_cluster/dMML/'
for(fn in dir(dir_in,pattern='csv')){
  
  tissue=gsub('.csv','',fn)
  csv_in=fread(paste0(dir_in,fn))
  dMML_out=csv_in
  dMML_out$cluster=dMML_cluster[[tissue]][dMML_out$region]
  dNME_out=csv_in
  dNME_out$cluster=dNME_cluster[[tissue]][dNME_out$region]
  write.csv(dMML_out,paste0(dir_out_dMML,tissue,'.csv'),row.names = F)
  write.csv(dNME_out,paste0(dir_out_dNME,tissue,'.csv'),row.names = F)
}

# recover GO run for all chromHMM for heart and forebrain -----------------
GO_out_all=list()
for(tissue in c('heart','forebrain')){
  dir_in='../downstream/output/mm10_result/zero_cutoff_mm10/chromHMM_enhancer_NA/cluster_GO/'
  GO_out_tissue=vector(mode = "list", length = 10)
  for(i in 1:10){
    GO_out_tissue[[i]]$GO_out_cluster_all=fread(paste0(dir_in,tissue,'-',i,'_cluster_GO.csv'))
    
  }
  GO_out_all[[tissue]]=GO_out_tissue
}
GO_out_all=list()
for(tissue in c('heart','forebrain')){
  dir_in='../downstream/output/mm10_result/zero_cutoff_mm10/permutechromHMM_enhancer_10000/cluster_GO/'
  GO_out_tissue=vector(mode = "list", length = 10)
  for(i in 1:10){
    GO_out_tissue[[i]]$GO_out_cluster_all=fread(paste0(dir_in,tissue,'-',i,'_cluster_GO.csv'))
    
  }
  GO_out_all[[tissue]]=GO_out_tissue
}

# Check permuted GO -------------------------------------------------------
dir_in='../downstream/output/'
fn= dir(dir_in,pattern='chromHMM_enhancer_-1')
fn= dir(dir_in,pattern='GO_FANTOM5_-1_permuted')
permute_heart_test=mclapply(fn,function(x){
  GO_in=readRDS(paste0(dir_in,x))
  filter_pass=do.call(rbind,lapply(GO_in,function(x){
    return(data.table(sig_num=nrow(x$GO_out_cluster_all[FC>=2&FDR<=0.1]),
               cluster=unique(x$GO_out_cluster_all$cluster)))}))
  filter_pass$perm_number=gsub('.rds','',gsub('.*_','',x))
  return(list(filter_pass,GO_in))
  
},mc.cores=18)
#Build ecdf for each cluster
permute_heart_test_ecdf=fastDoCall('rbind',lapply(permute_heart_test,function(x){
  
  return(fastDoCall('rbind',lapply(x[[2]],function(y){
    return(data.table(FC=y$GO_out_cluster$FC, cluster=y$GO_out_cluster$cluster))
    
    
  })))
} ))
permute_heart_test_ecdf_ls=lapply(1:10,function(i) {
  
  return(ecdf(permute_heart_test_ecdf[cluster==i]$FC))
  
})
heart_FC_pval=lapply(heart_FC,function(x){
  x$GO_out_cluster$p_perm=permute_heart_test[[unique(x$GO_out_cluster$cluster)]](x$GO_out_cluster$FC)
  x$GO_out_cluster$FDR_perm=p.adjust(x$GO_out_cluster$p_perm,method='BH')
return(x$GO_out_cluster)
  }
  )

permute_heart_test=do.call(rbind,lapply(permute_heart_test,function(x) x[[1]]))
dcast.data.table(permute_heart_test,cluster~perm_number,value.var = 'sig_num')
# Check overlap between all data and permuted data ------------------------
heart_original=fread('../downstream/input/mm10_cluster_all/heart.csv')
limb_original=fread('../downstream/input/mm10_cluster_all/limb.csv')
heart_perm=fread('../downstream/input/ucrandom_all/heart.csv')
limb_perm=fread('../downstream/input/ucrandom_all/limb.csv')
#Find overlapped regions for all
original=heart_original$region
perm=heart_perm$region
region_dt=data.table(region_type=factor(c('original only','both','perm only'),levels=c('original only','both','perm only')),
                     number_regions=c(sum(!original%in%perm),sum(original%in%perm),sum(!perm%in%original)))
region_dt$number_regions=region_dt$number_regions/sum(region_dt$number_regions)

ggplot(region_dt,aes(x='all',y=number_regions,group=region_type,fill=region_type))+geom_bar(stat='identity')+ylab('percent region')

region_dt_clu=lapply(unique(heart_original$cluster),function(x){
  
  original=heart_original[cluster==x]$region
  perm=heart_perm[cluster==x]$region
  region_dt_x=data.table(region_type=factor(c('original only','both','perm only'),levels=c('original only','both','perm only')),
                       number_regions=c(sum(!original%in%perm),sum(original%in%perm),sum(!perm%in%original)))
  region_dt_x$number_regions=region_dt_x$number_regions/sum(region_dt_x$number_regions)
  region_dt_x$clu=x
  return(region_dt_x)
})
region_dt_clu=do.call(rbind,region_dt_clu)
region_dt_clu$clu=factor(region_dt_clu$clu,levels=c(1:10))
ggplot(region_dt_clu,aes(x=clu,y=number_regions,group=region_type,fill=region_type))+geom_bar(stat='identity')+ylab('percent region')

#Check chromHMM regions

original=heart_original[chromHMM_enhancer ==TRUE]$region
perm=heart_perm[chromHMM_enhancer ==TRUE]$region
region_dt=data.table(region_type=factor(c('original only','both','perm only'),levels=c('original only','both','perm only')),
                     number_regions=c(sum(!original%in%perm),sum(original%in%perm),sum(!perm%in%original)))
region_dt$number_regions=region_dt$number_regions/sum(region_dt$number_regions)
ggplot(region_dt,aes(x='all',y=number_regions,group=region_type,fill=region_type))+geom_bar(stat='identity')+ylab('percent region')
#Check genes

original=unique(heart_original[chromHMM_enhancer ==TRUE]$gene)
perm=unique(heart_perm[chromHMM_enhancer ==TRUE]$gene)
gene_dt=data.table(region_type=factor(c('original only','both','perm only'),levels=c('original only','both','perm only')),
                     number_regions=c(sum(!original%in%perm),sum(original%in%perm),sum(!perm%in%original)))
gene_dt$number_regions=gene_dt$number_regions/sum(gene_dt$number_regions)
ggplot(gene_dt,aes(x='all',y=number_regions,group=region_type,fill=region_type))+geom_bar(stat='identity')+ylab('percent region')

gene_dt_clu=lapply(unique(heart_original$cluster),function(x){
  
  original=heart_original[cluster==x&chromHMM_enhancer ==TRUE]$gene
  perm=heart_perm[cluster==x&chromHMM_enhancer ==TRUE]$gene
  gene_dt_x=data.table(region_type=factor(c('original only','both','perm only'),levels=c('original only','both','perm only')),
                         number_regions=c(sum(!original%in%perm),sum(original%in%perm),sum(!perm%in%original)))
  #gene_dt_x$number_regions=gene_dt_x$number_regions/sum(gene_dt_x$number_regions)
  gene_dt_x$clu=x
  return(gene_dt_x)
})
gene_dt_clu=do.call(rbind,gene_dt_clu)
gene_dt_clu$clu=factor(gene_dt_clu$clu,levels=c(1:10))
ggplot(gene_dt_clu,aes(x=clu,y=number_regions,group=region_type,fill=region_type))+geom_bar(stat='identity')+ylab('percent region')

tissue=c('heart','forebrain','limb')
GO_out_all_Top30=GO_run_tissue_perm(tissue,'Top_30_cutoff',"chromHMM_enhancer")
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all_Top30)

GO_out_all_Top_20=GO_run_tissue_perm(tissue,'Top_20_cutoff',"chromHMM_enhancer")
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all_Top_20)

GO_out_all_Top_10=GO_run_tissue_perm(tissue,'Top_10_cutoff',"chromHMM_enhancer")
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all_Top_10)

GO_out_all_Top_1=GO_run_tissue_perm(tissue,'Top_1_cutoff',"chromHMM_enhancer")
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all_Top_1)


plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all,'FANTOM5')
#add bin enhancer if necessary
#Plotting GO annotations
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME_only",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_MML_only",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME_MML",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_non_NME_non_MML",GO_out_all)


plot_GO_heatmap_perm(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all_Top30)
# dir_in='mm10_cluster_chromHMM'
# GO_out_all=GO_run_tissue(tissue,dir_in,nme_cor,mml_cor,"chromHMM_enhancer")
# saveRDS(GO_out_all,'../downstream/output/GO_out_chromHMM_chromHMM_all_stat.rds')

# pdf('../downstream/output/graphs/FigureS6/correlation.pdf',width=5,height=5)
# 
# dNME_cor_plot=ggplot(data.frame(dNME_cor=unlist(dNME_cor)),aes(x=dNME_cor))+geom_density(fill="light blue",color='light blue')+
#   xlab("dNME UC correlation")+theme_glob
# dMML_cor_plot=ggplot(data.frame(dMML_cor=unlist(dMML_cor)),aes(x=dMML_cor))+geom_density(fill="pink",color='pink')+
#   xlab("dMML UC correlation")+theme_glob
# ggarrange(dNME_cor_plot,dMML_cor_plot,nrow=2,ncol=1)
# dev.off()





