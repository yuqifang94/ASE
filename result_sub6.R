rm(list=ls())
source("mainFunctions_sub.R")
#read in JSD data
UC_raw=readRDS('../downstream/output/uc_matrix_DNase.rds')
theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
# violin plot -----------------------------------------------------------
#From mm10_reformating.R
UC_in=readRDS('../downstream/output/uc_matrix_DNase.rds')
dnme=readRDS('../downstream/output/dnme_matrix_DNase.rds')
dmml=readRDS('../downstream/output/dnme_matrix_DNase.rds')
matrix_all=list()
for(sp in names(UC)){
  dmml_in=dmml[[sp]]
  dnme_in=dnme[[sp]]
  UC_in=UC[[sp]]
  colnames(dmml_in)=paste0(colnames(dmml_in),'_','dMML')
  colnames(dnme_in)=paste0(colnames(dnme_in),'_','dNME')
  colnames(UC_in)=paste0(colnames(UC_in),'_','UC')
  matrix_tissue=as.data.table(cbind(UC_in,dmml_in[rownames(UC_in),],dnme_in[rownames(UC_in),]),keep.rownames = T)
  colnames(matrix_tissue)[1]="regions"
  olap_chromHMM=findOverlaps(convert_GR(matrix_tissue$regions),chromHMM[chromHMM$tissue==sp])
  olap_promoter=findOverlaps(convert_GR(matrix_tissue$regions),promoters)
  matrix_tissue$promoter=FALSE
  matrix_tissue$enhancer=FALSE
  matrix_tissue$tissue=sp
  matrix_tissue$promoter[queryHits(olap_promoter)]=TRUE
  matrix_tissue$enhancer[queryHits(olap_chromHMM)]=TRUE
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
matrix_all_max_change=matrix_all[,list(value=max(value)),by=list(stat,regions,stage,tissue,region_type)]
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

pdf('../downstream/output/graphs/Figure5/FigureS7_dNME_dMML_enhancer_boxplot_max.pdf',width=3.5,height=3.5)
ggplot(matrix_all_max_change[region_type!="Non-regulatory"&stat!="UC"],aes(x=stat,y=value,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
   xlab('')+ylab('')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.4))
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
#add bin enhancer if necessary
#Plotting GO annotations
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_all",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME_only",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_MML_only",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_NME_MML",GO_out_all)
plot_GO_heatmap(c('forebrain','heart','limb'),"GO_out_cluster_non_NME_non_MML",GO_out_all)
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





