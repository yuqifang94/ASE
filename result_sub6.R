rm(list=ls())
source("mainFunctions_sub.R")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Mus.musculus)
library(Gmisc)
library(ggplot2)
library(data.table)
#read in JSD data
UC=readRDS('../downstream/output/uc_matrix_DNase.rds')
chromHMM=readRDS('../downstream/output/chromHMM_enhancer.rds')
mml <- readRDS('../downstream/output/mml_matrix_DNase.rds')
nme <- readRDS('../downstream/output/nme_matrix_DNase.rds')
convert_GR<-function(x){
  gr=GRanges(seqnames=sub(':.*','',x),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',x))),
                     end=as.numeric(sub('.*-','',x))))
  return(gr)
  
}
theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- GenomicFeatures::genes(txdb)
promoters <- promoters(genes,upstream=2000,downstream=1000)
promoters$gene_name=AnnotationDbi::select(Mus.musculus,key=as.character(promoters$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))

UC <- lapply(UC,function(i) {
  i <- i[rowSums(i > 0.1) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  i <- i[,order(match(colnames(i),timeorder))]
  i <- i[complete.cases(i),]
  return(i)
})

#dMML and dNME
dmml=list()
dnme=list()
for (n in names(UC)){
  dmml[[n]] <- sapply(colnames(UC[[n]]),function(i) {
    time <- strsplit(i,'-')
    sapply(time,function(time_in)
      {abs(mml[rownames(UC[[n]]),paste0(n,'-',time_in[1],'-all')]-mml[rownames(UC[[n]]),paste0(n,'-',time_in[2],'-all')])})
  })
  colnames(dmml[[n]])=colnames(UC[[n]])
  rownames(dmml[[n]])=rownames(UC[[n]])

  
  dnme[[n]] <- sapply(colnames(UC[[n]]),function(i) {
    time <- strsplit(i,'-')
  sapply(time,function(time_in)
    {abs(nme[rownames(UC[[n]]),paste0(n,'-',time_in[1],'-all')]-nme[rownames(UC[[n]]),paste0(n,'-',time_in[2],'-all')])})

  })
  colnames(dnme[[n]])=colnames(UC[[n]])
  rownames(dnme[[n]])=rownames(UC[[n]])
}


# quantile plot -----------------------------------------------------------
matrix_all=list()
for(sp in names(UC)){
  dmml_in=dmml[[sp]]
  dnme_in=dnme[[sp]]
  UC_in=UC[[sp]]
  colnames(dmml_in)=paste0(colnames(dmml_in),'_','dMML')
  colnames(dnme_in)=paste0(colnames(dnme_in),'_','dNME')
  colnames(UC_in)=paste0(colnames(UC_in),'_','UC')
  matrix_tissue=as.data.table(cbind(UC_in,dmml_in,dnme_in),keep.rownames = T)
  colnames(matrix_tissue)[1]="regions"
  olap_chromHMM=findOverlaps(convert_GR(matrix_tissue$regions),chromHMM)
  olap_promoter=findOverlaps(convert_GR(matrix_tissue$regions),promoters)
  matrix_tissue$promoter=FALSE
  matrix_tissue$enhancer=FALSE
  matrix_tissue$tissue=sp
  matrix_tissue$promoter[queryHits(olap_promoter)]=TRUE
  matrix_tissue$enhancer[queryHits(olap_chromHMM)]=TRUE
  matrix_tissue=melt.data.table(matrix_tissue,id.vars = c("regions","promoter","enhancer","tissue"))
  matrix_tissue$stage=sub('_.*','',matrix_tissue$variable)
  matrix_tissue$stat=sub('.*_','',matrix_tissue$variable)
  matrix_tissue=dcast.data.table(matrix_tissue,tissue+enhancer+promoter+stage+regions~stat)


  matrix_all[[sp]]=matrix_tissue
}
matrix_all=fastDoCall('rbind',matrix_all)
#matrix_all=matrix_all[matrix_all$UC!=1]
matrix_all$UC_quant=findInterval(matrix_all$UC,quantile(matrix_all$UC,prob=c(0,0.25,0.5,0.75,1),na.rm=T))
# matrix_all$dNME_quant=findInterval(matrix_all$dNME,quantile(matrix_all$dNME,prob=c(0,0.25,0.75,1),na.rm=T))
# matrix_all$dMML_quant=findInterval(matrix_all$dMML,quantile(matrix_all$dMML,prob=c(0,0.25,0.75,1),na.rm=T))
matrix_all$UC_quant[matrix_all$UC_quant==5]=4#5th quantile is the maximum number, move to 4th
quant_conv=c("Q1","Q2","Q3","Q4")
matrix_all$UC_quant=quant_conv[matrix_all$UC_quant]
matrix_all$region_type="Non-regulatory"
matrix_all$region_type[matrix_all$promoter]="Promoter"
matrix_all$region_type[matrix_all$enhancer]="Enhancer"
# matrix_all_agg=matrix_all[,list(dMML=median(dMML_quant),dNME=median(dNME_quant),UC=median(UC),
#                                 dMML_top25=quantile(dMML,prob=0.75),dNME_top25=quantile(dNME,prob=0.75),
#                                 dMML_bottom25=quantile(dMML,prob=0.25),dNME_bottom25=quantile(dNME,prob=0.25)),by=list(UC_quant,region_type)]
pdf('../downstream/output/graphs/Figure5/Figure5B_dNME_dMML_enhancer_UC_quantile.pdf',width=3.5,height=3.5)
dNME_plot=ggplot(matrix_all[!is.na(matrix_all$region_type)],aes(x=UC_quant,y=dNME,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.75))
dMML_plot=ggplot(matrix_all[!is.na(matrix_all$region_type)],aes(x=UC_quant,y=dMML,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.3))
ggarrange(dNME_plot,dMML_plot,nrow=2,ncol=1,common.legend=T)
dev.off()

#Plotting GO terms as examples
#Heat c9, C1
plot_GO_example<-function(GO_in,fill){
  GO_in=GO_in[order(GO_in$classicFisher)][1:10]
  GO_in=GO_in[order(GO_in$classicFisher,decreasing = T)]
  GO_in$Term=factor(GO_in$Term,levels=GO_in$Term)
  print(ggplot(GO_in,aes(x=Term,y=-log10(classicFisher)))+geom_bar(stat="identity",fill=fill)+theme_glob+ coord_flip()+
    theme(legend.position = "",axis.text.y = element_blank(),axis.title.y = element_blank())+
    geom_text(label=GO_in$Term,aes(y=0.1), hjust = 0,size=2))+xlab('-log10(P-value)')
  cat(unique(unlist(lapply(strsplit(GO_in$genes,";"),function(x) x[1:5]))),sep=',')
}
cluster_col<- brewer.pal(10,'Set3')
heartC9=fread('../downstream/output/mm10_result/all_regions/chromHMM_enhancer/cluster_GO/heart-9_cluster_GO.csv')
pdf('../downstream/output/graphs/Figure5/GO_heart_C9.pdf',height = 2,width=2)
plot_GO_example(heartC9,cluster_col[9])
dev.off()
heartC1=fread('../downstream/output/mm10_result/all_regions/chromHMM_enhancer/cluster_GO/heart-1_cluster_GO.csv')
pdf('../downstream/output/graphs/Figure5/GO_heart_C1.pdf',height = 2,width=2)
plot_GO_example(heartC1,cluster_col[1])
dev.off()


# Find top dNME genes for each cluster ------------------------------------
for (fn in dir('../downstream/output/mm10_result/chromHMM_enhancer/cluster_GO/')){
  GO_in=fread(paste0('../downstream/output/mm10_result/chromHMM_enhancer/cluster_GO/',fn))
  sp=sub('_.*','',fn)
  sp=sub("-",'_',sp)
  if(nrow(GO_in)>0){
  gene_GO=unique(unlist(strsplit(GO_in$genes,';')))
  gene_anno=fread(paste0('../downstream/input/mm10_cluster/',sp,'.csv'))
  gene_anno=gene_anno[gene_anno$chromHMM_enhancer&gene_anno$gene%in%gene_GO]
  gene_anno=gene_anno[,list(region=region[which.max(dNME_maxJSD)],dNME_maxJSD=dNME_maxJSD[which.max(dNME_maxJSD)],
                  distance=distance[which.max(dNME_maxJSD)],dMML_maxJSD=dMML_maxJSD[which.max(dNME_maxJSD)],
                  dMML_maxJSD_rank=dMML_maxJSD_rank[which.max(dNME_maxJSD)],
                  dNME_maxJSD_rank=dNME_maxJSD_rank[which.max(dNME_maxJSD)],
                  dNME_maxpair=dNME_maxpair[which.max(dNME_maxJSD)],
                  dMML_maxpair=dMML_maxpair[which.max(dNME_maxJSD)]),by=list(gene)]
  gene_anno=gene_anno[order(gene_anno$dNME_maxJSD,decreasing = T)]
  write.csv(gene_anno,paste0('../downstream/output/mm10_result/chromHMM_enhancer/cluster_gene/',sp,'_gene.csv'))
  }
}
#Double check genes with interesting motif
heartC5_CTCF=fread('../downstream/input/mouse motif cluster/heart_cluster_5_MA0139.1_CTCF_target_gene.csv')
heartC5=fread('../downstream/output/mm10_result/chromHMM_enhancer/cluster_gene/heart_5_gene.csv')
heartC5=heartC5[order(dNME_maxJSD_rank,decreasing = F)]
heartC5_CTCF_motif=heartC5[region %in% paste0(heartC5_CTCF$seqnames,':',heartC5_CTCF$start,'-',heartC5_CTCF$end)]
write.csv(heartC5_CTCF_motif[dNME_maxJSD>=0.3],'../downstream/output/mm10_result/chromHMM_enhancer/heartC5_CTCF.csv')

heartC4=fread('../downstream/output/mm10_result/chromHMM_enhancer/cluster_gene/heart_4_gene.csv')
heartC4=heartC4[order(dNME_maxJSD_rank,decreasing = F)]

heartC4_EGR1=fread('../downstream/input/mouse motif cluster/heart_cluster_4_MA0162.4_EGR1_target_gene.csv')
heartC4_EGR1_motif=heartC4[region %in% paste0(heartC4_EGR1$seqnames,':',heartC4_EGR1$start,'-',heartC4_EGR1$end)]
write.csv(heartC4_EGR1_motif[dNME_maxJSD>=0.3],'../downstream/output/mm10_result/chromHMM_enhancer/heartC4_EGR1.csv')

heartC4_KLF15=fread('../downstream/input/mouse motif cluster/heart_cluster_4_MA1513.1_KLF15_target_gene.csv')
heartC4_KLF15_motif=heartC4[region %in% paste0(heartC4_KLF15$seqnames,':',heartC4_KLF15$start,'-',heartC4_KLF15$end)]
write.csv(heartC4_KLF15_motif[dNME_maxJSD>=0.3],'../downstream/output/mm10_result/chromHMM_enhancer/heartC4_KLF15.csv')

heartC4_KLF3=fread('../downstream/input/mouse motif cluster/heart_cluster_4_MA1516.1_KLF3_target_gene.csv')
heartC4_KLF3_motif=heartC4[region %in% paste0(heartC4_KLF3$seqnames,':',heartC4_KLF3$start,'-',heartC4_KLF3$end)]
write.csv(heartC4_KLF3_motif[dNME_maxJSD>=0.3],'../downstream/output/mm10_result/chromHMM_enhancer/heartC4_KLF3.csv')

heartC4_KLF2=fread('../downstream/input/mouse motif cluster/heart_cluster_4_MA1515.1_KLF2_target_gene.csv')
heartC4_KLF2_motif=heartC4[region %in% paste0(heartC4_KLF2$seqnames,':',heartC4_KLF2$start,'-',heartC4_KLF2$end)]
write.csv(heartC4_KLF3_motif[dNME_maxJSD>=0.3],'../downstream/output/mm10_result/chromHMM_enhancer/heartC4_KLF2.csv')


#Compare old UC and fixed UC
UC_matrix_ls_old=readRDS('../UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')
UC_matrix_ls=readRDS('../UC_run_before_MDS/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
DNase=readRDS('../mm10_DNase.rds')
UC_matrix_ls_old=lapply(UC_matrix_ls_old,function(x) subsetByOverlaps(x,DNase[DNase$region_type=="DNase"],type='equal'))
UC_in_matrix_ls=lapply(UC_in_matrix_ls,function(x) subsetByOverlaps(x,DNase[DNase$region_type=="DNase"],type='equal'))
comp_UC_fix=lapply(names(UC_matrix_ls_old),function(x){
  total_regions=length(UC_matrix_ls_old[[x]])
  UC_matrix_old=as.matrix(mcols(UC_matrix_ls_old[[x]]))
  rownames(UC_matrix_old)=paste0(seqnames(UC_matrix_ls_old[[x]]),':',start(UC_matrix_ls_old[[x]]),'-',end(UC_matrix_ls_old[[x]]))
  cat(x,'\n')
  total_regions_01=sum(rowSums(UC_matrix_old>0.1,na.rm=T)>0)
  UC_matrix_old=UC_matrix_old[rowSums(UC_matrix_old==1,na.rm = T)>0,]
  UC_matrix_fix=as.matrix(mcols(UC_in_matrix_ls[[x]]))
  rownames(UC_matrix_fix)=paste0(seqnames(UC_in_matrix_ls[[x]]),':',start(UC_in_matrix_ls[[x]]),'-',end(UC_in_matrix_ls[[x]]))
  UC_matrix_fix=UC_matrix_fix[rownames(UC_matrix_old),]
  #Number of regions
  diff=rowSums(UC_matrix_old==1,na.rm = T)-rowSums(UC_matrix_fix==1,na.rm = T)
  N_region=sum(diff!=0)
  #Number of sample x regions
  N_sample_region=sum(diff)
  return(list(data.table(tissue=x,N_region=N_region,N_sample_region=N_sample_region,total_regions=total_regions,
                         total_regions_01=total_regions_01),
              regions=rownames(UC_matrix_old)))
})
saveRDS(comp_UC_fix,'../downstream/output/UC_jsd_fix.rds')
#Bug Fix
UC_fix=readRDS('../downstream/output/UC_jsd_fix.rds')
tissue=unlist(lapply(UC_fix,function(x) x[[1]]$tissue))
names(UC_fix)=tissue
tissue=sub('_.*','',dir('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/'))
csv_error=lapply(tissue,function(x) {
  csv_in=fread(paste0('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/',x,'_all.csv'))
  UC_region=UC_fix[[x]][[2]]
  return(csv_in[csv_in$region%in%UC_region])
})
names(csv_error)=tissue
csv_error_dNME=do.call(rbind,lapply(names(csv_error),function(x) return(data.table(csv_error[[x]]$dNME_maxJSD,tissue=x))))
ggplot(csv_error_dNME,aes(x=tissue,y=V1))+geom_boxplot()+theme_glob+ylab('dNME')
lapply(UC_fix,function(x) x[[1]]$N_region/x[[1]]$total_regions_01)

# density plot ------------------------------------------------------------
UC_gr=lapply(UC,function(x) {
  gr=GRanges(seqnames=sub(':.*','',rownames(x)),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',rownames(x)))),
                     end=as.numeric(sub('.*-','',rownames(x)))))
  gr$maxUC=rowMax(x)
  gr$maxUC_loc=apply(x,1,which.max)
  return(gr)
})

dnme_gr=lapply(dnme,function(x) {
  gr=GRanges(seqnames=sub(':.*','',rownames(x)),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',rownames(x)))),
                     end=as.numeric(sub('.*-','',rownames(x)))))
  mcols(gr)=x
  return(gr)
})

dmml_gr=lapply(dmml,function(x) {
  gr=GRanges(seqnames=sub(':.*','',rownames(x)),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',rownames(x)))),
                     end=as.numeric(sub('.*-','',rownames(x)))))
  mcols(gr)=x
  return(gr)
})
for(sp in names(dnme_gr)){
  dnme_gr[[sp]]$maxUC_loc=UC_gr[[sp]]$maxUC_loc
  dnme_gr[[sp]]$dNME_max_UC= apply(mcols(dnme_gr[[sp]]),1,function(x) x[x["maxUC_loc"]])
  dmml_gr[[sp]]$maxUC_loc=UC_gr[[sp]]$maxUC_loc
  dmml_gr[[sp]]$dMML_max_UC= apply(mcols(dmml_gr[[sp]]),1,function(x) x[x["maxUC_loc"]])
}

UC_enhancer=unlist(lapply(UC_gr,function(x) subsetByOverlaps(x,chromHMM)$maxUC))
UC_promoter=unlist(lapply(UC_gr,function(x) subsetByOverlaps(x,promoters)$maxUC))

dnme_enhancer=unlist(lapply(dnme_gr,function(x) subsetByOverlaps(x,chromHMM)$dNME_max_UC))
dnme_promoter=unlist(lapply(dnme_gr,function(x) subsetByOverlaps(x,promoters)$dNME_max_UC))


dmml_enhancer=unlist(lapply(dmml_gr,function(x) subsetByOverlaps(x,chromHMM)$dMML_max_UC))
dmml_promoter=unlist(lapply(dmml_gr,function(x) subsetByOverlaps(x,promoters)$dMML_max_UC))

# plot_density_dt=rbind(data.table(UC=UC_enhancer,region="enhancer"),
#                        data.table(UC=UC_promoter,region="promoter"))
# pdf('../downstream/output/graphs/Figure5/Figure5B_UC_enhancer.pdf',width=3.5,height=3.5)
# print(ggplot(plot_density_dt,aes(x=UC,fill=region))+geom_density(alpha=0.5)+theme_glob+theme(legend.position = "bottom"))
# dev.off()


pdf('../downstream/output/graphs/Figure5/Figure5B_dNME_dMML_density_enhancer.pdf',width=3.5,height=3.5)
plot_density_dt=rbind(data.table(dnme=dnme_enhancer,region="enhancer"),
                      data.table(dnme=dnme_promoter,region="promoter"))
dNME_density=ggplot(plot_density_dt,aes(x=dnme,fill=region))+geom_density(alpha=0.5)+theme_glob+theme(legend.position = "bottom")

plot_density_dt=rbind(data.table(dmml=dmml_enhancer,region="enhancer"),
                      data.table(dmml=dmml_promoter,region="promoter"))

dMML_density=ggplot(plot_density_dt,aes(x=dmml,fill=region))+geom_density(alpha=0.5)+theme_glob+theme(legend.position = "bottom")
ggarrange(dNME_density,dMML_density,nrow=2,ncol=1,common.legend=T)
dev.off()
