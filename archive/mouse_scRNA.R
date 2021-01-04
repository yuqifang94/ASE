rm(list=ls())
source("mainFunctions_sub.R")
library(Gmisc)
library(preprocessCore)
library(pheatmap)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}
corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}
theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
#Read in mouse NME
NME_in=readRDS('../downstream/input/NME_agnostic_mouse_all_merged.rds')
dir='../downstream/data/Mouse_C1/'

# Getting result from ENCODE3 ---------------------------------------------

NME_in=NME_in[NME_in$tissue=="limb"&NME_in$N>=2]
gtf <- fread('../downstream/input/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
genes <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
genes$gene_name <- gn
NME_in=dist_calc(NME_in,genes)
chromHMM=readRDS('../downstream/output/chromHMM_enhancer.rds')
NME_in$chromHMM=FALSE
for(sp in unique(NME_in$Sample)){
  olap=findOverlaps(NME_in[NME_in$Sample==sp],chromHMM[(chromHMM$tissue==unique(NME_in[NME_in$Sample==sp]$tissue))&
                                                         (chromHMM$stage==unique(NME_in[NME_in$Sample==sp]$stage))])
  NME_in$chromHMM[NME_in$Sample==sp][queryHits(olap)]=TRUE
}
#Bin enhancer
enhancer_bin=readRDS("../downstream/output/bin_enhancer.rds")
NME_in$bin_enhancer=NA
olap=findOverlaps(NME_in,enhancer_bin)
NME_in$bin_enhancer[queryHits(olap)]=enhancer_bin$`Target Gene`[subjectHits(olap)]
NME_in_dt=as.data.table(mcols(NME_in))
NME_in_dt$region=paste0(seqnames(NME_in),':',start(NME_in),'-',end(NME_in))
#Check it's correlation near TSS
NME_in_dt$hyper_var=-100
NME_in_dt$var=-100
NME_in_dt$mean=-100

# NME_in_dt$hyper_var_bin=-100
# NME_in_dt$mean_bin=-100
# NME_in_dt$var_bin=-100
for(st in unique(NME_in$stage)){
  tt1=proc.time()[[3]]
  if(file.exists(paste0(dir,sub('E','',st),'.rds'))){
  scRNA_in=readRDS(paste0(dir,sub('E','',st),'.rds'))
  scRNA_in=scRNA_in[rownames(scRNA_in)%in% unique(c(NME_in_dt[(stage==st)]$gene)),]
  if(nrow(scRNA_in)>0){
       #Add hypervar to TSS and chromHMM enhancer
    NME_in_dt[(stage==st)]$hyper_var=scRNA_in[NME_in_dt[(stage==st)]$gene,"hypervar_logvar"]
    NME_in_dt[(stage==st)]$var=scRNA_in[NME_in_dt[(stage==st)]$gene,"var"]
    NME_in_dt[(stage==st)]$mean=scRNA_in[NME_in_dt[(stage==st)]$gene,"mean"]
    # #Add hypervar to bin enhancer
    # 
    # NME_in_dt[(stage==st)]$hyper_var_bin=scRNA_in[NME_in_dt[(stage==st)]$bin_enhancer,"hypervar_logvar"]
    #NME_in_dt[(stage==st)]$var_bin=scRNA_in[NME_in_dt[(stage==st)]$bin_enhancer,"var"]
    #NME_in_dt[(stage==st)]$mean_bin=scRNA_in[NME_in_dt[(stage==st)]$bin_enhancer,"mean"]

  }
  }else{cat("File not exist for ",st,'\n')}
  cat('Finish processing ',sub('E','',st),'in: ',proc.time()[[3]]-tt1,'\n')

}
saveRDS(NME_in_dt,'../downstream/output/NME_in_dt_limb_ENCODE_C1_imputed.rds')
#NME_in_dt=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_C1.rds')
NME_in_dt=NME_in_dt[(!is.na(hyper_var)&hyper_var!=-100)]
#ggplot(NME_in_dt,aes(x=NME,y=hyper_var))+geom_bin2d(bins=10000)

# matrix and quantile normalization ---------------------------------------

matrix_conv<-function(dt_in,value.var){
  out_dc=dcast.data.table(dt_in,region~stage,value.var=value.var)
  rn=out_dc$region
  out_dc=as.matrix(out_dc[,-1])
  rownames(out_dc)=rn
  return(out_dc)
}
NME_dc=matrix_conv(NME_in_dt,"NME")
hyper_var_dc=matrix_conv(NME_in_dt,"hyper_var")

NME_dc=NME_dc[rowSums(is.na(NME_dc))==0,]
hyper_var_dc=hyper_var_dc[rowSums(is.na(hyper_var_dc))==0,]
rn=intersect(rownames(NME_dc),rownames(hyper_var_dc))
NME_dc=NME_dc[rn,]
hyper_var_dc=hyper_var_dc[rn,]
#Test quantile normalization
hyper_var_dc_nm=normalize.quantiles(hyper_var_dc)
rownames(hyper_var_dc_nm)=rownames(hyper_var_dc)
colnames(hyper_var_dc_nm)=colnames(hyper_var_dc)
#After quantile normalization, check plot
hyper_var_dc_nm_dt=data.table(region=rownames(hyper_var_dc_nm))
hyper_var_dc_nm_dt=cbind(hyper_var_dc_nm_dt,as.data.table(hyper_var_dc_nm))
hyper_var_dc_nm_dt=melt.data.table(hyper_var_dc_nm_dt,id.vars="region",variable.name = "stage",value.name = "hyper_var")
# ggplot(hyper_var_dc_nm_dt,aes(x=hyper_var,color=stage))+geom_density(size=1)+theme(legend.position = "bottom")
#assign to orignal values
NME_in_dt$hyper_var=NULL
NME_in_dt$hyper_var=hyper_var_dc_nm_dt[match(paste0(NME_in_dt$region,NME_in_dt$stage),paste0(region,stage))]$hyper_var
NME_in_dt=NME_in_dt[!is.na(hyper_var)]
NME_in_dt_qt=NME_in_dt[,list(quantile_scRNA=findInterval(hyper_var,quantile(hyper_var,prob=c(0,0.25,0.5,0.75)),rightmost.closed=F),
                             quantile_scRNA_100=findInterval(hyper_var,quantile(hyper_var,prob=seq(0,0.99,0.01)),rightmost.closed=F)),
                    by=list(stage)]

NME_in_dt=cbind(NME_in_dt,NME_in_dt_qt[,.(quantile_scRNA,quantile_scRNA_100)])

quantile_vec=c("0%-25%","25%-50%","50%-75%","75%-100%")
NME_in_dt$quantile_scRNA=quantile_vec[NME_in_dt$quantile_scRNA]
saveRDS(NME_in_dt,'../downstream/output/NME_in_dt_limb_ENCODE_C1_nrom.rds')
NME_in_dt=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_C1_nrom.rds')
# Compare 10x and C1 ------------------------------------------------------
cor.test(C1$hypervar_logvar,tenx$hypervar_logvar)
#Normalized
NME_in_dt_C1=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_C1_norm.rds')
NME_in_dt_C1_E135=NME_in_dt_C1[stage=="E13.5"]
NME_in_dt_C1_E135_np5=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_C1_np5_norm.rds')
NME_in_dt_C1_E135_np5=NME_in_dt_C1_E135_np5[stage=="E13.5"]
NME_in_dt_10x_E135_np5=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_10x_norm.rds')
NME_in_dt_10x_E135_np5=NME_in_dt_10x_E135_np5[stage=="E13.5"]
NME_in_dt_10x_E135_p5=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_10x_p5_norm.rds')
NME_in_dt_10x_E135_p5=NME_in_dt_10x_E135_p5[stage=="E13.5"]
shared_region=intersect(NME_in_dt_C1_E135$region,NME_in_dt_C1_E135_np5$region)
shared_region=intersect(shared_region,NME_in_dt_10x_E135_np5$region)
shared_region=intersect(shared_region,NME_in_dt_10x_E135_p5$region)

NME_in_dt_C1_E135=unique(NME_in_dt_C1_E135[region%in%shared_region,list(gene,hyper_var)])
NME_in_dt_C1_E135_np5=unique(NME_in_dt_C1_E135_np5[region%in%shared_region,list(gene,hyper_var)])
NME_in_dt_10x_E135_np5=unique(NME_in_dt_10x_E135_np5[region%in%shared_region,list(gene,hyper_var)])
NME_in_dt_10x_E135_p5=unique(NME_in_dt_10x_E135_p5[region%in%shared_region,list(gene,hyper_var)])
hyper_var_df=data.frame(E13_C1=NME_in_dt_C1_E135[order(gene  )]$hyper_var,
                        E13_5_C1=NME_in_dt_C1_E135_np5[order(gene)]$hyper_var,
                        E13_10x=NME_in_dt_10x_E135_np5[order(gene)]$hyper_var,
                        E13_5_10x=NME_in_dt_10x_E135_p5[order(gene)]$hyper_var)
ggplot(hyper_var_df,aes(x=E13_C1,y=E13_5_C1))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_C1,hyper_var_df$E13_5_C1)
ggplot(hyper_var_df,aes(x=E13_10x,y=E13_5_10x))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_10x,hyper_var_df$E13_5_10x)
ggplot(hyper_var_df,aes(x=E13_C1,y=E13_10x))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_C1,hyper_var_df$E13_10x)
ggplot(hyper_var_df,aes(x=E13_5_C1,y=E13_5_10x))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_5_C1,hyper_var_df$E13_5_10x)
#

E105_10x=readRDS('../downstream/data/10x_ENCODE_limb_point5/10.5.rds')
E155_10x=readRDS('../downstream/data/10x_ENCODE_limb/15.5.rds')
shared_gene=intersect(rownames(E105_10x),rownames(E155_10x))
hyper_var_df=data.table(E105_10x=E105_10x[shared_gene,"hypervar_logvar"],E155_10x=E155_10x[shared_gene,"hypervar_logvar"])
ggplot(hyper_var_df,aes(x=E105_10x,y=E155_10x))+geom_bin2d(bins=100)
ggplot(as.data.table(E105_10x),aes(x=mean,y=var))+geom_bin2d(bins=100)

NME_in_dt_C1_E135=readRDS('../downstream/data/c1_ENCODE_limb/12.5.rds')
NME_in_dt_C1_E135_np5=readRDS('../downstream/data/c1_ENCODE_limb_nop5/13.5.rds')
NME_in_dt_10x_E135_np5=readRDS('../downstream/data/10x_ENCODE_limb/13.5.rds')
NME_in_dt_10x_E135_p5=readRDS('../downstream/data/10x_ENCODE_limb_point5/13.5.rds')
#Filter by 0.1 mean?
NME_in_dt_C1_E135=NME_in_dt_C1_E135[NME_in_dt_C1_E135$mean>=1,]
NME_in_dt_C1_E135_np5=NME_in_dt_C1_E135_np5[NME_in_dt_C1_E135_np5$mean>=1,]
NME_in_dt_10x_E135_np5=NME_in_dt_10x_E135_np5[NME_in_dt_10x_E135_np5$mean>=1,]
NME_in_dt_10x_E135_p5=NME_in_dt_10x_E135_p5[NME_in_dt_10x_E135_p5>=1,]

shared_gene=intersect(rownames(NME_in_dt_C1_E135),rownames(NME_in_dt_C1_E135_np5))
shared_gene=intersect(shared_gene,rownames(NME_in_dt_10x_E135_np5))
shared_gene=intersect(shared_gene,rownames(NME_in_dt_10x_E135_p5))
hyper_var_df=data.frame(E13_C1=NME_in_dt_C1_E135[shared_gene,'hypervar_logvar'],
                        E13_5_C1=NME_in_dt_C1_E135_np5[shared_gene,'hypervar_logvar'],
                        E13_10x=NME_in_dt_10x_E135_np5[shared_gene,'hypervar_logvar'],
                        E13_5_10x=NME_in_dt_10x_E135_p5[shared_gene,'hypervar_logvar'])
mean_df=data.frame(E13_C1=NME_in_dt_C1_E135[shared_gene,'mean'],
                        E13_5_C1=NME_in_dt_C1_E135_np5[shared_gene,'mean'],
                        E13_10x=NME_in_dt_10x_E135_np5[shared_gene,'mean'],
                        E13_5_10x=NME_in_dt_10x_E135_p5[shared_gene,'mean'])
var_df=data.frame(E13_C1=NME_in_dt_C1_E135[shared_gene,'var'],
                   E13_5_C1=NME_in_dt_C1_E135_np5[shared_gene,'var'],
                   E13_10x=NME_in_dt_10x_E135_np5[shared_gene,'var'],
                   E13_5_10x=NME_in_dt_10x_E135_p5[shared_gene,'var'])

ggplot(hyper_var_df,aes(x=E13_C1,y=E13_5_C1))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_C1,hyper_var_df$E13_5_C1)
ggplot(hyper_var_df,aes(x=E13_10x,y=E13_5_10x))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_10x,hyper_var_df$E13_5_10x)
ggplot(hyper_var_df,aes(x=E13_C1,y=E13_10x))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_C1,hyper_var_df$E13_10x)
ggplot(hyper_var_df,aes(x=E13_5_C1,y=E13_5_10x))+geom_bin2d(bins=100)
cor.test(hyper_var_df$E13_5_C1,hyper_var_df$E13_5_10x)
ggplot(mean_df,aes(x=E13_C1,y=E13_10x))+geom_bin2d(bins=100)
cor.test(mean_df$E13_C1,mean_df$E13_10x)
ggplot(mean_df,aes(x=E13_5_C1,y=E13_5_10x))+geom_bin2d(bins=100)
cor.test(mean_df$E13_5_C1,mean_df$E13_5_10x)

ggplot(var_df,aes(x=E13_C1,y=E13_10x))+geom_bin2d(bins=100)
cor.test(var_df$E13_C1,var_df$E13_10x)
ggplot(var_df,aes(x=E13_5_C1,y=E13_5_10x))+geom_bin2d(bins=100)
cor.test(var_df$E13_5_C1,var_df$E13_5_10x)


#Compare old 10x and new one


# Plot line plot and heatmap ----------------------------------------------
pdf('../downstream/output/graphs/Figure6/scRNA_NME_C1_imputed_limb_rmRPLRPS.pdf',width=3.5,height=3.5)
print(ggplot(NME_in_dt[abs(dist)<=3000 &!is.na(quantile_scRNA)],aes(x=dist,y=NME,color=quantile_scRNA ))+geom_smooth(size=1)+xlab('distance to TSS')+
        theme_glob+theme(legend.position = "bottom",legend.title = element_blank()))

dev.off()

pdf('../downstream/output/graphs/Figure6/scRNA_NME_limb_C1_imputed_sample_rmRPLRPS.pdf',width=3.5,height=3.5)
for( sp in unique(NME_in_dt$Sample)){
  print(ggplot(NME_in_dt[abs(dist)<=3000 &!is.na(quantile_scRNA)&Sample==sp],aes(x=dist,y=NME,color=quantile_scRNA ))+geom_smooth(size=1)+xlab('distance to TSS')+
          theme_glob+theme(legend.position = "bottom")+ggtitle(sp))+theme(legend.title = element_blank())
  }

dev.off()

pdf('../downstream/output/graphs/Figure6/scRNA_NME_500_heatmap_limb.pdf',width=7,height=7)
NME_in_dt_reshape=NME_in_dt[dist>0&dist<=500&!is.na(quantile_scRNA),list(NME=median(NME),min_hypervar=min(hyper_var)),by=list(quantile_scRNA_100,Sample)]
print(ggplot(NME_in_dt_reshape,aes(x=quantile_scRNA_100,y=Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
        xlab('Hypervaribility quantile')+ylab('Sample')+theme_glob+theme(legend.position = 'bottom'))
dev.off()

#NME vs hypervaribility
ggplot(NME_in_dt,aes(x=NME,y=hyper_var))+geom_bin2d(bins=100)
ggplot(NME_in_dt[abs(dist)<=3000],aes(x=NME,y=hyper_var))+geom_bin2d(bins=100)
ggplot(NME_in_dt[chromHMM==TRUE],aes(x=NME,y=hyper_var))+geom_bin2d(bins=100)

#calculate correlation between NME vs hypervaribility
# NME_in_dt_mean_region_N=NME_in_dt[chromHMM==TRUE,list(N=length(Sample)),by=list(region)]
# NME_in_dt_cor_region=NME_in_dt[region %in% NME_in_dt_mean_region_N[N==5]$region,
#                                           list(cor=cor(NME,hyper_var),gene=unique(gene)),by=list(region)]
# NME_in_dt_cor=NME_in_dt_cor_region[,list(cor=max(cor)),by=list(gene)]

NME_in_dt_mean_NME=NME_in_dt[chromHMM==TRUE,list(NME=mean(NME),hyper_var=mean(hyper_var)),by=list(gene,Sample)]
NME_in_dt_mean_gene_N=NME_in_dt_mean_NME[,list(N=length(Sample)),by=list(gene)]
NME_in_dt_cor=NME_in_dt_mean_NME[gene %in% NME_in_dt_mean_gene_N[N>=5]$gene,
                                 list(cor=cor(NME,hyper_var)),by=list(gene)]
NME_in_dt_cor=NME_in_dt_cor[order(cor,decreasing=T)]
write(unique(NME_in_dt_cor[cor>=0.5]$gene),
      '../downstream/output/mm10_high_cor_limb_NME_hyperver_4th_quant_C1.txt')
write(unique(NME_in_dt_cor[cor>=0&cor<0.5]$gene),
      '../downstream/output/mm10_high_cor_limb_NME_hyperver_3rd_quant_C1.txt')
write(unique(NME_in_dt_cor[cor<0&cor>=-0.5]$gene),
      '../downstream/output/mm10_high_cor_limb_NME_hyperver_2rd_quant_C1.txt')
write(unique(NME_in_dt_cor[cor<(-0.5)]$gene),
      '../downstream/output/mm10_high_cor_limb_NME_hyperver_1st_quant_C1.txt')
write(unique(NME_in_dt_cor[cor>=quantile(NME_in_dt_cor$cor,prob=c(0.90),na.rm=T)]$gene),
      '../downstream/output/mm10_high_cor_limb_NME_hyperver_90_quant_C1_max.txt')
write(unique(NME_in_dt_cor[cor<=quantile(NME_in_dt_cor$cor,prob=c(0.10),na.rm=T)]$gene),
      '../downstream/output/mm10_high_cor_limb_NME_hyperver_10_quant_C1_max.txt')

# ggplot(NME_in_dt_mean_NME[gene %in% unique(NME_in_dt_cor[cor>=0.9]$gene)],aes(x=NME,y=hyper_var))+geom_point(aes(color=gene))+theme(legend.position = 'none')
# 
# NME_in_dt_gene_N=NME_in_dt[chromHMM==TRUE,list(N=length(Sample)),by=list(gene)]
# NME_in_dt_cor=NME_in_dt[chromHMM==TRUE&gene %in% NME_in_dt_gene_N[N>=5]$gene,
#                         list(cor=cor.test(NME,hyper_var)$estimate,pvalue=cor.test(NME,hyper_var)$p.value),by=list(gene)]
# NME_in_dt_cor$FDR=p.adjust(NME_in_dt_cor$pvalue,method="BH")
# hist(NME_in_dt_cor$cor)
# write(unique(NME_in_dt_cor[FDR<=0.1&cor>0]$gene),
#       '../downstream/output/mm10_high_cor_limb_NME_hyperver_pos_pval_C1.txt')
# write(unique(NME_in_dt_cor[FDR<=0.1&cor<0]$gene),
#       '../downstream/output/mm10_high_cor_limb_NME_hyperver_neg_pval_C1.txt')
#Plot heatmap of high correlated genes
pos_cor=unique(NME_in_dt_cor[cor>=quantile(NME_in_dt_cor$cor,prob=c(0.75),na.rm=T)]$gene)
NME_in_dt_mean_NME_high_cor=NME_in_dt_mean_NME[gene%in%pos_cor]
matrix_conv_gene<-function(dt_in,value.var){
  out_dc=dcast.data.table(dt_in,gene~Sample,value.var=value.var)
  rn=out_dc$gene
  out_dc=as.matrix(out_dc[,-1])
  rownames(out_dc)=rn
  return(out_dc)
}
NME_dc_mean=matrix_conv_gene(NME_in_dt_mean_NME_high_cor,"NME")
hyper_var_dc_mean=matrix_conv_gene(NME_in_dt_mean_NME_high_cor,"hyper_var")
hyper_var_dc_mean=hyper_var_dc_mean[rownames(NME_dc_mean),]
colnames(hyper_var_dc_mean)=paste0("hypervar-",sub('limb-','',sub('-all','',colnames(hyper_var_dc_mean))))
colnames(NME_dc_mean)=paste0("NME-",sub('limb-','',sub('-all','',colnames(NME_dc_mean))))
mean_all=cbind(scalematrix(hyper_var_dc_mean),scalematrix(NME_dc_mean))
set.seed(123456)
clu <- kmeans(mean_all,10,iter.max = 10000)$cluster
clu=clu[names(clu)%in% rownames(mean_all)]
clu=sort(clu)
mean_all=mean_all[names(clu),]
rowann <- data.frame(cluster=as.character(clu),cor=corfunc(hyper_var_dc_mean[names(clu),],NME_dc_mean[names(clu),]),stringsAsFactors = F)
rownames(rowann)=rownames(mean_all)
colann <- data.frame(statics=sub('-.*','',colnames(mean_all)),stringsAsFactors = F)
rownames(colann)=colnames(mean_all)
c1=factor(c("red",'green'))
names(c1)=c("hypervar","NME")
pdf(paste0('../downstream/output/graphs/Figure6/scRNA_hypervar_NME_cluster_chromHMM_gene_limb_C1_cor_075_quant.pdf'),heigh=50,width=10)
print(pheatmap(mean_all,cluster_rows = F,cluster_cols = F,
               annotation_row = rowann,annotation_col = colann,
               show_colnames = T,show_rownames = T,
               gaps_row = cumsum(table(clu)),gaps_col=cumsum(rle(colann[,1])$lengths),
               annotation_colors = list(statics=c1)))
dev.off()
# clustering using NME and hypervar ---------------------------------------
NME_in_dt=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_C1_nrom.rds')
NME_dc=matrix_conv(NME_in_dt,"NME")
hyper_var_dc=matrix_conv(NME_in_dt,"hyper_var")

NME_dc=NME_dc[rowSums(is.na(NME_dc))==0,]
hyper_var_dc=hyper_var_dc[rowSums(is.na(hyper_var_dc))==0,]
rn=intersect(rownames(NME_dc),rownames(hyper_var_dc))
NME_dc=NME_dc[rn,]
hyper_var_dc=hyper_var_dc[rn,]
colnames(hyper_var_dc)=paste0("hypervar-",colnames(hyper_var_dc))
colnames(NME_dc)=paste0("NME-",colnames(NME_dc))
NME_hypervar_mt=cbind(scalematrix(NME_dc),scalematrix(hyper_var_dc))
NME_hypervar_mt=NME_hypervar_mt[rowSums(is.na(NME_hypervar_mt))==0,]
set.seed(123456)
clu <- kmeans(NME_hypervar_mt,10,iter.max = 10000)$cluster
clu=clu[names(clu)%in% rownames(NME_hypervar_mt)]
clu=sort(clu)
NME_hypervar_mt=NME_hypervar_mt[names(clu),]

rowann <- data.frame(cluster=as.character(clu),cor=corfunc(hyper_var_dc[names(clu),],NME_dc[names(clu),]),stringsAsFactors = F)
rownames(rowann)=rownames(NME_hypervar_mt)
colann <- data.frame(statics=sub('-.*','',colnames(NME_hypervar_mt)),stringsAsFactors = F)
rownames(colann)=colnames(NME_hypervar_mt)
c1=factor(c("red",'green'))
names(c1)=c("hypervar","NME")
png(paste0('../downstream/output/graphs/Figure6/scRNA_hypervar_NME_cluster_all_limb_C1.png'))
print(pheatmap(NME_hypervar_mt,cluster_rows = F,cluster_cols = F,
               annotation_row = rowann,annotation_col = colann,
               show_colnames = T,show_rownames = F,
               gaps_row = cumsum(table(clu)),gaps_col=cumsum(rle(colann[,1])$lengths),
               annotation_colors = list(statics=c1)))
dev.off()
#region annotation
chromHMM_regions=unique(NME_in_dt[chromHMM==TRUE]$region)
chromHMM_regions=chromHMM_regions[chromHMM_regions%in% rownames(NME_hypervar_mt)]
TSS_regions=unique(NME_in_dt[abs(dist)<=500]$region)
TSS_regions=TSS_regions[TSS_regions%in% rownames(NME_hypervar_mt)]

png(paste0('../downstream/output/graphs/Figure6/scRNA_hypervar_NME_cluster_all_limb_chromHMM_C1.png'))
print(pheatmap(NME_hypervar_mt[rownames(NME_hypervar_mt)%in%chromHMM_regions,],cluster_rows = F,cluster_cols = F,
               annotation_row = rowann,annotation_col = colann,
               show_colnames = T,show_rownames = F,
               gaps_row = cumsum(table(clu[names(clu)%in% chromHMM_regions])),gaps_col=cumsum(rle(colann[,1])$lengths),
               annotation_colors = list(statics=c1)))
dev.off()


png(paste0('../downstream/output/graphs/Figure6/scRNA_hypervar_NME_cluster_all_limb_TSS_C1.png'))
print(pheatmap(NME_hypervar_mt[rownames(NME_hypervar_mt)%in%TSS_regions,],cluster_rows = F,cluster_cols = F,
               annotation_row = rowann,annotation_col = colann,
               show_colnames = T,show_rownames = F,
               gaps_row = cumsum(table(clu[names(clu)%in% TSS_regions])),gaps_col=cumsum(rle(colann[,1])$lengths),
               annotation_colors = list(statics=c1)))
dev.off()
# change in NME -----------------------------------------------------------
NME_dc=matrix_conv(NME_in_dt,"NME")
hyper_var_dc=matrix_conv(NME_in_dt,"hyper_var")
time_difference=combn(colnames(NME_dc),2)
time_difference=apply(time_difference,2,function(x) paste0(x[1],'-',x[2]))
NME_dc_diff=NULL
hyper_var_dc_diff=NULL
for(diff in time_difference){
  time=strsplit(diff,'-')[[1]]
  NME_diff=as.matrix(abs(NME_dc[,time[1]]-NME_dc[,time[2]]))
  colnames(NME_diff)=diff
  NME_dc_diff=cbind(NME_dc_diff,NME_diff)
  hyper_var_diff=as.matrix(abs(hyper_var_dc[,time[1]]-hyper_var_dc[,time[2]]))
  colnames(hyper_var_diff)=diff
  hyper_var_dc_diff=cbind(hyper_var_dc_diff,hyper_var_diff)
  
}

#filter regions
NME_dc_diff=NME_dc_diff[rowSums(NME_dc_diff>=quantile(NME_dc_diff,prob=0.9))>0,]
hyper_var_dc_diff=hyper_var_dc_diff[rowSums(hyper_var_dc_diff>=quantile(hyper_var_dc_diff,prob=0.9))>0,]
rn=intersect(rownames(NME_dc_diff),rownames(hyper_var_dc_diff))
NME_dc_diff=NME_dc_diff[rn,]
hyper_var_dc_diff=hyper_var_dc_diff[rn,]
#Use only adjacent changes
time_series=paste0("E",10:14,'.5-E',11:15,'.5')
NME_dc_diff=NME_dc_diff[,colnames(NME_dc_diff)%in%time_series]
hyper_var_dc_diff=hyper_var_dc_diff[,colnames(hyper_var_dc_diff)%in%time_series]
colnames(NME_dc_diff) = paste0("NME-",colnames(NME_dc_diff))
colnames(hyper_var_dc_diff) = paste0("hypervar-",colnames(hyper_var_dc_diff))
# cor=corfunc(NME_dc_diff,hyper_var_dc_diff)
# chromHMM_regions=unique(NME_in_dt[chromHMM==TRUE]$region)
# chromHMM_regions=chromHMM_regions[chromHMM_regions%in% rn]
#cor=corfunc(NME_dc_diff[chromHMM_regions,],hyper_var_dc_diff[chromHMM_regions,])
# TSS_regions=unique(NME_in_dt[abs(dist)<=500]$region)
# TSS_regions=TSS_regions[TSS_regions%in% rn]
# cor=corfunc(NME_dc_diff[TSS_regions,],hyper_var_dc_diff[TSS_regions,])


# HOX gene check ----------------------------------------------------------
NME_in=readRDS('../downstream/input/NME_agnostic_mouse_all_merged.rds')
UC_in=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
gtf <- fread('../downstream/input/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
genes <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
genes$gene_name <- gn
NME_in=dist_calc(NME_in[NME_in$N>=2],genes)

MML_in=readRDS('../downstream/input/MML_agnostic_mouse_all_merged.rds')
MML_in=dist_calc(MML_in[MML_in$N>=2],genes)
saveRDS(MML_in,'../downstream/input/MML_agnostic_mouse_all_merged_gene.rds')
NME_in_dt=as.data.table(mcols(NME_in))
NME_in_dt$regions=paste0(seqnames(NME_in),':',start(NME_in),'-',end(NME_in))
MML_in_dt=as.data.table(mcols(MML_in))
MML_in_dt$regions=paste0(seqnames(MML_in),':',start(MML_in),'-',end(MML_in))
saveRDS(MML_in_dt,'../downstream/input/MML_agnostic_mouse_all_merged_gene.rds')
saveRDS(NME_in_dt,'../downstream/input/NME_agnostic_mouse_all_merged_gene.rds')
#Same starting point
Hox_region_NME=NME_in_dt[grep("Hox",gene)]
Hox_region_MML=MML_in_dt[grep("Hox",gene)]
UC_in_dt=mclapply(UC_in,function(x){
  x_out=as.data.table(mcols(x))
  colnames(x_out)=colnames(mcols(x))
  x_out$regions=paste0(seqnames(x),':',start(x),'-',end(x))
  x_out=melt.data.table(x_out,id.var="regions",variable.name="Sample")
  x_out$tissue=gsub('-.*','',x_out$Sample)
  x_out$stage=gsub('-all','',gsub(paste0(unique(x_out$tissue),'-'),'',x_out$Sample))
  return(x_out)
},mc.cores=12)

UC_in_dt=fastDocall('rbind',UC_in_dt)
saveRDS(UC_in_dt,'../downstream/input/UC_agnostic_mouse_all_merged_dt.rds')
colnames(UC_in_dt)[3]="score"
UC_in_dt$stat="UC"
NME_in_dt=NME_in_dt[,.(score,tissue,stage,Sample,dist,gene,regions)]
NME_in_dt$stat="NME"
MML_in_dt=MML_in_dt[,.(score,tissue,stage,Sample,dist,gene,regions)]
MML_in_dt$stat="MML"
UC_in_dt$dist=MML_in_dt$dist[match(UC_in_dt$regions,MML_in_dt$regions)]
UC_in_dt$gene=MML_in_dt$gene[match(UC_in_dt$regions,MML_in_dt$regions)]
UC_in_dt=UC_in_dt[, (match(colnames(MML_in_dt),colnames(UC_in_dt))),with=FALSE]
all_dt=fastDoCall('rbind',list(UC_in_dt,NME_in_dt,MML_in_dt))
all_dt_hox=all_dt[grepl("Hox",gene)]
NME_dt_hox=all_dt_hox[stat=="NME"]
MML_dt_hox=all_dt_hox[stat=="MML"]
UC_dt_hox=all_dt_hox[stat=="UC"]
saveRDS(NME_dt_hox,'../downstream/output/NME_hox_all.rds')
saveRDS(MML_dt_hox,'../downstream/output/MML_hox_all.rds')
saveRDS(UC_dt_hox,'../downstream/output/UC_hox_all.rds')

# Ploting HOX based on regions --------------------------------------------
#NME
plot_heatmap_Hox<-function(Hox_dt_in,N_cutoff=6){
  Hox_dt_in_mean=Hox_dt_in[,list(score=mean(score)),by=list(tissue,stage,gene)]
  Hox_dt_in_mean$gene_tissue=paste0(Hox_dt_in_mean$gene,'-',Hox_dt_in_mean$tissue)
  Hox_dt_in_tissue_N=Hox_dt_in_mean[,list(N=length(stage)),by=list(tissue,gene)]
  Hox_dt_in_mean=Hox_dt_in_mean[tissue%in%unique(Hox_dt_in_tissue_N[N>=N_cutoff]$tissue)]
  HOX_order=c("Hoxa1","Hoxb1","Hoxd1","Hoxa2","Hoxb2","Hoxa3","Hoxb3","Hoxd3","Hoxa4","Hoxb4","Hoxc4",
              "Hoxd4","Hoxa5","Hoxb5","Hoxc5","Hoxa6","Hoxb6","Hoxc6","Hoxa7","Hoxb7","Hoxb8","Hoxc8",
              "Hoxd8","Hoxa9","Hoxb9","Hoxc9","Hoxd9","Hoxa10","Hoxc10","Hoxd10","Hoxa11","Hoxc11",
              "Hoxd11","Hoxc12","Hoxd12","Hoxa13","Hoxb13","Hoxc13","Hoxd13")
  Hox_dt_in_mean_ordered=data.table()
  for(ts in unique(Hox_dt_in_mean$tissue)){
    for(st in unique(Hox_dt_in_mean[tissue==ts]$stage)){
      Hox_dt_in_mean_ordered=rbind(Hox_dt_in_mean_ordered,Hox_dt_in_mean[tissue==ts&stage==st][match(HOX_order,gene)])
    }
  }
  Hox_dt_in_mean_ordered=Hox_dt_in_mean_ordered[!is.na(gene)&!is.na(stage)]
  Hox_dt_in_mean_ordered$gene_tissue=factor(Hox_dt_in_mean_ordered$gene_tissue,levels=unique(Hox_dt_in_mean_ordered$gene_tissue))
  
  print(ggplot(Hox_dt_in_mean_ordered)+geom_tile(aes(x=stage,y=gene_tissue,fill=score))+ scale_fill_distiller(palette = "RdBu"))
  
}
enhancer_subset<-function(dt_hox_in,chromHMM){
  dt_hox_in_enhancer=data.table()
  for(ts in unique(dt_hox_in$tissue)){
    dt_hox_in_ts=dt_hox_in[tissue==ts]
    dt_hox_in_ts_gr=convert_GR(dt_hox_in_ts$regions)
    dt_hox_in_ts_gr=subsetByOverlaps(dt_hox_in_ts_gr,chromHMM[chromHMM$tissue==ts])
    dt_hox_in_ts_gr=paste0(seqnames(dt_hox_in_ts_gr),':',start(dt_hox_in_ts_gr),'-',end(dt_hox_in_ts_gr))
    dt_hox_in_ts=dt_hox_in_ts[regions%in%dt_hox_in_ts_gr]
    dt_hox_in_enhancer=rbind(dt_hox_in_enhancer,dt_hox_in_ts)
  }
  return(dt_hox_in_enhancer)
}
NME_dt_hox=readRDS('../downstream/output/NME_hox_all.rds')
NME_dt_hox=NME_dt_hox[stage!="P0"]
NME_dt_hox_TSS=NME_dt_hox[abs(dist)<=3000]
NME_dt_hox_TSS=NME_dt_hox_TSS[order(gene)]
pdf('../downstream/output/HOX_tissue_TSS_NME.pdf',width=5,height=30)
plot_heatmap_Hox(NME_dt_hox_TSS)
dev.off()

chromHMM=readRDS('../downstream/output/chromHMM_enhancer.rds')
chromHMM=chromHMM[chromHMM$stage!="P0"]
NME_dt_hox_enhancer=enhancer_subset(NME_dt_hox,chromHMM)
pdf('../downstream/output/HOX_tissue_enhancer_NME.pdf',width=5,height=30)
plot_heatmap_Hox(NME_dt_hox_enhancer)
dev.off()

MML_dt_hox=readRDS('../downstream/output/MML_hox_all.rds')
MML_dt_hox=MML_dt_hox[stage!="P0"]
MML_dt_hox_TSS=MML_dt_hox[abs(dist)<=3000]
MML_dt_hox_TSS=MML_dt_hox_TSS[order(gene)]
pdf('../downstream/output/HOX_tissue_TSS_MML.pdf',width=5,height=30)
plot_heatmap_Hox(MML_dt_hox_TSS)
dev.off()

chromHMM=readRDS('../downstream/output/chromHMM_enhancer.rds')
chromHMM=chromHMM[chromHMM$stage!="P0"]
MML_dt_hox_enhancer=enhancer_subset(MML_dt_hox,chromHMM)
pdf('../downstream/output/HOX_tissue_enhancer_MML.pdf',width=5,height=30)
plot_heatmap_Hox(MML_dt_hox_enhancer)
dev.off()

UC_dt_hox=readRDS('../downstream/output/UC_hox_all.rds')
UC_dt_hox=UC_dt_hox[score>=0.1]
same_comp=paste0("E10.5-E",10:16,".5")
UC_dt_hox=UC_dt_hox[stage%in%same_comp]
UC_dt_hox_TSS=UC_dt_hox[abs(dist)<=3000]
UC_dt_hox_TSS=UC_dt_hox_TSS[order(gene)]
pdf('../downstream/output/HOX_tissue_TSS_UC.pdf',width=5,height=30)
plot_heatmap_Hox(UC_dt_hox_TSS,N_cutoff=4)+theme(axis.text.x = element_text(angle = 90))
dev.off()

chromHMM=readRDS('../downstream/output/chromHMM_enhancer.rds')
chromHMM=chromHMM[chromHMM$stage!="P0"]
UC_dt_hox_enhancer=enhancer_subset(UC_dt_hox,chromHMM)
pdf('../downstream/output/HOX_tissue_enhancer_UC.pdf',width=5,height=30)
plot_heatmap_Hox(UC_dt_hox_enhancer,N_cutoff=4)+theme(axis.text.x = element_text(angle = 90))
dev.off()




UC_in_mt=UC_in_mt[,colnames(UC_in_mt)%in%time_difference]
UC_in_mt=UC_in_mt[rowSums(UC_in_mt>=0.1)>0,]


UC_in_mt_Hox=as.data.table(UC_in_mt[rownames(UC_in_mt)%in%Hox_region$region,])
UC_in_mt_Hox$region=rownames(UC_in_mt)[rownames(UC_in_mt)%in%Hox_region$region]
UC_in_mt_Hox$gene=Hox_region$gene[match(UC_in_mt_Hox$region,Hox_region$region)]
UC_in_mt_Hox= melt.data.table(UC_in_mt_Hox,id.vars =c("region","gene"),value.name="UC",variable.name = "stage")
UC_in_mt_Hox=UC_in_mt_Hox[region%in%TSS_regions]
UC_in_mt_Hox_mean=UC_in_mt_Hox[,list(UC=mean(UC)),by=list(gene,stage)]
UC_in_mt_Hox_mean$gene=factor(UC_in_mt_Hox_mean$gene)
UC_in_mt_Hox_mean_max=data.table()
for(gene_in in unique(UC_in_mt_Hox_mean$gene)){
  
  UC_in_mt_Hox_mean_max=rbind(UC_in_mt_Hox_mean_max,
                              data.table(gene=gene_in,UC_max=max(UC_in_mt_Hox[gene==gene_in]$UC),
                                         stage_max=UC_in_mt_Hox[gene==gene_in]$stage[which.max(UC_in_mt_Hox[gene==gene_in]$UC)]))
  
}
UC_in_mt_Hox_mean_max=UC_in_mt_Hox_mean_max[order(stage_max)]
ggplot(UC_in_mt_Hox_mean,aes(x=stage,y=UC,color=gene,group=gene))+geom_point()+geom_line()
# clustering --------------------------------------------------------------
#adding UC
NME_hyper_diff_mt=cbind(scalematrix(NME_dc_diff),scalematrix(hyper_var_dc_diff))
NME_hyper_diff_mt=NME_hyper_diff_mt[rowSums(is.na(NME_hyper_diff_mt))==0,]
UC_in=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
UC_in=UC_in$limb
UC_in_mt=as.matrix(mcols(UC_in))
rownames(UC_in_mt)=paste0(seqnames(UC_in),':',start(UC_in),'-',end(UC_in))
colnames(UC_in_mt)=sub('-all','',sub('limb-','',colnames(UC_in_mt)))
UC_in_mt=UC_in_mt[rownames(NME_hyper_diff_mt),colnames(UC_in_mt)%in%time_series]
colnames(UC_in_mt)=paste0('UC-',colnames(UC_in_mt))
NME_hyper_diff_mt=cbind(NME_hyper_diff_mt,scalematrix(UC_in_mt))
# #Jason cluster
cluster=readRDS('../downstream/input/jsd.rds')
clu=cluster$limb
set.seed(123456)
#clu <- kmeans(NME_hyper_diff_mt[,sub('-.*','',colnames(NME_hyper_diff_mt))%in%c('NME','hypervar')],5,iter.max = 10000)$cluster


clu=clu[names(clu)%in% rownames(NME_hyper_diff_mt)]
clu=sort(clu)
NME_hyper_diff_mt=NME_hyper_diff_mt[names(clu),]

rowann <- data.frame(cluster=as.character(clu),cor=corfunc(NME_dc_diff[names(clu),],hyper_var_dc_diff[names(clu),]),stringsAsFactors = F)
rownames(rowann)=rownames(NME_hyper_diff_mt)
colann <- data.frame(statics=sub('-.*','',colnames(NME_hyper_diff_mt)),stringsAsFactors = F)
rownames(colann)=colnames(NME_hyper_diff_mt)
c1=factor(c("red","blue",'green'))
names(c1)=c("hypervar","NME",'UC')

png(paste0('../downstream/output/graphs/Figure6/scRNA_hypervar_cluster_all_diff_limb_JS_10x.png'))
print(pheatmap(NME_hyper_diff_mt,cluster_rows = F,cluster_cols = F,
               annotation_row = rowann,annotation_col = colann,
               show_colnames = T,show_rownames = F,
               gaps_row = cumsum(table(clu)),gaps_col=cumsum(rle(colann[,1])$lengths),
               annotation_colors = list(statics=c1)))
dev.off()

#region annotation
chromHMM_regions=unique(NME_in_dt[chromHMM==TRUE]$region)
chromHMM_regions=chromHMM_regions[chromHMM_regions%in% rownames(NME_hyper_diff_mt)]
TSS_regions=unique(NME_in_dt[abs(dist)<=500]$region)
TSS_regions=TSS_regions[TSS_regions%in% rownames(NME_hyper_diff_mt)]

png(paste0('../downstream/output/graphs/Figure6/scRNA_hypervar_cluster_all_diff_chromHMM_limb_10x_JS.png'))
print(pheatmap(NME_hyper_diff_mt[rownames(NME_hyper_diff_mt)%in%chromHMM_regions,],cluster_rows = F,cluster_cols = F,
               annotation_row = rowann,annotation_col = colann,
               show_colnames = T,show_rownames = F,
               gaps_row = cumsum(table(clu[names(clu)%in% chromHMM_regions])),gaps_col=cumsum(rle(colann[,1])$lengths),
               annotation_colors = list(statics=c1)))
dev.off()

png(paste0('../downstream/output/graphs/Figure6/scRNA_hypervar_cluster_all_diff_TSS_limb_10x_JS.png'))
print(pheatmap(NME_hyper_diff_mt[rownames(NME_hyper_diff_mt)%in%TSS_regions,],cluster_rows = F,cluster_cols = F,
               annotation_row = rowann,annotation_col = colann,
               show_colnames = T,show_rownames = F,
               gaps_row = cumsum(table(clu[names(clu)%in% TSS_regions])),gaps_col=cumsum(rle(colann[,1])$lengths),
               annotation_colors = list(statics=c1)))
dev.off()

# pdf('../downstream/output/graphs/Figure6/scRNA_NME_heatmap_limb.pdf',width=7,height=7)
# NME_in_dt_reshape=NME_in_dt[abs(dist)<=3000&!is.na(quantile_scRNA),list(NME=median(NME),min_hypervar=min(hyper_var)),by=list(quantile_scRNA_100,Sample)]
# print(ggplot(NME_in_dt_reshape,aes(x=quantile_scRNA_100,y=Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
#         xlab('Hypervaribility quantile')+ylab('Sample')+theme_glob+theme(legend.position = 'bottom'))
# dev.off()
# 
# 
# pdf('../downstream/output/graphs/Figure6/scRNA_NME_chromHMM_heatmap_limb.pdf',width=7,height=7)
# NME_in_dt_reshape=NME_in_dt[chromHMM&!is.na(quantile_scRNA),list(NME=median(NME),min_hypervar=min(hyper_var)),by=list(quantile_scRNA_100,Sample)]
# print(ggplot(NME_in_dt_reshape,aes(x=quantile_scRNA_100,y=Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
#         xlab('Hypervaribility quantile')+ylab('Sample')+theme_glob+theme(legend.position = 'bottom'))
# dev.off()
# 
# pdf('../downstream/output/graphs/Figure6/scRNA_NME_Bin_heatmap.pdf',width=7,height=7)
# print(ggplot(NME_in_dt[!is.na(bin_enhancer)],aes(x=quantile_scRNA_100_bin,y=Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
#         xlab('Hypervaribility quantile')+ylab('Sample')+theme_glob+theme(legend.position = 'bottom'))
# dev.off()

#Read scRNA data
Convert("200120_10x.h5ad", dest = "h5seurat", overwrite = TRUE)
LoadH5Seurat('200120_10x.h5seurat')
#Read scRNA data
Convert("200315_C1_categorical.h5ad", dest = "h5seurat", overwrite = TRUE)
C1=LoadH5Seurat('200315_C1_categorical.h5seurat')
saveRDS(C1,"C1_all.rds")

