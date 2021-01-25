rm(list=ls())
source('mainFunctions_sub.R')
tissue='heart'
stage1="day10_5"
stage2="day16_5"
dir_in='../downstream/data/rep_analysis/'
DNase=readRDS('../downstream/input/DNase_mm10_peak_merge_250bp.rds')
DNase=convert_GR(DNase,direction="DT")
UC_betweeen=c(
  read.agnostic.mouse.uc(paste0(dir_in,'mm10_',tissue,'_',stage1,'_merged1-vs-mm10_',tissue,'_',stage2,'_merged1_uc.bedGraph')),
  read.agnostic.mouse.uc(paste0(dir_in,'mm10_',tissue,'_',stage1,'_merged2-vs-mm10_',tissue,'_',stage2,'_merged2_uc.bedGraph'))
  
)
UC_betweeen$comparison=paste0('between_sample',UC_betweeen$replicate)
UC_betweeen=convert_GR(UC_betweeen,direction="DT")
UC_within=c(
    read.agnostic.mouse.uc(paste0(dir_in,'mm10_',tissue,'_',stage1,'_merged1-vs-mm10_',tissue,'_',stage1,'_merged2_uc.bedGraph')),
    read.agnostic.mouse.uc(paste0(dir_in,'mm10_',tissue,'_',stage2,'_merged1-vs-mm10_',tissue,'_',stage2,'_merged2_uc.bedGraph'))
)
UC_within$replicate=match(UC_within$Sample,unique(UC_within$Sample))
UC_within$comparison=paste0('within_sample',UC_within$replicate)
UC_within=convert_GR(UC_within,direction="DT")

UC_all=rbind(UC_betweeen,UC_within)
UC_all=UC_all[N>=2]
UC_all=dcast.data.table(UC_all,region~comparison,value.var = 'score')
UC_all=UC_all[region%in%DNase$region]
UC_all=UC_all[rowSums(is.na(UC_all[,-1]))==0]

# sanity check ------------------------------------------------------------
UC_all$between_mean=(UC_all$between_sample1+UC_all$between_sample2)/2
UC_all$within_mean=(UC_all$within_sample1+UC_all$within_sample2)/2
#Look for points moved from upper cornor to lower corner
ggplot(UC_all,aes(x=between_mean,y=within_mean))+geom_bin2d(bins=100)+theme_glob+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
    scale_fill_gradient2(high = "darkblue",low="lightblue")


# ggplot(UC_all,aes(x=between_sample1,y=within_sample1))+geom_bin2d(bins=100)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
#   #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
#   geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
#   geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")
# ggplot(UC_all,aes(x=between_sample2,y=within_sample2))+geom_bin2d(bins=100)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
#   #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
#   geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
#   geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")
# ggplot(UC_all,aes(x=between_mean,y=within_mean ))+geom_bin2d(bins=100)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
#   #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
#   geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
#   geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

# ggplot(UC_all,aes(x=between_mean,y=within_mean ))+geom_bin2d(bins=100)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
#   #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
#   geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
#   geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")
#   scale_fill_gradient2(high = "darkblue",low="lightblue")
  
#Filtered
jpeg('sanity_check.jpg')
ggplot(UC_all_10_11,aes(x=within_sample1,y=within_sample2 ))+geom_point(alpha=0.01)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

dev.off()
jpeg('sanity_check2.jpg')
ggplot(UC_all,aes(x=between_sample1 ,y=within_sample1   ))+geom_point(alpha=0.01)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

dev.off()
jpeg('sanity_check3.jpg')
ggplot(UC_all,aes(x=between_sample2 ,y=within_sample2   ))+geom_point(alpha=0.01)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

dev.off()

jpeg('sanity_check_high_cutoff.jpg')
ggplot(UC_all[between_mean>=0.1],aes(x=between_mean,y=within_mean ))+geom_point(alpha=0.01)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

dev.off()
#Find regions above fitted line
UC_all$within_sample1_est1=0.6-UC_all$between_sample1
UC_all$within_sample1_est2=UC_all$between_sample1+0.05
UC_all$within_sample1_est3=UC_all$between_sample1-0.05
UC_all_sub1=UC_all[((within_sample1>within_sample1_est1) & (within_sample1>within_sample1_est2))|
         ((within_sample1>within_sample1_est1) & (within_sample1<within_sample1_est3))]

ggplot(UC_all_sub1,aes(x=between_sample1,y=within_sample1))+geom_point(alpha=0.5)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
 geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

UC_all$within_sample2_est1=0.6-UC_all$between_sample2
UC_all$within_sample2_est2=UC_all$between_sample2+0.05
UC_all$within_sample2_est3=UC_all$between_sample2-0.05
UC_all_sub2=UC_all[((within_sample2>within_sample2_est1) & (within_sample2>within_sample2_est2))|
                    ((within_sample2>within_sample2_est1) & (within_sample2<within_sample2_est3))]

ggplot(UC_all_sub2,aes(x=between_sample2,y=within_sample2))+geom_point(alpha=0.5)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

ggplot(UC_all_sub2,aes(x=between_sample2,y=within_sample2))+geom_point(alpha=0.5)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
 geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")


ggplot(UC_all_sub2,aes(x=between_sample1,y=within_sample1))+geom_point(alpha=0.5)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

ggplot(UC_all_sub2,aes(x=between_mean,y=within_mean))+geom_point(alpha=0.5)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+
  #geom_abline(slope=-1,intercept = 0.5,color='red',size=1)+geom_abline(slope=1,intercept = 0,color='red',size=1)+
  geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")

jpeg('sanity_check4.jpg')
ggplot(UC_all[!region %in% c(UC_all_sub1$region,UC_all_sub2$region)],aes(x=between_mean ,y=within_mean))+geom_point(alpha=0.01)+
  theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")+ geom_abline(slope=-1,intercept = 0.6,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = 0.05,color='red',size=1,linetype = "dashed")+
  geom_abline(slope=1,intercept = -0.05,color='red',size=1,linetype = "dashed")
dev.off()


UC_all$within_mean_est1=0.6-UC_all$between_mean
UC_all$within_mean_est2=UC_all$between_mean+0.05
UC_all$within_mean_est3=UC_all$between_mean-0.05
UC_all_exclude_region=UC_all[!region %in% c(UC_all_sub1$region,UC_all_sub2$region)]
UC_all_exclude_region[within_mean>within_mean_est3&within_mean<=within_mean_est2&within_mean>within_mean_est1]
# compare between mean UC and merged UC -----------------------------------

UC_raw=readRDS('../downstream/output/uc_matrix_DNase.rds')
UC_raw=UC_raw[[tissue]]
UC_raw=UC_raw[rownames(UC_raw)%in% UC_all$region,"E10.5-E11.5"]

UC_all$merged_UC=UC_raw[UC_all$region]
ggplot(UC_all[between_mean>=0.05&merged_UC>=0.05],aes(x=between_mean,y=merged_UC))+geom_bin2d(bins=100)+theme_glob+xlim(c(0,0.6))+ylim(c(0,0.6))+theme(legend.position = "bottom")

#Subset the ones with all data

UC_all$between_mean=(UC_all$between_sample1+UC_all$between_sample2)/2

UC_all$within_mean=(UC_all$within_sample1+UC_all$within_sample2)/2

ggplot(UC_all,aes(x=within_mean,y=between_mean))+geom_bin2d(bins=100)

blue_tri=UC_all[between_mean<=within_mean]

#+0.001 to avoid -INF?
UC_all_lg=as.matrix(log(UC_all[,2:5]+0.001))
rownames(UC_all_lg)=UC_all$region
#UC_all_lg_norm=normalizeVSN(UC_all_lg)
saveRDS(UC_all_lg,'../downstream/output/UC_all_lg_heart.rds')
#lmFit
design=data.frame(within=as.numeric(grepl('within',colnames(UC_all_lg))),
                  between=as.numeric(grepl('between',colnames(UC_all_lg))))
design=data.frame(between_within=design$between-design$within)
rownames(design)=colnames(UC_all_lg)
# fit <- lmFit(UC_all_lg, design,weights=NULL)
# fit2 <- eBayes(fit)
# topTable(fit2, coef=ncol(design))
# meanSdPlot(UC_all_lg_norm)
# meanSdPlot(UC_all_lg)

#Use dge list?
# 
# UC_all_expand=round(as.matrix(UC_all[,2:5]*10^6))
# rownames(UC_all_expand)=UC_all$region
# 
# dge=DGEList(counts=UC_all_expand)
# dge <- calcNormFactors(dge)
# logCPM <- cpm(dge, log=TRUE, prior.count=3)
UC_all_lg=UC_all[,-1][,lapply(.SD,qnorm)]
UC_all_lg=as.matrix(UC_all_lg)
rownames(UC_all_lg)=UC_all$region
#UC_all_lg=UC_all_lg[rownames(UC_all_lg) %in% GR_out_dc_cor[FDR<=0.25]$region,]
fit <- lmFit(UC_all_lg, design)
 fit <- eBayes(fit, trend=F)
topTable(fit, coef=ncol(design),adjust="BH")
UC_all_lg_dt=as.data.table(UC_all_lg)
cor.test(UC_all_lg_dt$between_sample1,UC_all_lg_dt$within_sample1)
jpeg("sanity_check_5.jpg")
ggplot(UC_all_lg_dt,aes(x=between_sample1,y=within_sample1))+geom_point(alpha=0.01)
dev.off()
#Correlation of multiple UC between two replicates
dir_in='../downstream/input/heart_test/'
GR_out=GRanges()
for(fn in dir(dir_in)){
  GR_out=c(GR_out,read.agnostic.mouse.uc(paste0(dir_in,fn)))
}
GR_out=GR_out[GR_out$N>=2]
GR_out=convert_GR(GR_out,direction="DT")
GR_out$tissue_stage=paste0(GR_out$tissue,'_',GR_out$stage)
GR_out$replicate=paste0("rep",GR_out$replicate)
GR_out_dc=dcast.data.table(GR_out,region+tissue_stage~replicate,value.var="score")#30395865
GR_out_dc=GR_out_dc[rowSums(is.na(GR_out_dc[,list(rep1,rep2)]))==0]#29651477
GR_out_dc=GR_out_dc[region%in% names(table(GR_out_dc$region))[table(GR_out_dc$region)>=4]]#29615236
GR_out_dc_cor=GR_out_dc[,list(cor=cor(rep1,rep2),p_val=cor.test(rep1,rep2)$p.value),by=list(region)]#1442515
GR_out_dc_cor=GR_out_dc_cor[!is.na(GR_out_dc_cor$p_val)]#1442258
GR_out_dc_cor$FDR=p.adjust(GR_out_dc_cor$p_val,method='BH')#70189 significant
#
UC_raw=readRDS('../downstream/output/uc_matrix_DNase.rds')
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
UC_raw=UC_raw[[tissue]]
UC_raw=UC_raw[rownames(UC_raw)%in%GR_out_dc_cor$region,]
UC_raw=UC_raw[rowSums(UC_raw>=0.1)>0,]
#14770 out of (167389 covered in GR_out_dc)184896 in it, 8%
#Total significant: 70189 FDR<=0.1, ~38% of number of previously significant 
#Plot 2 example
low_correlation=GR_out_dc[region=='chr10:100004501-100004750']
high_correlation=GR_out_dc[region=='chr10:100014215-100014464']
high_correlation=GR_out_dc[region=='chr10:10010257-10010506']
GR_out_dc$mean_UC=(GR_out_dc$rep1+GR_out_dc$rep2)/2
GR_out_dc_mt=dcast.data.table(GR_out_dc,region~tissue_stage,value.var = 'mean_UC')
hist(rowMaxs(as.matrix(GR_out_dc_mt[region%in%GR_out_dc_cor[FDR<=0.1]$region,-1],na.rm=T)))
ggplot(low_correlation,aes(x=rep1,y=rep2))+geom_point()+geom_abline(slope=1,size=1,color='blue')
ggplot(high_correlation,aes(x=rep1,y=rep2))+geom_point()+geom_abline(slope=1,size=1,color='blue')

#2d kernal, choose bandwidth -> cross validation: estimated density to explain leftout data, what's likelihood of generating the leftout data
#choose bandwidth max likelihood.

#, add psedocunt -> smoothing the matrix -> FDR cutoff-> how many observation pass cutoff
UC_all_10_16=UC_all
UC_all_10_16$within_mean=(UC_all_10_16$within_sample1+UC_all_10_16$within_sample2)/2
UC_all_10_16$between_mean=(UC_all_10_16$between_sample1+UC_all_10_16$between_sample2)/2

# Method 1: estimate upper and lower triangle density separately ----------

expected_UC=UC_all_10_16[within_mean>=between_mean]
observed_UC=UC_all_10_16[within_mean<between_mean]
jpeg('../downstream/output/raw_data.jpg')
ggplot(UC_all_10_16,aes(x=between_mean,y=within_mean))+geom_point(alpha=0.01)+xlim(c(0,0.6))+ylim(c(0,0.6))+geom_abline(slope=1,size=1,color='red')
dev.off()
jpeg('../downstream/output/density_est_method1_expected.jpg')
ggplot(expected_UC,aes(x=between_mean,y=within_mean))+geom_point(alpha=0.01)+xlim(c(0,0.6))+ylim(c(0,0.6))
dev.off()
jpeg('../downstream/output/density_est_method1_observed.jpg')
ggplot(observed_UC,aes(x=between_mean,y=within_mean))+geom_point(alpha=0.01)+xlim(c(0,0.6))+ylim(c(0,0.6))
dev.off()
jpeg('../downstream/output/density_est_method1_expected_rotated.jpg')
ggplot(expected_UC,aes(y=between_mean,x=within_mean))+geom_point(alpha=0.01)+xlim(c(0,0.6))+ylim(c(0,0.6))
dev.off()
expected_UC_mt=as.matrix(data.frame(var1=expected_UC$within_mean,var2=expected_UC$between_mean))
observed_UC_mt=as.matrix(data.frame(var1=observed_UC$between_mean,var2=observed_UC$within_mean))
pse_count=0.0001
#1. Z value is an estimated pdf which is not equivalent to cdf, cdf should sum to 1, but pdf can be larger than one
# cdf = sum(pdf*bandwidth), bandwidth can be small
#Local FDR = expected/observed: http://statweb.stanford.edu/~ckirby/brad/papers/2005LocalFDR.pdf
#http://statweb.stanford.edu/~ckirby/brad/LSI/chapter5.pdf
library(MASS)
# UC_all_10_16$uppder_filter=UC_all_10_16$between_mean+0.01
# UC_all_10_16$lower_filter=UC_all_10_16$between_mean-0.01
#try ks package with Hscv
#https://cran.r-project.org/web/packages/ks/vignettes/kde.pdf
library(ks)
#Cross validation
fdr_smooth<-function(fdr){
  #add pseducount
  fdr[fdr>1]=1
  library(oce)
  return(matrixSmooth(fdr,passes=100))
}
Hpi1<- Hscv(x=expected_UC_mt)
# fhat.exp_pse<- kde(x=expected_UC_mt, H=Hpi1)
# #adding psedocount?
# pse_countn=1
# pse_count=NULL
# for(i in 1:pse_countn){
#   pse_count=rbind(pse_count,as.matrix(expand.grid(fhat.exp_pse$eval.points[[1]],fhat.exp_pse$eval.points[[2]])))
#
# }
# expected_UC_mt_pse=rbind(expected_UC_mt,pse_count)
# observed_UC_mt_pse=rbind(observed_UC_mt,pse_count)
fhat.exp<- kde(x=expected_UC_mt, H=Hpi1)
fhat.obs<- kde(x=observed_UC_mt, H=t(Hpi1))
#check smooth function
library(plotly)
jpeg('../downstream/output/method_1_expected_ks_estimate.jpg')
filled.contour(x = fhat.exp$eval.points[[1]], y = fhat.exp$eval.points[[2]], z = fhat.exp$estimate) 
dev.off()
jpeg('../downstream/output/method_1_expected_ks_observed.jpg')
filled.contour(x = fhat.obs$eval.points[[1]], y = fhat.obs$eval.points[[2]], z = fhat.obs$estimate) 
dev.off()

fhat.obs$local.fdr=(fhat.exp$estimate+pse_count)/(fhat.obs$estimate+pse_count)
fhat.obs$local.fdr[ fhat.obs$local.fdr>1]=1
jpeg('../downstream/output/method_1_local_fdr_unsmoothed.jpg')
filled.contour(x = fhat.obs$eval.points[[1]], y = fhat.obs$eval.points[[2]], z = fhat.obs$local.fdr) 
dev.off()
fhat.obs$local.fdr_sm=fdr_smooth(fhat.obs$local.fdr)
jpeg('../downstream/output/method_1_local_fdr_smoothed_100.jpg')
filled.contour(x = fhat.obs$eval.points[[1]], y = fhat.obs$eval.points[[2]], z = fhat.obs$local.fdr_sm)
dev.off()
observed_UC$x_loc=unlist(lapply(observed_UC$between_mean,function(x) which.min(abs(fhat.obs$eval.points[[1]]-x))))
observed_UC$y_loc=unlist(lapply(observed_UC$within_mean,function(x)which.min(abs(fhat.obs$eval.points[[2]]-x))))
observed_UC$loc.fdr=unlist(lapply(1:nrow(observed_UC),function(x) return(fhat.obs$local.fdr[observed_UC$x_loc[x], observed_UC$y_loc[x]])))

jpeg('../downstream/output/method_1_local_fdr_unsmoothed_01.jpg')

ggplot(observed_UC[loc.fdr<=0.1],aes(x=between_mean,y=within_mean))+geom_point(alpha=0.1)+xlim(c(0,0.6))+ylim(c(0,0.6))

dev.off()
#K-fold CV for bandwidth using ucv:
var1_bw=ucv(expected_UC_mt[,'var1'])
var2_bw=ucv(expected_UC_mt[,'var2'])
saveRDS(var1_bw,'../downstream/var1_bw_m1.rds')
saveRDS(var2_bw,'../downstream/var2_bw_m1.rds')
fhat.exp=kde2d(expected_UC_mt[,'var1'], expected_UC_mt[,'var2'],h=c(var1_bw,var2_bw),lims=c(range(expected_UC_mt[,'var1']),range(expected_UC_mt[,'var2'])))
fhat.obs=kde2d(observed_UC_mt[,'var1'], observed_UC_mt[,'var2'],h=c(var2_bw,var1_bw),lims=c(range(expected_UC_mt[,'var2']),range(expected_UC_mt[,'var1'])))
fhat.obs$local.fdr=(fhat.exp$z+pse_count)/(fhat.obs$z+pse_count)
fhat.obs$local.fdr_sm=fdr_smooth(fhat.obs$local.fdr)

jpeg('../downstream/output/method_1_kde2d_estimate.jpg')
filled.contour(x = fhat.exp$x, y = fhat.exp$y, z = fhat.exp$z) 
dev.off()
jpeg('../downstream/output/method_1_kde2d_observed.jpg')
filled.contour(x = fhat.obs$x, y = fhat.obs$y, z = fhat.obs$z) 
dev.off()

fhat.obs$local.fdr[ fhat.obs$local.fdr>1]=1
jpeg('../downstream/output/method_1_local_fdr_unsmoothed_kde2d.jpg')
filled.contour(x = fhat.obs$x, y = fhat.obs$y, z = fhat.obs$local.fdr) 
dev.off()
fhat.obs$local.fdr_sm=fdr_smooth(fhat.obs$local.fdr)
jpeg('../downstream/output/method_1_local_fdr_smoothed_100_kde2d.jpg')
filled.contour(x = fhat.obs$x, y = fhat.obs$y, z = fhat.obs$local.fdr_sm)
dev.off()

observed_UC$x_loc=unlist(lapply(observed_UC$between_mean,function(x) which.min(abs(fhat.obs$x-x))))
observed_UC$y_loc=unlist(lapply(observed_UC$within_mean,function(x)which.min(abs(fhat.obs$y-x))))
observed_UC$loc.fdr=unlist(lapply(1:nrow(observed_UC),function(x) return(fhat.obs$local.fdr_sm[observed_UC$x_loc[x], observed_UC$y_loc[x]])))
jpeg('../downstream/output/method_1_local_fdr_smoothed_100_01_kde2d.jpg')

ggplot(observed_UC[loc.fdr<=0.25],aes(x=between_mean,y=within_mean))+geom_point(alpha=0.1)+xlim(c(0,0.6))+ylim(c(0,0.6))

dev.off()
# Method 2: estimate upper and lower triangle density together and use upper as null ----------
UC_mt=as.matrix(data.frame(var1=UC_all_10_16$within_mean,var2=UC_all_10_16$between_mean))
Hpi1<- Hscv(x=UC_mt)
fhat.all<- kde(x=UC_mt, H=Hpi1)

jpeg('../downstream/output/method_2_raw_density.jpg')
filled.contour(x = fhat.all$eval.points[[1]], y = fhat.all$eval.points[[2]], z = fhat.all$estimate) 
dev.off()
jpeg('../downstream/output/method_2_raw_density_rotate.jpg')
filled.contour(x = fhat.all$eval.points[[1]], y = fhat.all$eval.points[[2]], z = fhat.all$estimate) 
dev.off()
jpeg('../downstream/output/method_2_raw_density_expected_obs.jpg')
filled.contour(x = fhat.all$eval.points[[1]], y = fhat.all$eval.points[[2]], z = (fhat.all$estimate+pse_count)/(t(fhat.all$estimate)+pse_count)) 
dev.off()
fhat.all$loc.fdr=(fhat.all$estimate+pse_count)/(t(fhat.all$estimate)+pse_count)
fhat.all$loc.fdr[fhat.all$loc.fdr>1]=1
fhat.all$loc.fdr=(fhat.all$estimate+pse_count)/(t(fhat.all$estimate)+pse_count)
jpeg('../downstream/output/method_2_raw_loc_fdr_unsmoothed.jpg')
filled.contour(x = fhat.all$eval.points[[1]], y = fhat.all$eval.points[[2]], z = fhat.all$loc.fdr) 
dev.off()
UC_all_10_16$x_loc=unlist(lapply(UC_all_10_16$between_mean,function(x) which.min(abs(fhat.all$x-x))))
UC_all_10_16$y_loc=unlist(lapply(UC_all_10_16$within_mean,function(x)which.min(abs(fhat.all$y-x))))
UC_all_10_16$loc.fdr=unlist(lapply(1:nrow(UC_all_10_16),function(x) return(fhat.all$local.fdr_sm[UC_all_10_16$x_loc[x], UC_all_10_16$y_loc[x]])))

jpeg('../downstream/output/method_2_local_fdr_smoothed_100_01_ks.jpg')

ggplot(observed_UC[loc.fdr<=0.25],aes(x=between_mean,y=within_mean))+geom_point(alpha=0.1)+xlim(c(0,0.6))+ylim(c(0,0.6))

dev.off()

#151*151 matrix
var1_bw=ucv(UC_mt[,'var1'],nb=10000)
var2_bw=ucv(UC_mt[,'var2'],nb=10000)

# UC_all_10_16_sub$exp_x=unlist(lapply(UC_all_10_16_sub$within_mean,function(x) min(which(f1$x>=x))))
# UC_all_10_16_sub$exp_y=unlist(lapply(UC_all_10_16_sub$between_mean,function(x) min(which(f1$y>=x))))
# UC_all_10_16_sub$exp_z=unlist(lapply(1:nrow(UC_all_10_16_sub),function(x) return(f1$z[UC_all_10_16_sub$exp_x[x], UC_all_10_16_sub$exp_y[x]]+pse_count)))
# 
# UC_all_10_16_sub$obs_x=unlist(lapply(UC_all_10_16_sub$within_mean,function(x) min(which(f2$x>=x))))
# UC_all_10_16_sub$obs_y=unlist(lapply(UC_all_10_16_sub$between_mean,function(x) min(which(f2$y>=x))))
# UC_all_10_16_sub$obs_z=unlist(lapply(1:nrow(UC_all_10_16_sub),function(x) return(f2$z[UC_all_10_16_sub$exp_x[x], UC_all_10_16_sub$exp_y[x]]+pse_count)))

