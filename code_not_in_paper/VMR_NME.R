source('mainFunctions_sub.R')
OR_VMR<-function(NME_dat,vmr,percent,NME_quant='quant_score'){
  NME_dat$VMR=FALSE
  olap=findOverlaps(NME_dat,vmr)
  NME_dat$VMR[unique(queryHits(olap))]=TRUE
  NME_VMR=sum((elementMetadata(NME_dat)[[NME_quant]]%in%percent)&NME_dat$VMR)
  nonNME_VMR=sum((!elementMetadata(NME_dat)[[NME_quant]]%in%percent)&NME_dat$VMR)
  nonNME_nonVMR=sum((!elementMetadata(NME_dat)[[NME_quant]]%in%percent)&!NME_dat$VMR)
  NME_nonVMR=sum((elementMetadata(NME_dat)[[NME_quant]]%in%percent)&!NME_dat$VMR)
  #print(matrix(c(NME_VMR,nonNME_VMR,NME_nonVMR,nonNME_nonVMR),nrow = 2))
  fisher.test(matrix(c(NME_VMR,nonNME_VMR,NME_nonVMR,nonNME_nonVMR),nrow = 2))
}
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.text.x=element_text(size=10),
                                 axis.text.y=element_text(size=10))
# NME_VMR -----------------------------------------------------------------
NME_in=readRDS(NME_agnostic_file)
#Brain
load("../downstream/input/human_analysis/vmrs_hg19_brain.rda")
vmr_HC2=vmrs_hg19$HC2
vmr_HC1=vmrs_hg19$HC1
names(vmr_HC2)=NULL
names(vmr_HC1)=NULL
#Do HC2
vmr=do.call(c,vmr_HC2)
saveRDS(vmr,'../downstream/output/human_analysis/VMR_NME/vmr_HC2.rds')
vmr=readRDS('../downstream/output/human_analysis/VMR_NME/vmr_HC2.rds')
NME_in_brain=NME_in[NME_in$Sample%in%c('Brain_Hippocampus_middle_paired - 149','Brain_Hippocampus_middle_paired - 150')]
olap=findOverlaps(NME_in_brain,vmr)
NME_in_brain$VMR=NA
NME_in_brain$VMR[queryHits(olap)]=vmr$meanSDS[subjectHits(olap)]
#cor.test(NME_in_brain$VMR,NME_in_brain$NME,method="spearman")#-0.01201045, 0.00784
NME_in_brain=NME_in_brain[order(NME_in_brain$NME,decreasing = T)]
NME_in_brain$quant_score_1p=findInterval(NME_in_brain$NME,quantile(NME_in_brain$NME,prob=seq(0,1,0.01)))
NME_in_brain$quant_score_1p=NME_in_brain$quant_score_1p-1
OR_quant_1p=data.frame()
for(percent in unique(NME_in_brain$quant_score_1p)){
  OR=OR_VMR(NME_in_brain,vmr,percent,NME_quant='quant_score_1p')
  OR_quant_1p=rbind(OR_quant_1p,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
  
}
pdf('../downstream/output/human_analysis/VMR_NME/brain_NME_quantile.pdf',height=3.5,width=3.5)
ggplot(OR_quant_1p,aes(x=quant,y=OR))+geom_point(size=0.5)+geom_smooth()+xlab('NME quantile')+ylab('OR of VMR')+theme_glob
dev.off()

NME_in_brain$quant_score=findInterval(NME_in_brain$NME,quantile(NME_in_brain$NME,prob=seq(0.25,0.75,0.25)))
OR_quant=data.frame()
for(percent in unique(NME_in_brain$quant_score)){
  OR=OR_VMR(NME_in_brain,vmr,percent,NME_quant='quant_score')
  OR_quant=rbind(OR_quant,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
  
}
OR_quant$quant=c("0-25%","25%-50%","50%-75%","75%-100%")[OR_quant$quant+1]
write.csv(OR_quant,'../downstream/output/human_analysis/VMR_NME/brain_NME_quantile.csv')
#Lung
load("../downstream/input/human_analysis/List_of_VMRs_lung.rda")
VMR_Lung=list_of_VMRs$Lung_all
NME_in_lung=NME_in[NME_in$Sample%in%c("Lung_single - STL001","Lung_single - STL002")]
NME_in_lung$quant_score=findInterval(NME_in_brain$NME,quantile(NME_in_brain$NME,prob=seq(0.25,0.75,0.25)))
OR_quant=data.frame()
for(percent in unique(NME_in_lung$quant_score)){
  OR=OR_VMR(NME_in_lung,VMR_Lung,percent,NME_quant='quant_score')
  OR_quant=rbind(OR_quant,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
  
}
OR_quant$quant=c("0-25%","25%-50%","50%-75%","75%-100%")[OR_quant$quant+1]
write.csv(OR_quant,'../downstream/output/human_analysis/VMR_NME/lung_NME_quantile.csv')

NME_in_lung$quant_score_1p=findInterval(NME_in_lung$NME,quantile(NME_in_lung$NME,prob=seq(0,1,0.01)))
NME_in_lung$quant_score_1p=NME_in_lung$quant_score_1p-1
OR_quant_1p=data.frame()
for(percent in unique(NME_in_lung$quant_score_1p)){
  OR=OR_VMR(NME_in_lung,vmr,percent,NME_quant='quant_score_1p')
  OR_quant_1p=rbind(OR_quant_1p,data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2]))
  
}
pdf('../downstream/output/human_analysis/VMR_NME/lung_NME_quantile.pdf',height=3.5,width=3.5)
ggplot(OR_quant_1p,aes(x=quant,y=OR))+geom_point(size=0.5)+geom_smooth()+xlab('NME quantile')+ylab('OR of VMR')+theme_glob
dev.off()

##corSIV analysis: no significant result
Waterland_CorSIV <- as.data.frame(read_excel("../downstream/input/human_analysis/Waterland_CorSIV.xls", sheet = "S3"),stringsAsFactors=F)
CorSIV_loc=strsplit(Waterland_CorSIV$USCS_Coordinates_CoRSIV,':')
CorSIV_loc_start=unlist(lapply(CorSIV_loc,function(x) x[1]))
CorSIV_loc_ranges=lapply(CorSIV_loc,function(x) IRanges(as.numeric(strsplit(x[2],'-')[[1]][1]),width=601))
CorSIV_gr=unique(GRanges(seqnames=CorSIV_loc_start,ranges=do.call('c',CorSIV_loc_ranges)))
NME_in$quant_score=as.data.table(mcols(NME_in))[,list(quant_score=findInterval(NME,quantile(NME,prob=seq(0.25,0.75,0.25)))),by=list(Sample)]$quant_score
OR_quant_corSIV_sp=list()
for(sp in unique(NME_in$Sample)){
  OR_quant_corSIV=data.frame()
  tt1=proc.time()[[3]]
  for(percent in unique(NME_in$quant_score)){
    OR=OR_VMR(NME_in[NME_in$Sample==sp],CorSIV_gr,percent,NME_quant='quant_score')
    OR_quant_corSIV=rbind(OR_quant_corSIV,
                          data.frame(quant=percent,OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.in[1],upperCI=OR$conf.in[2],Sample=sp))
    OR_quant_corSIV_sp[[sp]]=OR_quant_corSIV
  }
  cat('Finish processing',sp,'in',proc.time()[[3]]-tt1,'\n')
}
saveRDS(OR_quant_corSIV_sp,'../downstream/output/human_analysis/VMR_NME/OR_quant_corSIV_sp.rds')

