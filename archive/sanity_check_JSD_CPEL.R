rm(list=ls())
source('mainFunctions_sub.R')
#Read in JSD files from informME
dir_JSD='../downstream/data/JSD/'
bed_out=GRanges()
for(fn in dir(dir_JSD,pattern="(JSD|dNME|dMML|NME|MML)")){
  cat("Processing:",fn,'\n')
  tt=proc.time()[[3]]
  bed_in=import.bed(paste0(dir_JSD,fn))
  fn=gsub("_all","",fn)
  fn=sub('.bed','',fn)
  fn=gsub('mm10_','',fn)
  stat_in=sub('-.*','',fn)
  fn=sub(paste0(stat_in,'-'),'',fn)
  tissue=sub('_.*','',fn)
  stage=gsub(paste0(tissue,'_'),'',fn)
  Sample=paste0(tissue,'-',stage)
  names(mcols(bed_in))="score"
  bed_in$statistics=stat_in
  bed_in$Sample=Sample
  bed_out=c(bed_out,bed_in)
  cat("Finished in:",proc.time()[[3]]-tt,'\n')
}
bed_out_tb=as.data.table(mcols(bed_out))
bed_out_tb$regions=paste0(seqnames(bed_out),':',start(bed_out),'-',end(bed_out))
rm(bed_in)
saveRDS(bed_out,'../downstream/output/bed_out_JSD.rds')
saveRDS(bed_out_tb,'../downstream/output/bed_out_tb_JSD.rds')
bed_out_tb=dcast.data.table(bed_out_tb,regions+Sample~statistics,value.var = 'score')
bed_out_tb=readRDS('../downstream/output/bed_out_tb_JSD.rds')

#Compare NME and MML

NME_in=readRDS('../downstream/input/NME_agnostic_mouse_all_merged.rds')
NME_in=NME_in[NME_in$tissue=="heart"]
NME_in$Sample=gsub("-all","",NME_in$Sample)
NME_in$NME_informME=NA
for(sp in unique(NME_in$Sample)){
  bed_out_tb_sp=bed_out_tb[Sample==sp&!is.na(NME)]
  olap=findOverlaps(NME_in[NME_in$Sample==sp],convert_GR(bed_out_tb_sp$regions),minoverlap=100)
  bed_out_tb_sp=data.table(qt=queryHits(olap),NME_informME=bed_out_tb_sp[subjectHits(olap)]$NME)
  bed_out_tb_sp=bed_out_tb_sp[,list(NME=mean(as.numeric(NME_informME))),by=list(qt)]
  NME_in[NME_in$Sample==sp][bed_out_tb_sp$qt]$NME_informME=bed_out_tb_sp$NME
}

NME_in=NME_in[!is.na(NME_in$NME_informME)]
saveRDS(NME_in,'../downstream/input/NME_agnostic_mouse_heart_merged.rds')
NME_in_dt=as.data.table(mcols(NME_in))
NME_in_dt$regions=paste0(seqnames(NME_in),':',start(NME_in),'-',end(NME_in))
heart_in=fread('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/heart_all.csv')
NME_in_dt_ft=NME_in_dt[N>=4]
NME_in_dt_ft=NME_in_dt_ft[regions %in% heart_in[GO_result!=""]$region]
NME_in_dt_ft=dcast(NME_in_dt_ft,regions~Sample,value.var = "NME")
NME_in_dt_ft_mt=as.matrix(NME_in_dt_ft[,2:8])
rownames(NME_in_dt_ft_mt)=NME_in_dt_ft$regions

heart_in=heart_in[region%in%NME_in_dt_ft$regions]
NME_in_dt_ft_mt=NME_in_dt_ft_mt[heart_in[GO_result!=""][order(dNME_maxJSD,decreasing=T)]$region,]
NME_in_dt_ft_mt[which(rowSums(NME_in_dt_ft_mt<=0.15)>=1),][1:20,]clu=readRDS('../downstream/input/jsd.rds')
clu=clu$heart
pdf('../downstream/output/sanity_check/NME_comparison_used_region.pdf')
ggplot(NME_in_dt[regions%in%names(clu)],aes(x=NME,y=NME_informME))+geom_bin2d(bins=100)+  scale_fill_continuous(type = "viridis") +
  xlab("CPEL")+ylab("infromME")+geom_abline(slope=1,color='red')+geom_smooth(color="white")
dev.off()
rm(NME_in)
#MML

MML_in=readRDS('../downstream/input/MML_agnostic_mouse_all_merged.rds')
MML_in=MML_in[MML_in$tissue=="heart"]
MML_in$Sample=gsub("-all","",MML_in$Sample)
MML_in$MML_informME=NA
for(sp in unique(MML_in$Sample)){
  bed_out_tb_sp=bed_out_tb[Sample==sp&!is.na(MML)]
  olap=findOverlaps(MML_in[MML_in$Sample==sp],convert_GR(bed_out_tb_sp$regions),minoverlap = 100)
  bed_out_tb_sp=data.table(qt=queryHits(olap),MML_informME=bed_out_tb_sp[subjectHits(olap)]$MML)
  bed_out_tb_sp=bed_out_tb_sp[,list(MML=mean(as.numeric(MML_informME))),by=list(qt)]
  MML_in[MML_in$Sample==sp][bed_out_tb_sp$qt]$MML_informME=bed_out_tb_sp$MML
}
MML_in=MML_in[!is.na(MML_in$MML_informME)]
saveRDS(MML_in,'../downstream/input/MML_agnostic_mouse_heart_merged_125.rds')
MML_in_dt=as.data.table(mcols(MML_in))
MML_in_dt$regions=paste0(seqnames(MML_in),':',start(MML_in),'-',end(MML_in))
pdf('../downstream/output/sanity_check/MML_comparison2.pdf')
ggplot(MML_in_dt[regions%in%names(clu)],aes(x=MML,y=MML_informME))+geom_bin2d(bins=100)+  scale_fill_continuous(type = "viridis") +
  xlab("CPEL")+ylab("infromME")+geom_abline(slope=1,color='red')+geom_smooth(color="white")
dev.off()
rm(MML_in)

#JSD: plot the same heatmap at same region
clu=readRDS('../downstream/input/jsd.rds')
JSD=bed_out_tb[!is.na(JSD)]
JSD=dcast.data.table(JSD,regions~Sample,value.var = 'JSD')
clu=clu$heart
olap=findOverlaps(convert_GR(JSD$regions),convert_GR(names(clu)),minoverlap = 100)
JSD_olap=JSD[queryHits(olap)]
JSD_olap$regions_clu=names(clu)[subjectHits(olap)]
JSD_olap$clu=clu[subjectHits(olap)]

JSD_olap=JSD_olap[,list(`heart-E10.5-VS-E11.5`=mean(as.numeric(`heart-E10.5-VS-E11.5`)),
              `heart-E11.5-VS-E12.5`=mean(as.numeric(`heart-E11.5-VS-E12.5`)),
              `heart-E12.5-VS-E13.5`=mean(as.numeric(`heart-E12.5-VS-E13.5`)),
              `heart-E13.5-VS-E14.5`=mean(as.numeric(`heart-E13.5-VS-E14.5`)),
              `heart-E14.5-VS-E15.5`=mean(as.numeric(`heart-E14.5-VS-E15.5`)),
              `heart-E15.5-VS-E16.5`=mean(as.numeric(`heart-E15.5-VS-E16.5`))),
         by=list(regions_clu,clu)]
JSD_olap=JSD_olap[order(clu)]
JSD_olap=JSD_olap[rowSums(is.na(JSD_olap[,3:8]))==0]#86% of region covered
JSD_olap_mt=as.matrix(JSD_olap[,3:8])
rownames(JSD_olap_mt)=JSD_olap$regions_clu
rowann <- data.frame(cluster=JSD_olap$clu,stringsAsFactors = F)
rownames(rowann) <- JSD_olap$regions_clu
c2 <- brewer.pal(10,'Set3')
names(c2) <- 1:10
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}
pdf('../downstream/output/sanity_check/UC_comparison_heatmap.pdf')
pheatmap(scalematrix(JSD_olap_mt),cluster_rows = F,annotation_row = rowann,cluster_cols = F,
         show_colnames = F,show_rownames = F,gaps_row = cumsum(table(JSD_olap$clu)),
         annotation_colors = list(cluster=c2))
dev.off()
UC_in=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
UC_in=UC_in$heart
UC_in_mt=as.matrix(mcols(UC_in))
UC_in_mt=UC_in_mt[,colnames(UC_in_mt)%in%paste0(gsub('VS-','',colnames(JSD_olap_mt)),'-all')]
rownames(UC_in_mt)=paste0(seqnames(UC_in),':',start(UC_in),'-',end(UC_in))
UC_in_mt=UC_in_mt[rownames(JSD_olap_mt),]
UC_in_mt=as.data.table(cbind(UC_in_mt,JSD_olap_mt))
UC_in_mt=melt.data.table(UC_in_mt,variable.name = 'sample',value.name = 'score')
UC_in_mt$stat_type="UC"
UC_in_mt[grepl('VS',sample)]$stat_type="JSD"
UC_in_mt=as.data.table(unstack(UC_in_mt,score~stat_type))
UC_in_mt$JSD_rank=rank(UC_in_mt$JSD)
UC_in_mt$UC_rank=rank(UC_in_mt$UC)
pdf('../downstream/output/sanity_check/JSD_UC_comparison.pdf')
ggplot(UC_in_mt,aes(x=UC,y=JSD))+geom_bin2d(bins=500)+  scale_fill_continuous(type = "viridis") +
  xlab("CPEL")+ylab("infromME")+geom_abline(slope=1,color='red')+geom_smooth(color="white")
dev.off()
pdf('../downstream/output/sanity_check/JSD_UC_comparison_ranked.pdf')
ggplot(UC_in_mt,aes(x=JSD_rank,y=UC_rank))+geom_bin2d(bins=500)+  scale_fill_continuous(type = "viridis") +
  xlab("CPEL")+ylab("infromME")+geom_abline(slope=1,color='red')+geom_smooth(color="white")
dev.off()
#Max UC vs max JSD
UC_in_mt=as.matrix(mcols(UC_in))
rownames(UC_in_mt)=paste0(seqnames(UC_in),':',start(UC_in),'-',end(UC_in))
UC_in_mt=UC_in_mt[,colnames(UC_in_mt)%in%paste0(gsub('VS-','',colnames(JSD_olap_mt)),'-all')]
UC_in_mt=UC_in_mt[rownames(JSD_olap_mt),]
UC_JSD_max=data.table(maxUC=rowMaxs(UC_in_mt),max_which_UC=max.col(UC_in_mt),
                      maxJSD=rowMaxs(JSD_olap_mt),max_which_JSD=max.col(JSD_olap_mt))

pdf('../downstream/output/sanity_check/JSD_UC_comparison_max_JSD.pdf')
ggplot(UC_JSD_max,aes(x=max_which_UC,y=max_which_JSD))+geom_point(alpha=0.1)+  
  xlab("CPEL")+ylab("infromME")+geom_abline(slope=1,color='red')
dev.off()

pdf('../downstream/output/sanity_check/JSD_UC_comparison_max_JSD_value.pdf')
ggplot(UC_JSD_max,aes(x=maxUC,y=maxJSD))+geom_bin2d(bins=500)+  scale_fill_continuous(type = "viridis") +
  xlab("CPEL")+ylab("infromME")+geom_abline(slope=1,color='red')+geom_smooth(color="white")
dev.off()


sum(UC_JSD_max$max_which_UC==UC_JSD_max$max_which_JSD)/nrow(UC_JSD_max)#0.1868383



#Simulate CPEL and informME
p=c(0.25,0.5,0.25)#infromME
p=c(0.25,0.25,0.25,0.25)#CPEL
MML=sum(p*c(0,1,2))/2
NME=-1/2*sum(p*log2(p))
barplot(p,col='blue',names.arg=c(0,1,2),ylab='p(x)')

p=c(0.99,0.01/3,0.01/3,0.01/3)#CPEL
MML=sum(p*c(0,1,1,2))/2
NME=-1/2*sum(p*log2(p))
p=c(0.99,0.02/3,0.01/3)#informME
MML=sum(p*c(0,1,2))/2
NME=-1/2*sum(p*log2(p))
barplot(p,col='blue',names.arg=c(0,1,2),ylab='p(x)')
