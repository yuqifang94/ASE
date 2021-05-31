#Compare JSD and UC
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=48),
                                 axis.title.x=element_text(hjust=0.5,size=72,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=72,face="bold"),
                                 axis.text.x=element_text(size=70),
                                 axis.text.y=element_text(size=70))


#Use EFP10.5vs11.5 as example
UC_in=read.agnostic.mouse.uc('../downstream/input/UC/mm10_EFP_day10_5_all-vs-mm10_EFP_day11_5_all_uc.bedGraph')
UC_in=UC_in[UC_in$N>=2]
jsd_in=read.agnostic.mouse.uc('../downstream/input/UC/mm10_EFP_day10_5_all-vs-mm10_EFP_day11_5_all_jsd.bedGraph')
jsd_in=jsd_in[jsd_in$N>=2]
identical(paste0(seqnames(UC_in),':',start(UC_in),'-',end(UC_in)),paste0(seqnames(jsd_in),':',start(jsd_in),'-',end(jsd_in)))
cor(UC_in$score,jsd_in$score,method='spearman')# 0.8455,  0.8352@ DNase region
UC_jsd=data.table(regions=paste0(seqnames(UC_in),':',start(UC_in),'-',end(UC_in)),UC=UC_in$score,JSD=jsd_in$score,N=jsd_in$N)


#Find example
UC_in=import.bedGraph('../downstream/input/UC/mm10_EFP_day10_5_all-vs-mm10_EFP_day12_5_all_uc.bedGraph')
UC_in$N=UC_in$NA.
UC_in=UC_in[which(UC_in$N>=2)]
jsd_in=import.bedGraph('../downstream/input/jsd_output/mm10_EFP_day10_5_all-vs-mm10_EFP_day12_5_all_jsd.bedGraph')
jsd_in$N=jsd_in$NA.
jsd_in=jsd_in[which(jsd_in$N>=2)]
olap=findOverlaps(UC_in,jsd_in)
UC_in$jsd=NA
UC_in$jsd[queryHits(olap)]=jsd_in$score[subjectHits(olap)]
#dNME and dMML
#dNME
NME_10_5=read.agnostic.mouse.uc('mm10_EFP_day10_5_all_allele_agnostic_nme.bedGraph')
NME_11_5=read.agnostic.mouse.uc('mm10_EFP_day11_5_all_allele_agnostic_nme.bedGraph')
NME_10_5_mt=as.matrix(mcols(NME_10_5))
NME_11_5_mt=as.matrix(mcols(NME_11_5))
rownames(NME_10_5_mt)=paste0(seqnames(NME_10_5),':',start(NME_10_5),'-',end(NME_10_5))
rownames(NME_11_5_mt)=paste0(seqnames(NME_11_5),':',start(NME_11_5),'-',end(NME_11_5))
UC_jsd$dNME=abs(as.numeric(NME_10_5_mt[UC_jsd$regions,"score"])-as.numeric(NME_11_5_mt[UC_jsd$regions,"score"]))

#dMML
MML_10_5=read.agnostic.mouse.uc('mm10_EFP_day10_5_all_allele_agnostic_mml.bedGraph')
MML_11_5=read.agnostic.mouse.uc('mm10_EFP_day11_5_all_allele_agnostic_mml.bedGraph')
MML_10_5_mt=as.matrix(mcols(MML_10_5))
MML_11_5_mt=as.matrix(mcols(MML_11_5))
rownames(MML_10_5_mt)=paste0(seqnames(MML_10_5),':',start(MML_10_5),'-',end(MML_10_5))
rownames(MML_11_5_mt)=paste0(seqnames(MML_11_5),':',start(MML_11_5),'-',end(MML_11_5))
UC_jsd$dMML=abs(as.numeric(MML_10_5_mt[UC_jsd$regions,"score"])-as.numeric(MML_11_5_mt[UC_jsd$regions,"score"]))
#plotting
UC_jsd=readRDS('../downstream/output/UC_jsd_EFP.rds')
png('../downstream/output/UC_JSD/UC_JSD_cor01.png',width=1000,height=1000)
ggplot( UC_jsd[UC>=0.1&JSD>=0.1],aes(x=JSD,y=UC))+geom_point(size=0.1,alpha=0.2)+theme_glob
dev.off()

png('../downstream/output/UC_JSD/UC_JSD_all.png',width=1000,height=1000)
ggplot( UC_jsd,aes(x=JSD,y=UC))+geom_point(size=0.1,alpha=0.2)+theme_glob
dev.off()
png("../downstream/output/UC_JSD/UC_dNME.png",width=1000,height=1000)
ggplot(UC_jsd[UC>=0.1&JSD>=0.1],aes(x=UC,y=dNME))+geom_point(size=0.1,alpha=0.1)+theme_glob
dev.off()

png("../downstream/output/UC_JSD/UC_dMML.png",width=1000,height=1000)
ggplot(UC_jsd[UC>=0.1&JSD>=0.1],aes(x=UC,y=dMML))+geom_point(size=0.1,alpha=0.1)+theme_glob
dev.off()

png("../downstream/output/UC_JSD/JSD_dNME.png",width=1000,height=1000)
ggplot(UC_jsd[UC>=0.1&JSD>=0.1],aes(x=JSD,y=dNME))+geom_point(size=0.1,alpha=0.1)+theme_glob
dev.off()

png("../downstream/output/UC_JSD/JSD_dMML.png",width=1000,height=1000)
ggplot(UC_jsd[UC>=0.1&JSD>=0.1],aes(x=JSD,y=dMML))+geom_point(size=0.1,alpha=0.1)+theme_glob
dev.off()
#correlation
DNase=readRDS('../downstream/input/mm10_DNase.rds')
DNase=DNase[DNase$region_type=="DNase"]
UC_jsd=UC_jsd[UC>=0.1&JSD>=0.1]
cor(UC_jsd$UC,UC_jsd$dNME,method="spearman")
cor(UC_jsd$JSD,UC_jsd$dNME,method="spearman")
cor(UC_jsd$UC,UC_jsd$dMML,method="spearman")
cor(UC_jsd$JSD,UC_jsd$dMML,method="spearman")
UC_jsd_DNase=UC_jsd[UC_jsd$regions%in%paste0(seqnames(DNase),':',start(DNase),'-',end(DNase))]
cor(UC_jsd_DNase$UC,UC_jsd_DNase$dNME,method="spearman")
cor(UC_jsd_DNase$JSD,UC_jsd_DNase$dNME,method="spearman")
cor(UC_jsd_DNase$UC,UC_jsd_DNase$dMML,method="spearman")
cor(UC_jsd_DNase$JSD,UC_jsd_DNase$dMML,method="spearman")
#Correlation with N
UC_jsd$N_plot=UC_jsd$N
UC_jsd$N_plot[UC_jsd$N_plot>=10]=">=10"
UC_jsd$N_plot=factor(UC_jsd$N_plot,levels=c(1:10,">=10"))
jpeg('../downstream/output/UC_JSD/cor_with_N.jpeg')
UC_N_plot=ggplot(UC_jsd,aes(x=N_plot,y=UC))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('N CpG')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,1))
JSD_N_plot=ggplot(UC_jsd,aes(x=N_plot,y=JSD))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('N CpG')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,1))
ggarrange(UC_N_plot,JSD_N_plot,nrow=2,ncol=1,common.legend=T)
dev.off()
#Enhancer vs promoter

enhancer=readRDS("../downstream/output/chromHMM_enhancer.rds")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- GenomicFeatures::genes(txdb)
promoters <- promoters(genes,upstream=2000,downstream=1000)

enhancer=subsetByOverlaps(DNase,enhancer)
promoter=subsetByOverlaps(DNase,promoters)
enhancer=paste0(seqnames(enhancer),':',start(enhancer),'-',end(enhancer))
promoter=paste0(seqnames(promoter),':',start(promoter),'-',end(promoter))
UC_jsd$region_type=NA
UC_jsd$region_type[UC_jsd$regions%in%promoter]="promoter"
UC_jsd$region_type[UC_jsd$regions%in%enhancer]="enhancer"
quant_conv=c("Q1","Q2","Q3","Q4")
UC_jsd$UC_quant=findInterval(UC_jsd$UC,quantile(UC_jsd$UC,prob=c(0,0.25,0.5,0.75,1),na.rm=T))
UC_jsd$UC_quant[UC_jsd$UC_quant==5]=4
UC_jsd$UC_quant=quant_conv[UC_jsd$UC_quant]
UC_jsd$JSD_quant=findInterval(UC_jsd$JSD,quantile(UC_jsd$JSD,prob=c(0,0.25,0.5,0.75,1),na.rm=T))
UC_jsd$JSD_quant[UC_jsd$JSD_quant==5]=4
UC_jsd$JSD_quant=quant_conv[UC_jsd$JSD_quant]
pdf('../downstream/output/UC_JSD/Figure5B_dNME_dMML_enhancer_UC_quantile.pdf',width=3.5,height=3.5)

dNME_plot=ggplot(UC_jsd[!is.na(UC_jsd$region_type)],aes(x=UC_quant,y=dNME,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.75))
dMML_plot=ggplot(UC_jsd[!is.na(UC_jsd$region_type)],aes(x=UC_quant,y=dMML,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.3))
ggarrange(dNME_plot,dMML_plot,nrow=2,ncol=1,common.legend=T)
dev.off()

pdf('../downstream/output/UC_JSD/Figure5B_dNME_dMML_enhancer_JSD_quantile.pdf',width=3.5,height=3.5)

dNME_plot=ggplot(UC_jsd[!is.na(UC_jsd$region_type)],aes(x=JSD_quant,y=dNME,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('JSD quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.75))
dMML_plot=ggplot(UC_jsd[!is.na(UC_jsd$region_type)],aes(x=JSD_quant,y=dMML,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('JSD quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.3))
ggarrange(dNME_plot,dMML_plot,nrow=2,ncol=1,common.legend=T)
dev.off()
