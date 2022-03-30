source('mainFunctions_sub.R')
NME_in=readRDS(NME_agnostic_ASM_file)
MML_in=readRDS(MML_agnostic_ASM_file)
NME_in_tb=as.data.table(mcols(NME_in))
NME_in_tb$region=paste0(seqnames(NME_in),':',start(NME_in),'-',end(NME_in))
MML_in_tb=as.data.table(mcols(MML_in))
MML_in_tb$region=paste0(seqnames(MML_in),':',start(MML_in),'-',end(MML_in))
GR_merge_tb=rbind(NME_in_tb,MML_in_tb)
rm(NME_in_tb)
rm(MML_in_tb)
rm(NME_in)
rm(MML_in)
GR_merge_tb$K=NULL
GR_merge_tb$N=NULL
GR_merge=readRDS(GR_merge_file)
GR_merge_tb_asm=rbind(
  data.table(score=GR_merge$dNME,Sample=GR_merge$Sample,statistics='dNME',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dNME_pval,Sample=GR_merge$Sample,statistics='dNME_pval',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dMML,Sample=GR_merge$Sample,statistics='dMML',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dMML_pval,Sample=GR_merge$Sample,statistics='dMML_pval',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=(GR_merge$MML1+GR_merge$MML2)/2,Sample=GR_merge$Sample,statistics='MML_ASM',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=(GR_merge$NME1+GR_merge$NME2)/2,Sample=GR_merge$Sample,statistics='NME_ASM',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge)))
)
rm(GR_merge)

GR_merge_tb=rbind(GR_merge_tb,GR_merge_tb_asm)
GR_merge_tb=dcast.data.table(GR_merge_tb,Sample+region~statistics,value.var = "score")
saveRDS(GR_merge_tb,'../downstream/output/human_analysis/QC/GR_merge_ASM_comp.rds')
# # check NME vs NME from allele specific way -------------------------------
GR_merge_tb=readRDS('../downstream/output/human_analysis/QC/GR_merge_ASM_comp.rds')
GR_merge_tb=GR_merge_tb[!is.na(MML_ASM)&!is.na(MML)&!is.na(NME)]#3284912/52263042,total ASM regions=3332744
pdf('../downstream/output/graphs/Figure3/Fig-S5A-dNME_NME_all_pt.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb,aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb$NME_ASM,GR_merge_tb$NME)
pdf('../downstream/output/graphs/Figure3/Fig-S5B-dNME_NME_all_dMML_ASM_pt.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb[dMML_pval<=pval_cutoff],aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb[dMML_pval<=pval_cutoff]$NME_ASM,GR_merge_tb[dMML_pval<=pval_cutoff]$NME)
pdf('../downstream/output/graphs/Figure3/Fig-S5C-dNME_NME_non_dMML.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb[dMML_pval>pval_cutoff],aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+ xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+  
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb[dMML_pval>pval_cutoff]$NME_ASM,GR_merge_tb[dMML_pval>pval_cutoff]$NME)