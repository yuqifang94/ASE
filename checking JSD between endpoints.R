#Checking correlation between 1st and last time points

MML_in=readRDS('MML_matrix_mouse_DNAase_all_DNase_dedup_N2.rds')
MML_in_sub=granges(MML_in)
mcols(MML_in_sub)=mcols(MML_in)[,c("EFP-day10_5-1","EFP-day10_5-2","EFP-day15_5-1","EFP-day15_5-2","EFP-day11_5-1","EFP-day11_5-2")]
MML_in_sub=MML_in_sub[which(apply(mcols(MML_in_sub),1,function(x) sum(is.na(x))==0))]
MML_in_sub$diff_rep1=abs(MML_in_sub$`EFP-day10_5-1`-MML_in_sub$`EFP-day15_5-1`)
MML_in_sub$diff_rep2=abs(MML_in_sub$`EFP-day10_5-2`-MML_in_sub$`EFP-day15_5-2`)
MML_in_sub$diff_rep1_close=abs(MML_in_sub$`EFP-day10_5-1`-MML_in_sub$`EFP-day11_5-1`)
MML_in_sub$diff_rep2_close=abs(MML_in_sub$`EFP-day10_5-2`-MML_in_sub$`EFP-day11_5-2`)

MML_in_sub=granges(MML_in)
mcols(MML_in_sub)=mcols(MML_in)[,c("kidney-day14_5-1","kidney-day14_5-2","kidney-day15_5-1","kidney-day15_5-2","kidney-day16_5-1","kidney-day16_5-2")]
MML_in_sub=MML_in_sub[which(apply(mcols(MML_in_sub),1,function(x) sum(is.na(x))==0))]
MML_in_sub$diff_rep1=abs(MML_in_sub$`kidney-day14_5-1`-MML_in_sub$`kidney-day16_5-1`)
MML_in_sub$diff_rep2=abs(MML_in_sub$`kidney-day14_5-2`-MML_in_sub$`kidney-day16_5-2`)
MML_in_sub$diff_rep1_close=abs(MML_in_sub$`kidney-day14_5-1`-MML_in_sub$`kidney-day15_5-1`)
MML_in_sub$diff_rep2_close=abs(MML_in_sub$`kidney-day14_5-2`-MML_in_sub$`kidney-day15_5-2`)

MML_in_sub=granges(MML_in)
mcols(MML_in_sub)=mcols(MML_in)[,c("hindbrain-day10_5-1","hindbrain-day10_5-2","hindbrain-day15_5-1","hindbrain-day15_5-2","hindbrain-day11_5-1","hindbrain-day11_5-2")]
MML_in_sub=MML_in_sub[which(apply(mcols(MML_in_sub),1,function(x) sum(is.na(x))==0))]
MML_in_sub$diff_rep1=abs(MML_in_sub$`hindbrain-day10_5-1`-MML_in_sub$`hindbrain-day15_5-1`)
MML_in_sub$diff_rep2=abs(MML_in_sub$`hindbrain-day10_5-2`-MML_in_sub$`hindbrain-day15_5-2`)
MML_in_sub$diff_rep1_close=abs(MML_in_sub$`hindbrain-day10_5-1`-MML_in_sub$`hindbrain-day11_5-1`)
MML_in_sub$diff_rep2_close=abs(MML_in_sub$`hindbrain-day10_5-2`-MML_in_sub$`hindbrain-day11_5-2`)


UC_in_matrix=readRDS('../downstream/output/UC_agnostic_mouse_DNase_matrix_dedup_DNase_N2_endpoints.rds')
UC_in_matrix_sub=UC_in_matrix[which(apply(mcols(UC_in_matrix),1,function(x) sum(is.na(x))==0))]
#PCA plots
#UC_in_matrix_sub=UC_in_matrix_sub[-which(apply(mcols(UC_in_matrix_sub),1,function(x) sum(x>0.1)==0))]
UC_in_matrix_sub_df=PCA_df_prep(UC_in_matrix_sub)
pdf('../downstream/output/PCA_DNase_ft01.pdf')
autoplot(UC_in_matrix_sub_PCA,data=UC_in_matrix_sub_df,colour='sample')+geom_text(label=UC_in_matrix_sub_df$sample)
autoplot(UC_in_matrix_sub_PCA,data=UC_in_matrix_sub_df,colour='sample',x=3,y=4)+geom_text(label=UC_in_matrix_sub_df$sample)
autoplot(UC_in_matrix_sub_PCA,data=UC_in_matrix_sub_df,colour='sample',x=5,y=6)+geom_text(label=UC_in_matrix_sub_df$sample)
autoplot(UC_in_matrix_sub_PCA,data=UC_in_matrix_sub_df,colour='sample',x=7,y=8)+geom_text(label=UC_in_matrix_sub_df$sample)
dev.off()
#select regions
sample_all=gsub('-1','',colnames(mcols(UC_in_matrix)))
sample_all=gsub('-2','',sample_all)
UC_matrix_sub_out=list()
UC_enrich_out=data.frame()
gene_3kb=list()
mm10_DNase <- readRDS("../downstream/output/mm10_DNase.rds")
DNase_enrich<-function(UC_in,UC_sig,DNase_in){
  UC_sig_DNase_in=subsetByOverlaps(DNase_in,UC_in[UC_sig],type='equal')
  UC_non_sig_DNase_in=subsetByOverlaps(DNase_in,UC_in[-UC_sig],type='equal')
  UC_sig_DNase=sum(UC_sig_DNase_in$region_type=='DNase')
  UC_sig_control=sum(UC_sig_DNase_in$region_type=='control')
  UC_non_sig_DNase=sum(UC_non_sig_DNase_in$region_type=='DNase')
  UC_non_sig_control=sum(UC_non_sig_DNase_in$region_type=='control')
  return(fisher.test(matrix(c(UC_sig_DNase,UC_sig_control,UC_non_sig_DNase,UC_non_sig_control),nrow=2)))
}
for(sp in unique(sample_all)){
  rep1=paste(sp,'1',sep='-')
  rep2=paste(sp,'2',sep='-')
  UC_matrix_mt=granges(UC_in_matrix_sub)
  mcols(UC_matrix_mt)=mcols(UC_in_matrix_sub)[,c(rep1,rep2)]
  cutoff=quantile(unlist(mcols(UC_matrix_mt)),prob=0.99)
  UC_sig=which(apply(mcols(UC_matrix_mt),1,function(x) x[1]>=cutoff&x[2]>=cutoff))
  UC_matrix_mt=UC_matrix_mt[UC_sig]
  UC_matrix_sub_out[[sp]]=UC_matrix_mt
  UC_enrich=DNase_enrich(granges(UC_in_matrix_sub),UC_sig,mm10_DNase)
  UC_enrich_out=rbind(UC_enrich_out,data.frame(sample=sp,OR=UC_enrich$estimate,
                                               lowerCI=UC_enrich$conf.int[1],upperCI=UC_enrich$conf.int[2]))
  gene_3kb[[sp]]=subsetByOverlaps(mm10_DNase[abs(mm10_DNase$dist)<=3000&mm10_DNase$region_type=='DNase'],UC_in_matrix_sub[UC_sig])
  #write(unique(gene_3kb[[sp]])$gene,paste('../downstream/output/gene_',sp,'.txt',sep = ''))
  
}
write(unique( gene_3kb$`stomach-day14_5-day16_5`$gene),'../downstream/output/gene_stomach.txt')
write(unique( mm10_DNase$gene[abs(mm10_DNase$dist)<=3000&mm10_DNase$region_type=='DNase']),'../downstream/output/gene_mm10.txt')
#names(UC_matrix_sub_out)=NULL
UC_in_matrix_sig=unique(do.call('c',UC_matrix_sub_out))
UC_in_matrix_sig_all_dat=subsetByOverlaps(UC_in_matrix_sub,UC_in_matrix_sig,type='equal')
UC_in_matrix_sig_mt=t(as.matrix(mcols(UC_in_matrix_sig_all_dat)))
meta_mt=unlist(lapply(strsplit(rownames(UC_in_matrix_sig_mt),'-'),function(x) paste(x[1],x[2],x[3],sep = '-')))
comp_color=data.frame(comp=unique(meta_mt),
                      color= rainbow(length(unique(meta_mt))),stringsAsFactors = F)

comp_color=comp_color$color[match(meta_mt,comp_color$comp)]
library(gplots)
pdf('../downstream/output/JSD_endpoints_heatmap.pdf',width=15)
pheatmap(UC_in_matrix_sig_mt,col = col_fun(1024),dendrogram='none',trace='none',margins=c(5,20),
         Rowv=FALSE,RowSideColors = comp_color,na.rm = T,scale="none",symm=F)
dev.off()

#heatmap of JSD in between

heatmap_NME_MML<-function(UC_matrix_sub_out,sp,matrix_in){
  UC_matrix_sub_out_sp=UC_matrix_sub_out[[sp]]
  UC_matrix_all_dat_sig=subsetByOverlaps(matrix_in,UC_matrix_sub_out_sp,type='equal')
  
  UC_matrix_all_dat_sig_mt=as.matrix(mcols(UC_matrix_all_dat_sig))
  
  tissue=strsplit(sp,'-')[[1]][1]
  
  UC_matrix_all_dat_sig_mt=UC_matrix_all_dat_sig_mt[,which(unlist(lapply(strsplit(colnames(UC_matrix_all_dat_sig_mt),'-'),
                                                                         function(x) x[1]))==tissue)]
  meta_mt=unlist(lapply(strsplit(colnames(UC_matrix_all_dat_sig_mt),'-'),function(x) paste(x[1],x[2],x[3],sep = '-')))
  comp_color=data.frame(comp=unique(meta_mt),
                        color= heat.colors(length(unique(meta_mt))),stringsAsFactors = F)
  comp_color=comp_color$color[match(meta_mt,comp_color$comp)]
  heatmap.2(UC_matrix_all_dat_sig_mt,col = col_fun(1024),dendrogram='none',trace='none',margins=c(20,5),
            Colv=FALSE,ColSideColors  = comp_color,na.rm = T,scale="none",symm=F)
}
NME_matrix_in=readRDS('../downstream/output/NME_matrix_mouse_DNAase_all_DNase_dedup_N2.rds')
MML_matrix_in=readRDS('../downstream/output/MML_matrix_mouse_DNAase_all_DNase_dedup_N2.rds')
heatmap_NME_MML(UC_matrix_sub_out,'forebrain-day10_5-day16_5',NME_matrix_in)
heatmap_NME_MML(UC_matrix_sub_out,'forebrain-day10_5-day16_5',MML_matrix_in)
UC_in_matrix_sub_df_ft=as.data.frame(t(as.matrix(mcols(UC_in_matrix_sub_ft))))
UC_in_matrix_sub_PCA_ft=prcomp(UC_in_matrix_sub_df_ft, scale. = TRUE)
UC_in_matrix_sub_df_ft$sample=unlist(lapply(strsplit(rownames(UC_in_matrix_sub_df_ft),'-'),function(x) x[1]))
pdf('../downstream/output/PCA_DNase_ft.pdf')
autoplot(UC_in_matrix_sub_PCA_ft,data=UC_in_matrix_sub_df_ft,colour='sample')
dev.off()

cor_out=data.frame()
UC_matrix_sub_out=list()
for(sp in unique(sample_all)){
  rep1=paste(sp,'1',sep='-')
  rep2=paste(sp,'2',sep='-')
  UC_matrix_sub=granges(UC_in_matrix)
  mcols(UC_matrix_sub)=mcols(UC_in_matrix)[,c(rep1,rep2)]
  #UC_matrix_sub=UC_matrix_sub[which(apply(mcols(UC_matrix_sub),1,function(x) sum(is.na(x))==0))]
  UC_matrix_sub=UC_matrix_sub[which(apply(mcols(UC_matrix_sub),1,function(x) abs(x[1]-x[2])<=0.1))]
  #UC_matrix_sub_out[[sp]]=UC_matrix_sub
  cor_test_out=cor.test(mcols(UC_matrix_sub)[,rep1]/mean(mcols(UC_matrix_sub)[,rep1],na.rm = T),
                        mcols(UC_matrix_sub)[,rep2]/mean(mcols(UC_matrix_sub)[,rep2],na.rm = T))
  cor_out=rbind(cor_out,data.frame(sample=sp,cor=cor_test_out$estimate,regions=length(UC_matrix_sub)))
  
}
