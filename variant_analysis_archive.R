#Heatmap:

density_df=data.frame(density=factor(round(NME_ASM_het_sub$density,digits = 3)),diff=factor(round(NME_ASM_het_sub$diff,digits=2)))
density_df=count(density_df,c('density','diff'))
density_df_dc=dcast(density_df,diff~density,value.var = 'freq')
density_df_dc[is.na(density_df_dc)]=0

density_heatmap=density_df_dc[,2:ncol(density_df_dc)]
rownames(density_heatmap)=density_df_dc$diff
density_heatmap=scale(density_heatmap,center=FALSE)
heatmap.2(density_heatmap,scale='none',col = bluered(100),trace = "none", density.info = "none",dendrogram='none', Rowv=FALSE, Colv=FALSE)


#NME vs density
NME_ASM_allele=NME_allele_ASM_calc[[2]]
NME_ASM_allele_het=NME_ASM_allele[NME_ASM_allele$HetCpG]

#Higher density indicate higer NME
density_df=data.frame(density=round(log10(NME_ASM_allele_het$density),digits = 2),NME=NME_ASM_allele_het$Value)
density_agg=aggregate(density_df,by=list(density_df$density),FUN=median)
plot(density_agg$Group.1,density_agg$NME,xlab='log10(density)',ylab='allelic NME')
####Calculating enrichment in het CpG vs ASM####
NME_all_gr$density_bin=round(NME_all$density*2,digits = 2)/2
NME_all_gr_enrich = data.frame(density=unique(NME_all_gr$density_bin),enrichment=NA)
for (log_den in unique(NME_all_gr$density_bin)){
  NME_all_gr_enrich$enrichment[NME_all_gr_enrich$density==log_den]=ASM_het_enrichment(NME_all_gr[NME_all_gr$density_bin==log_den])$estimate
}
NME_all_gr_enrich=NME_all_gr_enrich[NME_all_gr_enrich$density<=0.1,]
plot(NME_all_gr_enrich$density,NME_all_gr_enrich$enrichment,xlab='CpG density',ylab='Odds Ratio')


#############################Het CpG enrichment analysis##################################

for(sp in unique(variant_HetCpG_NME$Sample)){
  variant_sub_dNME=c(variant_sub_dNME,
                     subsetByOverlaps(variant_HetCpG_NME[variant_HetCpG_NME$Sample==sp],GR_merge_exclude_GM[GR_merge_exclude_GM$Sample==sp]))
}

ggplot(ASM_enrich_all$dMML,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 8)+
  ggtitle('dMML Het CpG enrichment')+xlab('Sample name')+ylab('Odds Ratio')+theme_bar

ggplot(ASM_enrich_all$UC,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 8) +
  ggtitle('UC Het CpG enrichment')+xlab('Sample name')+ylab('Odds Ratio')+theme_bar

#variant_HetCpG_meta=variant_HetCpG_meta[variant_HetCpG_meta$Sample %in% c( unique(variant_HetCpG_meta$Sample)[1:43])]
ASM_het_enrichment(variant_HetCpG_meta)
ASM_enrich_all=list()


p2 <- ggplot(data.frame(NME=c(GR_merge$NME1,GR_merge$NME2),
                        MML=c(GR_merge$MML1,GR_merge$MML2)), aes(NME, MML)) +
  geom_pointdensity() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(axis.text = element_text(size=1), axis.title = element_text(size=1))

################################Density analysis#########################################
gr_allele_CpG$TpA_CpG=gr_allele_CpG$TpA_count_extend/gr_allele_CpG$CG_hg19_extend
#Find ways to calculate CpG density, here use number of CG/size of the region
#Plot density distribution
#hg19 CG density
cpg_hg19_density=getCpgdensH19()
hist(log10(unlist(cpg_hg19_density[[2]])),xlab='log10(distance)',main='Distance between CpG site across hg19')

#####################################Global observation on CpG density etc #################################
#Put something into NME all
#CpG content distribution
NME_all$CG_cont_type[NME_all$ASM=='No' & NME_all$CGcount_diff!=0]='Imbalanced CG Non ASM'
NME_all$CG_cont_type[NME_all$ASM=='No' & NME_all$CGcount_diff==0]='Balanced CG Non ASM'
NME_all$CG_cont_type[NME_all$ASM=='Yes' & NME_all$CGcount_diff!=0]='Imbalanced CG ASM'
NME_all$CG_cont_type[NME_all$ASM=='Yes' & NME_all$CGcount_diff==0]='Balanced CG ASM'
#Choose from CG content, CG density, TpA/CpG, TA/CG
plot_density(data.frame(density=NME_all$TpA_count_extend/NME_all$CpG_count_extend,CG_type=NME_all$CG_cont_type),
             ylab='TpA count/CpG count',ylim=c(0,30),title='')
plot_density(data.frame(density=NME_all$TpA_count_extend,CG_type=NME_all$CG_cont_type),
             ylab='TpA count',ylim=c(0,100),title='')
plot_density(data.frame(density=NME_all$CGcount_hg19_extend,CG_type=NME_all$CG_cont_type),
             ylab='Oberseved/Expected CpG',ylim=c(0,0.6),title='')
plot_density(data.frame(density=NME_all$CpG_count_extend/NME_all$gff_size_extend,CG_type=NME_all$CG_cont_type),
             ylab='CpG density',ylim=c(0,0.03),title='')
plot_density(data.frame(density=NME_all$AT_count_extend,CG_type=NME_all$CG_cont_type),
             ylab='AT count',ylim=c(0,500),title='')
plot_density(data.frame(density=NME_all$AT_count_extend/NME_all$CG_count_extend,CG_type=NME_all$CG_cont_type),
             ylab='AT/CG count',ylim=c(0,0.3),title='')

############################## MML analysis #######################################################
MML_allele=gr_allele_CpG[gr_allele_CpG$Statistic=='MML']
MML_allele=add_ASM(MML_allele,GR[GR$Statistic=='dMML'])
MML_allele_ASM=MML_allele[which(MML_allele$ASM=='Yes')]
MML_allele_calc=allele_diff(MML_allele)
MML_allele_ASM_calc=allele_calc_plot(MML_allele_ASM,'MML','../downstream/output/dMML/','ASM Region')
allele_plot(MML_allele_ASM_calc,'MML','../downstream/output/dMML/','ASM Region')
MML_allele_all_calc=allele_calc_plot(MML_allele,'MML','../downstream/output/dMML/','all')
saveRDS(MML_allele_ASM_calc,'../downstream/output/dMML/MML_ASM_het_calc.rds')
#Generate bed files for genome browser
outdir='../downstream/output/dMML/'
MML_ASM=MML_allele_ASM_calc[[1]]
for (sp in unique(MML_ASM_het$Sample)){ASM_bed_gen_sp(MML_ASM_het,gr_allele_CpG,sp,'../downstream/output/dMML/bedallele/')}
#####Find the hist gram of MML
MML_allele_ASM=MML_allele[which(MML_allele$pval<=pval_cutoff)]
hist(MML_allele_ASM$Value,xlab='Mean methylation level',main='Mean methylation value at each allele at ASM')

##############################Het Cpg Analysis: Density effect #######################################################
############################## NME analysis #######################################################
#CG count enrichment
#Using all NME calculation, define ASM vs non-ASM, het CG vs non-het CG 
NME_all=NME_allele_calc[[1]] #Granges for each allele
NME_all=NME_all[!is.na(NME_all$pval)]
NME_all$CG_type=NA
NME_all$CG_type[NME_all$CGcount_diff!=0] ='Imbalanced CG'
NME_all$CG_type[NME_all$CGcount_diff==0] ='Balanced CG'
NME_all$ASM=NA
NME_all$ASM[NME_all$pval<=pval_cutoff] ='Yes'
NME_all$ASM[NME_all$pval>pval_cutoff] ='No'
NME_ASM=NME_all[NME_all$ASM=='Yes']
NME_het=NME_all[NME_all$CG_type=='Imbalanced CG']
NME_all$HetCpG=NME_all$CGcount_diff!=0
NME_ASM_het=NME_het[which(NME_het$ASM=='Yes')]
subject_old=unique(NME_all$Sample)[1:43]
subject_new=unique(NME_all$Sample)[44:49]
NME_all=NME_all[NME_all$Sample%in%c(subject_old,'merged - H1','merged - GM12878')]
#Subset for SNP-containing ranges
NME_hap=SNP_conmtaining_hap(NME_all,variant_HetCpG)
NME_SNP=NME_hap[[1]]
NME_non_SNP=NME_hap[[2]]
#Enrichement test for SNP-overlapping region
ASM_het_enrichment(NME_SNP)
ASM_het_enrich(NME_SNP,'dNME HetCpG enrichment for SNP overlapping haplotypes')

#Enrichement test for non-SNP containing het CpG
ASM_het_enrichment(NME_non_SNP)
ASM_het_enrich(NME_non_SNP,'dNME HetCpG enrichment for non SNP overlapping haplotypes')

#Enrichement test for all het CpG
ASM_het_enrichment(NME_all)
ASM_het_enrich(NME_all,'dNME HetCpG enrichment for all')
##############################Science paper: find average NME at 2 alleles at dMML region #################################################
NME_all=NME_allele_calc[[2]]
NME_all=NME_all[width(NME_all)>2]
NME_mean_dMML=c()
#Looking at dMML ASM
for (sp in unique(GR$Sample)){NME_mean_dMML=c(NME_mean_dMML,
                                              subsetByOverlaps(NME_all[NME_all$Sample==sp],GR[GR$Sample==sp & GR$Statistic=='dMML' &GR$pvalue<=pval_cutoff])$Value)
}
NME_ASM=data.frame(NME_mean=NME_mean,NME_type='dMML ASM',line_type='solid',stringsAsFactors = F)
#Looking at imprinted region dMML >0.9
NME_mean_imp=c()
for (sp in unique(GR$Sample)){NME_mean_imp=c(NME_mean_imp,subsetByOverlaps(NME_all[NME_all$Sample==sp],
                                                                           GR[GR$Sample==sp & GR$Statistic=='dMML' &GR$pvalue<=pval_cutoff & GR$Value>=0.9])$Value)
}
NME_ASM_imp=data.frame(NME_mean=NME_mean_imp,NME_type='dMML ASM bimodal',line_type='solid',stringsAsFactors = F)
NME_mean_non_ASM=c()
#Looking at dMML nonASM
for (sp in unique(GR$Sample)){NME_mean_non_ASM=c(NME_mean_non_ASM,subsetByOverlaps(NME_all[NME_all$Sample==sp],
                                                                                   GR[GR$Sample==sp & GR$Statistic=='dMML' &GR$pvalue>pval_cutoff])$Value)}
NME_non_ASM=data.frame(NME_mean=NME_mean_non_ASM,NME_type='Non dMML ASM',line_type='dashed',stringsAsFactors = F)
dNME_ASM=data.frame(NME_mean=NME_allele_ASM_calc[[2]]$Value,NME_type='dNME ASM',line_type='soild')
NME_mean_df=rbind(NME_ASM,NME_ASM_imp,NME_non_ASM,dNME_ASM)

ggplot(NME_mean_df,aes(x=NME_mean,color=NME_type))+
  geom_density(aes(linetype=line_type),alpha=0.6,size=1)+xlab('Mean NME')+ggtitle('Average NME at 2 alleles')+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+ylim(0,4)+scale_color_manual(values=c("blue","red","black",'purple'))+
  geom_hline(yintercept=0, colour="white", size=1)


ggplot(NME_mean_df,aes(x=NME_type,y=NME_mean,fill=NME_type))+
  geom_boxplot()+xlab('Type of ASM')+ggtitle('NME at 2 alleles')+ylab('NME')+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+ylim(0,1)



#It looks ASM have higher AT/CG count redo analysis
plot_density(data.frame(density=NME_all$TpA_count_extend/NME_all$CpG_count_extend,CG_type=NME_all$ASM),
             ylab='TpA count/CpG count',ylim=c(0,30),title='',xlab='ASM')
plot_density(NME_CG_df=data.frame(density=NME_all$AT_count_extend/NME_all$CG_count_extend,CG_type=NME_all$ASM),
             ylab='AT/CG count',ylim=c(0,0.3),title='',xlab='ASM')
plot_density(data.frame(density=NME_all$CGcount_hg19_extend,CG_type=NME_all$ASM),
             ylab='Oberseved/Expected CpG',ylim=c(0,0.6),title='',xlab='ASM')
#Checking those regions that have no CpG number difference but with Het CpG

cor(NME_ASM_het$CGcount_diff/NME_ASM_het$CGcount_hg19_extend,NME_ASM_het$diff)
cor(NME_ASM_het$CpGdiff/width(NME_ASM_het),NME_ASM_het$diff)
cor(NME_ASM_het$CGcount_diff,NME_ASM_het$diff)




############################## MML analysis #######################################################
#CG count enrichment
#Using all MML calculation, define ASM vs non-ASM, het CG vs non-het CG 
MML_all=MML_allele_calc[[1]] #Granges for each allele
MML_all=MML_all[!is.na(MML_all$pval)]
MML_all$CG_type=NA
MML_all$CG_type[MML_all$CGcount_diff!=0] ='Imbalanced CG'
MML_all$CG_type[MML_all$CGcount_diff==0] ='Balanced CG'
MML_all$ASM=NA
MML_all$ASM[MML_all$pval<=pval_cutoff] ='Yes'
MML_all$ASM[MML_all$pval>pval_cutoff] ='No'
MML_ASM=MML_all[MML_all$ASM=='Yes']
MML_het=MML_all[MML_all$CG_type=='Imbalanced CG']
MML_all$HetCpG=MML_all$CGcount_diff!=0
MML_ASM_het=MML_het[which(MML_het$ASM=='Yes')]
subject_old=unique(MML_all$Sample)[1:43]
subject_new=unique(MML_all$Sample)[44:49]
MML_all=MML_all[MML_all$Sample%in%c(subject_old,'merged - H1','merged - GM12878')]
#Subset for SNP-containing ranges
MML_hap=SNP_conmtaining_hap(MML_all,variant_HetCpG)
MML_SNP=MML_hap[[1]]
MML_non_SNP=MML_hap[[2]]
#Enrichement test for SNP-overlapping region: Het CpG only enriched in dNME, likely MML difference are caused by dNME
ASM_het_enrichment(MML_SNP)
ASM_het_enrich(MML_SNP,'dMML HetCpG enrichment for SNP overlapping haplotypes')

#Enrichement test for non-SNP containing het CpG
ASM_het_enrichment(MML_non_SNP)
ASM_het_enrich(MML_non_SNP,'dMML HetCpG enrichment for non SNP overlapping haplotypes')

#Enrichement test for all het CpG
ASM_het_enrichment(MML_all)
ASM_het_enrich(MML_all,'dMML HetCpG enrichment for all')


###########################Plot enrichment at each genomic feature####################
MML_all=MML_allele_calc[[1]]
NME_all=NME_allele_calc[[1]]
genomic_features=readRDS('../downstream/input/genomic_features.rds')
#Islands
genome_feature_plot(MML_all,genomic_features$`CpG island`,'dMML','CpG Island enrichment',ylim=c(0,10))
genome_feature_plot(NME_all,genomic_features$`CpG island`,'dNME','CpG Island enrichment',ylim=c(0,10))
#Open seas
genome_feature_plot(MML_all,genomic_features$`CpG open sea`,'dMML','CpG open sea enrichment',ylim=c(0,10))
#Promoters
genome_feature_plot(MML_all,genomic_features$promoter,'dMML','dMML Promoter enrichment',ylim=c(0,6))
genome_feature_plot(NME_all,genomic_features$promoter,'dNME','dNME Promoter enrichment',ylim=c(0,6))

#enhancer
genome_feature_plot(MML_all,genomic_features$enhancer,'dMML','enhancer enrichment dMML',ylim=c(0,2))
genome_feature_plot(NME_all,genomic_features$enhancer,'dNME','enhancer enrichment dNME',ylim=c(0,2))

###############Multiple analysis############################
#dNME containing motif
GR_merge_multiple_motif=subsetByOverlaps(GR_merge_multiple,motif_gene_strong)
write(unique(GR_merge_multiple_motif$genes_promoter),"../downstream/output/dNME_promoter_motif_sig.txt")
#Checking motif in GR_merge_multiple
motif_sig_strong_olap=subsetByOverlaps(motif_gene_strong,GR_merge_gr_multiple,maxgap = 50)
motif_gene_count=data.frame(gene=c(),sig_motif_n=c(),sig_non_motif_n=c(), non_sig_motif_n=c(),non_sig_non_motif_n=c(),OR_lowCI=c(),
                            OR_highCI=c(),OR=c())
for (motif in unique(motif_gene_strong$geneSymbol)){
  
  GR_merge_uq=unique(subsetByOverlaps(GR_merge,genomic_features$promoter,maxgap = 1000))
  # GR_merge_uq=unique(subsetByOverlaps(GR_merge_uq,motif_gene_strong[motif_gene_strong$geneSymbol==motif],maxgap=1000))
  all_olap=findOverlaps(GR_merge_uq,motif_gene_strong[motif_gene_strong$geneSymbol==motif],maxgap = 100)
  motif_all_n=length(unique(queryHits(all_olap)))
  non_motif_all_n=length(GR_merge_uq)- motif_all_n
  
  GR_merge_gr_multiple_sub=subsetByOverlaps(GR_merge_gr_multiple,genomic_features$promoter,maxgap = 1000)
  #GR_merge_gr_multiple_sub=subsetByOverlaps(GR_merge_gr_multiple_sub,motif_gene_strong[motif_gene_strong$geneSymbol==motif],maxgap =1000)
  sig_motif_olap=findOverlaps(GR_merge_gr_multiple_sub,motif_gene_strong[motif_gene_strong$geneSymbol==motif],maxgap = 100)
  sig_motif_n=length(unique(queryHits(sig_motif_olap)))
  sig_non_motif_n=length(GR_merge_gr_multiple_sub)-sig_motif_n
  non_sig_motif_n=motif_all_n-sig_motif_n
  non_sig_non_motif_n=non_motif_all_n-sig_non_motif_n
  OR_out=fisher.test(matrix(c(sig_motif_n,sig_non_motif_n,non_sig_motif_n,non_sig_non_motif_n),nrow=2))
  
  motif_gene_count=rbind(motif_gene_count,data.frame(gene=motif,sig_motif_n=sig_motif_n,sig_non_motif_n=sig_non_motif_n,
                                                     non_sig_motif_n=non_sig_motif_n,non_sig_non_motif_n=non_sig_non_motif_n,OR_lowCI=OR_out$conf.int[1],
                                                     OR_highCI=OR_out$conf.int[2],OR=OR_out$estimate[[1]]))
  
  
  
}
motif_gene_count[order(motif_gene_count$OR_lowCI,decreasing = T),]
saveRDS(motif_gene_count,"../downstream/output/motif_gene_count_100.rds")
saveRDS(motif_gene_count,"../downstream/output/motif_gene_count_promoter_1k.rds")
motif_gene_count_sig=motif_gene_count[motif_gene_count$n_sig>0,]
motif_gene_count_sig$freq=motif_gene_count_sig$n_sig/motif_gene_count_sig$n_all
motif_gene_count_sig=motif_gene_count_sig[order(motif_gene_count_sig$freq,decreasing = T),]



# H1_rep1_g1=read.table("../downstream/input/RNA_seq/H1_rep1.genome1.genes.out",sep='\t',header=T)
# H1_rep2_g1=read.table("../downstream/input/RNA_seq/H1_rep2.genome1.genes.out",sep='\t',header=T)
# H1_rep3_g1=read.table("../downstream/input/RNA_seq/H1_rep3.genome1.genes.out",sep='\t',header=T)
# 
# H1_rep1_g2=read.table("../downstream/input/RNA_seq/H1_rep1.genome2.genes.out",sep='\t',header=T)
# H1_rep2_g2=read.table("../downstream/input/RNA_seq/H1_rep2.genome2.genes.out",sep='\t',header=T)
# H1_rep3_g2=read.table("../downstream/input/RNA_seq/H1_rep3.genome2.genes.out",sep='\t',header=T)
# 
# RNA_genome= data.frame(genes=H1_rep1_g1$Gene.Name,rep1.g1=H1_rep1_g1$Coverage,
#                        rep2.g1=H1_rep2_g1$Coverage[match(H1_rep1_g1$Gene.Name,H1_rep2_g1$Gene.Name)],
#                        rep3.g1=H1_rep3_g1$Coverage[match(H1_rep1_g1$Gene.Name,H1_rep3_g1$Gene.Name)],
#                        rep1.g2=H1_rep1_g2$Coverage[match(H1_rep1_g1$Gene.Name,H1_rep1_g2$Gene.Name)],
#                        rep2.g2=H1_rep2_g2$Coverage[match(H1_rep1_g1$Gene.Name,H1_rep2_g2$Gene.Name)],
#                        rep3.g2=H1_rep3_g2$Coverage[match(H1_rep1_g1$Gene.Name,H1_rep3_g2$Gene.Name)])
# 
# 
# 
# RNA_genome$mean.g1=rowMeans(RNA_genome[,c("rep1.g1","rep2.g1","rep3.g1")])
# RNA_genome$mean.g2=rowMeans(RNA_genome[,c("rep1.g2","rep2.g2","rep3.g2")])
# 
# genome=factor(c(1,1,1,2,2,2))
# design=model.matrix(~0+genome)
# RNA_dge <- DGEList(counts=as.matrix(RNA_genome[,2:7]),group=c(1,1,1,2,2,2),genes=RNA_genome$genes)
# keep <- filterByExpr(RNA_dge)
# RNA_dge=RNA_dge[keep,,keep.lib.sizes=FALSE]
# ###if the genes are identical?
# which(!H1_rep1_g1$Gene.Name %in% H1_rep3_g1$Gene.Name)


#########H1###################
chromHMM_H1_dNME<-chromHMM_OR(GR_merge,"AH46858","merged - H1")
chromHMM_H1_dMML<-chromHMM_OR(GR_merge,"AH46858","merged - H1",stat="dMML_pval")
#################H9################
chromHMM_H9_dNME<-chromHMM_OR(GR_merge,"AH46863","42_embryonic_stem_cell_single - H9")
chromHMM_H9_dMML<-chromHMM_OR(GR_merge,"AH46863","42_embryonic_stem_cell_single - H9",stat="dMML_pval")
#################HUES64################
chromHMM_HUES64_dNME<-chromHMM_OR(GR_merge,"AH46871","stem_27_undifferentiated_paired - HUES64")
chromHMM_HUES64_dMML<-chromHMM_OR(GR_merge,"AH46871","stem_27_undifferentiated_paired - HUES64",stat="dMML_pval")


#Part of import_ASM

if (!isEmpty(gff_in)){
  gff_in=gff_in[gff_in$Subject==subjects]
  #resize to gff ranges
  olap_gff=findOverlaps(GR,gff_in,type='within')
  GR_out=gff_in[subjectHits(olap_gff)]
  olap_gff=findOverlaps(GR,GR_out,type='within')
}else{
  GR_out=GRanges(GR)
  olap_gff=findOverlaps(GR,GR_out)
}


#Resize for each subject, subject to change
GR_resize_sub<-function(sub,gr_allele,CpG,vcf_in){
  cat(paste('Analyzing',sub,'\n',sep=' '))
  tt1=proc.time()[3]
  sub_vcf=vcf_in[[sub]]
  het_vcf=sub_vcf[sub_vcf$HetCpg]
  gr_out=gr_allele[gr_allele$Subject==sub]
  gr_out=GR_resize(gr_out,CpG,het_vcf,sub)
  #saveRDS(gr_out,paste('../downstream/temp/gr_resize_sub2',sub,'.rds',sep=''))
  cat(paste('Finishing analyzing',sub,'in',proc.time()[3]-tt1,'\n',sep=' '))
  return(gr_out)
}


allele_diff<-function(allele_gr_in){
  out_gr=GRanges()
  out_allele=GRanges()
  for(sample in unique(allele_gr_in$Sample)){
    cat('Processing:',sample,'\n')
    #Sort genome
    genome1=sort(allele_gr_in[allele_gr_in$Genome==1 & allele_gr_in$Sample==sample])
    genome2=sort(allele_gr_in[allele_gr_in$Genome==2 & allele_gr_in$Sample==sample])
    #Mroe CpG or less CpG
    genome2$CpGstat=NA
    genome1$CpGstat=NA
    cat('Identical check:',all(genome1$HetCpG==genome2$HetCpG),'\n')
    gr=granges(genome1)
    gr$N=genome1$N
    #always ref - alt, but need to be more CpG-less, if ref > alt, need ref-ale, then g1-g2 >0
    # else if ref<alt, need alt-ref, then g1-g2<0
    sign=(genome1$CpGallele-genome2$CpGallele)
    #Assign allele stat based on sign
    #For genome 2=alt, sign >0 mean ref > alt, then the alt is less
    genome2$CpGstat[sign!=0]='More'
    genome2$CpGstat[sign>0]='Less'
    genome2$alleleCpG=sign
    #Assign allele stat based on sign
    #For genome 1=ref, sign >0 mean ref > alt, then the alt is more
    genome1$alleleCpG=sign
    genome1$CpGstat[sign!=0]='Less'
    genome1$CpGstat[sign>0]='More'
    #Put sign here
    gr$CpGdiff=sign
    #Normalize sign
    sign[sign!=0]=sign[sign!=0]/abs(sign[sign!=0])
    sign[sign==0]=1
    #Use normalized sign to apply sign to the genome value
    #all is ref - alt, if sign ==1, then nCpG ref > alt, keep same
    #if sign ==-1, then nCpg alt > ref, then apply negative sign here because it should be alt -ref
    gr$diff=(genome1$Value-genome2$Value)*sign
    gr$mean=(genome1$Value+genome2$Value)/2
    #add meta information
    gr$Sample=sample
    gr$Statistic=paste('d',unique(genome1$Statistic),sep='')
    gr$HetCpG=genome1$HetCpG
    gr$pval=genome1$pval
    gr$ASM=genome1$ASM
    gr$Subject =genome1$Subject 
    #Calculating density
    gr$density =genome1$CGcount_hg19_extend
    #gr$density_diff=(genome1$density-genome2$density)*sign
    #gr$TpA_CpG_diff=(genome1$TpA_CpG-genome2$TpA_CpG)*sign
    #gr$CGcount_diff=(genome1$CGcount_allele_extend-genome2$CGcount_allele_extend)*sign
    gr$CG_nonhet_extend=genome1$CG_nonhet_extend
    #gr$CG_count_extend=genome1$CG_count_extend
    #gr$AT_count_extend=genome1$AT_count_extend
    #gr$TpA_count_extend=genome1$TpA_count_extend
    #gr$CpG_count_extend=genome1$CG_hg19_extend
    gr$gff_size_extend=genome1$gff_size_extend
    
    #Put allele information and gr information
    out_allele=c(out_allele,genome1,genome2)
    out_gr=c(out_gr,gr)
  }
  return(list(out_gr,out_allele))
}


#Add ASM information to allele values:

add_ASM<-function(allele_GR,ASM_GR){
  allele_GR_out=GRanges()
  for (sp in unique(allele_GR$Sample)){
    ASM_sub=ASM_GR[ASM_GR$Sample==sp]
    allele_sub=allele_GR[allele_GR$Sample==sp]
    allele_sub$ASM=NA
    allele_sub$pval=NA
    olap=findOverlaps(allele_sub,ASM_sub,type='equal')
    allele_sub$ASM[queryHits(olap)]=ASM_sub$ASM[subjectHits(olap)]
    allele_sub$pval[queryHits(olap)]=ASM_sub$pvalue[subjectHits(olap)]
    #allele_sub$N[queryHits(olap)]=ASM_sub$N[subjectHits(olap)]
    allele_GR_out=c(allele_GR_out,allele_sub)
  }
  return(allele_GR_out) 
}

GR_resize<-function(GR.in,CpG_sites,hetCpG,sub,gene_size=500){
  ##From definitino of CpG island, use 200 bp regions
  GR.extend=resize(GR.in,width=width(GR.in)+gene_size,fix='center')
  GR.in$CG_hg19_extend=countOverlaps(GR.extend,CpG_sites)
  #Change here for ref
  GR.in$CG_nonhet_extend=GR.in$CG_hg19_extend-countOverlaps(GR.extend,hetCpG[hetCpG$genome1_plus=='CG'|hetCpG$genome1_minus=='CG'])
  #Count CpG in genome 1
  GR.in$CG_allele_extend=NA
  GR.in$CG_allele_extend[GR.in$Genome==1]=GR.in$CG_nonhet_extend[GR.in$Genome==1]+
    countOverlaps(GR.extend[GR.extend$Genome==1],hetCpG[hetCpG$genome1_plus=='CG'|hetCpG$genome1_minus=='CG'])
  GR.in$CG_allele_extend[GR.in$Genome==2]=GR.in$CG_nonhet_extend[GR.in$Genome==2]+
    countOverlaps(GR.extend[GR.extend$Genome==2],hetCpG[hetCpG$genome2_plus=='CG'|hetCpG$genome2_minus=='CG'])
  GR.in$gff_size_extend=width(GR.extend)
  #calculate Odds ratio for expected CG vs actual CG
  #Expected CG number C * number G/total length
  # Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes". Journal of Molecular Biology.
  #Wiki, actual: ((number of C + Number of G)/2)^2/length of genomics Normalized CpG content, whole genome ~25%
  # #Check this command, no longer need TpA or AT, save processing time
  # gr_seq=getSeq(Hsapiens,GR.extend,as.character=T)
  # GR.in$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
  # GR.in$TpA_count_extend=do.call('c',lapply(gr_seq,function(x){countPattern('TA',x)}))
  # GR.in$AT_count_extend=do.call('c',lapply(gr_seq,function(x){countPattern('T',x)}))+do.call('c',lapply(gr_seq,function(x){countPattern('A',x)}))
  # GR.in$CG_count_extend=do.call('c',lapply(gr_seq,function(x){countPattern('C',x)}))+do.call('c',lapply(gr_seq,function(x){countPattern('G',x)}))
  # GR.in$CGcount_hg19_extend=GR.in$CG_hg19_extend/GR.in$CGcont_exp
  # GR.in$CGcount_allele_extend=GR.in$CG_allele_extend/GR.in$CGcont_exp
  #saveRDS(GR.in,paste('../downstream/temp/gr_resize_sub2',sub,'.rds',sep=''))
  return(GR.in)
}

#This is mainly modifed to fit the output of CPEL ASM
hetCGallele_sub<-function(sub,gr_allele,gff,CpG,vcf_in){
  cat('Analyzing',sub,'\n')
  t1=proc.time()[[3]]
  #Import vcf file
  sub_vcf=vcf_in[[sub]]
  het_vcf=sub_vcf[sub_vcf$HetCpg]
  sub_allele=gr_allele[gr_allele$Subject==sub]
  sub_het=gff[[sub]]
  #gr_allele got resized with start and end +1, use +2 to include equal start & end, no longer needed for new output?
  #sub_het=resize(sub_het, width(sub_het) + 4, fix="center")
  sub_ref=sub_allele[sub_allele$Genome==1]
  sub_alt=sub_allele[sub_allele$Genome==2]
  #Add cg number to the subject
  sub_ref=allele_hetCG(sub_ref,sub_het,'g1CG')
  sub_alt=allele_hetCG(sub_alt,sub_het,'g2CG')
  gr_out=c(sub_ref,sub_alt)
  gr_out=GR_resize(gr_out,CpG,het_vcf,gene_size=500,sub)
  cat('Finish analyzing',sub,proc.time()[[3]]-t1,'\n')
  return(gr_out)
}

#for each type of CG find gff overlap, add width and N, no longer need this complicated method 
allele_hetCG<-function(allele,gff,cgtype){
  #Find overlap region, use equal for new output
  olap=findOverlaps(allele,gff,type='equal',select="all")
  #NAs to be solved: check size of na
  olap_na=which(is.na(olap))
  cat("size of NA in olap:",length(olap_na),'\n')
  if (length(olap_na)!=0){
    allele=allele[-olap_na]
    olap=olap[-olap_na]
  }
  #separate 1 overlap and more than 1 overlap
  qt=queryHits(olap)
  #qt_table=table(qt)
  olap_sub=subjectHits(olap)
  #Find overlap regions =1
  #olap_1=as.numeric(names(qt_table[qt_table==1]))
  #for unique data, gff may have more than 1 times subsetted
  #olap_gff_1=olap_sub[which(qt%in%olap_1)]
  #allele=allele_addMeta(allele,olap_1,olap_gff_1,gff,cgtype)
  #This is for new output
  allele=allele_addMeta(allele,qt,olap_sub,gff,cgtype)
  #find overlap regions >1
  #olap_more=as.numeric(names(qt_table[qt_table>1]))
  #for data more than 1 hits, use larger region in gff
  #gff_out=unlist(lapply(olap_more,gff_min,olap_sub=olap_sub,qt=qt,gffwid=width(gff)))
  #merge unique ranges
  #allele=allele_addMeta(allele,olap_more,gff_out,gff,cgtype)
  
  return(allele)
}
#add metadata to allele GR
allele_addMeta<-function(allele,qh,sh,gff,cgtype){
  allele$CpGallele[qh]=elementMetadata(gff)[sh,cgtype]
  allele$N[qh]=gff$N[sh]
  allele$HetCpG[qh]=gff$HetCpG[sh]
  allele$N_hg19[qh]=gff$N_hg19[sh]
  allele$gff_size[qh]=width(gff[sh])
  return(allele)
}


#Generate bed files for genome browser
NME_ASM=NME_allele_ASM_calc[[1]]
NME_ASM$density_diff=NME_ASM$CpGdiff/width(NME_ASM)

#NME_allele_all_calc=readRDS('../downstream/output/dNME/NME_all_het_calc_new.rds')


#density vs dNME at HetCpG, make graph
density_df=data.frame(density=round(log10(NME_all_calc$density[NME_all_calc$density_diff!=0]),digits=3),
                      dNME=NME_all_calc$diff[NME_all_calc$density_diff!=0])
density_df_agg=aggregate(abs(density_df$dNME),by=list(density_df$density),FUN=median)
ggplot(density_df_agg,aes(x=Group.1, y=x))+
  ylim(c(0,0.25))+ggtitle("All HetCpG regions")+geom_smooth(method="lm")+
  ylab("dNME")+xlab("log10(CpG density)")+geom_point(alpha=0.5)+theme(plot.title = element_text(hjust = 0.5))
summary(lm(density_df_agg$Group.1~density_df_agg$x))#slope=-9.61189
cor.test(density_df_agg$Group.1,density_df_agg$x)#-0.72
##density vs dNME at ASM
density_df=data.frame(density=round(NME_all_calc$density[NME_all_calc$density_diff!=0 & NME_all_calc$pva<=pval_cutoff],digits=2),
                      dNME=NME_all_calc$diff[NME_all_calc$density_diff!=0& NME_all_calc$pva<=pval_cutoff])
density_df_agg=aggregate(abs(density_df$dNME),by=list(density_df$density),FUN=median)
ggplot(density_df_agg,aes(x=Group.1, y=x))+
  ylim(c(0,1))+ggtitle("All HetCpG regions")+geom_smooth(method="lm")+
  ylab("dNME")+xlab("local CpG density")+geom_point()

#########################################GO analysis for hetCpG regions#################################
#reference: within 5k of SNP
names(variant_HetCpG)=NULL
GO_ref=subsetByOverlaps(genomic_features$promoter,do.call('c',variant_HetCpG),maxgap = 2000)
write(unlist(GO_ref$gene_name),'../downstream/output/dNME/GO_ref_SNP.txt')
NME_ASM_het_sub=NME_ASM_het[NME_ASM_het$Sample %in% unique(NME_ASM_het$Sample)[1:43]]
NME_ASM_het_sub_promoter=subsetByOverlaps(genomic_features$promoter,NME_ASM_het_sub,maxgap = 1000)
write(unlist(NME_ASM_het_sub_promoter$gene_name),'../downstream/output/dNME/GO_het_ASM_all.txt')
#ESC
NME_ASM_het_sub=NME_ASM_het[NME_ASM_het$Sample %in% unique(NME_ASM_het$Sample)[1:5]]
NME_ASM_het_sub_promoter=subsetByOverlaps(genome_features$promoter,NME_ASM_het_sub,maxgap = 500)
write(unlist(NME_ASM_het_sub_promoter$gene_name),'../downstream/output/dNME/GO_het_ASM_ESC.txt')
#Differentiated
NME_ASM_het_sub=NME_ASM_het[NME_ASM_het$Sample %in% unique(NME_ASM_het$Sample)[6:43]]
NME_ASM_het_sub_promoter=subsetByOverlaps(genome_features$promoter,NME_ASM_het_sub,maxgap = 500)
write(unlist(NME_ASM_het_sub_promoter$gene_name),'../downstream/output/dNME/GO_het_ASM_diff.txt')

#ASM only
NME_ASM_sub=NME_ASM[NME_ASM$Sample %in% unique(NME_ASM$Sample)[1:43]]
NME_ASM_sub_promoter=subsetByOverlaps(genome_features$promoter,NME_ASM_sub,maxgap = 5000)
write(unlist(NME_ASM_sub_promoter$gene_name),'../downstream/output/dNME/GO_ASM_all_5k.txt')


###########################Plot enrichment at each genomic feature####################
GR_merge=readRDS("../downstream/output/gr_merge_new_promter_enhancer_density.rds")
genomic_features=readRDS('../downstream/input/genomic_features2020.rds')
GR_merge$Subject=unlist(lapply(GR_merge$Sample, function(x) strsplit(x," - ")[[1]][2]))
####Label imprinted regions
GR_merge$Group=NA
GR_merge$Group[GR_merge$Subject %in% c("H1","H9","HUES64")]="Undifferentiated"
GR_merge$Group[GR_merge$Subject %in% c("skin03","STL001","STL002","STL003","STL011")]="Differentiated"
#Find promoter regions
Imprinted_Genes <- read_excel("../downstream/input/Imprinted Genes.xlsx")
olap_promoter=findOverlaps(GR_merge,genomic_features$promoter,maxgap = 1000)
GR_merge$genes_promoter=NA
GR_merge$genes_promoter[queryHits(olap_promoter)]=genomic_features$promoter$gene_name[subjectHits(olap_promoter)]
GR_merge_promoter=GR_merge[!is.na(GR_merge$genes_promoter)]
#########Plotting odds ratio
#Islands
genome_feature_plot(GR_merge,genomic_features$`CpG island`,'dMML','CpG Island enrichment',ylim=c(0,15))
genome_feature_plot(GR_merge_promoter,genomic_features$`CpG island`,'dMML','CpG Island enrichment',ylim=c(0,20))
genome_feature_plot( GR_merge_promoter[!GR_merge_promoter$genes_promoter %in% Imprinted_Genes$Gene & GR_merge_promoter$N>1],
                     genomic_features$`CpG island`,'dMML','CpG Island enrichment',ylim=c(0,20))
#Open seas
genome_feature_plot(NME_SNP,genomic_features$`CpG open sea`,'dNME','CpG open sea enrichment')
#Test for HUES64 run 2 including boundary conditions
gr_distance(HUES64_new[which(HUES64_new$dNME_pval<=0.1)],genomic_features$`CpG island`,
            'distance to nearest CpG island (kb)','dNME ASM distance to nearest CpG island',ylim=c(0,1))
gr_distance(HUES64_new_dmml[which(HUES64_new_dmml$dMML_pval<=0.1)],genomic_features$`CpG island`,
            'distance to nearest CpG island (kb)','dMML ASM distance to nearest CpG island',ylim=c(0,1))

gr_distance(HUES64_new,genomic_features$`CpG island`,
            'distance to nearest CpG island (kb)','dNME ASM distance to nearest CpG island',ylim=c(0,1))

jpeg("../downstream/output/dMML_NME_allele_at_all_region_t1.jpeg")
ggplot(GR_merge_NME, aes(x=MML, y=NME)) + geom_smooth()+
  xlim(c(0,1))+ylim(c(0,1))+ggtitle("dMML and dNME relationship at all region")+geom_point(alpha = 0.005)
dev.off()


#####Heatmap trial###########
p2 <- ggplot(data.frame(NME=c(GR_merge$NME1,GR_merge$NME2),
                        MML=c(GR_merge$MML1,GR_merge$MML2)), aes(NME, MML)) +
  geom_pointdensity() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(axis.text = element_text(size=1), axis.title = element_text(size=1))

GR_merge_df_heat=data.frame(dMML=round(GR_merge$dMML*2,digits=2),dNME=round(GR_merge$dNME*2,digits=2))
GR_merge_agg_heat=as.data.frame(aggregate(GR_merge_df_heat,by=list(dMML=GR_merge_df_heat$dMML,dNME=GR_merge_df_heat$dNME),FUN=length))
names(GR_merge_agg_heat)=c("dMML","dNME","count1","count_normal")
GR_merge_agg_heat$dMML_total=NA
for(i in unique(GR_merge_agg_heat$dMML)){
  GR_merge_agg_heat$dMML_total[which(GR_merge_agg_heat$dMML==i)]=sum(GR_merge_agg_heat$count1[which(GR_merge_agg_heat$dMML==i)])
}
GR_merge_agg_heat$count_nromal=log(GR_merge_agg_heat$count1)/log(GR_merge_agg_heat$dMML_total)
ggplot(GR_merge_agg_heat , aes(x = dMML, y = dNME)) +
  geom_raster(aes(fill = count_nromal), interpolate=T) +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(0,1))


###########################
GR_ASM=GR[GR$pvalue<=pval_cutoff]
olap_sample=data.frame()
for(sp in unique(GR_ASM$Sample)){olap_sample=rbind(olap_sample,olap_ASM(GR_ASM[GR_ASM$Sample==sp]))}
summary_olap=colSums(olap_sample[1:43,1:10])
#dMML_dNME      dMML_UC      dNME_UC        olap3         dMML         dNME           UC dMML_nonolap dNME_nonolap   UC_nonolap 
#494         3527          494          517         9621       322840         8146         4049       320301         2574 

###############################E14 analysis############################################
E14=read.table('../Data/E14_ASM.vcf',skip=30,header = T)
colnames(E14)=c('chr','start','filter','Ref','Alt','QUAL','Het','Info','E14')
E14$end=E14$start
E14=makeGRangesFromDataFrame(E14)


#################################################Test enrichment of imprinted region in cell lines##########################################
GR_merge=readRDS("../downstream/output/gr_merge.H1.GM12878.promoter1k.rds")

genomic_features=getGeneralFeats_CpG("../downstream/input/")
saveRDS(genomic_features,"D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/genomic_features2020.rds")
GR_merge=GR_merge[!GR_merge$Sample %in% c("1 - GM12878","2 - GM12878","merged - GM12878","rep1 - H1","rep2 - H1")]
#GR_merge=GR_merge[GR_merge$N>1]

GR_mrege_islands=subsetByOverlaps(GR_merge,genomic_features$`CpG island`)
GR_merge$genes_body=NA
olap_body=findOverlaps(GR_merge,genomic_features$`gene body`)
GR_merge$genes_body[queryHits(olap_body)]=genomic_features$promoter$gene_name[subjectHits(olap_body)]

ESC_Sample=c("42_embryonic_stem_cell_single - H9","stem_27_undifferentiated_paired - HUES64",
             "ectoderm_paired - HUES64","endoerm_27_paired - HUES64","mesoderm_23_paired - HUES64",
             "rep1 - H1","rep2 - H1","merged - H1")
GR_merge_ESC=GR_merge[GR_merge$Sample %in% ESC_Sample]
GR_merge_non_ESC=GR_merge[!GR_merge$Sample %in% ESC_Sample]

GR_merge_ESC_sig=GR_merge_ESC[GR_merge_ESC$dMML_pval<=pval_cutoff]
GR_merge_ESC_sig[GR_merge_ESC_sig$genes_promoter %in% Imprinted_Genes$Gene]#23
GR_merge_non_ESC_sig=GR_merge_non_ESC[GR_merge_non_ESC$dMML_pval<=pval_cutoff]
GR_merge_non_ESC_sig[GR_merge_non_ESC_sig$genes_promoter %in% Imprinted_Genes$Gene] #68
#########OR of imprinting
cont_table_imprinting=matrix(c(23,1346-23,68,8412-68),ncol=2)
fisher.test(cont_table_imprinting)
olap_ESC=findOverlaps(GR_merge_ESC_sig,GR_merge_non_ESC_sig,type="equal")
GR_merge_ESC_sig_only=GR_merge_ESC_sig[-queryHits(olap_ESC)]
GR_merge_ESC_sig_both=GR_merge_ESC_sig[queryHits(olap_ESC)]
GR_merge_non_ESC_sig_only=GR_merge_non_ESC[subjectHits(olap_ESC)]
genes_ESC_islands=unique(GR_merge_ESC_sig_only$genes_promoter[order(GR_merge_ESC_sig_only$dMML,decreasing=T)])

genes_non_ESC_islands=unique(GR_merge_non_ESC$genes_promoter[GR_merge_non_ESC$dMML_pval<=pval_cutoff])

genes_unique_ESC=genes_ESC_islands[!genes_ESC_islands %in% genes_non_ESC_islands]
genes_unique_ESC=genes_unique_ESC[!genes_unique_ESC%in%Imprinted_Genes$Gene]
gene_promoter_islands=unique(subsetByOverlaps(GR_merge,genomic_features$`CpG island`)$genes_promoter)
write(gene_promoter_islands,"../downstream/output/gene_promoter_island.txt")

Imprinted_Genes <- read_excel("../code/downstream/input/Imprinted Genes.xlsx")
#check number of dMML-ASM in imprinted or not
sum(GR_merge_ESC$genes_promoter[GR_merge_ESC$dMML_pval<=pval_cutoff] %in% Imprinted_Genes$Gene)#69
sum(GR_merge_ESC$genes_promoter[GR_merge_ESC$dMML_pval>pval_cutoff] %in% Imprinted_Genes$Gene)#564
sum(!GR_merge_ESC$genes_promoter[GR_merge_ESC$dMML_pval<=pval_cutoff] %in% Imprinted_Genes$Gene)#3772
sum(!GR_merge_ESC$genes_promoter[GR_merge_ESC$dMML_pval>pval_cutoff] %in% Imprinted_Genes$Gene)#1867869
cont_imprinting=matrix(c(69,564,3772,1867869),nrow=2)
fisher.test(cont_imprinting)#60.58

cont_imprinting=matrix(c(69,564,204,1500),nrow=2)
fisher.test(cont_imprinting)#60.58
sum(GR_merge_non_ESC$genes_promoter[GR_merge_non_ESC$dMML_pval<=pval_cutoff] %in% Imprinted_Genes$Gene)#204
sum(GR_merge_non_ESC$genes_promoter[GR_merge_non_ESC$dMML_pval>pval_cutoff] %in% Imprinted_Genes$Gene)#1500
sum(!GR_merge_non_ESC$genes_promoter[GR_merge_non_ESC$dMML_pval<=pval_cutoff] %in% Imprinted_Genes$Gene)#15993
sum(!GR_merge_non_ESC$genes_promoter[GR_merge_non_ESC$dMML_pval>pval_cutoff] %in% Imprinted_Genes$Gene)#7200264
cont_imprinting=matrix(c(204,1500,15993,7200264),nrow=2)
fisher.test(cont_imprinting)#55
#Redo OR
ESC_dMML_gene=unique(GR_merge_ESC$genes_promoter[GR_merge_ESC$dMML_pval<=pval_cutoff])#148
ESC_non_dMML_gene=unique(GR_merge_ESC$genes_promoter[GR_merge_ESC$dMML_pval>pval_cutoff])#148
non_ESC_dMML_gene=unique(GR_merge_non_ESC$genes_promoter[GR_merge_non_ESC$dMML_pval<=pval_cutoff])#260
non_ESC_non_dMML_gene=unique(GR_merge_non_ESC$genes_promoter[GR_merge_non_ESC$dMML_pval>pval_cutoff])#260
sum(Imprinted_Genes$Gene %in% ESC_dMML_gene) #8 with nonASM region
sum(Imprinted_Genes$Gene %in% non_ESC_dMML_gene)#2 with nonASM region
sum(Imprinted_Genes$Gene %in% ESC_non_dMML_gene) #8,45 with nonASM region
sum(Imprinted_Genes$Gene %in% non_ESC_non_dMML_gene)#2,47 with nonASM region
#Ask Jordi where he got his many regions?
cont_imprinting=matrix(c(8,43,13,56),nrow=2)
fisher.test(cont_imprinting)

write(unique(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff]),"../downstream/output/dNME_sig_gene_promoter.txt")
write(unique(GR_merge$genes_promoter),"../downstream/output/all_gene_promoter.txt")

write(unique(GR_merge$genes_body[GR_merge$dNME_pval<=pval_cutoff]),"../downstream/output/dNME_sig_gene_body.txt")
write(unique(GR_merge$genes_body),"../downstream/output/all_gene_body.txt")

write(unique(GR_merge$genes_body[GR_merge$dMML_pval<=pval_cutoff ]),"../downstream/output/dMML_sig_gene_body.txt")


write(unique(GR_merge$genes_promoter[GR_merge$dMML_pval<=pval_cutoff& GR_merge$Sample %in% ESC_Sample]),"../downstream/output/dMML_sig_gene_promoter_ESC.txt")


#Find GR overlap the MAE:dMML gene promoter
MAE_GR_dMML=GR_merge[GR_merge$dMML_pval<=pval_cutoff][GR_merge$genes_promoter[GR_merge$dMML_pval<=pval_cutoff] %in% MAE]
write(unique(MAE_GR_dMML$genes_promoter),"../downstream/output/dMML_sig_gene_body_MAE.txt")
#Enrichment in gene body: gene body is enriched
MAE=MAE_BAE_data_Gimelbrant$Gene[ MAE_BAE_data_Gimelbrant$`MAE=1_BAE=0`==1]
dMML_MAE=sum(unique(GR_merge$genes_body[GR_merge$dMML_pval<=pval_cutoff]) %in% MAE)#528
nondMML_MAE=sum(unique(GR_merge$genes_body[GR_merge$dMML_pval>pval_cutoff]) %in% MAE)#3443
nondMML_nonMAE=sum(!(unique(GR_merge$genes_body[GR_merge$dMML_pval>pval_cutoff]) %in% MAE))#14766
dMML_nonMAE=sum(!(unique(GR_merge$genes_body[GR_merge$dMML_pval<=pval_cutoff])) %in% MAE)#1609
fisher.test(matrix(c(dMML_MAE,dMML_nonMAE,nondMML_MAE,nondMML_nonMAE),nrow=2))

#Enrichment in gene body: gene body is enriched
MAE=MAE_BAE_data_Gimelbrant$Gene[ MAE_BAE_data_Gimelbrant$`MAE=1_BAE=0`==1]
sum(unique(GR_merge$genes_body[GR_merge$dNME_pval<=pval_cutoff]) %in% MAE)#2150
sum(unique(GR_merge$genes_body[GR_merge$dNME_pval>pval_cutoff]) %in% MAE)#3446
sum(!unique(GR_merge$genes_body[GR_merge$dNME_pval>pval_cutoff]) %in% MAE)#14571
sum(!unique(GR_merge$genes_body[GR_merge$dNME_pval<=pval_cutoff]) %in% MAE)#8968
fisher.test(matrix(c(2150,8968,3446,14571),nrow=2))


####################Enrichment analysis using motif_break result#############################
motif_gene <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/motif_break_all_unique_SNP.rds")
motif_gene$subject=unlist(lapply(names(motif_gene),function(x) strsplit(x,"-")[[1]][1]))
#Motif break using nearby regions: within 10kb
#For a given gene, check enrichment in any ASM: about 300
#Subset by strong effect genes/SNPs
motif_gene_strong=motif_gene
length(unique(motif_gene_strong$geneSymbol))
#For a given gene,filter SNPs within the range
genomic_features=getGeneralFeats_CpG("../downstream/input/")

#Looking for which one have motif in range
# motif_gene_in_range=GRanges()
# for (i in unique(motif_gene_strong$geneSymbol)){
#   gene_in=motif_gene_strong[motif_gene_strong$geneSymbol==i]
#   motif_gene_in_range=c(motif_gene_in_range,subsetByOverlaps(gene_in,c(genomic_features$promoter[which(genomic_features$promoter$gene_name==i)],
#                              genomic_features$`gene body`[which(genomic_features$`gene body`$gene_name==i)]),maxgap = 10000))
# }
#unique(motif_gene_in_range$geneSymbol) #4 genes
subsetByOverlaps(motif_gene_strong,GR_merge[GR_merge$dNME_pval<=pval_cutoff],maxgap = 100)$geneSymbol
subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],motif_gene_strong)
subsetByOverlaps(motif_gene_strong,GR_merge[GR_merge$dMML_pval<=pval_cutoff],maxgap = 100) #HNF4A
subsetByOverlaps(GR_merge[GR_merge$dMML_pval<=pval_cutoff],motif_gene_strong)
variant_HetCpG_meta=readRDS('../downstream/output/ASM_enrich_meta.rds')
variant_HetCpG_meta=subsetByOverlaps(variant_HetCpG_meta,GR_merge,maxgap = 100)

#Find genes that are affected
GR_in_motif=subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],motif_gene_strong,maxgap = 100)
write(unique(GR_in_motif$genes_promoter),"../downstream/output/motif_gene_promoter_gene_body.txt")
#Subset Motifs that are sensitive to methylation
library(readr)
methyl_sensitive_gene=as.data.frame(read_csv("../downstream/input/Methyl_call_SELEX.csv"),stringsASfactors=F)
methyl_plus=methyl_sensitive_gene[methyl_sensitive_gene$Call=="MethylPlus",]
methyl_minus=methyl_sensitive_gene[methyl_sensitive_gene$Call=="MethylMinus",]
motif_gene_strong_Methyl_plus=motif_gene_strong[motif_gene_strong$geneSymbol %in% methyl_plus$`TF name`]
motif_gene_strong_Methyl_minus=motif_gene_strong[motif_gene_strong$geneSymbol %in% methyl_minus$`TF name`]

####################Motif enriched and annotation#########################################################################
#Find the motif that are enriched in dNME ASM
motif_enriched=readRDS("../downstream/motif_enrich_all.rds")
cat(motif_enriched$gene_name[motif_enriched$OR>1 & motif_enriched$qval<=0.1 & motif_enriched$gene_name %in% methyl_plus$`TF name`],sep=',')
cat(motif_enriched$gene_name[motif_enriched$OR>1 & motif_enriched$qval<=0.1 & motif_enriched$gene_name %in% methyl_minus$`TF name`],sep=',')

cat(motif_enriched$gene_name[motif_enriched$OR>1 & motif_enriched$qval<=0.1 & motif_enriched$gene_name %in% methyl_plus$`TF name`],sep=',')
cat(motif_enriched$gene_name[motif_enriched$OR>1 & motif_enriched$qval<=0.1 & motif_enriched$gene_name %in% methyl_minus$`TF name`],sep=',')
###########################Example regions#############################
#GR_in_motif_methyl_plus=subsetByOverlaps(variant_HetCpG_meta[which(variant_HetCpG_meta$dNME_pval<=pval_cutoff)],
#                                         motif_gene_strong_Methyl_plus[motif_gene_strong_Methyl_plus$geneSymbol %in% motif_enriched$gene_name[motif_enriched$OR>1 & motif_enriched$qval<=0.1]])
variant_HetCpG_meta_NME_sig_pos=variant_HetCpG_meta[which(variant_HetCpG_meta$dNME_pval<=pval_cutoff)]
GR_in_motif_methyl_plus_olap=findOverlaps(variant_HetCpG_meta_NME_sig_pos,motif_gene_strong_Methyl_plus)
GR_in_motif_methyl_plus_pos=variant_HetCpG_meta_NME_sig_pos[queryHits(GR_in_motif_methyl_plus_olap)]
GR_in_motif_methyl_plus_pos$alleleDiff=NA
GR_in_motif_methyl_plus_pos$alleleDiff=motif_gene_strong_Methyl_plus$alleleDiff[subjectHits(GR_in_motif_methyl_plus_olap)]
GR_in_motif_methyl_plus_pos$pctRef=NA
GR_in_motif_methyl_plus_pos$pctRef=motif_gene_strong_Methyl_plus$pctRef[subjectHits(GR_in_motif_methyl_plus_olap)]
GR_in_motif_methyl_plus_pos_sub=GR_in_motif_methyl_plus_pos[GR_in_motif_methyl_plus_pos$N>5]
#Binomial test
GR_in_motif_methyl_plus_pos=GR_in_motif_methyl_plus
binom.test( sum(GR_in_motif_methyl_plus_pos$refNME<GR_in_motif_methyl_plus_pos$altNME),
            length(GR_in_motif_methyl_plus_pos),p=sum(variant_HetCpG_meta_NME_sig$refNME<variant_HetCpG_meta_NME_sig$altNME)/length(variant_HetCpG_meta_NME_sig))#55%

GR_in_motif_methyl_plus_promoter=GR_in_motif_methyl_plus_pos[!is.na(GR_in_motif_methyl_plus_pos$genes_promoter) & GR_in_motif_methyl_plus_pos$refNME<GR_in_motif_methyl_plus_pos$altNME]
GO_promoter_plus=GO_anno(unique(GR_in_motif_methyl_plus$genes_promoter),unique(variant_HetCpG_meta$genes_promoter))
GO_promoter_plus[GO_promoter_plus$Significant>=5 & GO_promoter_plus$classicFisher<=0.05,]

GO_genebody_plus=GO_anno(unique(GR_in_motif_methyl_plus$genes_body),unique(variant_HetCpG_meta$genes_body))
GO_genebody_plus[GO_genebody_plus$Significant>=5 & GO_genebody_plus$classicFisher<=0.05,]
#######selecting regions
GR_in_motif_methyl_plus_promoter=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff & !is.na(GR_merge$genes_promoter)],motif_gene_strong_Methyl_plus)
GR_in_methyl_plus_sub=GR_in_motif_methyl_plus_promoter[GR_in_motif_methyl_plus_promoter$N>10 & 
                                                         (GR_in_motif_methyl_plus_promoter$MML1+GR_in_motif_methyl_plus_promoter$MML2)>1]
GR_in_methyl_plus_sub[order(GR_in_methyl_plus_sub$dNME,decreasing=T)]
subsetByOverlaps(motif_gene_strong_Methyl_plus,GR_in_methyl_plus_sub,maxgap = 100)$geneSymbol %in% 
  motif_enriched$gene_name[motif_enriched$OR>1 & motif_enriched$qval<=0.1 & motif_enriched$gene_name %in% methyl_plus$`TF name`]

plotMB(subsetByOverlaps(motif_gene_strong_Methyl_plus,GR_in_methyl_plus_sub[8]),
       names(subsetByOverlaps(motif_gene_strong_Methyl_plus,GR_in_methyl_plus_sub[8])))
#check motif in it
plotMB(subsetByOverlaps(motif_gene_strong_Methyl_plus,GR_in_methyl_plus_sub[2]),
       rsid=unique(names(subsetByOverlaps(motif_gene_strong_Methyl_plus,GR_in_methyl_plus_sub[2]))))

GR_in_motif_methyl_plus_pos_sub=subsetByOverlaps(GR_in_motif_methyl_plus_pos,GR_in_methyl_plus_sub)

write(unique(GR_in_motif_methyl_plus$genes_promoter),"../downstream/output/motif_gene_promoter_methyl_plus.txt")
write(unique(GR_in_motif_methyl_plus$genes_body),"../downstream/output/motif_gene_body_methyl_plus.txt")
#################Methylminus analysis##########################################################################
variant_HetCpG_meta_NME_sig_minus=variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff]
GR_in_motif_methyl_minus_olap=findOverlaps(variant_HetCpG_meta_NME_sig_minus,motif_gene_strong_Methyl_minus)
GR_in_motif_methyl_minus=variant_HetCpG_meta_NME_sig_minus[queryHits(GR_in_motif_methyl_minus_olap)]
GR_in_motif_methyl_minus$alleleDiff=NA
GR_in_motif_methyl_minus$alleleDiff=motif_gene_strong_Methyl_plus$alleleDiff[subjectHits(GR_in_motif_methyl_minus_olap)]
GR_in_motif_methyl_minus$pctRef=NA
GR_in_motif_methyl_minus$pctRef=motif_gene_strong_Methyl_plus$pctRef[subjectHits(GR_in_motif_methyl_minus_olap)]

#GR_in_motif_methyl_minus_promoter=GR_in_motif_methyl_minus[!is.na(GR_in_motif_methyl_minus$genes_promoter)]
binom.test( sum(GR_in_motif_methyl_minus$refNME<GR_in_motif_methyl_minus$altNME),
            length(GR_in_motif_methyl_minus),p=sum(variant_HetCpG_meta_NME_sig_minus$refNME<variant_HetCpG_meta_NME_sig_minus$altNME)/length(variant_HetCpG_meta_NME_sig_minus))#44%

GO_promoter_minus_neg=GO_anno(unique(GR_in_motif_methyl_minus$genes_promoter),unique(variant_HetCpG_meta$genes_promoter))
GO_promoter_minus_neg[GO_promoter_minus_neg$Significant>=10 & GO_promoter_minus_neg$classicFisher<=0.05,]

GO_body_minus_neg=GO_anno(unique(GR_in_motif_methyl_minus$genes_body),unique(variant_HetCpG_meta$genes_body))
GO_body_minus_neg[GO_body_minus_neg$Significant>=10 & GO_body_minus_neg$classicFisher<=0.05,]
#############################GO together########################
motif_enriched_sig=motif_enriched[motif_enriched$OR>1 &motif_enriched$qval<=0.1,"gene_name" ]
motif_gene_sig=motif_gene[motif_gene$geneSymbol %in% motif_enriched_sig]
GO_promoter_all=GO_anno(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff],unique(GR_merge$genes_promoter))
GO_promoter_all_motif=GO_anno(subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],motif_gene_sig)$genes_promoter,unique(GR_merge$genes_promoter))
sensory_gene=unique(select(org.Hs.eg.db, unlist(as.list(org.Hs.egGO2EG)[GO_promoter_all$GO.ID[1]]), c("SYMBOL"), "ENTREZID")$SYMBOL)
GR_merge_sensory=GR_merge[GR_merge$genes_promoter %in% sensory_gene & GR_merge$dNME_pval<=pval_cutoff]
GR_merge_sensory_N1=GR_merge_sensory[order(GR_merge_sensory$N,decreasing=T)][1:48]
GR_merge_sensory_N1[order(GR_merge_sensory_N1$Subject)]
######Checking GO in each sample#################
GO_promoter_all_sample_NME=list()

for(sp in unique(GR_merge$Sample)){
  GO_promoter_all_sample_NME[[sp]]=GO_anno(GR_merge[GR_merge$Sample==sp & GR_merge$dNME_pval<=pval_cutoff]$genes_promoter,unique(GR_merge$genes_promoter))
  
}
GO_promoter_all_sample_NME_ft=lapply(GO_promoter_all_sample_NME,function(x) x[x$Significant >=10 & x$classicFisher<=0.01,])
GO_promoter_all_sample_NME_ft_bind=do.call(rbind,GO_promoter_all_sample_NME_ft)
GO_promoter_all_sample_NME_ft_bind[GO_promoter_all_sample_NME_ft_bind$classicFisher<=0.05 &GO_promoter_all_sample_NME_ft_bind$elimFisher<=0.05 &GO_promoter_all_sample_NME_ft_bind$weightFisher<=0.05,]
########genes in each sample are natually enriched in olfactory genes#####################

GO_promoter_all_sample_null=list()

for(sp in unique(GR_merge$Sample)){
  GO_promoter_all_sample_NME[[sp]]=GO_anno(GR_merge[GR_merge$Sample==sp]$genes_promoter,unique(GR_merge[GR_merge$Sample==sp]$genes_promoter))
  
}

GO_promoter_all[GO_promoter_all$Significant>=10 & GO_promoter_all$classicFisher<=0.05,]
GR_in_motif=c(GR_in_motif_methyl_minus,GR_in_motif_methyl_plus)

GO_body_all=GO_anno(unique(c(GR_in_motif_methyl_minus$genes_body,GR_in_motif_methyl_plus$genes_body)),unique(variant_HetCpG_meta$genes_body))
GO_body_all[GO_body_all$Significant>=10 & GO_body_all$classicFisher<=0.05,]

GO_nme_all=GO_anno(unique(variant_HetCpG_meta_NME_sig$genes_promoter),unique(variant_HetCpG_meta$genes_promoter))
GO_nme_all[GO_nme_all$Significant>=10 & GO_nme_all$classicFisher<=0.05,]


#############Finding regions#################################
GR_in_motif_methyl_minus_promoter=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff & !is.na(GR_merge$genes_promoter)],motif_gene_strong_Methyl_minus)
GR_in_methyl_minus_sub=GR_in_motif_methyl_minus_promoter[GR_in_motif_methyl_minus_promoter$N>10 & 
                                                           (GR_in_motif_methyl_minus_promoter$MML1+GR_in_motif_methyl_minus_promoter$MML2)<1]
GR_in_methyl_minus_sub[order(GR_in_methyl_minus_sub$dNME,decreasing=T)]
plotMB(subsetByOverlaps(motif_gene_strong_Methyl_minus,GR_in_methyl_minus_sub[4]),
       names(subsetByOverlaps(motif_gene_strong_Methyl_minus,GR_in_methyl_minus_sub[4])))

GR_in_motif_methyl_minus_sub=subsetByOverlaps(GR_in_motif_methyl_minus,GR_in_methyl_minus_sub)
sum(sign(GR_in_motif_methyl_minus_sub$refNME-GR_in_motif_methyl_minus_sub$altNME)!=sign(GR_in_motif_methyl_minus_sub$alleleDiff))

write(unique(GR_in_motif_methyl_minus$genes_promoter),"../downstream/output/motif_gene_promoter_methyl_minus.txt")
write(unique(c(GR_in_motif_methyl_minus$genes_promoter,GR_in_motif_methyl_plus$genes_promoter)),"../downstream/output/motif_gene_promoter_all.txt")
###############################Go analysis#########################################################

write(unique(GR_in_motif_methyl_plus$genes_body),"../downstream/output/motif_gene_body_methyl_plus.txt")
write(unique(GR_in_motif_methyl_plus$genes_promoter),"../downstream/output/motif_gene_promoter_methyl_plus.txt")

write(unique(GR_in_motif_methyl_minus$genes_body),"../downstream/output/motif_gene_body_methyl_minus.txt")
write(unique(GR_in_motif_methyl_minus$genes_promoter),"../downstream/output/motif_gene_promoter_methyl_minus.txt")

#TopGO analysis

methyl_plus_promoter_GO=GO_anno(unique(GR_in_motif_methyl_plus$genes_body),unique(GR_merge$genes_promoter))

GR_merge_motif=subsetByOverlaps(GR_merge,motif_gene_strong,maxgap = 100)
write(unique(GR_merge_motif$genes_body),"../downstream/output/motif_gene_body_all.txt")
write(unique(GR_merge_motif$genes_promoter),"../downstream/output/motif_gene_promoter_all.txt")
#Motif analysis using variant HetCpG meta only

motif_gene <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/motif_break_all_unique_SNP.rds")
variant_HetCpG_meta=readRDS('../downstream/output/ASM_enrich_meta_new.rds')
#variant_out=motif_enrich(motif_gene,variant_HetCpG_meta = variant_HetCpG_meta,pval_cutoff = 0.1,dist=500)
saveRDS(variant_out,"../downstream/motif_enrich_all_SNP.rds")
variant_out_methyl_plus=variant_out$gene_name[variant_out$gene_name %in% methyl_plus$`TF name`]
variant_out_methyl_minus=variant_out$gene_name[variant_out$gene_name %in% methyl_minus$`TF name`]
sig_SNPs=variant_sig_out[variant_sig_out$OR>=2 & variant_sig_out$qval<=0.1]#182
unique(variant_sig_out$gene_name[variant_sig_out$qval>0.1])#139
unique(variant_sig_out$gene_name[variant_sig_out$OR<1 & variant_sig_out$qval<=0.1])#1
#Look for which motif that have dirctionality
#variant_out=readRDS("../downstream/motif_enrich_all_SNP.rds")

#motif_sig_df_sample=list()
#for (sp in unique(variant_HetCpG_meta$Sample)){
#motif_sig_df_sample[[sp]]=motif_enrich(motif_gene,variant_HetCpG_meta[variant_HetCpG_meta$Sample==sp],
#                                 pval_cutoff =pval_cutoff,dist=500)
#}
#saveRDS(motif_sig_df_sample,"../downstream/output/motif_sig_df_sample.rds")
##Check number of enriched motif
#sum(unlist(lapply(motif_sig_df_sample,function(x) sum(x$OR>1 & x$qval<=0.1))))#2947
#sum(unlist(lapply(motif_sig_df_sample,function(x) sum(x$OR<1 & x$qval<=0.1))))#33
#apply new CMH mmethod

###############Sample specific DNase analysis
#Match each sample to it's ENCODE
ENCODE_number_to_sample=ENCODE_to_sample(unique(variant_HetCpG_meta$Sample))
variant_HetCpG_meta$ENCODE=NA
variant_HetCpG_meta$ENCODE=ENCODE_number_to_sample$ENCODE[match(variant_HetCpG_meta$Sample,ENCODE_number_to_sample$sample)]
#Check if all ENOCDE name are in colum name
ENCODE_number_to_sample$ENCODE %in% colnames(elementMetadata(DNase_all))
#Filter regioins based on DNase region
variant_HetCpG_meta_dnase=GRanges()
for(sp in unique(variant_HetCpG_meta$Sample)){
  ENCODE=ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample==sp]
  if(!is.na(ENCODE)){
    dnase_sp=DNase_all[unlist(elementMetadata(DNase_all)[ENCODE])!=0]
    variant_HetCpG_meta_dnase=c(variant_HetCpG_meta_dnase,subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$Sample==sp],dnase_sp,maxgap=500))
  }
}
#Only 5% regions left


#Comparing different method

sum(motif_sig_df_DNAase$gene_name[motif_sig_df_DNAase$qval<=0.1&motif_sig_df_DNAase$OR>1.5] %in% motif_sig_df_all_DNaseI$gene_name[motif_sig_df_all_DNaseI$qval<=0.1&motif_sig_df_all_DNaseI$OR>1.5])

sum( motif_sig_df_500$gene_name[ motif_sig_df_500$qval<=0.1& motif_sig_df_500$OR>1] %in% motif_sig_df_all_DNaseI$gene_name[motif_sig_df_all_DNaseI$qval<=0.1&motif_sig_df_all_DNaseI$OR>1])

sum( motif_sig_df_all_DNaseI$gene_name[ motif_sig_df_all_DNaseI$qval<=0.1& motif_sig_df_all_DNaseI$OR>1] %in%motif_sig_df_500$gene_name[motif_sig_df_500$qval<=0.1&motif_sig_df_500$OR>1])

############################Looking at genes at those motif
motif_gene_ent=motif_gene[motif_gene$geneSymbol %in% unique(motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1 & motif_dir_sig$prob>0.5])]
motif_gene_less_ent=motif_gene[motif_gene$geneSymbol %in% unique(motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1 & motif_dir_sig$prob<0.5])]
#####Find example regions
motif_gene_ent_exp=subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$sub=='H1' & 
                                                          variant_HetCpG_meta$dNME_pval<=pval_cutoff &
                                                          variant_HetCpG_meta$HetCpG&
                                                          variant_HetCpG_meta$N>1 &
                                                          variant_HetCpG_meta$dMML_pval>0.1],motif_gene_less_ent)
motif_gene_ent_exp2=subsetByOverlaps(motif_gene_ent_exp[order(motif_gene_ent_exp$N* motif_gene_ent_exp$dNME,decreasing=T)],genomic_features$TSS,maxgap=5000)
subsetByOverlaps(motif_gene_ent_exp,DNase_H1)

GR_merge_H1=GR_merge[GR_merge$Subject=='H1']
gene_dNME_ent=subsetByOverlaps(GR_merge_H1[GR_merge_H1$dNME_pval<=pval_cutoff],motif_gene_ent)
#Subset by DNase I hypersensitive region, H1 E003
load("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/state_calls_enh.RData")
enhancer_E003=max_probs[,"E003"]
enhancer_E003_all=enchancer_DNase_gr[which(!is.na(enhancer_E003))]
enhancer_E003_all$score=enhancer_E003[which(!is.na(enhancer_E003))]
export.bedGraph(enhancer_E003_all,'../downstream/output/H1_enhancer.bedGraph')
load("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/state_calls_prom.RData")
prom_E003=max_probs[,"E003"]
prom_E003_all=prom_DNase_gr[which(!is.na(prom_E003))]
prom_E003_all$score=prom_E003[which(!is.na(prom_E003))]
export.bedGraph(prom_E003_all,'../downstream/output/H1_promoter.bedGraph')
DNase_H1=c(enhancer_E003_all,prom_E003_all)

gene_dNME_ent_DNase=subsetByOverlaps(gene_dNME_ent,DNase_H1)
gene_dNME_ent_DNase_select=gene_dNME_ent_DNase[gene_dNME_ent_DNase$N>1]
motif_gene_ent_exp=subsetByOverlaps(motif_gene,gene_dNME_ent_DNase_select[!is.na(gene_dNME_ent_DNase_select$genes_promoter)][8])
plotMB(motif_gene_ent_exp[3],'H1-510385')
#2,
gene_dNME_less_ent=subsetByOverlaps(GR_merge_H1[GR_merge_H1$dNME_pval<=pval_cutoff],motif_gene_less_ent)
gene_dNME_less_ent_DNase=subsetByOverlaps(gene_dNME_less_ent,DNase_H1)
gene_dNME_less_ent_DNase_select=subsetByOverlaps(gene_dNME_less_ent_DNase[gene_dNME_less_ent_DNase$N>1],variant_HetCpG_meta[variant_HetCpG_meta$HetCpg & variant_HetCpG_meta$sub=='H1'])
gene_dNME_less_ent_DNase_select_exp=gene_dNME_less_ent_DNase_select[!is.na(gene_dNME_less_ent_DNase_select$genes_promoter)][4]
motif_gene_less_ent_exp=subsetByOverlaps(motif_gene,gene_dNME_less_ent_DNase_select_exp)
plotMB(motif_gene_less_ent_exp,'HUES64-1265540')

genes_hypervar_ent=genes_hypervar[genes_hypervar$gene_name%in%gene_dNME_ent$genes_promoter,]
genes_hypervar_less_ent=genes_hypervar[genes_hypervar$gene_name%in%gene_dNME_less_ent$genes_promoter,]

plot(ecdf(genes_hypervar_ent$hypervar_log2),do.points=F,col="red", verticals=TRUE)
plot(ecdf(genes_hypervar_less_ent$hypervar_log2),add=TRUE,do.points=F, verticals=TRUE)
####Enrichment analysis
human_mono_motif_TSV=as.data.frame(read.table('../downstream/input/HUMAN_mono_motifs.tsv',sep='\t',header=T,stringsAsFactors = F))
human_mono_motif_TSV$Transcription.factor=gsub("HUMAN:","",human_mono_motif_TSV$Transcription.factor)
which(!unique(motif_gene$geneSymbol) %in% human_mono_motif_TSV$Transcription.factor)
motif_family_enrich(unique(motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1]),unique(motif_gene$geneSymbol),human_mono_motif_TSV)
more_ent_enrich=motif_family_enrich(unique(motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1 & motif_dir_sig$prob>0.5]),unique(motif_gene$geneSymbol),human_mono_motif_TSV)
less_ent_enrich=motif_family_enrich(unique(motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1& motif_dir_sig$prob<0.5]),unique(motif_gene$geneSymbol),human_mono_motif_TSV)
#All factors
SNP_all=motif_gene[motif_gene$geneSymbol %in% motif_dir_sig$TF]
SNP_all_sig=SNP_all[abs(SNP_all$alleleDiff)>=1.5 & (SNP_all$pctAlt>=0.9 | SNP_all$pctRef>=0.9)]
dNME_ASM_all=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],SNP_all_sig,maxgap = 0)
dNME_ASM_all_hyper_var=dNME_ASM_all[dNME_ASM_all$genes_promoter %in% genes_hypervar_top75]
SNP_example_all=subsetByOverlaps(SNP_all_sig,dNME_ASM_all_hyper_var[dNME_ASM_all_hyper_var$Subject=='H1'][c(2,4)])
plotMB(SNP_example_all,"STL001-705344")
plotMB(SNP_example_all[2],"H1-748882")
#EST-related factors
motif_dir_sig_ETS=motif_dir_sig[motif_dir_sig$TF %in% human_mono_motif_TSV$Transcription.factor[human_mono_motif_TSV$TF.family=='Ets-related factors{3.5.2}'],]
SNP_ETS=motif_gene[motif_gene$geneSymbol %in% motif_dir_sig_ETS$TF]
SNP_ETS_sig=SNP_ETS[abs(SNP_ETS$alleleDiff)>=1.5 & (SNP_ETS$pctAlt>=0.9 | SNP_ETS$pctRef>=0.9)]
dNME_ASM_ETS=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],SNP_ETS_sig,maxgap = 0)
dNME_ASM_ETS_prom=dNME_ASM_ETS[!is.na(dNME_ASM_ETS$genes_promoter)]
genes_hypervar=as.data.frame(read.csv("../downstream/input/H1 data/H1_hypervar_result.csv"))

genes_hypervar_top75=as.character(genes_hypervar$gene_name[genes_hypervar$hypervar_log2>=quantile(genes_hypervar$hypervar_log2,prob=0.75)])
dNME_ASM_ETS_prom_hyper=dNME_ASM_ETS_prom[dNME_ASM_ETS_prom$genes_promoter%in%genes_hypervar_top75]
dNME_SNP_ETS_sig=dNME_ASM_ETS_prom_hyper[dNME_ASM_ETS_prom_hyper$dNME>0.5 & dNME_ASM_ETS_prom_hyper$dMML_pval>pval_cutoff]
dNME_ASM_ETS_prom_ATAC=subsetByOverlaps(dNME_ASM_ETS,gr_CPM[which(gr_CPM$padj<=0.1)],maxgap = 1000)#None
dNME_ASM_ETS_prom_scATAC=subsetByOverlaps(dNME_ASM_ETS,scATAC_hyper[scATAC_hyper$hypervar_log2>=quantile(scATAC_hyper$hypervar_log2,prob=0.75)])#NoneATAC, 74 in H1

#Hyper var gene example
hyper_sig=subsetByOverlaps(SNP_ETS_sig,dNME_SNP_ETS_sig[5])
plotMB(hyper_sig[4],'STL001-1058665')
GO_ETS=GO_anno(unique(dNME_ASM_ETS$genes_promoter),unique(GR_merge$genes_promoter),topNodes = 15371)
GO_ETS_table=GO_ETS[[1]]
hyper_sig=subsetByOverlaps(SNP_ETS,dNME_ASM_ETS_prom_scATAC[!is.na(dNME_ASM_ETS_prom_scATAC$genes_promoter)][2])
hyper_sig=subsetByOverlaps(SNP_ETS_sig,dNME_ASM_ETS_prom[order(dNME_ASM_ETS_prom$N,decreasing=T)][2])
plotMB(hyper_sig[18],'HUES64-417160')
NME_inforME=import.bed('../downstream/input/H1_inforME/NME-H1_merged_all_st.bed')
#look for genomic features: 

subsetByOverlaps(NME_inforME,dNME_ASM_ETS_prom[order(dNME_ASM_ETS_prom$N,decreasing=T)][2])
quantile(as.numeric(NME_inforME$name),prob=0.90)
#in gene body only 1 that have 2 for H1
dNME_ASM_ETS_body=dNME_ASM_ETS[!is.na(dNME_ASM_ETS$genes_body)]
dNME_ASM_ETS_body_hyper=dNME_ASM_ETS_body[dNME_ASM_ETS_body$genes_body%in%genes_hypervar_top75]
dNME_SNP_ETS_sig_body=dNME_ASM_ETS_body_hyper[dNME_ASM_ETS_body_hyper$dNME>0.5 & dNME_ASM_ETS_body_hyper$dMML_pval>pval_cutoff]

#TAL-related factors:only 3
GO_obj=GO_ETS[[2]]
ETS_gene=genesInTerm(GO_obj,GO_ETS_ft$GO.ID[1])#get term related to regulation of cell death
GR_merge_EST=GR_merge[GR_merge$dNME_pval<=pval_cutoff &GR_merge$genes_promoter %in% ETS_gene$`GO:0010941` & GR_merge$dMML<0.3]
GR_merge_EST_subj=GR_merge[GR_merge$dNME_pval<=pval_cutoff &GR_merge$genes_promoter %in% ETS_gene$`GO:0010941` & GR_merge$dMML<0.3 &
                             GR_merge$Sample=='mesoderm_23_paired - HUES64']
subsetByOverlaps(motif_gene,GR_merge_EST_subj[order(GR_merge_EST_subj$N *GR_merge_EST_subj$dNME,decreasing = T)][1])


write.csv(as.data.frame(GO_ETS[GO_ETS$Significant>=5 & GO_ETS$classicFisher<=0.05,]),'../downstream/output/ETS_GO.csv')
EST_death_gene=unique(select(org.Hs.eg.db, unlist(as.list(org.Hs.egGO2EG)[GO_ETS_ft$GO.ID[1]]), c("SYMBOL"), "ENTREZID")$SYMBOL)
GR_merge[GR_merge$genes_promoter %in% EST_death_gene]
#This is not working, use all motif and try again
dNME_ASM_ETS_H1=dNME_ASM_ETS[dNME_ASM_ETS$Subject=='H1']
genes_hypervar=as.data.frame(read.csv("../downstream/input/H1 data/H1_hypervar_result.csv"))
genes_hypervar_top25=genes_hypervar$gene_name[genes_hypervar$hypervar_log2>=0.7]
dNME_ASM_ETS_H1_hyper_var=dNME_ASM_ETS_H1[dNME_ASM_ETS_H1$genes_promoter %in% genes_hypervar_top25]
#Looking for TF that have methylation preference


#############################Find regions that have ASM events in multiple samples#################################################
GR_merge=readRDS("../downstream/output/gr_merge.H1.GM12878.promoter1k.rds")
#GR_merge$Subject=unlist(lapply(GR_merge$Sample,function(x) strsplit(x,' - ')[[1]][2]))
#GR_merge$tissue=unlist(lapply(GR_merge$Sample,function(x) strsplit(x,' - ')[[1]][1]))
saveRDS(GR_merge,"../downstream/output/gr_merge.H1.GM12878.promoter1k.rds")
GR_merge=GR_merge[!GR_merge$Sample %in% c("rep1 - H1","rep2 - H1","1 - GM12878","2 - GM12878","merged - GM12878")]
#GR_merge=GR_merge[GR_merge$N>1]
#Subset by common tissue
GR_merge_STL=GR_merge[GR_merge$Subject %in% c("STL001","STL002","STL003")]
Sample_STL=unique(GR_merge_STL$Sample)
tissue_STL=unlist(lapply(Sample_STL,function (x) strsplit(x,' - ')[[1]][1]))
tissue_shared=names(table(tissue_STL)[table(tissue_STL)==3]) #Gastric, Psoas_Muscle,Small_Intestine,Spleen
#Find share tissue between samples
GR_merge_STL_tissue=GR_merge_STL[GR_merge_STL$tissue %in% tissue_shared]
GR_merge_STL_gr=unique(granges(GR_merge_STL_tissue))
GR_merge_STL_gr$gastric_all=countOverlaps(GR_merge_STL_gr,GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Gastric_single"],type="equal")
GR_merge_STL_gr$gastric_dNME_ASM=countOverlaps(GR_merge_STL_gr,
                                               GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Gastric_single" & GR_merge_STL_tissue$dNME_pval<=pval_cutoff],
                                               type="equal")

GR_merge_STL_gr$Psoas_Muscle_all=countOverlaps(GR_merge_STL_gr,GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Psoas_Muscle_single"],type="equal")
GR_merge_STL_gr$Psoas_Muscle_dNME_ASM=countOverlaps(GR_merge_STL_gr,
                                                    GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Psoas_Muscle_single" & GR_merge_STL_tissue$dNME_pval<=pval_cutoff],
                                                    type="equal")

GR_merge_STL_gr$Small_Intestine_all=countOverlaps(GR_merge_STL_gr,GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Small_Intestine_single"],type="equal")
GR_merge_STL_gr$Small_Intestine_dNME_ASM=countOverlaps(GR_merge_STL_gr,
                                                       GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Small_Intestine_single" & GR_merge_STL_tissue$dNME_pval<=pval_cutoff],
                                                       type="equal")

GR_merge_STL_gr$Spleen_all=countOverlaps(GR_merge_STL_gr,GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Spleen_single"],type="equal")
GR_merge_STL_gr$Spleen_dNME_ASM=countOverlaps(GR_merge_STL_gr,
                                              GR_merge_STL_tissue[GR_merge_STL_tissue$tissue=="Spleen_single" & GR_merge_STL_tissue$dNME_pval<=pval_cutoff],
                                              type="equal")

gastric_GR_merge=GR_merge_STL_gr[which(GR_merge_STL_gr$gastric_dNME_ASM==3)]
PM_GR_merge=GR_merge_STL_gr[which(GR_merge_STL_gr$Psoas_Muscle_dNME_ASM==3)]
SI_GR_merge=GR_merge_STL_gr[which(GR_merge_STL_gr$Small_Intestine_dNME_ASM==3)]
SP_GR_merge=GR_merge_STL_gr[which(GR_merge_STL_gr$Spleen_dNME_ASM==3)]

#Subset by sample frequency
#GR_merge_gr=granges(unique(GR_merge))
# GR_merge_gr$olap_all=countOverlaps(GR_merge_gr,GR_merge)
# GR_merge_gr$olap_dMML_ASM=countOverlaps(GR_merge_gr,GR_merge[GR_merge$dMML_pval<=pval_cutoff],type="equal")
# GR_merge_gr$olap_dNME_ASM=countOverlaps(GR_merge_gr,GR_merge[GR_merge$dNME_pval<=pval_cutoff],type="equal")
# GR_merge_gr$olap_dNME_ASM_freq=GR_merge_gr$olap_dNME_ASM/GR_merge_gr$olap_all
# GR_merge_gr_multiple=GR_merge_gr[GR_merge_gr$olap_dNME_ASM_freq>=0.5 &GR_merge_gr$olap_dNME_ASM>2]
# GR_merge_multiple=subsetByOverlaps(GR_merge,GR_merge_gr_multiple)
# GR_merge_multiple[GR_merge_multiple$dNME_pval<=pval_cutoff]
#NME_all_calc_raw= readRDS('../downstream/output/dNME/NME_all_het_calc.rds')
NME_all_calc=readRDS('../downstream/output/dNME/NME_all_het_calc.rds')[[1]]
#NME_all_calc=subsetByOverlaps(NME_all_calc_raw[[1]],GR_merge_multiple[GR_merge_multiple$dNME_pval<=pval_cutoff])
#log10(density diff) vs dNME
density_df=data.frame(density_diff=log10(round(NME_all_calc$density_diff[NME_all_calc$density_diff!=0],digits=4)),
                      dNME=NME_all_calc$diff[NME_all_calc$density_diff!=0])
density_df_agg=aggregate(abs(density_df$dNME),by=list(density_df$density_diff),FUN=median)
ggplot(density_df_agg,aes(x=Group.1, y=x))+
  ylim(c(0,0.3))+ggtitle("dNME-ASM regions")+geom_smooth(method="lm")+
  ylab("dNME")+xlab("log10(density difference)")+geom_point()



# write(unique(GR_merge_multiple$genes_promoter[GR_merge_multiple$dNME_pval<=pval_cutoff]),"../downstream/output/mutlipe_gr_sig_gene_promoter.txt")
# 
# 
# GR_merge_multiple$genes_promoter=NA
# olap_promoter=findOverlaps(GR_merge_multiple,genomic_features$promoter,maxgap = 0)
# GR_merge_multiple$genes_promoter[queryHits(olap_promoter)]=genomic_features$promoter$gene_name[subjectHits(olap_promoter)]
# write(unique(GR_merge_multiple$genes_promoter[GR_merge_multiple$dNME_pval<=pval_cutoff]),"../downstream/output/mutlipe_gr_sig_gene_promoter.txt")
# 
# GR_merge_gr_multiple$genes_body=NA
# olap_body=findOverlaps(GR_merge_gr_multiple,genomic_features$`gene body`)
# GR_merge_gr_multiple$genes_body[queryHits(olap_body)]=genomic_features$`gene body`$gene_name[subjectHits(olap_body)]
# write(unique(GR_merge_multiple$genes_body[GR_merge_multiple$dNME_pval<=pval_cutoff]),"../downstream/output/mutlipe_gr_sig_gene_body.txt")
# 
# GR_merge$genes_promoter=NA
# olap_promoter=findOverlaps(GR_merge,genomic_features$promoter,maxgap = 0)
# GR_merge$genes_promoter[queryHits(olap_promoter)]=genomic_features$promoter$gene_name[subjectHits(olap_promoter)]
# 
# GR_merge$genes_body=NA
# olap_body=findOverlaps(GR_merge,genomic_features$`gene body`)
# GR_merge$genes_body[queryHits(olap_body)]=genomic_features$`gene body`$gene_name[subjectHits(olap_body)]
###################################################################
GR_merge=readRDS("../downstream/output/gr_merge.H1.GM12878.promoter1k.rds")
GR_merge=GR_merge[!GR_merge$Sample %in% c("rep1 - H1","rep2 - H1","1 - GM12878","2 - GM12878","merged - GM12878")]
dMML_ASM_unique=unique(GR_merge[GR_merge$dMML_pval<=pval_cutoff])#7836-3272, total 9758 regions
dNME_ASM_unique=unique(GR_merge[GR_merge$dNME_pval<=pval_cutoff])#215184-3272, total 344502
subsetByOverlaps(dMML_ASM_unique,dNME_ASM_unique,type="equal") #3272
subsetByOverlaps(dMML_ASM_unique,GR_merge_gr_multiple,type="equal") #198

#repetitive element: have a lot of overlap between repeat masker
#Without het CpG
NME_all=readRDS('../downstream/output/dNME/NME_all_het_calc.rds')
NME_all=NME_all[[1]]
NME_all$ASM=NA
NME_all$ASM[NME_all$pva<=pval_cutoff]="Yes"
NME_all$ASM[NME_all$pva>pval_cutoff]="No"
NME_all=NME_all[!is.na(NME_all$ASM) &!is.na(NME_all$HetCpG)]
rep=import.bed('../downstream/input/simple_repeats.bed')
rep=import.bed('../downstream/input/hg19_repeat_masker.bed')
subjects=unique(NME_all$Sample)
OR_repeats=data.frame(subject=subjects,OR=0,CpG_type='simple_repeats')
for(subj in subjects){
  OR_repeats$OR[OR_repeats$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Sample==subj],rep,'dNME',maxgap = 1000)$estimate
} 
########mouse vs human#######################################
############Looking at enrichment in those? HUMAN OLFACTORY RECEPTORS: NOVEL CELLULAR FUNCTIONS OUTSIDE OF THE NOSE
gene_shared_tissue=unlist(read.table("../downstream/input/olfactory_gene.txt",stringsAsFactors = F))
sum(unlist(gene_shared_tissue) %in% GR_merge_sensory$genes_promoter) #18/40
length(unique(GR_merge_sensory$genes_promoter)) #175, OR=1.3 not sig
GR_merge_sensory_shared_tissue=GR_merge_sensory[GR_merge_sensory$genes_promoter %in% gene_shared_tissue]

###################################Looking for wheather gene body enrichment are in motif######################
GR_merge=readRDS("../downstream/output/gr_merge.H1.GM12878.promoter1k.rds")
motif_out=readRDS("../downstream/output/motif_break_all_unique_SNP.rds")
GR_merge_gene_body=GR_merge[!is.na(GR_merge$genes_body) & GR_merge$dNME_pval<=pval_cutoff]
dNME_GO_gene_body=GO_anno(unique(GR_merge$genes_body[GR_merge$dNME_pval<=pval_cutoff]),unique(GR_merge$genes_body))
dNME_GO_gene_body[order(dNME_GO_gene_body$elimFisher),]
dNME_GO_gene_body_sig_all=dNME_GO_gene_body[as.numeric(dNME_GO_gene_body$classicFisher)<=0.01 &
                                              as.numeric(dNME_GO_gene_body$elimFisher<=0.01) &as.numeric(dNME_GO_gene_body$weightFisher)<=0.01,]
GO_sig_gene_all=unique(select(org.Hs.eg.db, unlist(as.list(org.Hs.egGO2EG)[dNME_GO_gene_body_sig_all$GO.ID]), c("SYMBOL"), "ENTREZID")$SYMBOL)
GR_merge_gene_body #82107,128044
GR_merge[GR_merge$dNME_pval<=pval_cutoff] #215184,344502
#non gene body ASM = 344502-128044=216458
#gene body non ASM:1996578
#non gene body non ASM: 2857342
fisher.test(matrix(c(128044,1996578,216458,2857342),nrow=2)) #OR = 0.84
GR_merge_sig_gene_body=GR_merge[GR_merge$genes_body %in% GO_sig_gene_all &GR_merge$dNME_pval<=pval_cutoff]
GR_merge_sig_gene_body_motif=subsetByOverlaps(GR_merge_sig_gene_body,motif_out,maxgap=100)
#Motif_analysis
library(readr)
methyl_sensitive_gene=as.data.frame(read_csv("../downstream/input/Methyl_call_SELEX.csv"),stringsASfactors=F)
methyl_plus=methyl_sensitive_gene[methyl_sensitive_gene$Call=="MethylPlus",]
methyl_minus=methyl_sensitive_gene[methyl_sensitive_gene$Call=="MethylMinus",]
motif_out_methyl_plus=motif_out[motif_out$geneSymbol %in% methyl_plus$`TF name`]
motif_out_methyl_minus=motif_out[motif_out$geneSymbol %in% methyl_minus$`TF name`]
GR_merge_sig_gene_body_motif=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],motif_out,maxgap=100)
GR_merge_sig_gene_body_motif_plus=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],motif_out_methyl_plus,maxgap=100)
GR_merge_sig_gene_body_motif_minus=subsetByOverlaps(GR_merge[GR_merge$dNME_pval<=pval_cutoff],motif_out_methyl_minus,maxgap=100)
dNME_GO_gene_body_motif=GO_anno(unique(GR_merge_sig_gene_body_motif$genes_body),unique(GR_merge$genes_body))
dNME_GO_gene_body_motif[as.numeric(dNME_GO_gene_body_motif$elimFisher)<=0.01 & 
                          as.numeric(dNME_GO_gene_body_motif$classicFisher)<=0.01 & as.numeric(dNME_GO_gene_body_motif$weightFisher)<=0.01,]


dNME_GO_gene_body_plus=GO_anno(unique(GR_merge_sig_gene_body_motif_plus$genes_body),unique(GR_merge$genes_body))
dNME_GO_gene_body_plus[as.numeric(dNME_GO_gene_body_plus$elimFisher)<=0.01 & 
                         as.numeric(dNME_GO_gene_body_plus$classicFisher)<=0.01 & as.numeric(dNME_GO_gene_body_plus$weightFisher)<=0.01,][1:10,]

dNME_GO_gene_body_minus=GO_anno(unique(GR_merge_sig_gene_body_motif_minus$genes_body),unique(GR_merge$genes_body))
dNME_GO_gene_body_minus[as.numeric(dNME_GO_gene_body_minus$elimFisher)<=0.01 & 
                          as.numeric(dNME_GO_gene_body_minus$classicFisher)<=0.01 & as.numeric(dNME_GO_gene_body_minus$weightFisher)<=0.01,][1:10,]


################distance between dMML and dNME########################################################### pretty flat
dNME_dist_dMML=GRanges()
for (subj in unique(GR_merge$Sample)){
  dNME_dist_dMML=c(dNME_dist_dMML,dist_calc(GR_merge[GR_merge$dNME_pval<=pval_cutoff & GR_merge$Sample ==  subj],
                                            GR_merge[GR_merge$dMML_pval<=pval_cutoff & GR_merge$Sample ==  subj]))
}
dNME_dist_dMML[dNME_dist_dMML$dist<=1000] #within 5kb 171317
GR_merge[GR_merge$dNME_pval<=pval_cutoff] #344502
non_dNME_dist_dMML=GRanges()
for (subj in unique(GR_merge$Sample)){
  non_dNME_dist_dMML=c(non_dNME_dist_dMML,dist_calc(GR_merge[GR_merge$dNME_pval>pval_cutoff & GR_merge$Sample ==  subj],
                                                    GR_merge[GR_merge$dMML_pval<=pval_cutoff & GR_merge$Sample ==  subj]))
}

non_dNME_dist_dMML[non_dNME_dist_dMML$dist<=1000] #within 5kb 2410393
GR_merge[GR_merge$dNME_pval>pval_cutoff] #4853920

fisher.test(matrix(c(171317,344502-171317,2410393,4853920-2410393),nrow=2))#OR=1.002

dNME_dist_tb=table(dNME_dist_dMML$dist_round)
dNME_plot_df=data.frame(dist=as.numeric(names(dNME_dist_tb)),percent_ASM=dNME_dist_tb/length(dNME_dist_dMML))
plot(dNME_plot_df$dist,dNME_plot_df$percent_ASM.Freq,pch=1,cex=0.8,ylab='Proportion of ASM',ylim=c(0,0.001),xlim=c(-10,10),xlab="distance to nearest dMML")

#########CTCF binding, motif analysis should be per subject basis, need to subset from variant HetCpG and do analysis
CTCF=motif_gene[motif_gene$geneSymbol=='CTCF']
variant_HetCpG_meta_CTCF=subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$Statistic=='dNME'],CTCF)
GR_merge_CTCF=GRanges()
for (subj in unique(variant_HetCpG_meta_CTCF$Subject)){
  variant_HetCpG_meta_CTCF_subj=variant_HetCpG_meta_CTCF[variant_HetCpG_meta_CTCF$Subject==subj & variant_HetCpG_meta_CTCF$Statistic=='dNME']
  olap_CTCF_subj=findOverlaps(variant_HetCpG_meta_CTCF_subj,CTCF)
  variant_HetCpG_meta_CTCF_subj$alleleDiff[queryHits(olap_CTCF_subj)]=unlist(CTCF$alleleDiff[subjectHits(olap_CTCF_subj)])
  olap_GR=findOverlaps(GR_merge[GR_merge$Subject==subj],variant_HetCpG_meta_CTCF_subj,maxgap = 100)
  GR_merge_CTCF_subj=GR_merge[GR_merge$Subject==subj][queryHits(olap_GR)]
  GR_merge_CTCF_subj$alleleDiff=variant_HetCpG_meta_CTCF_subj$alleleDiff[subjectHits(olap_GR)]
  GR_merge_CTCF=c(GR_merge_CTCF,GR_merge_CTCF_subj)
  
}

#####################Simulate 300 bp up and downstream of TSS



##############Checking output from informME and CPEL##################
informME_NME=import.bed('../downstream/input/H1_inforME/NME-H1_merged_all_st.bed')
CPEL_NME=import.bedGraph('../downstream/input/H1_inforME/H1_allele_agnostic_nme.bedGraph')
olap=findOverlaps(CPEL_NME,informME_NME)
CPEL_NME$informME_NME=NA
CPEL_NME$informME_MML=NA
CPEL_NME[queryHits(olap)]$informME_NME=as.numeric(informME_NME$name[subjectHits(olap)])
CPEL_NME=CPEL_NME[!is.na(CPEL_NME$informME_NME) & !is.na(CPEL_NME$score)]

informME_MML=import.bed('../downstream/input/H1_inforME/MML-H1_merged_all_st.bed')
GR_merge_H1=GR_merge[GR_merge$Subject=='H1']
olap=findOverlaps(GR_merge_H1,informME_NME)
GR_merge_H1$informME_NME=NA
GR_merge_H1$informME_MML=NA
GR_merge_H1[queryHits(olap)]$informME_NME=as.numeric(informME_NME$name[subjectHits(olap)])
olap=findOverlaps(GR_merge_H1,informME_MML)
GR_merge_H1[queryHits(olap)]$informME_MML=as.numeric(informME_MML$name[subjectHits(olap)])
GR_merge_H1$NME_mean=(GR_merge_H1$NME1+GR_merge_H1$NME2)/2
GR_merge_H1$MML_mean=(GR_merge_H1$MML1+GR_merge_H1$MML2)/2
GR_merge_H1_sub=GR_merge_H1[!is.na(GR_merge_H1$informME_NME) &!is.na(GR_merge_H1$informME_MML)]
cor(GR_merge_H1_sub$informME_NME,GR_merge_H1_sub$NME_mean)#0.488
cor(GR_merge_H1_sub$informME_NME[GR_merge_H1_sub$dNME_pval<=pval_cutoff],GR_merge_H1_sub$NME_mean[GR_merge_H1_sub$dNME_pval<=pval_cutoff])#0.2
cor(GR_merge_H1_sub$informME_NME[GR_merge_H1_sub$dNME_pval>pval_cutoff],GR_merge_H1_sub$NME_mean[GR_merge_H1_sub$dNME_pval>pval_cutoff])#0.496
#MML
cor(GR_merge_H1_sub$informME_MML,GR_merge_H1_sub$MML_mean)#0.74
cor(GR_merge_H1_sub$informME_MML[GR_merge_H1_sub$dMML_pval<=pval_cutoff],GR_merge_H1_sub$MML_mean[GR_merge_H1_sub$dMML_pval<=pval_cutoff])#0.585
cor(GR_merge_H1_sub$informME_MML[GR_merge_H1_sub$dMML_pval>pval_cutoff],GR_merge_H1_sub$MML_mean[GR_merge_H1_sub$dMML_pval>pval_cutoff])#0.740
GR_merge_H1_ggplot=data.frame(NME=GR_merge_H1_sub$NME_mean,MML=GR_merge_H1_sub$MML_mean,dNME_pval=GR_merge_H1_sub$dNME_pval,
                              dMML_pval=GR_merge_H1_sub$dMML_pval,inf_NME=GR_merge_H1_sub$informME_NME,inf_MML=GR_merge_H1_sub$informME_MML)
GR_merge_H1_ggplot$dNME_ASM="non-ASM"
GR_merge_H1_ggplot$dNME_ASM[GR_merge_H1_ggplot$dNME_pval<=pval_cutoff]="ASM"
GR_merge_H1_ggplot$dNME_ASM=as.factor(GR_merge_H1_ggplot$dNME_ASM)
GR_merge_H1_ggplot$dMML_ASM="non-ASM"
GR_merge_H1_ggplot$dMML_ASM[GR_merge_H1_ggplot$dMML_pval<=pval_cutoff]="ASM"
GR_merge_H1_ggplot$dMML_ASM=as.factor(GR_merge_H1_ggplot$dMML_ASM)
ggplot(GR_merge_H1_ggplot,aes(x=NME,y=inf_NME,color=dNME_ASM))+geom_point(alpha=0.01)+xlab('CPEL')+ylab('InformME')+
  ggtitle("NME comparison between pipeline")+theme(legend.position="bottom")+scale_color_manual(values=c("red", "blue"))
ggplot(GR_merge_H1_ggplot,aes(x=MML,y=inf_MML,color=dMML_ASM))+geom_point(alpha=0.01)+xlab('CPEL')+ylab('InformME')+
  ggtitle("MML comparison between pipeline")+theme(legend.position="bottom")+scale_color_manual(values=c("red", "blue"))


#######Hypervaribility analysis draft###################
GR_merge_H1=GR_merge[GR_merge$Sample=='merged - H1']
GR_merge_H1$mean_NME=(GR_merge_H1$NME1+GR_merge_H1$NME2)/2
GR_merge_H1$mean_MML=(GR_merge_H1$MML1+GR_merge_H1$MML2)/2

#quantiles within 1k of promoter
genomic_features=readRDS("../downstream/input/genomic_features2020.rds")
#Assign gene body and promoter
GR_merge_H1$genes_promoter=genomic_features$promoter$gene_name[findOverlaps(GR_merge_H1,genomic_features$promoter,maxgap = 1000,select="first")]
GR_merge_H1$genes_body=genomic_features$`gene body`$gene_name[findOverlaps(GR_merge_H1,genomic_features$`gene body`,maxgap = 3000,select="first")]
GR_merge_H1$hyper_var_promoter=genes_hypervar$hypervar_log2[match(GR_merge_H1$genes_promoter,genes_hypervar$gene_name)]
GR_merge_H1$hyper_var_gene_body=genes_hypervar$hypervar_log2[match(GR_merge_H1$genes_body,genes_hypervar$gene_name)]
#Null distribution
plot(density(genes_hypervar$hypervar_log2),xlab='log2(hypervaribility)',main="null distribution of hypervaribility")
plot(density(GR_merge_H1$mean_NME),xlab='mean NME',main="null distribution of mean NME")
abline(v=quantile(GR_merge_H1$mean_NME,prob=c(0.9,0.95,0.975)),col=c("black","red","blue"))

density_plot_hyper(GR_merge_H1,genes_hypervar,genomic_features)
density_plot_hyper(subsetByOverlaps(GR_merge_H1,genomic_features$`CpG open sea`),genes_hypervar,genomic_features)
density_plot_hyper(subsetByOverlaps(GR_merge_H1,genomic_features$`CpG shelf`),genes_hypervar,genomic_features)
density_plot_hyper(subsetByOverlaps(GR_merge_H1,genomic_features$`CpG shore`),genes_hypervar,genomic_features)
density_plot_hyper(subsetByOverlaps(GR_merge_H1,genomic_features$`CpG island`),genes_hypervar,genomic_features)
#density_plot_hyper(GR_merge_H1[GR_merge_H1$mean_MML<=0.75 & GR_merge_H1$mean_MML>=0.25],genes_hypervar,genomic_features)
density_plot_hyper(GR_merge_H1[GR_merge_H1$dNME_pval<=pval_cutoff],genes_hypervar,genomic_features)

#Hyper varibility distribution at different level of NME
GR_merge_H1=dist_calc(GR_merge_H1,genomic_features$promoter,k_round=100)
NME_quant=quantile(GR_merge$mean_NME,prob=c(0.25,0.5,0.75))
GR_merge_H1_025=GR_merge_H1[GR_merge_H1$gene %in% genes_hypervar$gene_name[genes_hypervar$hypervar_log2<=hyper_varibility_quant[1]]] #322
NME_plot_promoter=data.frame(NME=GR_merge_H1_025$mean_NME,dNME=GR_merge_H1_025$dNME,quant="0-25%",dist=GR_merge_H1_025$dist_round)
GR_merge_H1_2550=GR_merge_H1[
  GR_merge_H1$genes_promoter %in% genes_hypervar$gene_name[(genes_hypervar$hypervar_log2>hyper_varibility_quant[1]& genes_hypervar$hypervar_log2<=hyper_varibility_quant[2])]] #322
NME_plot_promoter=rbind(NME_plot_promoter,data.frame(NME=GR_merge_H1_2550$mean_NME,dNME=GR_merge_H1_2550$dNME,quant="25%-50%",dist=GR_merge_H1_2550$dist))
GR_merge_H1_5075=GR_merge_H1[GR_merge_H1$genes_promoter %in% genes_hypervar$gene_name[genes_hypervar$hypervar_log2>hyper_varibility_quant[2] & 
                                                                                        genes_hypervar$hypervar_log2<=hyper_varibility_quant[3]]]#119
NME_plot_promoter=rbind(NME_plot_promoter,data.frame(NME=GR_merge_H1_5075$mean_NME,dNME=GR_merge_H1_5075$dNME,quant="50% -75%",dist=GR_merge_H1_5075$dist))
GR_merge_H1_75=GR_merge_H1[GR_merge_H1$genes_promoter %in% genes_hypervar$gene_name[genes_hypervar$hypervar_log2>hyper_varibility_quant[3]]] #44
NME_plot_promoter=rbind(NME_plot_promoter,data.frame(NME=GR_merge_H1_75$mean_NME,dNME=GR_merge_H1_75$dNME,quant="75%-100%",dist=GR_merge_H1_75$dist))


#quantiles within gene body
NME_plot_body=data.frame(NME=c(GR_merge_H1$NME1,GR_merge_H1$NME2),quant="NULL")
GR_merge_H1_09=GR_merge_H1[GR_merge_H1$genes_body %in% genes_hypervar$gene_name[genes_hypervar$hypervar_log2>=hyper_varibility_quant[1]]]
NME_plot_body=rbind(NME_plot_body,data.frame(NME=c(GR_merge_H1_09$NME1,GR_merge_H1_09$NME2),quant="90% quantile"))
GR_merge_H1_095=GR_merge_H1[GR_merge_H1$genes_body %in% genes_hypervar$gene_name[genes_hypervar$hypervar_log2>=hyper_varibility_quant[2]]]
NME_plot_body=rbind(NME_plot_body,data.frame(NME=c(GR_merge_H1_095$NME1,GR_merge_H1_095$NME2),quant="95% quantile"))
GR_merge_H1_0975=GR_merge_H1[GR_merge_H1$genes_body %in% genes_hypervar$gene_name[genes_hypervar$hypervar_log2>=hyper_varibility_quant[3]]]
NME_plot_body=rbind(NME_plot_body,data.frame(NME=c(GR_merge_H1_0975$NME1,GR_merge_H1_0975$NME2),quant="97.5% quantile"))
ggplot(NME_plot_body,aes(x=NME,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+ggtitle("NME distribution within gene body")

#Finding examples
#High NME and high varibility
GR_merge_H1_var_sorted=GR_merge_H1[order(GR_merge_H1$hyper_var_promoter,decreasing=T)]
hist(GR_merge_H1_var_sorted$NME_avg[1:500])
hist(GR_merge_H1_var_sorted$dNME[1:100])
hist(GR_merge_H1_var_sorted$dNME)
GR_merge_H1_var_sorted[GR_merge_H1_var_sorted$dNME_pval<=pval_cutoff][1:100]
####Finding olfactory genes
dNME_promo_GO=GO_anno(unique(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff]),unique(GR_merge$genes_promoter))
sensory_gene=unique(select(org.Hs.eg.db, unlist(as.list(org.Hs.egGO2EG)[dNME_promo_GO$GO.ID[1]]), c("SYMBOL"), "ENTREZID")$SYMBOL)
GR_merge[GR_merge$genes_promoter %in% sensory_gene & GR_merge$dNME_pval<=pval_cutoff & !is.na(GR_merge$hyper_var_promoter)]
############Genes have both high NME and high variance?
GR_merge_sig_var=GR_merge_H1[which(GR_merge_H1$dNME_pval<=pval_cutoff & GR_merge_H1$N>=3)]
GR_merge_sig_var=GR_merge_sig_var[order(GR_merge_sig_var$hyper_var_promoter,decreasing = T)]
subsetByOverlaps(GR_merge_sig_var,motif_gene)
subsetByOverlaps(motif_gene,GR_merge_sig_var[which(GR_merge_sig_var$genes_promoter=='SEC22B')])
GR_merge_hyper_var=GO_anno(GR_merge_H1$genes_promoter[which(GR_merge_H1$dNME_pval<=pval_cutoff & GR_merge_H1$hyper_var_promoter>=2)],GR_merge_H1$genes_promoter)
GR_merge_hyper_var[GR_merge_hyper_var$Significant>5,]


###########Temporarily check HUES64######################
HUES64_new_diff=read.diffGR('HUES64','stem_27_undifferentiated_paired','../downstream/data/bedGraph_diff_new/')
HUES64_new_allele=read.alleleGR('HUES64','stem_27_undifferentiated_paired','../downstream/data/bedGraph_diff_new/')
HUES64_merge=stat_merge(HUES64_new_diff,HUES64_new_allele)
HUES64_merge=hetCGallele_merged("HUES64",HUES64_merge,hetCpG_gff,CpG_hg19,variant_HetCpG)
HUES64_new_df=elementMetadata(HUES64_merge)
GR_merge_HUES64_old=GR_merge[GR_merge$Sample=='stem_27_undifferentiated_paired - HUES64']
HUES64_old_df=elementMetadata(GR_merge_HUES64_old)
cutoff_N_min=data.frame()
for (i in unique(HUES64_new_df$N)){
  cutoff_out=data.frame(dNME_new=min(HUES64_new_df$dNME[HUES64_new_df$dNME_pval<=pval_cutoff & HUES64_new_df$N==i]),
                        dNME_old=min(HUES64_old_df$dNME[HUES64_old_df$dNME_pval<=pval_cutoff & HUES64_old_df$N==i]),
                        dMML_new=min(HUES64_new_df$dMML[HUES64_new_df$dMML_pval<=pval_cutoff & HUES64_new_df$N==i]),
                        dMML_old=min(HUES64_old_df$dMML[HUES64_old_df$dMML_pval<=pval_cutoff & HUES64_old_df$N==i]),
                        N=i)
  cutoff_N_min=rbind(cutoff_N_min,cutoff_out)
}
cutoff_N_min[order(cutoff_N_min$N),]
cutoff_N_min_mt_dNME=melt(cutoff_N_min[,c('dNME_old','dNME_new','N')],id='N',value.name = 'dNME_cutoff',variable.name = 'Stat')
ggplot(cutoff_N_min_mt_dNME,aes(x=N,y=dNME_cutoff,fill=Stat))+ geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=round(dNME_cutoff,digits=3)), position=position_dodge(1))+theme(legend.position = 'bottom')

cutoff_N_min_mt_dMML=melt(cutoff_N_min[,c('dMML_old','dMML_new','N')],id='N',value.name = 'dMML_cutoff',variable.name = 'Stat')
ggplot(cutoff_N_min_mt_dMML,aes(x=N,y=dMML_cutoff,fill=Stat))+ylim(0,1)+
  geom_bar(stat="identity", position=position_dodge())+theme(legend.position = 'bottom')
hist(HUES64_new_df$N,breaks=20,xlab='N')


cutoff_N_max=data.frame()
for (i in unique(HUES64_new_df$N)){
  cutoff_out=data.frame(dNME_cutoff_N_new=max(HUES64_new_df$dNME[HUES64_new_df$dNME_pval>pval_cutoff & HUES64_new_df$N==i]),
                        dMML_cutoff_N_new=max(HUES64_new_df$dMML[HUES64_new_df$dMML_pval>pval_cutoff & HUES64_new_df$N==i]),
                        dNME_cutoff_N_old=max(HUES64_old_df$dNME[HUES64_old_df$dNME_pval>pval_cutoff & HUES64_old_df$N==i]),
                        dMML_cutoff_N_old=max(HUES64_old_df$dMML[HUES64_old_df$dMML_pval>pval_cutoff & HUES64_old_df$N==i]),
                        N=i)
  cutoff_N_max=rbind(cutoff_N,cutoff_out)
}
cutoff_N_max[order(cutoff_N_max$N),]

HUES64_new_df$ASM_dNME='No'
HUES64_new_df$ASM_dNME[HUES64_new_df$dNME_pval<=pval_cutoff]='Yes'
HUES64_new_df=as.data.frame(HUES64_new_df)
ggplot(HUES64_new_df[HUES64_new_df$N==3,],aes(x=dNME,color=ASM_dNME))+geom_density(size=1)+theme(legend.position = 'bottom')

HUES64_old_df$ASM_dNME='No'
HUES64_old_df$ASM_dNME[HUES64_old_df$dNME_pval<=pval_cutoff]='Yes'
HUES64_old_df=as.data.frame(HUES64_old_df)
ggplot(HUES64_old_df[HUES64_old_df$N==1,],aes(x=dNME,color=ASM_dNME))+geom_density(size=1)+theme(legend.position = 'bottom')
ah=AnnotationHub()
chromHMM=ah[['AH46871']]
sp= "stem_27_undifferentiated_paired - HUES64"
dNME_HUES64_new=chromHMM_OR(HUES64_merge, chromHMM,sp)
dMML_HUES64_new=chromHMM_OR(HUES64_merge, chromHMM,sp,stat="dMML_pval")
ggplot(dNME_HUES64_new[[1]],aes(x=state,y=OR,fill=state))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("dNME-ASM states annotation")+geom_text(aes(label=round(OR,digits = 2)), vjust=3, color="black", size=3.5)

ggplot(dMML_HUES64_new[[1]],aes(x=state,y=OR,fill=state))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("dMML-ASM states annotation")+geom_text(aes(label=round(OR,digits = 2)), vjust=3, color="black", size=3.5)

ggplot(chromHMM_dNME_all_ls$`stem_27_undifferentiated_paired - HUES64`[[1]],aes(x=state,y=OR,fill=state))+geom_bar(stat="identity",color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),width=0.2,position=position_dodge(0.9))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")+xlab("chromHMM states")+
  ggtitle("dNME-ASM states annotation")+geom_text(aes(label=round(OR,digits = 2)), vjust=3, color="black", size=3.5)

#######Checking het CpG location
HetCG_loc=which(unlist(lapply(gff_in$hetCpGg1,function(x) length(x)!=0))|unlist(lapply(gff_in$hetCpGg2,function(x) length(x)!=0)))
gff_Het=gff_in[HetCG_loc]
HetCG_out=GRanges()
for(i in 1:length(gff_Het)){
  chr=seqnames(gff_Het[i])
  hetCpG=unique(c(gffsplit(gff_Het$hetCpGg1[i]),
                  gffsplit(gff_Het$hetCpGg2[i])))
  gff_het_out=GRanges(seqnames = rep(chr,length(hetCpG)),ranges = IRanges(start=hetCpG,end=hetCpG),strand="*")
  gff_het_out$Subject=gff_Het$Subject[i]
  HetCG_out=c(HetCG_out,gff_het_out)
}
gffsplit<-function(x){
  x=unlist(x)
  if(length(x)!=0){
    x=gsub("\\[",'',x)
    x=gsub("\\]",'',x)
    x=as.numeric(unlist(strsplit(x,", ")))
    return(unlist(x))
  }
}
HetCG_out_H9=HetCG_out[HetCG_out$Subject=='H9']
vcf_H9=variant_HetCpG$H9
vcf_H9=vcf_H9[vcf_H9$HetCpg]
olap_H9=findOverlaps(HetCG_out_H9,vcf_H9,maxgap=1)
unique(paste(vcf_H9$REF[subjectHits(olap_H9)],vcf_H9$ALT[subjectHits(olap_H9)],sep='-'))
unique(paste(vcf_H9$REF[-subjectHits(olap_H9)],vcf_H9$ALT[-subjectHits(olap_H9)],sep='-'))



############N=1 issue########################
genomic_features=readRDS("../downstream/input/genomic_features2020.rds")
GR_merge_N_df=data.frame(N1=rep("N>1",length(GR_merge)),ASM=rep("Not ASM",length(GR_merge)),
                         open_sea=rep("not open sea",length(GR_merge)),stringsAsFactors = F)
GR_merge_N_df$N1[GR_merge$N==1]="N=1"
GR_merge_N_df$ASM[GR_merge$dNME_pval<=0.1]="ASM"
olap_open_sea=findOverlaps(GR_merge,genomic_features$`CpG open sea`)
GR_merge_N_df$open_sea[unique(queryHits(olap_open_sea))]="open sea"
barplot(table(GR_merge_N_df$N1))
barplot(table(GR_merge_N_df$N1[GR_merge_N_df$ASM=='ASM']))
barplot(table(GR_merge_N_df$open_sea[GR_merge_N_df$N1=='N=1']))


###############chekcing the enrichment of genomic features stratifying by N###################
GR_merge_1=GR_merge[GR_merge$N==1]
GR_merge_1$ASM='No'
GR_merge_1$ASM[GR_merge_1$dMML_pval<=0.05 & GR_merge_1$Subject=='STL003']='Yes'
genomic_features=readRDS("../downstream/input/genomic_features_new.rds")
feature_dNME_enrich=list()
for(ft in names(genomic_features)){
  feature_dNME_enrich[[ft]]=testEnrichmentFeature_stat(GR_merge_1,genomic_features[[ft]])[[2]]
}

lapply(feature_dNME_enrich,function(x) x$estimate)
GR_merge_2=GR_merge[GR_merge$N==2 &GR_merge$Subject=='STL003']
GR_merge_2$ASM='No'
GR_merge_2$ASM[GR_merge_2$dNME_pval<=0.05]='Yes'
feature_dNME_enrich=list()
for(ft in names(genomic_features)){
  feature_dNME_enrich[[ft]]=testEnrichmentFeature_stat(GR_merge_2,genomic_features[[ft]])[[2]]
}
lapply(feature_dNME_enrich,function(x) x$p.value)

GR_merge$density_round=round(GR_merge$density,digits = 1)
density_frequency=data.frame(log10_density=unique(GR_merge$density_round))
GR_merge_island_olap=findOverlaps(GR_merge,genomic_features$`CpG open sea`)
GR_merge_island=GR_merge[queryHits(GR_merge_island_olap)]
GR_merge_not_island=GR_merge[-queryHits(GR_merge_island_olap)]
prob_island=c()
prob_not_island=c()
N=2
##In feature and out feature
for (den in density_frequency$log10_density){
  prob_island=c(prob_island,sum(GR_merge_island$dNME_pval<=0.1 & GR_merge_island$density_round==den )/
                  sum(GR_merge_island$density_round==den))
  prob_not_island=c(prob_not_island,sum(GR_merge_not_island$dNME_pval<=0.1 & GR_merge_not_island$density_round==den)/
                      sum(GR_merge_not_island$density_round==den))
  
  
}
GR_out=data.frame(denstiy=density_frequency$log10_density,prob_island=prob_island,prob_not_island=prob_not_island)
GR_out[order(GR_out$denstiy),]
GR_2_out=data.frame(denstiy=density_frequency$log10_density,prob_island=prob_island,prob_not_island=prob_not_island)
GR_1_out=data.frame(denstiy=density_frequency$log10_density,prob_island=prob_island,prob_not_island=prob_not_island)


dist_NME_plot('../downstream/input/H1_inforME/MML-H1_merged_all_st.bed',"var")+ggtitle("MML distribution vs distance to TSS")
dist_NME_plot('../downstream/input/H1_inforME/NME-H1_merged_all_st.bed',"var")+ggtitle("NME distribution vs distance to TSS")
dist_NME_plot('../downstream/input/H1_inforME/MML-H1_merged_all_st.bed',"mean")+ggtitle("MML distribution vs distance to TSS")
dist_NME_plot('../downstream/input/H1_inforME/NME-H1_merged_all_st.bed',"mean")+ggtitle("NME distribution vs distance to TSS")


#Read in CpG location and read in data for each allele
gr_allele=readRDS('../downstream/output/GRs_allele2.rds')
#Resize and calculate number of CGs
gr_allele_CpG_resize=GRanges()
for (subj in subjects){gr_allele_CpG_resize=c(gr_allele_CpG_resize,GR_resize_sub(subj,gr_allele_CpG,CpG_hg19,variant_HetCpG))}
#This one should not have density in it
saveRDS(gr_allele_CpG,'../downstream/input/gr_allele_CpG_500_2.rds')
##Checking if difference comes from all C-T?
#Reading in GR
GR=readRDS('../downstream/output/GR.all.diff.H1.GM12878.rds')#Use ASM cutoff=0.05
gr_allele_CpG=readRDS('../downstream/output/gr_allele_CpG_new.rds')#dropping one CpG region in future

#Use GR_merge to do it?
NME_allele=gr_allele_CpG[gr_allele_CpG$Statistic=='NME']
#Assign ASM information and find the ASM region
GR$ASM[GR$pvalue<=pval_cutoff]="Yes"
NME_allele=add_ASM(NME_allele,GR[GR$Statistic=='dNME'])
NME_allele_ASM=NME_allele[which(NME_allele$ASM=='Yes')]


#dMML-ASM
gr_distance(GR_merge[GR_merge$dMML_pval<=pval_cutoff],genomic_features$`CpG island`,
            'distance to nearest CpG island (kb)','dMML ASM distance to nearest CpG island',ylim=c(0,0.5))
#all regions
gr_distance(GR_merge,genomic_features$`CpG island`,
            'distance to nearest CpG island (kb)','all regions distance to CpG island',ylim=c(0,0.5))
#dNME-ASM
gr_distance(GR_merge[GR_merge$dNME_pval<=pval_cutoff],genomic_features$`CpG island`,
            'distance to nearest CpG island (kb)','dNME ASM distance to CpG island',ylim=c(0,0.5))


# Checking distribution of chromHMM with each N ---------------------------


cutoff_N_min=data.frame()
for (i in unique(GR_merge$N)){
  cutoff_out=data.frame(dNME_new=min(GR_merge$dNME[GR_merge$dNME_pval<=pval_cutoff & GR_merge$N==i]),
                        dMML_new=min(GR_merge$dMML[GR_merge$dMML_pval<=pval_cutoff & GR_merge$N==i]),
                        N=i)
  cutoff_N_min=rbind(cutoff_N_min,cutoff_out)
}

ggplot(cutoff_N_min,aes(x=N,y=dNME_new))+ylim(0,1)+ylab('dNME cutoff')+
  geom_bar(stat="identity", position=position_dodge())+theme(legend.position = 'bottom')

ggplot(cutoff_N_min,aes(x=N,y=dMML_new))+ylim(0,1)+ylab('dMML cutoff')+
  geom_bar(stat="identity", position=position_dodge())+theme(legend.position = 'bottom')

for (i in unique(GR_merge$N)){
  cutoff_out=data.frame(dNME_new=max(GR_merge$dNME[GR_merge$dNME_pval>pval_cutoff & GR_merge$N==i]),
                        dMML_new=max(GR_merge$dMML[GR_merge$dMML_pval>pval_cutoff & GR_merge$N==i]),
                        N=i)
  cutoff_N_min=rbind(cutoff_N_min,cutoff_out)
}


#############Getting Chromm regions for each available sample####################
chromHMM_all=list()

for(state in unique(chromHMM_region_all$`rep1 - H1`$name)){
  chromHMM_state=GRanges()
  for(sp in unique(GR_merge$Sample)){
    if (!is.null(chromHMM_region_all[[sp]])){
      chromHMM_state=c(chromHMM_state,
                       subsetByOverlaps(GR_merge[GR_merge$Sample==sp],
                                        chromHMM_region_all[[sp]][chromHMM_region_all[[sp]]$name==state])
      )
    }
  }
  chromHMM_all[[state]]=chromHMM_state
  
}
#Flanking Bivalent TSS/Enh: not enriched in any dNME<0.1 or >0.1
GO_out_sig=list()
GO_out_non_sig=list()
for(state in names(chromHMM_all)){
  GO_in=chromHMM_all[[state]]
  GO_out_sig[[state]]=GO_anno(unique(GO_in$genes_promoter[GO_in$dNME_pval<=pval_cutoff]),unique(GO_in$genes_promoter))[[1]]
  GO_out_non_sig[[state]]=GO_anno(unique(GO_in$genes_promoter[GO_in$dNME_pval>pval_cutoff]),unique(GO_in$genes_promoter))[[1]]
  
}



import.subject.test<-function(inDir,calc='diff'){
  #for calc: diff -> dMML etc, allele -> NME etc
 
  # H9
  H9_subject <- rep("H9",1)
  H9_subject_labels <- rep("CL1",1)
  H9_tissues <- c(
    "42_embryonic_stem_cell_single"
  )             
  H9_tissue_labels <- c(
    "Embryonic Stem Cell"
  )
  H9_gtex_labels <- c(
    ""
  )
  H9_diff_labels <- c(
    "Undifferentiated"
  )
  
  # HUES64
  HUES64_subject <- rep("HUES64",1)
  HUES64_subject_labels <- rep("CL2",1)
  HUES64_tissues <- c(
    "ectoderm_paired"

  )
  HUES64_tissue_labels <- c(

    "Ectoderm"

  )
  HUES64_gtex_labels <- c("")
  HUES64_diff_labels <- c(

    "Semidifferentiated"

  )
  
 
  # STL003
  stl003_subject <- rep("STL003",1)
  stl003_subject_labels <- rep("D7",1)
  stl003_tissues <- c(
    
    "Adrenal_Gland_single"
  )
  stl003_tissue_labels <- c(

    "Adrenal Gland"
  )
  stl003_gtex_labels <- c(

    "Adrenal_Gland"
  )
  stl003_diff_labels <- rep("Differentiated",length(stl003_subject))
  
 
  # Create single vectors
  subjects <- c(H9_subject,HUES64_subject,stl003_subject)
  tissues <- c(H9_tissues,HUES64_tissues,stl003_tissues)
  subject_labels <- c(H9_subject_labels,HUES64_subject_labels,stl003_subject_labels)
  tissue_labels <- c(H9_tissue_labels,HUES64_tissue_labels,stl003_tissue_labels)
  gtex_labels <- c(H9_gtex_labels,HUES64_gtex_labels,stl003_gtex_labels)
  diff_labels <- c(H9_diff_labels,HUES64_diff_labels,stl003_diff_labels)
  GRs=GRanges()
  for (i in 1:length(subjects)) {
    # Print sample being loaded
    print(paste("Loading sample:",subjects[i],tissues[i]))
    if (calc=='diff'){
      GR.in=read.diffGR(subjects[i],tissues[i],inDir,cutoff=0.05)
    }else if(calc=='allele'){
      GR.in=read.alleleGR(subjects[i],tissues[i],inDir)
    }else {cat('Wrong calc \n')}
    GR.in$SubjectLabel <- subject_labels[i]
    GR.in$Tissue <- tissue_labels[i]
    GR.in$GTEx <- gtex_labels[i]
    GR.in$State <- diff_labels[i]
    GRs=append(GRs,GR.in)
  }
  return(GRs)
  
}
cov8_bound= import.subject.test('../downstream/data/cov8_real/',calc='allele')#1432880
GR_allele=import.subject.test('../downstream/data/Run_version2/bedGraph_allele/',calc='allele')
length(GR_allele[GR_allele$Sample %in% cov8_bound$Sample])#1448604
for (sp in unique(cov8_bound$Sample)){
  
  cat(sp,'percent regions left at coverage 8:',length(cov8_bound[cov8_bound$Sample==sp])/length(GR_allele[GR_allele$Sample==sp]),'\n')
  sp2=paste(unique(cov8_bound$Subject[cov8_bound$Sample==sp]),strsplit(sp,' - ')[[1]][1],'phased',sep='_')
  cat(sp,'theoritical percent regions left at coverage 8:',
      coverage_cutoff_agg$N[coverage_cutoff_agg$Sample==sp2&coverage_cutoff_agg$cutoff==12]/
        coverage_cutoff_agg$N[coverage_cutoff_agg$Sample==sp2&coverage_cutoff_agg$cutoff==5],'\n')
}
coverage_cutoff_percent=data.frame()
for(sp in unique(coverage_cutoff_agg$Sample)){
  coverage_cutoff_percent=rbind(coverage_cutoff_percent,
      data.frame(sample=sp,percent=coverage_cutoff_agg$N[coverage_cutoff_agg$Sample==sp&coverage_cutoff_agg$cutoff==12]/
        coverage_cutoff_agg$N[coverage_cutoff_agg$Sample==sp&coverage_cutoff_agg$cutoff==5]))
}

#library(org.Mm.eg.db)
sensory_gene_mouse=unique(select(org.Mm.eg.db, unlist(as.list(org.Mm.egGO2EG)[dNME_GO$GO.ID[1:8]]), c("SYMBOL"), "ENTREZID")$SYMBOL)
GR_merge_sensory=GR_merge[GR_merge$genes_promoter %in% sensory_gene &GR_merge$dNME_pval<=pval_cutoff]
GR_merge_sensory[order(GR_merge_sensory$N,decreasing = T)]
library(readxl)
# mmc1 <-as.data.frame(read_excel("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/mmc1.xls",col_names = FALSE))
# colnames(mmc1)=c("genes","transcription")
# mmc1_GO_genes=mmc1$genes[mmc1$genes %in% sensory_gene_mouse]

# Motif_GO ----------------------------------------------------------------

###GO analysis and more on significant ones, you don't have to do GO analysis because we're not doing whole genome survey any way
motif_gene_sig=motif_gene[motif_gene$geneSymbol%in% motif_dir_sig_qval$TF]
motif_gene_sig=motif_gene[motif_gene$geneSymbol=='CTCF']
#Exclude sensory genes?
motif_GO=GO_anno(subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],motif_gene_sig)$genes_promoter,
                 variant_HetCpG_meta$genes_promoter)
sensory_gene <- genesInTerm(GOdata, sel.terms)
motif_GO_genes=motif_GO[[1]]
motif_GO_genes[motif_GO_genes$Significant>=5,]

# Looking at the chromHMM states separately -------------------------------

#Checking regions with same direction or opposite direction
GR_merge_opposite=(sign(GR_merge_original$g1CG-GR_merge_original$g2CG)!=0) &
  (sign(GR_merge_original$NME1-GR_merge_original$NME2)!=sign(GR_merge_original$g1CG-GR_merge_original$g2CG))
GR_merge_same=(sign(GR_merge_original$NME1-GR_merge_original$NME2)==sign(GR_merge_original$g1CG-GR_merge_original$g2CG))
GR_merge_non=(sign(GR_merge_original$g1CG-GR_merge_original$g2CG)==0)
chromHMM_dNME_opposite_ls=list()
chromHMM_dNME_same_ls=list()
chromHMM_dNME_non_ls=list()

chromHMM_dNME_opposite_ls[[sp]]=chromHMM_OR(GR_merge[GR_merge_opposite], chromHMM,sp)
chromHMM_dNME_same_ls[[sp]]=chromHMM_OR(GR_merge[GR_merge_same], chromHMM,sp)
chromHMM_dNME_non_ls[[sp]]=chromHMM_OR(GR_merge[GR_merge_non], chromHMM,sp)

chromHMM_dNME_opposite_all=chromHMM_combine(chromHMM_dNME_opposite_ls)
chromHMM_dNME_same_all=chromHMM_combine(chromHMM_dNME_same_ls)
chromHMM_dNME_non_all=chromHMM_combine(chromHMM_dNME_non_ls)



#######Covert entrezid to gene symbol##############
gene_all$gene_name= AnnotationDbi::select(Homo.sapiens,key=as.character(gene_all$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
diff_genes=rownames(res_RNA[which(res_RNA$pvalue<=0.1),])
###Check unassigned genes with RNA-seq change
diff_genes=gene_all[gene_all$gene_name %in% diff_genes]

#######Check if there's region overlapping with GR merge######TODO
GR_merge_H1_diff_exp=subsetByOverlaps(GR_merge_H1,diff_genes,maxgap = 0)
GR_merge_H1_diff_exp[GR_merge_H1_diff_exp$dMML_pval<=pval_cutoff]

######Getting ASE with high hypervaribility>=1 overlapping with allele-specific RNA expression to find examples
res_hyper=res[which(res$padj<=0.1),][rownames(res[which(res$padj<=0.1),]) %in% genes_hypervar$gene_name[genes_hypervar$hypervar_log2>=1],]
gene_hyper=gene_all[gene_all$gene_name %in% rownames(res_hyper)] #MIR6723
gene_MIR6723=gene_all[which(gene_all$gene_name =="MIR6723")]
# MIR4461 Chromosome 5: 134,263,720-134,264,016
#RP11-600F24.2 chr14:103,878,271-103,879,282
#RP5-857K21.7: imbalance in coverage, Chromosome 1: 568,137-568,818
gene_not_in_txdb = makeGRangesFromDataFrame(data.frame(seqnames=c("chr5","chr14","chr1"),start=c(134263720,103878271,568137),end=c(134264016,103879282,568818),
                                                       gene=c('MIR4461','RP11-600F24','RP5-857K21.7'),stringsAsFactors = F),keep.extra.columns = T)


######ATAC
dds_ATAC<-DESeq(dds_ATAC,betaPrior=TRUE)
res_ATAC<-results(dds_ATAC,contrast=c("condition","genome2","genome1"))

#reformat to granges
ATAC_CPM=do.call(rbind,strsplit(rownames(res_ATAC_lfc),"_"))
colnames(ATAC_CPM)=c("seqnames","start","end")
ATAC_CPM=makeGRangesFromDataFrame(as.data.frame(gr_CPM))
elementMetadata(ATAC_CPM)=res_ATAC_lfc


# ATAC-dMML ---------------------------------------------------------------

#######Plotting
density_plot_hyper(GR_merge_H1_ATAC,genes_hypervar,genomic_features,"ATAC_FC_CPM")
########Checking if the sign of ATAC and NME different agree, sign is correct
GR_merge_H1_ATAC$dNME_sig=sign(GR_merge_H1_ATAC$NME2-GR_merge_H1_ATAC$NME1)
GR_merge_H1_ATAC$FC_sig=sign(GR_merge_H1_ATAC$ATAC_FC)
sum(GR_merge_H1_ATAC$dNME_sig==GR_merge_H1_ATAC$FC_sig)#9287
sum(GR_merge_H1_ATAC$dNME_sig!=GR_merge_H1_ATAC$FC_sig)#9298
#Look at dNME ASM region
sum(GR_merge_H1_ATAC$dNME_sig[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff]==GR_merge_H1_ATAC$FC_sig[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff])#417
sum(GR_merge_H1_ATAC$dNME_sig[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff]!=GR_merge_H1_ATAC$FC_sig[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff])#397

#######From data from with important and matched values 
ATAC_df=data.frame(log2FC_CPM=GR_merge_H1_ATAC$ATAC_FC_CPM,
                   dNME=GR_merge_H1_ATAC$dNME,
                   dNME_pval=GR_merge_H1_ATAC$dNME_pval,
                   log2FC_pval=GR_merge_H1_ATAC$ATAC_FC_pval,
                   dMML=GR_merge_H1_ATAC$dMML,
                   dMML_pval=GR_merge_H1_ATAC$dMML_pval,
                   promoter=!is.na(GR_merge_H1_ATAC$genes_promoter))
ATAC_df$log2FC_pval_adj=p.adjust(ATAC_df$log2FC_pval,method="BH")
ATAC_df$log2FC_CPM_abs=abs(ATAC_df$log2FC_CPM)
########plot density plot of dNME and ATAC fold change, also plot volcano plot
plot(ATAC_df$log2FC_CPM,-log10(ATAC_df$log2FC_pval),ylim=c(0,10),xlab="log2FC",ylab="-log10(pval)")
plot(ATAC_df$log2FC_CPM,-log10(ATAC_df$log2FC_pval_adj),ylim=c(0,10),xlab="log2FC",ylab="-log10(qval)")
ATAC_df$dNME_sig=NA
ATAC_df$dNME_sig[ATAC_df$dNME_pval<=pval_cutoff]="dNME-ASM"
ATAC_df$dNME_sig[ATAC_df$dNME_pval>pval_cutoff]="non dNME-ASM"
sum(ATAC_df$dNME_pval<=pval_cutoff)#814
sum(ATAC_df$dNME_pval>pval_cutoff) #17771
ATAC_df$ATAC_sig=NA
ATAC_df$ATAC_sig[ATAC_df$log2FC_pval<=0.05 & ATAC_df$log2FC_CPM>=1]="ATAC change"
ATAC_df$ATAC_sig[ATAC_df$log2FC_pval>0.05 | ATAC_df$log2FC_CPM<1]="no ATAC change"
sum(ATAC_df$log2FC_pval<=0.05 & ATAC_df$log2FC_CPM>=1) #107
sum(ATAC_df$log2FC_pval>0.05 | ATAC_df$log2FC_CPM<1)#17604
ggplot(ATAC_df,aes(x=dNME,color=ATAC_sig))+geom_density()+theme(legend.position = "bottom")
ggplot(ATAC_df,aes(x=dMML,color=ATAC_sig))+geom_density()+theme(legend.position = "bottom")
ggplot(ATAC_df,aes(x=log2FC_CPM_abs,color=dNME_sig))+geom_density()+theme(legend.position = "bottom")+xlab("abs(Log2(FC))")
median(ATAC_df$dNME[ATAC_df$log2FC_pval<=0.05 & ATAC_df$log2FC_CPM>=1])#0.180
median(ATAC_df$dNME[ATAC_df$log2FC_pval>0.05 | ATAC_df$log2FC_CPM<1])#0.131
#median(ATAC_df$dNME[abs(ATAC_df$log2FC_CPM)>3])#0.146
#median(ATAC_df$dNME[abs(ATAC_df$log2FC_CPM)<=3])#0.129
#Checking difference in distribution using K-S method
ks.test(ATAC_df$dNME[ATAC_df$log2FC_pval<=0.05 & ATAC_df$log2FC_CPM>=1],ATAC_df$dNME[ATAC_df$log2FC_pval>0.05 | ATAC_df$log2FC_CPM<1]) #107 vs 17604
wilcox.test(ATAC_df$dNME[ATAC_df$log2FC_pval<=0.05 & ATAC_df$log2FC_CPM>=1],ATAC_df$dNME[ATAC_df$log2FC_pval>0.05 | ATAC_df$log2FC_CPM<1], alternative = "two.sided")

ks.test(ATAC_df$dMML[ATAC_df$log2FC_pval<=0.05 & ATAC_df$log2FC_CPM>=1],ATAC_df$dMML[ATAC_df$log2FC_pval>0.05 | ATAC_df$log2FC_CPM<1]) #107 vs17604
wilcox.test(ATAC_df$dMML[ATAC_df$log2FC_pval<=0.05 & ATAC_df$log2FC_CPM>=1],ATAC_df$dMML[ATAC_df$log2FC_pval>0.05 | ATAC_df$log2FC_CPM<1], alternative = "two.sided")



##################Find the correlation between ATAC-seq values CPM vs NME and FC vs dNME
#####################FC vs dNME and dMML###################
###Raw comparison
##Log2FC_CPM have uneven distribution
hist(abs(ATAC_df$log2FC_CPM),breaks=100,xlab="abs(log2FC)",main="")
hist(ATAC_df$log2FC_pval[abs(ATAC_df$log2FC_CPM)>=3],xlab="raw pval",main="")
#At all regions
ggplot(data=ATAC_df,aes(x=abs(log2FC_CPM),y=dNME))+geom_smooth()+geom_point(alpha = 0.005)+ylim(c(0,1))
ggplot(data=ATAC_df,aes(x=abs(log2FC_CPM),y=dMML))+geom_smooth()+geom_point(alpha = 0.005)
summary(lm(abs(ATAC_df$log2FC_CPM)~ATAC_df$dNME)) #0.056, pval=3.418e-14
#at regions with dNME pval <=0.2
ggplot(data=ATAC_df[ATAC_df$dNME_pval<=0.2,],aes(x=abs(log2FC_CPM),y=dNME))+geom_smooth()+geom_point(alpha = 0.05)
ggplot(data=ATAC_df[ATAC_df$dMML_pval<=0.2,],aes(x=abs(log2FC_CPM),y=dMML))+geom_smooth()+geom_point(alpha = 0.5)+ylim(0,1)
ggplot(data=ATAC_df[abs(ATAC_df$log2FC_CPM)>=3,],aes(x=abs(log2FC_CPM),y=dNME))+geom_smooth()+geom_point(alpha = 0.05)
summary(lm(abs(ATAC_df$log2FC_CPM[ATAC_df$dNME_pval<=0.2])~ATAC_df$dNME[ATAC_df$dNME_pval<=0.2])) #-0.0166594, pval=0.6351
summary(lm(abs(ATAC_df$log2FC_CPM[abs(ATAC_df$log2FC_CPM)>=3])~ATAC_df$dNME[abs(ATAC_df$log2FC_CPM)>=3])) #-0.0166594, pval=0.6351
###Bining FC?
#####################CPM vs NME or MML###################
ATAC_df_allele=data.frame(CPM=c(GR_merge_H1_ATAC$g1_CPM,GR_merge_H1_ATAC$g2_CPM),NME=c(GR_merge_H1_ATAC$NME1,GR_merge_H1_ATAC$NME2),
                          dNME_pval=c(GR_merge_H1_ATAC$dMML_pval,GR_merge_H1_ATAC$dMML_pval),MML=c(GR_merge_H1_ATAC$MML1,GR_merge_H1_ATAC$MML2))
ATAC_df_allele$log2CPM=log2(ATAC_df_allele$CPM+1)
ggplot(data=ATAC_df_allele,aes(x=log2CPM,y=NME))+geom_smooth()+geom_point(alpha = 0.05)
ggplot(data=ATAC_df_allele,aes(x=CPM,y=MML))+geom_smooth()+geom_point(alpha = 0.05)+xlim(0,200)
ggplot(data=ATAC_df_allele[ATAC_df_allele$log2CPM<=5.5,],aes(x=log2CPM,y=NME))+geom_smooth()+geom_point(alpha = 0.05)
summary(lm(ATAC_df_allele$NME[ATAC_df_allele$log2CPM<=5.5]~ATAC_df_allele$log2CPM[ATAC_df_allele$log2CPM<=5.5]))
###########Looking at ATAC-seq vs dNME using contengency table
pval_cutoff=0.2
ATAC_cutoff=0.1
ATAC_dNME=sum(GR_merge_H1_ATAC$dNME_pval<=pval_cutoff & (GR_merge_H1_ATAC$ATAC_FC_pval<=ATAC_cutoff & GR_merge_H1_ATAC$ATAC_FC_CPM>=1))
ATAC_nondNME=sum(GR_merge_H1_ATAC$dNME_pval>pval_cutoff & (GR_merge_H1_ATAC$ATAC_FC_pval<=ATAC_cutoff & GR_merge_H1_ATAC$ATAC_FC_CPM>=1))
nonATAC_dNME=sum(GR_merge_H1_ATAC$dNME_pval<=pval_cutoff & (GR_merge_H1_ATAC$ATAC_FC_pval>ATAC_cutoff | GR_merge_H1_ATAC$ATAC_FC_CPM<1))
nonATAC_nondNME=sum(GR_merge_H1_ATAC$dNME_pval>pval_cutoff & (GR_merge_H1_ATAC$ATAC_FC_pval>ATAC_cutoff | GR_merge_H1_ATAC$ATAC_FC_CPM<1))
cont_table=matrix(c(ATAC_dNME,ATAC_nondNME,nonATAC_dNME,nonATAC_nondNME),nrow=2)
fisher.test(cont_table)#1.8121, 1.032484-2.983621 
##########Different pvalue cutoff for ATAC-seq################################
pval_cutoff=0.2
ATAC_dNME=sum(GR_merge_H1_ATAC$dNME_pval<=pval_cutoff & abs(GR_merge_H1_ATAC$ATAC_FC_CPM)>3)
ATAC_nondNME=sum(GR_merge_H1_ATAC$dNME_pval>pval_cutoff & abs(GR_merge_H1_ATAC$ATAC_FC_CPM)>3)
nonATAC_dNME=sum(GR_merge_H1_ATAC$dNME_pval<=pval_cutoff & abs(GR_merge_H1_ATAC$ATAC_FC_CPM)<=3)
nonATAC_nondNME=sum(GR_merge_H1_ATAC$dNME_pval>pval_cutoff & abs(GR_merge_H1_ATAC$ATAC_FC_CPM)<=3)
cont_table=matrix(c(ATAC_dNME,ATAC_nondNME,nonATAC_dNME,nonATAC_nondNME),nrow=2)
fisher.test(cont_table)#1.524, 5.275e-06,  1.272081 1.818785

#########Looking at example regions including promoter and gene body#####
GR_merge_H1_ATAC_sig_dNME=GR_merge_H1_ATAC[GR_merge_H1_ATAC$dNME_pval<=pval_cutoff & abs(GR_merge_H1_ATAC$ATAC_FC_CPM)>2]#14 CRADD, #16 DCAF5 promoter?, #12 MYPN
GR_merge_H1_ATAC_sig_dMML=GR_merge_H1_ATAC[GR_merge_H1_ATAC$dMML_pval<=pval_cutoff & abs(GR_merge_H1_ATAC$ATAC_FC_CPM)>2]#14 CRADD, #16 DCAF5 promoter?, #12 MYPN
#dMML
olap_promoter_dMML=findOverlaps(GR_merge_H1_ATAC_sig_dMML,genomic_features$promoter,maxgap = 2000)
GR_merge_H1_ATAC_sig_promoter_dMML=GR_merge_H1_ATAC_sig_dMML[queryHits(olap_promoter_dMML)]#cohesin for separating chromatin, 
GR_merge_H1_ATAC_sig_promoter_dMML[GR_merge_H1_ATAC_sig_promoter_dMML$N>1]
subsetByOverlaps(GR_merge_H1,GR_merge_H1_ATAC_sig_promoter_dMML[14],maxgap = 20000)
#dNME
olap_promoter_dNME=findOverlaps(GR_merge_H1_ATAC_sig_dNME,genomic_features$promoter,maxgap = 2000)
GR_merge_H1_ATAC_sig_promoter_dNME=GR_merge_H1_ATAC_sig_dNME[queryHits(olap_promoter_dNME)]#cohesin for separating chromatin, 
GR_merge_H1_ATAC_sig_promoter_dNME[GR_merge_H1_ATAC_sig_promoter_dNME$N>1]
subsetByOverlaps(GR_merge_H1,GR_merge_H1_ATAC_sig_promoter_dNME[14],maxgap = 20000)
#dNME accessibility promoter: similar but with promoter, lower dNME
ATAC_df$log2FC_CPM_abs=abs(ATAC_df$log2FC_CPM)
ggplot(ATAC_df,aes(y=dNME,x=log2FC_CPM_abs,color=promoter))+geom_smooth()+geom_point(alpha=0.05)
#dMML:DUSP22, log2 fold change expression: -0.553096320051179 
#PDS5A, hypervar=0.02887, var expected log2 = 9.077, no SNP at transcript
#ATP8B4,hypervar=1.636,no SNP there
#ZNF229, FC change:-2.17183282995201, hypervarbile with dNME and ASE
olap_body=findOverlaps(GR_merge_H1_ATAC,genomic_features$`gene body`)
#######too few to find genes in promoter####################################
# GR_merge_H1_ATAC_promoter=GR_merge_H1_ATAC[queryHits(olap_promoter)]
# GR_merge_H1_ATAC_promoter$genes_promoter=genomic_features$promoter$gene_name[subjectHits(olap_promoter)]
# GR_merge_H1_ATAC_promoter_sig=GR_merge_H1_ATAC_promoter[GR_merge_H1_ATAC_promoter$dNME_pval<=pval_cutoff & GR_merge_H1_ATAC_promoter$ATAC_FC_pval<=0.05]
#######Look at gene body####################################
GR_merge_H1_ATAC_body=GR_merge_H1_ATAC[queryHits(olap_body)]
GR_merge_H1_ATAC_body$genes_body=genomic_features$`gene body`$gene_name[subjectHits(olap_body)]

GR_merge_H1_ATAC_body_sig=GR_merge_H1_ATAC_body[GR_merge_H1_ATAC_body$dNME_pval<=pval_cutoff & abs(GR_merge_H1_ATAC_body$ATAC_FC_CPM)>3&
                                                  GR_merge_H1_ATAC_body$genes_body%in% genes_hypervar$gene_name]
GR_merge_H1_ATAC_body_sig[order(GR_merge_H1_ATAC_body_sig$N,decreasing=T)][6:16]#same example:1,3,6 counter: 4,
###################Check whether ATAC data near any gene############
gr_promoter_olap=findOverlaps(gr,genomic_features$promoter)
gr_promoter=gr[queryHits(gr_promoter_olap)]
gr_promoter$genes_promoter=genomic_features$promoter$gene_name[subjectHits(gr_promoter_olap)]
##################Check whether ATAC data agree with RNA-seq data?


# scATAC-seq data ---------------------------------------------------------

library(readr)
GR_merge=readRDS(GR_merge_file)
GR_merge_H1=GR_merge[GR_merge$Subject=='H1']
#GR_merge_H1=GR_merge_H1[GR_merge_H1$N>=3]
scATAC_hyper=makeGRangesFromDataFrame(as.data.frame(read_csv("../downstream/input/H1_scATAC_hypervar_result.csv"),
                                                    stringsAsFactors=F),keep.extra.columns = T)
olap_GR_ATAC_Hyper=findOverlaps(GR_merge_H1,scATAC_hyper,maxgap = 500,select="first")
sum(!is.na(olap_GR_ATAC_Hyper))#9134 GR have hypervarble data,total 126027, 7% of the region, 2k 11516/60062

GR_merge_H1$ATAC_hyper=NA
GR_merge_H1$ATAC_hyper=scATAC_hyper$hypervar_log2[olap_GR_ATAC_Hyper]
GR_merge_H1=GR_merge_H1[!is.na(GR_merge_H1$ATAC_hyper)]
GR_merge_H1$mean_NME=(GR_merge_H1$NME1+GR_merge_H1$NME2)/2
GR_merge_H1$mean_MML=(GR_merge_H1$MML1+GR_merge_H1$MML2)/2
density_plot_hyper(GR_merge_H1,genes_hypervar,genomic_features,"ATAC_hyper")
GR_merge_H1$ATAC_hyper_interval=findInterval(GR_merge_H1$ATAC_hyper,quantile(GR_merge_H1$ATAC_hyper,prob=c(0,0.25,0.5,0.75)))
quant=c("0-25%","25%-50%","50%-75%","75%-100%")
GR_merge_H1$ATAC_hyper_interval=quant[GR_merge_H1$ATAC_hyper_interval]
GR_hyper_df=data.frame(ATAC_hyper=GR_merge_H1$ATAC_hyper,NME=GR_merge_H1$mean_NME,dNME=GR_merge_H1$dNME,
                       dMML=GR_merge_H1$dMML,MML=GR_merge_H1$mean_MML,promoter=!is.na(GR_merge_H1$genes_promoter),
                       interval=GR_merge_H1$ATAC_hyper_interval)
ggplot(data=GR_hyper_df,aes(x=ATAC_hyper,y=NME))+geom_smooth()+geom_point(alpha=0.05)
ggplot(data=GR_hyper_df,aes(x=ATAC_hyper,y=dNME))+geom_smooth()+geom_point(alpha=0.05)
ggplot(data=GR_hyper_df,aes(x=ATAC_hyper,y=MML))+geom_smooth()+geom_point(alpha=0.05)
ggplot(data=GR_hyper_df,aes(x=ATAC_hyper,y=dMML))+geom_smooth()+geom_point(alpha=0.05)
#dNME and NME density at different interval
ggplot(data=GR_hyper_df,aes(x=NME,color=interval))+geom_density(size=1)+theme(legend.position = "bottom")
ggplot(data=GR_hyper_df,aes(x=dNME,color=interval))+geom_density(size=1)+theme(legend.position = "bottom")
ggplot(data=GR_hyper_df,aes(x=MML,color=interval))+geom_density(size=1)+theme(legend.position = "bottom")
ggplot(data=GR_hyper_df,aes(x=dMML,color=interval))+geom_density(size=1)+theme(legend.position = "bottom")
#Only looking at island region: nothing
genomic_features=readRDS("../downstream/input/genomic_features2020.rds")
GR_merge_H1_island=subsetByOverlaps(GR_merge_H1,genomic_features$`CpG island`,maxgap = 1000)
GR_hyper_df_island=data.frame(ATAC_hyper=GR_merge_H1_island$ATAC_hyper,NME=GR_merge_H1_island$mean_NME,dNME=GR_merge_H1_island$dNME,
                              dMML=GR_merge_H1_island$dMML,MML=GR_merge_H1_island$mean_MML,promoter=!is.na(GR_merge_H1_island$genes_promoter),
                              interval=GR_merge_H1_island$ATAC_hyper_interval)
ggplot(data=GR_hyper_df_island,aes(x=ATAC_hyper,y=NME))+geom_smooth()+geom_point(alpha=0.05)
ggplot(data=GR_hyper_df_island,aes(x=ATAC_hyper,y=dNME))+geom_smooth()+geom_point(alpha=0.05)
ggplot(data=GR_hyper_df_island,aes(x=ATAC_hyper,y=MML))+geom_smooth()+geom_point(alpha=0.05)
ggplot(data=GR_hyper_df_island,aes(x=ATAC_hyper,y=dMML))+geom_smooth()+geom_point(alpha=0.05)



#######Getting motif analysis result from Jason
prop <- readRDS("D../downstream/input/prop.rds")
prop_enrich=lapply(prop,function(x) x[x$FDR<=0.1,])
prop_enrich_pos=lapply(prop_enrich,function(x) x[x$diffprop>0,])
prop_enrich_pos=lapply(prop_enrich_pos,function(x) {x$TF=rownames(x)
return(x)})
do.call(rbind,prop_enrich_pos)

# motif_analsysi_jason ----------------------------------------------------


do.call(rbind,prop_enrich_pos)#27 motif
prop_enrich_neg=lapply(prop_enrich,function(x) x[x$diffprop<0,])
prop_enrich_neg=lapply(prop_enrich_neg,function(x) {x$TF=rownames(x)
return(x)})
do.call(rbind,prop_enrich_neg) #28623 motif


sum(H1_hyper_var$gene_name %in% GR_merge_H1$genes_promoter) #2742
#genes_hypervar[!genes_hypervar$gene_name %in% GR_merge$genes_promoter,]
# dNME_promo_GO=GO_anno(unique(genes_hypervar$gene_name[which(genes_hypervar$gene_name %in% GR_merge_H1$genes_promoter[GR_merge_H1$dNME_pval<=pval_cutoff])]),
#                       unique(genes_hypervar$gene_name))


motif_sig_df=motif_enrich(subsetByOverlaps(motif_gene,GR_merge_sp_states),GR_merge_sp_states,pval_cutoff =pval_cutoff,dist=500)


#Old directionality test
H1_cont_table=NME_enrich(GR_merge_H1,H1_hyper_var,'H1')
HESC_cont_table=NME_enrich(GR_merge[GR_merge$Sample=='stem_27_undifferentiated_paired - HUES64'],HESC_hyper_var,'HESC')#2.3
adult_adipose=readRDS('../downstream/input/scRNA/AdultAdipose_1.rds')
Adipose_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Adipose_Tissue_single','Adipose_single')],adult_adipose,'Adipose')#2.3
adult_bladder=readRDS('../downstream/input/scRNA/AdultBladder_2.rds')
Bladder_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Bladder_single')],adult_bladder,'Bladder')#2.3
adult_liver=readRDS('../downstream/input/scRNA/AdultLiver_2.rds')
liver_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Liver_single')],adult_liver,'Liver')#2.3
adult_lung=readRDS('../downstream/input/scRNA/AdultLung_1.rds')
lung_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Lung_single')],adult_lung,'Lung')#2.3
adult_sigmoid_colon=readRDS('../downstream/input/scRNA/AdultSigmoidColon_1.rds')
sigmoid_colon_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Sigmoid_colon')],adult_sigmoid_colon,'Sigmoid Colon')#2.3
adult_spleen=readRDS('../downstream/input/scRNA/AdultSpleen_1.rds')
spleen_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Spleen_single')],adult_spleen,'Spleen')#2.3
adult_gastric=readRDS('../downstream/input/scRNA/AdultStomach_3.rds')
gastric_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Gastric_single')],adult_gastric,'Gastric')#2.3
adult_AG=readRDS('../downstream/input/scRNA/AdultAdrenalGland_2.rds')
AG_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Adrenal_Gland_single')],adult_AG,'Adrenal Gland')#2.3
adult_AT=readRDS('../downstream/input/scRNA/AdultArtery_1.rds')
AT_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Aorta_single')],adult_AT,'Aorta')#2.3
adult_ES=readRDS('../downstream/input/scRNA/AdultEsophagus_2.rds')
ES_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Esophagus_single')],adult_ES,'Esophagus')#2.3
adult_Pan=readRDS('../downstream/input/scRNA/AdultPancreas_1.rds')
Pan_cont_table=NME_enrich(GR_merge[GR_merge$tissue%in% c('Pancreas_single')],adult_Pan,'Pancreas')#2.3
CHM_total=rbind(H1_cont_table,
                HESC_cont_table,
                Adipose_cont_table,
                Bladder_cont_table,
                liver_cont_table,
                lung_cont_table,
                sigmoid_colon_cont_table,
                spleen_cont_table,
                gastric_cont_table,
                AG_cont_table,
                AT_cont_table,
                ES_cont_table,
                Pan_cont_table)
CMH_test(CHM_total)

#NME_enrich old
#Checking if from Jason read
# if(length(genes_hypervar$gene_name)==0){
#   colnames(genes_hypervar)[4]='hypervar_log2'
#   genes_hypervar$gene_name=rownames(genes_hypervar)
# }
# genes_hypervar$hypervar=genes_hypervar$hypervar_log2>=quantile(genes_hypervar$hypervar_log2,prob=0.75)
# GR_merge_sp$hyper_var_genes=GR_merge_sp$TSS%in%genes_hypervar$gene_name[genes_hypervar$hypervar]
# GR_merge_sp=GR_merge_sp[GR_merge_sp$TSS %in% genes_hypervar$gene_name]
# #print(GR_merge_sp)

#Density Heatmap trial
#Density heatmap trial
density_df$dNME=round(density_df$dNME,digits = 2)
density_df$density=round(density_df$density,digits = 2)
density_df_agg_heat=aggregate(density_df,by=list(density_df$density,density_df$dNME),FUN=length)
density_df_agg_heat=density_df_agg_heat[,c(1,2,3)]
colnames(density_df_agg_heat)=c('density','dNME','count')
ggplot(density_df_agg_heat , aes(x = density, y = dNME,fill=count)) +
  geom_tile()
###Allele-agnostic H1

NME_H1_calc=dist_plot_calc('../downstream/input/H1_inforME/H1_allele_agnostic_nme.bedGraph',
                           H1_hyper_var,genomic_features,GR_merge[GR_merge$Sample=='merged - H1'])
dist_plot_run(NME_H1_calc)
MML_H1_calc=dist_plot_calc('../downstream/input/H1_inforME/H1_allele_agnostic_mml.bedGraph',
                           H1_hyper_var,genomic_features,GR_merge[GR_merge$Sample=='merged - H1'])
dist_plot_run(MML_H1_calc)
#H1 RNA-seq
cor.test(H1_df$dMML_relative[H1_df$dMML_pval<=pval_cutoff],H1_df$ASE_log2FC[H1_df$dMML_pval<=pval_cutoff])
plot(H1_df$dMML_relative[H1_df$dMML_pval<=pval_cutoff],H1_df$ASE_log2FC[H1_df$dMML_pval<=pval_cutoff])
cor(H1_df$dNME[H1_df$dNME_pval<=pval_cutoff],abs(H1_df$ASE_log2FC[H1_df$dNME_pval<=pval_cutoff]))
plot(H1_df$dNME[H1_df$dNME_pval<=pval_cutoff],abs(H1_df$ASE_log2FC[H1_df$dNME_pval<=pval_cutoff]))

dNME_df=rbind(data.frame(age=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$NME1>NME_quant|variant_HetCpG_meta$NME2>NME_quant])$AgeMode_Jnt,
                         type='dNME-Hap'),
              data.frame(age=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$NME1<=NME_quant&variant_HetCpG_meta$NME2<=NME_quant])$AgeMode_Jnt,
                         type='non-dNME-Hap')
)

# Age analysis ------------------------------------------------------------


ggplot(dNME_df,aes(x=age,color=type))+geom_density(size=1)
wilcox.test(dNME_df$age[dNME_df$type=='dNME-Hap'],dNME_df$age[dNME_df$type=='non-dNME-Hap'])#1.2838e+13,
t.test(dNME_df$age[dNME_df$type=='dNME-Hap'],dNME_df$age[dNME_df$type=='non-dNME-Hap'])#1.2838e+13,
dMML_df=rbind(data.frame(age=log10(subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval<=pval_cutoff])$AgeMean_Jnt),
                         MAF=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval<=pval_cutoff])$MAF,
                         type='dMML-Hap'),
              data.frame(age=log10(subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval>pval_cutoff])$AgeMean_Jnt),
                         type='non-dMML-Hap')
)
ggplot(dMML_df,aes(x=age,color=type))+geom_density()
wilcox.test(dMML_df$age[dMML_df$type=='dMML-Hap'],dMML_df$age[dMML_df$type=='non-dMML-Hap'])#6.1552e+11
t.test(dMML_df$age[dMML_df$type=='dMML-Hap'],dMML_df$age[dMML_df$type=='non-dMML-Hap'])#6.1552e+11
all_df=rbind(data.frame(age=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff])$AgeMedian_Jnt,
                        DAF=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff])$MAF,
                        type='dNME-Hap'),
             data.frame(age=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval<=pval_cutoff])$AgeMedian_Jnt,
                        DAF=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval<=pval_cutoff])$MAF,
                        type='dMML-Hap'),
             data.frame(age=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval>pval_cutoff&
                                                                           variant_HetCpG_meta$dNME_pval>pval_cutoff])$AgeMedian_Jnt, 
                        DAF=subsetByOverlaps(age_out,variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval>pval_cutoff&
                                                                           variant_HetCpG_meta$dNME_pval>pval_cutoff])$MAF,
                        type='non-dMML-dNME-Hap')
)

#For each region in output, add heterogyzous CpG number in it, change het CpG count and output ranges
gr_allele_CpG=GRanges()
for (subj in subjects){gr_allele_CpG=c(gr_allele_CpG,hetCGallele_sub(subj,GR_allele,hetCpG_gff,CpG_hg19,variant_HetCpG))}
saveRDS(gr_allele_CpG,'../downstream/output/gr_allele_CpG_run4.rds')


# dNME vs local density and ecdf ------------------------------------------


NME_ASM=readRDS('../downstream/output/dNME/NME_ASM_het_calc_run3.rds')
#Plot cumulative distribution of NME with HetCpG and without HetCpG at ASM region
plot(ecdf(NME_ASM$diff_NME[NME_ASM$CpGdiff!=0]),main='cumulative distribution of dNME at ASM',xlab='dNME',xlim=c(-1,1),lwd=1)
lines(ecdf(NME_ASM$diff_NME[NME_ASM$CpGdiff==0]),col='red',lwd=3)
#Plot density difference vs dNME at ASM
density_df=data.frame(density_diff=round(log10(abs(NME_ASM$density_diff[NME_ASM$density_diff!=0])),digits=3),
                      dNME=NME_ASM$diff_NME[NME_ASM$density_diff!=0])
density_df_agg=aggregate(abs(density_df$dNME),by=list(density_df$density_diff),FUN=median)
ggplot(density_df_agg,aes(x=Group.1, y=x))+
  ylim(c(0,1))+ggtitle("dNME change as density difference change")+geom_smooth(method="lm")+
  ylab("dNME")+xlab("log10(density difference)")+geom_point()

ggplot(density_df[density_df$dNME_pval<=pval_cutoff],aes(x=density_diff, y=dNME))+
  ylim(c(0,1))+ggtitle("dNME change as density difference change")+geom_smooth(method="lm")+
  ylab("dNME")+xlab("log10(density difference)")+geom_point()
OR_VMR_perm<-function(NME_dat,vmr,nolap){
  NME_dat$VMR=FALSE
  NME_dat$VMR[sample(1:length(NME_dat),nolap,replace = F)]=TRUE
  NME_VMR=sum(NME_dat$quant_score%in%c('75%-100%')&NME_dat$VMR)
  nonNME_VMR=sum(!NME_dat$quant_score%in%c('75%-100%')&NME_dat$VMR)
  nonNME_nonVMR=sum(!NME_dat$quant_score%in%c('75%-100%')&NME_dat$VMR)
  NME_nonVMR=sum(NME_dat$quant_score%in%c('75%-100%')&NME_dat$VMR)
  fisher.test(matrix(c(NME_VMR,nonNME_VMR,NME_nonVMR,nonNME_nonVMR),nrow = 2))$estimate#OR=3.03
}
#Permutation test
olap_length=length(subsetByOverlaps(NME_in,vmr))
perm_OR=replicate(5,OR_VMR_perm(NME_in,vmr, olap_length))

#Motif break
#Check difference in number of variants
# enchancer_DNase_gr=readRDS("../downstream/input/enchancer_DNase.rds")
# prom_DNase_gr=readRDS("../downstream/input/prom_DNase.rds")
# DNase_all=c(prom_DNase_gr,enchancer_DNase_gr)
# variant_in_run1_from_2nd_run -------------------------------------------
# variant1=readRDS('../downstream/output/ASM_enrich_meta.rds')
# variant2=readRDS('../downstream/output/ASM_enrich_meta_new.rds')
# olap=findOverlaps(variant2,variant1)
# variant_in=unique(variant2[-queryHits(olap)])
# elementMetadata(variant_in)=elementMetadata(variant_in)[1:4]
#saveRDS(variant_in,'../downstream/output/variant_in_run1.rds')
#variant_in_Dnase=subsetByOverlaps(variant_in,DNase_all)
#saveRDS(variant_in_Dnase,'../downstream/output/variant_in_run1_Dnase.rds')
# variant_in_run3_from_3nd_run -------------------------------------------
# variant3=readRDS('../downstream/output/ASM_enrich_meta_run3.rds')
# variant2=readRDS('../downstream/output/ASM_enrich_meta_new.rds')
# olap=findOverlaps(variant3,variant2)
# variant_in=unique(variant3[-queryHits(olap)])
# variant_in_Dnase=subsetByOverlaps(variant_in,DNase_all)
# saveRDS(variant_in_Dnase,'../downstream/output/variant_in_run3_Dnase.rds')
# saveRDS(variant_in,'../downstream/output/variant_in_run3.rds')#Note this is using more regions than ASM_enrich_meta_new.rds


# Motif break example -----------------------------------------------------

###Finding example
H1_hyper_var=as.data.frame(read.csv("../downstream/input/H1 data/H1_hypervar_result.csv"))
adult_adipose=readRDS('../downstream/input/scRNA/AdultAdipose_1.rds')
adult_bladder=readRDS('../downstream/input/scRNA/AdultBladder_1.rds')
adult_AG=readRDS('../downstream/input/scRNA/AdultAdrenalGland_2.rds')#7,6
adult_AT=readRDS('../downstream/input/scRNA/AdultArtery_1.rds')#9
adult_ES=readRDS('../downstream/input/scRNA/AdultEsophagus_1.rds')#SPIN1,UNC13D
adult_GA=readRDS('../downstream/input/scRNA/AdultStomach_3.rds')#
adult_SC=readRDS('../downstream/input/scRNA/AdultSigmoidColon_1.rds')#AMPD3,JAML
adult_CE=readRDS('../downstream/input/scRNA/AdultCerebellum_1.rds')#AMPD3,JAML
fetal_brain=readRDS('../downstream/input/scRNA/FetalBrain_5.rds')
adult_LV=readRDS('../downstream/input/scRNA/AdultLiver_2.rds')#AMPD3,JAML
adult_PS=readRDS('../downstream/input/scRNA/AdultPancreas_1.rds')#AMPD3,JAML
#14
#Try all first if not try individuals
GR_merge_dir[which(GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$hyper_var_promoter>=GR_merge_dir$hyper_var_upper)]

GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue%in%c('Adipose_single','Adipose_Tissue_single')&
               GR_merge_dir$genes_promoter%in%rownames(adult_adipose[adult_adipose$hypervar_logvar>=quantile(adult_adipose$hypervar_logvar,prob=0.75),])]

GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue%in%c('Bladder_single')&
               GR_merge_dir$genes_promoter%in%rownames(adult_bladder[adult_bladder$hypervar_logvar>=quantile(adult_bladder$hypervar_logvar,prob=0.75),])]

GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue%in%c('Adrenal_Gland_single')&
               GR_merge_dir$genes_promoter%in%rownames(adult_AG[adult_AG$hypervar_logvar>=quantile(adult_AG$hypervar_logvar,prob=0.75),])]
#One hit (Best)
hit=GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue%in%c('Esophagus_single')&
                   GR_merge_dir$genes_promoter%in%rownames(adult_ES[adult_ES$hypervar_logvar>=quantile(adult_ES$hypervar_logvar,prob=0.75),])]
#One hit
GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue%in%c('Gastric_single')&
               GR_merge_dir$genes_promoter%in%rownames(adult_GA[adult_GA$hypervar_logvar>=quantile(adult_GA$hypervar_logvar,prob=0.75),])]
#one hit
GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue%in%c('Sigmoid_Colon_single')&
               GR_merge_dir$genes_promoter%in%rownames(adult_SC[adult_SC$hypervar_logvar>=quantile(adult_SC$hypervar_logvar,prob=0.75),])]

GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$Sample=='HuFGM02'&
               GR_merge_dir$genes_promoter%in%rownames(fetal_brain[fetal_brain$hypervar_logvar>=quantile(fetal_brain$hypervar_logvar,prob=0.75),])]

GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue=='Aorta_single'&
               GR_merge_dir$genes_promoter%in%rownames(adult_AT[adult_AT$hypervar_logvar>=quantile(adult_AT$hypervar_logvar,prob=0.75),])]

GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue=='Liver_single'&
               GR_merge_dir$genes_promoter%in%rownames(adult_LV[adult_LV$hypervar_logvar>=quantile(adult_LV$hypervar_logvar,prob=0.75),])]

GR_merge_dir[GR_merge_dir$N>=2&GR_merge_dir$dNME_pval<=0.1&GR_merge_dir$tissue=='Pancreas_single'&
               GR_merge_dir$genes_promoter%in%rownames(adult_PS[adult_PS$hypervar_logvar>=quantile(adult_PS$hypervar_logvar,prob=0.75),])]


plot_motif=subsetByOverlaps(motif_gene[abs(motif_gene$alleleDiff)>=1.5],hit)
plot_motif=plot_motif[plot_motif$geneSymbol%in%motif_dir_sig_qval$TF]
plotMB(plot_motif,'HUES64-1216110')


###dot plot of the trend
motif_gene_example='NFIB'
motif_gene=motif_gene[abs(motif_gene$alleleDiff)>1.5,]

motif_gene_TF=motif_gene[motif_gene$geneSymbol==motif_gene_example]
olap=findOverlaps(variant_HetCpG_meta,motif_gene_TF)
variant_HetCpG_meta_TF=variant_HetCpG_meta[queryHits(olap)]
variant_HetCpG_meta_TF$TF=motif_gene_TF$geneSymbol[subjectHits(olap)]
variant_HetCpG_meta_TF$binddiff=motif_gene_TF$alleleDiff[subjectHits(olap)]
variant_HetCpG_meta_TF$NME_diff=variant_HetCpG_meta_TF$altNME-variant_HetCpG_meta_TF$refNME
variant_HetCpG_meta_binding=variant_HetCpG_meta_TF[!is.na(variant_HetCpG_meta_TF$binddiff) & variant_HetCpG_meta_TF$dNME_pval<=pval_cutoff]
plot(variant_HetCpG_meta_binding$NME_diff,variant_HetCpG_meta_binding$binddiff,xlab='Entropy difference',ylab='Binding likelyhood difference')
abline(v=0)
abline(h=0)


# olap_table_promoter=table(queryHits(olap_promoter))
# unique_gene=names(olap_table_promoter[olap_table_promoter==1])
# olap_promoter_unique=olap_promoter[which(queryHits(olap_promoter)%in%unique_gene)]
# GR_merge$genes_promoter[queryHits(olap_promoter_unique)]=genomic_features$promoter$gene_name[subjectHits(olap_promoter_unique)]
# non_uique_gene=names(olap_table_promoter[olap_table_promoter>1])
# olap_promoter_non_unique=olap_promoter[which(queryHits(olap_promoter)%in%non_uique_gene)]
####Archive

motif_dir_IPA=data.frame(TF=motif_dir$TF,FC=motif_dir$prob/0.5,FDR=motif_dir$qvalue,stringsAsFactors = F )
motif_dir_IPA$TF=gsub('\\(var.2\\)','',motif_dir_IPA$TF)
motif_dir_IPA$TF=gsub('\\(var.3\\)','',motif_dir_IPA$TF)
motif_dir_IPA$TF=strsplit(motif_dir_IPA$TF,"::")
motif_dir_IPA=lapply(1:nrow(motif_dir_IPA),function(x,motif_dir_in){
  motif_dir_in=motif_dir_in[x,]
  motif_dir_out=data.frame()
  for(TF in motif_dir_in$TF){
    motif_dir_out=rbind(motif_dir_out,data.frame(TF=TF,FC=motif_dir_in$FC,FDR=motif_dir_in$FDR,stringsAsFactors = F))
    
  }
  return(motif_dir_out)
},motif_dir_IPA)
motif_dir_IPA_out=do.call(rbind,motif_dir_IPA)
write.csv(motif_dir_IPA_out,'../downstream/output/IPA/motif_IPA.csv')

motif_dir_ic_GO=GO_anno(motif_reformat(motif_dir_ic[motif_dir_ic$qvalue<=0.1,]),motif_reformat(motif_dir_ic),topNodes=4451)
motif_dir_ic_GO_ls=motif_dir_ic_GO[[1]]
motif_dir_ic_GO_ls=motif_dir_ic_GO_ls[motif_dir_ic_GO_ls$Expected>=5,]
motif_dir_ic_GO_ls$qvalue=p.adjust(motif_dir_ic_GO_ls$classicFisher,method='BH')

variant_pref_ic=motif_pref(variant_HetCpG_meta,readRDS('../downstream/output/motif_all_JASPAR_ic.rds'),motif_dir_ic)
variant_pref_ic_GO=GO_anno(unique(variant_pref_ic$genes_promoter),unique(variant_HetCpG_meta$genes_promoter))
variant_pref_ic_GO_ls=variant_pref_ic_GO[[1]]
variant_pref_ic_GO_ls=variant_pref_ic_GO_ls[variant_pref_log_GO_ls$Expected>=5,]

motif_dir_default=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_dir_default_GO=GO_anno(motif_reformat(motif_dir_default[motif_dir_default$qvalue<=0.1,]),motif_reformat(motif_dir_default),topNodes=4374)
motif_dir_default_GO_ls=motif_dir_default_GO[[1]]
motif_dir_default_GO_ls=motif_dir_default_GO_ls[motif_dir_default_GO_ls$Expected>=5,]
motif_dir_default_GO_ls$qvalue=p.adjust(motif_dir_default_GO_ls$classicFisher,method='BH')

variant_pref_default=motif_pref(variant_HetCpG_meta,readRDS('../downstream/output/motif_all_JASPAR_default.rds'),motif_dir_default)
variant_pref_default_GO=GO_anno(unique(variant_pref_default$genes_promoter),unique(variant_HetCpG_meta$genes_promoter))
variant_pref_default_GO_ls=variant_pref_default_GO[[1]]
variant_pref_default_GO_ls=variant_pref_default_GO_ls[variant_pref_default_GO_ls$Expected>=5,]

motif_dir_log=readRDS('../downstream/output/motif_dirction_all_JASPAR_log.rds')
motif_dir_log$qvalue=motif_dir_log$qval_binom
motif_dir_log_GO=GO_anno(motif_reformat(motif_dir_log[motif_dir_log$qvalue<=0.1,]),motif_reformat(motif_dir_log),topNodes=4839)
motif_dir_log_GO_ls=motif_dir_log_GO[[1]]
motif_dir_log_GO_ls=motif_dir_log_GO_ls[motif_dir_log_GO_ls$Expected>=5,]
motif_dir_log_GO_ls$qvalue=p.adjust(motif_dir_log_GO_ls$classicFisher,method='BH')
variant_pref_log=motif_pref(variant_HetCpG_meta,readRDS('../downstream/output/motif_all_JASPAR_log.rds'),motif_dir_log)
variant_pref_log_GO=GO_anno(unique(variant_pref_log$genes_promoter),unique(variant_HetCpG_meta$genes_promoter))
variant_pref_log_GO_ls=variant_pref_log_GO[[1]]
variant_pref_log_GO_ls=variant_pref_log_GO_ls[variant_pref_log_GO_ls$Expected>=5,]

##GMT custom GO
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
KEGG=read.gmt('../downstream/input/KEGG.gmt')
GO_KEGG=GO_custom(unique(unlist(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff])),KEGG,unique(unlist(GR_merge$genes_promoter)))
variant_pref_default=motif_pref(variant_HetCpG_meta,readRDS('../downstream/output/motif_all_JASPAR_default.rds'),motif_dir)
GO_KEGG=GO_custom(unique(unlist(variant_pref_default$genes_promoter[variant_pref_default$dNME_pval<=pval_cutoff])),
                  KEGG,unique(unlist(variant_HetCpG_meta$genes_promoter)))
GO_KEGG=GO_custom(unique(unlist(variant_pref_default$genes_promoter[variant_pref_default$dNME_pval<=pval_cutoff])),
                  KEGG,unique(unlist(variant_pref_default$genes_promoter)))
KEGG=read.gmt('../downstream/input/C2_all.gmt')
GO_KEGG=GO_custom(unique(motif_reformat(motif_dir_sig[motif_dir_sig$qval_binom<=0.1 & motif_dir_sig$prob<0.5,])),
                  KEGG,unique(motif_reformat(motif_dir_sig)))

CM=read.gmt('../downstream/input/cancer_modules.gmt')
GO_CM=GO_custom(unique(unlist(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff])),CM,unique(unlist(GR_merge$genes_promoter)))
GO_CM[GO_CM$FDR<=0.1,]#None
GO_CM=GO_custom(unique(unlist(variant_pref_default$genes_promoter[variant_pref_default$dNME_pval<=pval_cutoff])),
                CM,unique(unlist(variant_HetCpG_meta$genes_promoter)))
OS=read.gmt('../downstream/input/Oncogenic_signatures.gmt')
GO_OS=GO_custom(unique(unlist(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff])),OS,unique(unlist(GR_merge$genes_promoter)))
GO_OS[GO_OS$FDR<=0.1,]#None
GO_OS=GO_custom(unique(unlist(variant_pref_default$genes_promoter[variant_pref_default$dNME_pval<=pval_cutoff])),
                OS,unique(unlist(variant_HetCpG_meta$genes_promoter)))
C7=read.gmt('../downstream/input/C7_all.gmt')
GO_KEGG=GO_custom(unique(motif_reformat(motif_dir_sig[motif_dir_sig$qval_binom<=0.1 & motif_dir_sig$prob<0.5,])),
                  OS,unique(motif_reformat(motif_dir_sig)))
GO_C7=GO_custom(unique(unlist(GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff])),C7,unique(unlist(GR_merge$genes_promoter)))

#old =156
#no cutoff: olap = 77, pval perm <0.0001
#with cutoff olap=17, pval_perm =0.0026
motif_dir_sig$qval_binom=p.adjust(motif_dir_sig$binom.pval,method='BH')

motif_sig=unlist(strsplit(motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1&motif_dir_sig$prob>0.5],"::"))

write(motif_all,'../downstream/output/motif_all_JASPR.txt')
write(motif_sig,'../downstream/output/motif_sig_JASPAR.txt')
#dMML even have smaller number of motifs
motif_dir_sig_qval=motif_dir_sig[motif_dir_sig$qval_binom<=0.1,]
motif_short_list=motif_dir_sig_qval$TF[motif_dir_sig_qval$prob>0.5]
motif_short_list=c('GATA4','GATA5','NFIA','FOXP3','MEF2A','CDX2','CDX1','MEIS2','FOXB1','FOXL1','FOXC2','HOXA10','HOXB13',
                   'FOXC1','POU1F1','POU5F1P1','FOXO4','POU3F2','POU3F4','POU2F2','POU2F3','POU3F1','POU2F1','POU3F3','SOX2',
                   'SOX15','SOX8','SOX9')
write(motif_short_list,'../downstream/output/motif_sig_short.txt')
SNP_motif_high_ent=subsetByOverlaps(variant_HetCpG_meta,motif_gene[motif_gene$geneSymbol%in%motif_short_list])
SNP_motif=subsetByOverlaps(variant_HetCpG_meta,motif_gene)

SNP_motif_high_ent=SNP_motif_high_ent[which(sign(SNP_motif_high_ent$altNME-SNP_motif_high_ent$refNME)==sign(SNP_motif_high_ent$alleleDiff))]
SNP_motif_high_ent=SNP_motif_high_ent[abs(SNP_motif_high_ent$alleleDiff)>=1.5]
merge_motif_high_ent=GRanges()
for (subj in unique(GR_merge$Subject)){
  merge_motif_high_ent=c(merge_motif_high_ent,subsetByOverlaps(GR_merge[GR_merge$Subject==subj],SNP_motif_high_ent[SNP_motif_high_ent$sub==subj]))
  
}
merge_motif_all=GRanges()
for (subj in unique(GR_merge$Subject)){
  merge_motif_all=c(merge_motif_all,subsetByOverlaps(GR_merge[GR_merge$Subject==subj],SNP_motif[SNP_motif$sub==subj]))
  
}
write(unique(unlist(SNP_motif_high_ent$genes_promoter)),'../downstream/output/SNP_motif_high_ent_promoter.txt')
write( unique(unlist(variant_HetCpG_meta$genes_promoter)),'../downstream/output/SNP_all_promoter.txt')
write(unique(unlist(SNP_motif_high_ent$genes_promoter[SNP_motif_high_ent$dNME_pval<=pval_cutoff])),'../downstream/output/SNP_motif_high_ent_promoter_dNME.txt')
write(unique(unlist(variant_HetCpG_meta$genes_promoter[variant_HetCpG_meta$dNME_pval<=pval_cutoff])),'../downstream/output/SNP_dNME_promoter.txt')

export.bed(unique(merge_motif_high_ent),'../downstream/output/motif_region_selected.bed')
export.bed(unique(merge_motif_high_ent[merge_motif_high_ent$dNME_pval<=pval_cutoff]),'../downstream/output/motif_region_dNME.bed')
export.bed(unique(GR_merge[GR_merge$N>1]),'../downstream/output/all_regions.bed')
export.bed(unique(GR_merge[GR_merge$dNME_pval<=pval_cutoff]),'../downstream/output/all_regions_dNME.bed')
export.bed(unique(merge_motif_all),'../downstream/output/motif_region_all.bed')
job = submitGreatJob(as.data.frame(unique(merge_motif_high_ent[merge_motif_high_ent$dNME_pval<=pval_cutoff]))[,1:3],
                     as.data.frame(unique(merge_motif_all))[,1:3],species = 'hg19')

GO_high_ent_motif=GO_anno(unique(SNP_motif_high_ent$genes_promoter),
                          unique(variant_HetCpG_meta$genes_promoter))
GO_high_ent_motif_term=GO_high_ent_motif[[1]]
GO_high_ent_motif_term=GO_high_ent_motif_term[GO_high_ent_motif_term$Significant>=10&GO_high_ent_motif_term$classicFisher<=0.05,c(2,6)]
#Check overlap with Ken data
motif_site_high_NME_promoter <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/motif_site_high_NME_promoter.rds")
motif_Ken
motif_Ken=lapply(motif_site_high_NME_promoter,function(x) 
  unlist(lapply(strsplit(names(x),'_'),function(x) x[2])))


#check motif inside family
#motif_dir_sig$TF[motif_dir_sig$qval_binom<=0.1&motif_dir_sig$TF%in% human_mono_motif_TSV$Transcription.factor[human_mono_motif_TSV$TF.family=='Forkhead box (FOX) factors{3.3.1}']]

motif_gene_sig_dir=motif_gene[motif_gene$geneSymbol%in%motif_dir_sig_qval$TF]
#motif_gene_sig_dir=motif_gene[motif_gene$geneSymbol=='CTCF']
#variant_HetCpG_meta_dir=subsetByOverlaps(variant_HetCpG_meta,motif_gene_sig_dir)
GR_merge_dir=subsetByOverlaps(GR_merge,motif_gene_sig_dir[abs(motif_gene_sig_dir$alleleDiff)>=1.5])
GR_merge_dir=subsetByOverlaps(GR_merge_dir,DNase_all)
write.csv(motif_dir_sig_qval$TF[motif_dir_sig_qval$prob>0.5],'../downstream/output/motif_prefer_ent_no_OR.csv')
write.csv(motif_dir_sig_qval$TF[motif_dir_sig_qval$prob<0.5],'../downstream/output/motif_not_prefer_ent_no_OR.csv')

#trait_qvalue=lapply(unique(dNME_traits$trait),function(x) any(dNME_traits$qvalue[dNME_traits$trait==x]<=0.1))
#trait_sig=unique(dNME_traits$trait)[unlist(trait_qvalue)]
#dNME_traits_sig=dNME_traits[dNME_traits$trait%in%trait_sig,]
dNME_traits_sig$OR[dNME_traits_sig$qvalue>0.1]=NA
dNME_traits$OR[dNME_traits$qvalue>0.1]=NA
pdf('../downstream/output/dNME_traits_heatmap.pdf')
ggplot(dNME_traits[dNME_traits$qvalue<=0.1,],aes(subject,trait,fill=log(OR)))+geom_tile()+
  scale_fill_gradient2(low="navy", mid="white", high="red", midpoint=0,na.value = 'grey') +
  xlab('Sample')+ylab('trait')+
  theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=16))
dev.off()

pdf('../downstream/output/dNME_traits_heatmap.pdf')
ggplot(dNME_traits,aes(sp,trait,fill=p_value))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
  xlab('Sample')+ylab('trait')+
  theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=6))
dev.off()


pdf('../downstream/output/dMML_traits_heatmap.pdf')
ggplot(dMML_traits,aes(sp,trait,fill=log(OR)))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
  xlab('Sample')+ylab('trait')+
  theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size=6))
dev.off()
#Prepare for IPA analysis
genomic_features=readRDS(genomic_features_file)
olap=findOverlaps(genomic_features$promoter,GR_merge)
dNME_IPA=data.frame(genes=genomic_features$promoter$gene_name[queryHits(olap)],dNME=-GR_merge$dNME[subjectHits(olap)],
                    dNME_FDR=GR_merge$dNME_pval[subjectHits(olap)],NME_FC=-(GR_merge$NME2/GR_merge$NME1)[subjectHits(olap)],stringsAsFactors = F)
dNME_IPA_ag=aggregate(dNME_IPA[,c(2,3,4)],by=list(dNME_IPA$genes),FUN=min)
colnames(dNME_IPA_ag)=c('gene_name','dNME','FDR','FC')
dNME_IPA_ag$dNME=-dNME_IPA_ag$dNME
dNME_IPA_ag$FC=log2(-dNME_IPA_ag$FC)
write.csv(dNME_IPA_ag,'../downstream/output/dNME_IPA.csv')
#check SNPs inside the motif

enchancer_DNase_gr=readRDS("../downstream/input/enchancer_DNase.rds")
prom_DNase_gr=readRDS("../downstream/input/prom_DNase.rds")
DNase_all=c(prom_DNase_gr,enchancer_DNase_gr)
DNase_all=makeGRangesFromDataFrame(readRDS('../downstream/output/DNase_data_ken.rds'),keep.extra.columns = T)
DNase_all=DNase_all[which(abs(DNase_all$DH_gene_cor_spearman)>=0.1&abs(DNase_all$dist_TSS)<=500)]
motif_gene=readRDS('../downstream/output/motif_all_JASPAR_default.rds')
#motif_gene=motif_gene[abs(motif_gene$alleleDiff)>=1.5]
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
variant_HetCpG_meta=subsetByOverlaps(variant_HetCpG_meta,DNase_all)
motif_gene=subsetByOverlaps(motif_gene,variant_HetCpG_meta,type='equal')
olap=findOverlaps(variant_HetCpG_meta,motif_gene)
variant_HetCpG_meta=variant_HetCpG_meta[queryHits(olap)]
variant_HetCpG_meta$alleleDiff=motif_gene$alleleDiff[subjectHits(olap)]
variant_HetCpG_meta$NME_diff=variant_HetCpG_meta$altNME-variant_HetCpG_meta$refNME
variant_HetCpG_meta$MML_diff=variant_HetCpG_meta$altMML-variant_HetCpG_meta$refMML
variant_HetCpG_meta=variant_HetCpG_meta[!is.na(variant_HetCpG_meta$alleleDiff)]
sum(sign(variant_HetCpG_meta$NME_diff[variant_HetCpG_meta$dNME_pval<=pval_cutoff])==
      sign(variant_HetCpG_meta$alleleDiff[variant_HetCpG_meta$dNME_pval<=pval_cutoff]))/length(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff])
sum(sign(variant_HetCpG_meta$MML_diff[variant_HetCpG_meta$dMML_pval<=pval_cutoff])==
      sign(variant_HetCpG_meta$alleleDiff[variant_HetCpG_meta$dMML_pval<=pval_cutoff]))/length(variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval<=pval_cutoff])

N_cutoff=1
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
motif_dir$qval_binom=p.adjust(motif_dir$binom.pval,method='BH')
motif_sig=motif_gene[motif_gene$geneSymbol %in% motif_dir$TF[motif_dir$qval_binom<=0.1]]
variant_motif_sig=subsetByOverlaps(variant_HetCpG_meta,motif_sig)
variant_motif_sig_dNME=variant_motif_sig[variant_motif_sig$dNME_pval<=pval_cutoff]
variant_motif_sig_dNME_same_dir=variant_motif_sig_dNME[sign(variant_motif_sig_dNME$NME2-variant_motif_sig_dNME$NME1) == sign(variant_motif_sig_dNME$alleleDiff)]
variant_motif_sig_dNME_same_dir=variant_motif_sig_dNME_same_dir[variant_motif_sig_dNME_same_dir$N>=N_cutoff]

tb=table(unlist(variant_motif_sig_dNME_same_dir$trait))
top_traits=tb[tb>=1]
trait=as.list(names(top_traits))
perm_length=length(variant_motif_sig_dNME_same_dir)
sample_length=length(variant_HetCpG_meta[variant_HetCpG_meta$N>=N_cutoff])
trait_perm=variant_HetCpG_meta$trait[variant_HetCpG_meta$N>=N_cutoff]
trait_perm_out=data.frame()
nperm=10000
for(i in 1:nperm){
  trait_perm_name=unlist(trait_perm[sample(1:sample_length,perm_length,replace=F)])
  trait_perm_length=t(as.data.frame(unlist(lapply(trait,function(x) sum(trait_perm_name %in% x)))))
  rownames(trait_perm_length)=NA
  colnames(trait_perm_length)=trait
  trait_perm_out=rbind(trait_perm_out,trait_perm_length)
}
trait_pval_perm=data.frame()
for(tt in trait){
  trait_pval_perm=rbind(trait_pval_perm,data.frame(trait=tt,GWAS_number=top_traits[tt],pvalue=sum(trait_perm_out[,tt]>=top_traits[tt])/nperm,stringsAsFactors = F))
  
  
}
trait_pval_perm$qvalue=p.adjust(trait_pval_perm$pvalue,method='BH')
trait_pval_perm_sig=trait_pval_perm[trait_pval_perm$pvalue<=0.05,]
trait_pval_perm_disease=c('acute lymphoblastic leukemia','acute myeloid leukemia','adolescent idiopathic scoliosis',
                          'amyotrophic lateral sclerosis','anterior uveitis','anxiety','attention deficit hyperactivity disorder',
                          'autism spectrum disorder','bipolar disorder','bladder carcinoma','borderline personality disorder',
                          'breast carcinoma','bronchopulmonary dysplasia','cancer','cardiovascular disease','cataract','celiac disease',
                          'chronic kidney disease','chronic lymphocytic leukemia','chronic obstructive pulmonary disease','colorectal cancer',
                          'congenital heart disease','coronary artery disease','coronary heart disease','drug-induced agranulocytosis','erectile dysfunction',
                          'gallstones','HIV-1 infection','HIV infection','hypertension','hypothyroidism','insomnia measurement','kidney disease','keratinocyte carcinoma',
                          'longevity','lung carcinoma','lupus nephritis','lymphoma','mood instability measurement','multiple myeloma','multiple sclerosis','multisite chronic pain',
                          'myelodysplastic syndrome','obesity','ototoxicity','ovarian carcinoma','pancreatic carcinoma','prostate carcinoma','renal cell carcinoma',
                          'respiratory system disease','response to antidepressant','response to antineoplastic agent','response to antipsychotic drug','rheumatoid arthritis',
                          'schizoaffective disorder','schizophrenia','sensory perception of sweet taste','sickle cell anemia','squamous cell lung carcinoma','sporadic amyotrophic lateral sclerosis',
                          'stroke','sudden cardiac arrest','susceptibility to infectious disease measurement','systemic lupus erythematosus','testicular carcinoma','thyroid carcinoma',
                          'Trypanosoma cruzi seropositivity','unipolar depression')
trait_pval_perm_disease=trait_pval_perm[trait_pval_perm$trait%in%trait_pval_perm_disease,]
variant_motif_sig_dNME_same_dir_traits=variant_motif_sig_dNME_same_dir[-which(unlist(lapply(variant_motif_sig_dNME_same_dir$trait,is.null)))]
variant_motif_sig_dNME_same_dir_traits=subsetByOverlaps(variant_motif_sig_dNME_same_dir_traits,DNase_all,maxgap=500)
variant_motif_sig_dNME_same_dir_traits[order(variant_motif_sig_dNME_same_dir_traits$N,decreasing=T)]
variant_select=variant_motif_sig_dNME_same_dir_traits[which(unlist(lapply(variant_motif_sig_dNME_same_dir_traits$trait,function(x) any(trait_pval_perm_disease %in% x))))]
elementMetadata(variant_select)=elementMetadata(variant_select)[,c('dNME_pval','dNME','dMML_pval','dMML','genes_promoter','genes_body','TSS','trait','germlayer','alleleDiff','Sample','N')]
variant_select=variant_select[order(variant_select$N,decreasing = T)]
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)

vQTL_SNP[vQTL_SNP$SNP%in%unique(variant_HetCpG_meta$variant_id[variant_HetCpG_meta$dNME_pval<=pval_cutoff]),]

variant_HetCpG_meta_vQTL=variant_HetCpG_meta[variant_HetCpG_meta$variant_id%in%vQTL_SNP$SNP]
#Total 5 in 91 have dNME

variant_trait <- readRDS("../downstream/input/genome_1k_variant.rds")
elementMetadata(variant_trait)=elementMetadata(variant_trait)[c('Dbxref','Reference_seq','Variant_seq')]
variant_trait$rsid=unlist(lapply(strsplit(unlist(variant_trait$Dbxref),'\\:'),function(x) x[2]))
association_variant_df_all=readRDS('../downstream/output/association_variant_df_all.rds')
variant_trait$trait=association_variant_df_all$traits[match(variant_trait$rsid,association_variant_df_all$variant_id)]
variant_trait=variant_trait[unlist(lapply(variant_trait$trait,function(x) length(x)>0))]
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
elementMetadata(variant_HetCpG_meta)=elementMetadata(variant_HetCpG_meta)[,c('trait','dNME_pval','dMML_pval')]

T2D_trait=read.table('../downstream/output/T2D_ALL_Primary.txt.gz',header=T)
T2D_trait$start=T2D_trait$Pos
T2D_trait$end=T2D_trait$Pos
T2D_trait=makeGRangesFromDataFrame(T2D_trait,keep.extra.columns = T)
T2D_trait$qvalue=p.adjust(T2D_trait$P,method='BH')
seqlevels(T2D_trait)=paste('chr',seqlevels(T2D_trait),sep='')
variant_HetCpG_meta_sub=subsetByOverlaps(variant_HetCpG_meta,T2D_trait,maxgap = 1000)
T2D_trait_sub=T2D_trait[T2D_trait$qvalue<=0.1]
get_traits_GWAS_all_trait_single(variant_HetCpG_meta_sub,T2D_trait_sub,0.1,5,'dNME_pval',FALSE,maxgap=10000)

#breast cancer:GCST004988
# BC=read.csv('../downstream/output/breast_cancer_risk.csv')
# BC=makeGRangesFromDataFrame(BC,seqnames.field='chrom',start.field='position',end.field = 'position')
# seqlevels(BC)=paste('chr',seqlevels(BC),sep='')
BC=get_variants(study_id = 'GCST004988')
BC=as.data.frame(BC@variants)
BC=BC[!is.na(BC$chromosome_position),]
BC=makeGRangesFromDataFrame(BC, seqnames.field = 'chromosome_name',start.field = 'chromosome_position',end.field = 'chromosome_position')
seqlevels(BC)=paste('chr',seqlevels(BC),sep='')
get_traits_GWAS_all_trait_single(variant_HetCpG_meta,BC,0.1,5,'dNME_pval',FALSE,maxgap=10000)

#inflammatory bowl disease
BD=read.csv('../downstream/output/inflammatory bowl disease loci.csv')
#BD=BD[!is.na(BD$credible_end),]
BD=makeGRangesFromDataFrame(BD,seqnames.field='chr',start.field='credible_start',end.field = 'credible_end')
seqlevels(BD)=paste('chr',seqlevels(BD),sep='')
# BD=get_variants(study_id = 'GCST005837')
# BD=makeGRangesFromDataFrame(as.data.frame(BD@variants),
#                                            seqnames.field = 'chromosome_name',start.field = 'chromosome_position',end.field = 'chromosome_position')
# seqlevels(BD)=paste('chr',seqlevels(BD),sep='')


get_traits_GWAS_all_trait_single(variant_HetCpG_meta,BD,0.1,5,'dNME_pval',FALSE,maxgap=1000)

#pediatric autoimmue disease
PAD=get_variants(study_id = 'GCST003097')
variant_trait_PAD=makeGRangesFromDataFrame(as.data.frame(PAD@variants),
                                           seqnames.field = 'chromosome_name',start.field = 'chromosome_position',end.field = 'chromosome_position')
seqlevels(variant_trait_PAD)=paste('chr',seqlevels(variant_trait_PAD),sep='')
get_traits_GWAS_all_trait_single(variant_HetCpG_meta,variant_trait_PAD,0.1,5,'dNME_pval',FALSE,maxgap=10000)
#Birth weight
BW=get_variants(study_id = 'GCST005146')
variant_trait_BW=makeGRangesFromDataFrame(as.data.frame(BW@variants),
                                          seqnames.field = 'chromosome_name',start.field = 'chromosome_position',end.field = 'chromosome_position')
seqlevels(variant_trait_BW)=paste('chr',seqlevels(variant_trait_BW),sep='')
get_traits_GWAS_all_trait_single(variant_HetCpG_meta,variant_trait_BW,0.1,5,'dNME_pval',FALSE,maxgap=1000)


genome_freq=readRDS('../downstream/input/genome_1k_variant.rds')
association_db=readRDS('../downstream/input/association_db.rds')
genome_freq$variant_id=gsub('dbSNP_153\\:','',genome_freq$Dbxref)
genome_freq$variant_id=unlist(genome_freq$variant_id)
variant_all=granges(genome_freq)
variant_all$variant_id=genome_freq$variant_id
association_df=aggregate(association_db@risk_alleles$association_id,by=list(association_db@risk_alleles$variant_id),list)
association_df$x=lapply(association_df$x,as.character)
colnames(association_df)=c('variant_id','association_id')
variant_all$association_id=association_df$association_id[match(variant_all$variant_id,association_df$variant_id)]
rm(genome_freq)
variant_all_na=lapply(variant_all$association_id,is.null)
variant_all=variant_all[-which(unlist(variant_all_na))]
variant_all$traits=mclapply(variant_all$association_id,function(x) {
  traits=get_traits(association_id=x)@traits$trait
  return(traits)
},mc.cores=20)  
proc.time()[[3]] - tt1
GR_merge=readRDS(GR_merge_file)
variant_all=readRDS("../downstream/output/variant_all.rds")
variant_all$trait_length=unlist(lapply(variant_all$traits,length))
variant_all=variant_all[variant_all$trait_length>1]
olap_GWAS=findOverlaps(GR_merge,variant_all,maxgap=1000)
olap_GWAS_df=data.frame(qt=queryHits(olap_GWAS),st=subjectHits(olap_GWAS))
olap_GWAS_agg=aggregate(olap_GWAS_df$st,by=list(olap_GWAS_df$qt),c)
colnames(olap_GWAS_agg)=c('qt','st')
GR_merge$trait=NA
GR_merge$trait[olap_GWAS_agg$qt]=lapply(olap_GWAS_agg$st,function(x) unlist(variant_all$traits[x]))
GR_merge$trait_length=0
GR_merge$trait_length[olap_GWAS_agg$qt]=lapply(olap_GWAS_agg$st,function(x) unlist(variant_all$trait_length[x]))
GR_merge$rsid=NA
GR_merge$rsid[olap_GWAS_agg$qt]=lapply(olap_GWAS_agg$st,function(x) unlist(variant_all$variant_id[x]))
GR_merge=tissue_to_germlayer(GR_merge)
saveRDS(GR_merge,'../downstream/output/GR_merge_traits_1k_more_than_1_trait.rds')
#variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
olap=findOverlaps(unique(genome_freq),variant_HetCpG_meta,type='equal')
variant_HetCpG_meta$variant_id[subjectHits(olap)]=as.character(unique(genome_freq)$variant_id)[queryHits(olap)]

# Add GWAS annotation to each SNP -----------------------------------------


association_variant_df_all=data.frame(variant_id=unique(genome_freq$variant_id),
                                      stringsAsFactors = F)
association_variant_df_all$association_id=association_db@risk_alleles$association_id[match(association_variant_df_all$variant_id,
                                                                                           association_db@risk_alleles$variant_id)]
association_variant_df_all=association_variant_df_all[!is.na(association_variant_df_all$association_id),]
association_variant_df_all$traits=NA
association_variant_df_all$traits=lapply(association_variant_df_all$association_id,function(x) {
  traits=get_traits(association_id=x)@traits$trait
  return(traits)
})
saveRDS(association_variant_df_all,'../downstream/output/association_variant_df_all.rds')
association_variant_df_all=readRDS('../downstream/output/association_variant_df_all.rds')

# Plot NME vs MML heatmap---------------------------------------------------------
NME_df=data.frame(NME=c(GR_merge$NME1,GR_merge$NME2),
                  MML=c(GR_merge$MML1,GR_merge$MML2))
NME_df$NME=round(NME_df$NME,digits=4)
NME_df$MML=round(NME_df$MML,digits=4)
NME_df_agg=aggregate(NME_df,by=list(NME_df$NME,NME_df$MML),FUN=length)
colnames(NME_df_agg)=c('NME','MML','count','count2')
pdf('../downstream/output/MML_NME_heatmap_density2d.pdf')
ggplot(as.data.frame(NME_df), aes(MML,NME)) +
  #geom_tile() +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE)+
  scale_fill_gradient(low='lightblue',high='blue')+
  theme_bw() +theme(legend.position = 'bottom')+
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12))
dev.off()
pdf('../downstream/output/MML_NME_heatmap_sub_10k_11067.pdf')
p2 <- ggplot(NME_df[sample(1:nrow(NME_df),nrow(NME_df)/1000,replace = F),], aes(MML,NME)) +
  geom_pointdensity() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(axis.text = element_text(size=1), axis.title = element_text(size=1))
print(p2)
dev.off()

#GR_merge_df=GR_merge_df[(GR_merge_df$dMML_pval<=pval_cutoff |GR_merge_df$dNME_pval<=pval_cutoff),]
# ###Plotting without 95% range
# pdf('../downstream/output/graphs/Figure1/dNME_vs_dMML_no_quantile_differential.pdf')
# ggplot(GR_merge_df[GR_merge_df$dNME_pval<=pval_cutoff |GR_merge_df$dMML_pval<=pval_cutoff,],aes(x=dMML, y=dNME))+
#   xlim(c(0,1))+ggtitle("dMML and dNME relationship")+geom_smooth(size=1,se=TRUE)+ylim(c(0,1))+
#   ylab("dNME")+theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))#+geom_point()
# dev.off()
# GR_merge_df=GR_merge_df[GR_merge_df$dNME_pval<=pval_cutoff |GR_merge_df$dMML_pval<=pval_cutoff,]

# #Backup
# # ggplot(GR_df_agg,aes(x=dMML, y=median))+
# #   xlim(c(0,1))+ylim(c(0,1))+ggtitle("dMML and dNME relationship")+geom_smooth(method="loess",se=FALSE,aes(linetype=variable))+
# #   ylab("dNME")+theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
# #scale_linetype_manual(values=c("solid","twodash", "twodash"))+scale_color_manual(values=c("Blue","Blue","Blue"))
# #method = "loess", formula = y ~ x, size = 1

# GR_calc_agg=aggregate(GR_calc$dMML_ASM,by=list(GR_calc$sample),FUN=sum)
# GR_calc_nondMML=GR_calc[!GR_calc$dMML_ASM,]
# GR_calc_dMML_plot=rbind(GR_calc[GR_calc$dMML_ASM,],GR_calc_nondMML[sample(1:nrow(GR_calc_nondMML),sum(GR_calc$dMML_ASM)),])
# ggplot(GR_calc_dMML_plot,aes(x=NME,color=dMML_ASM))+geom_density(size=1)+labs(color='MML-Hap')+theme(legend.position = 'bottom')
# 
# 
# ggplot(GR_calc,aes(x=NME,color=dNME_ASM))+geom_density(size=1)+labs(color='NME-Hap')+theme(legend.position = 'bottom')
# ggplot(GR_calc,aes(x=NME))+geom_density(size=1)+theme(legend.position = 'bottom')
# GR_calc_nondMML_nondNME=GR_calc[!GR_calc$dNME_ASM_non_dMML,]
# GR_calc_nondMML_dNME_plot=rbind(GR_calc[GR_calc$dNME_ASM_non_dMML,],GR_calc_nondMML_nondNME[sample(1:nrow(GR_calc_nondMML_nondNME),sum(GR_calc$dNME_ASM_non_dMML)),])
# ggplot(GR_calc_nondMML_dNME_plot,aes(x=NME,color=dNME_ASM_non_dMML))+geom_density(size=1)+labs(color='dNME-ASM but not dMML-ASM')+theme(legend.position = 'bottom')
# ggplot(GR_calc,aes(x=dNME,y=NME))+geom_smooth()#+geom_point(alpha=0.1)
# ggplot(GR_calc,aes(x=dMML,y=NME))+geom_smooth()#+geom_point(alpha=0.1)


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

# heatmap for all tissue --------------------------------------------------
#clustering using specific columns
DNase=readRDS('../downstream/output/mm10_DNase.rds')
UC_in_Epiblast=subsetByOverlaps(UC_in_Epiblast,DNase[DNase$region_type=='DNase'],type='equal')
# mclapply(tissue_all,row_sd_tissue,UC_in_mt=UC_in_mt,mc.cores=11)
# row_sd_tissue<-function(tissue,UC_in_mt){
#   UC_in_mt_tissue=UC_in_mt[,unlist(lapply(strsplit(colnames(UC_in_mt),'-'),function(x) x[1]))==tissue]
#   return(quantile(rowSds(UC_in_mt_tissue),prob=0.95))
# }

heatmap_output=lapply(tissue,heatmap_tissue,UC_in=UC_in_Epiblast,sd_cutoff=0.05,JSD_cutoff=0.15,plot_scale=T,dist_scale=T)
heatmap_tissue<-function(tissue,UC_in,sd_cutoff=0.05,plot_scale=F,dist_scale=T,JSD_cutoff=0.1,
                         header_pic='../downstream/output/UC_tissue_'){
  cat('Start processing:',tissue,'\n')
  tt1=proc.time()[[3]]
  #Getting tissue information
  tissue_all=unlist(lapply(strsplit(colnames(mcols(UC_in)),'-'),function(x) x[1]))
  col_break=cumsum(table(tissue_all)[unique(tissue_all)])
  #filter regions with high row sd and JSD
  UC_in_tissue=granges(UC_in)
  mcols(UC_in_tissue)=mcols(UC_in)[,unlist(lapply(strsplit(colnames(mcols(UC_in)),'-'),function(x) x[1]))==tissue]
  
  UC_in_tissue=UC_in_tissue[which(apply(mcols(UC_in_tissue),1, function(x) (any(x>=JSD_cutoff)))&
                                    rowSds(as.matrix(mcols(UC_in_tissue)))>=sd_cutoff)]
  UC_in=subsetByOverlaps(UC_in,UC_in_tissue,type='equal')
  
  cat('total regions pass filter:',length(UC_in),'\n')
  UC_in_mt=as.matrix(mcols(UC_in))
  
  #scale is necessary
  if(plot_scale){
    UC_in_names=colnames(UC_in_mt)
    UC_in_mt=fastDoCall("cbind",lapply(unique(tissue_all),function(x) {
      out=t(scale(t(UC_in_mt[,which(unlist(lapply(strsplit(colnames(UC_in_mt),'-'),function(y) y[1]))==x)])))
      attributes(out)[c('scaled:center','scaled:scale')] <- NULL
      out=apply(out,2,function(y) {y[which(is.na(y))]=0
      return(y)})
      return(out)}))
    cat('finishing scaling in:',proc.time()[[3]]-tt1,'\n')
  }
  UC_in_tissue=as.matrix(mcols(UC_in_tissue))
  if(dist_scale){UC_in_tissue=t(scale(t(UC_in_tissue)))}
  
  #kmeans clustering
  kmeans_out=kmeans(UC_in_tissue,center=10,nstart=50,iter.max = 100,algorithm ="Hartigan-Wong")
  #plotting
  png(paste(header_pic,tissue,'.png',sep=''),height=600,width = 800)
  heatmap_out=pheatmap(UC_in_mt[order(kmeans_out$cluster,decreasing=F),],cluster_cols = FALSE,cluster_rows = FALSE,scale='none',
                       col=col_fun(1024), gaps_row=as.vector(cumsum(table(sort(kmeans_out$cluster,decreasing=F)))),
                       gaps_col=as.vector(col_break),main=tissue)
  print(heatmap_out)
  dev.off()
  cat('finishing in:',proc.time()[[3]]-tt1,'\n')
  return(list(granges(UC_in),kmeans_out,heatmap_out))
}

sink("../downstream/output/UC_submission_all_pair_sub")
cat('#!/bin/bash\n')
for(ts in unique(meta_df$tissue)){
  meta_df_ts=meta_df[meta_df$tissue==ts,]
  stage_num=nrow(meta_df_ts)
  meta_df_ts$sample=paste('mm10',meta_df_ts$tissue,meta_df_ts$stage,sep='_')
  for(i in 1:stage_num){
    #cat i+1 through stage num
    if(i<stage_num-1){
      for(j in (i+2):stage_num){
        cat('sbatch cpelasm_allele_agnostic_uc.slrm ',meta_df_ts$sample[i],'_merged1 ',meta_df_ts$sample[j],'_merged1\n',sep='')
        cat('sbatch cpelasm_allele_agnostic_uc.slrm ',meta_df_ts$sample[i],'_merged2 ',meta_df_ts$sample[j],'_merged2\n',sep='')
      }
    }
    
  }
  
}
sink()
sink("../downstream/output/UC_submission_all_BL6DBA")
cat('#!/bin/bash\n')
meta_df$sample=paste('mm10',meta_df$tissue,meta_df$stage,sep='_')
for(sp in meta_df$sample){
  cat('sbatch cpelasm_allele_agnostic_uc_BL6DBA.slrm ','BL6DBA_ICM_all ',sp,'_all\n',sep='')
  cat('sbatch cpelasm_allele_agnostic_uc_BL6DBA.slrm ','BL6DBA_Epiblast_all ',sp,'_all\n',sep='')
  
  
}
sink()


# loading the varibility data for each sample -----------------------------
# regions with high NME at one allele is more likely to have hyper vari or with preferential motif binding
GR_merge$hyper_var_promoter=NA
GR_merge$hyper_var_TSS=NA
GR_merge$hyper_var_body=NA
GR_merge$hyper_var_lower=NA
GR_merge$hyper_var_upper=NA
GR_merge$hyper_var_median=NA
GR_merge$tissue[GR_merge$tissue=='Adipose_Tissue_single']='Adipose_single'
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/HESC_1.rds'),
                       c('42_embryonic_stem_cell_single - H9','stem_27_undifferentiated_paired - HUES64','merged - H1'))
#GR_merge=add_hyper_var(GR_merge,'../downstream/output/H1_hypervar.rds',c('merged - H1'))
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultAdipose_1.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Adipose_single']))
#../downstream/input/scRNA/AdultBladder_2.rds 
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultBladder_1.rds','../downstream/input/scRNA/AdultBladder_2.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Bladder_single']))
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultIleum_2.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Small_Intestine_single']))
#../downstream/input/scRNA/AdultStomach_3.rds 
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultStomach_2.rds','../downstream/input/scRNA/AdultStomach_1.rds',
                                  '../downstream/input/scRNA/AdultStomach_3.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Gastric_single']))
#../downstream/input/scRNA/AdultHeart_2.rds
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultHeart_1.rds','../downstream/input/scRNA/AdultHeart_2.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Left_Ventricle_single']))
#../downstream/input/scRNA/AdultLung_3.rds 
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultLung_1.rds','../downstream/input/scRNA/AdultLung_2.rds',
                                  '../downstream/input/scRNA/AdultLung_3.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Lung_single']))
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultMuscle_1.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Psoas_Muscle_single']))
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultColon_1.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Sigmoid_Colon_single']))
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultSpleen_1.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Spleen_single']))
#../downstream/input/scRNA/AdultAdrenalGland_3.rds 
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultAdrenalGland_2.rds','../downstream/input/scRNA/AdultAdrenalGland_3.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Adrenal_Gland_single']))
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultArtery_1.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Aorta_single']))
#../downstream/input/scRNA/AdultEsophagus_2.rds 
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultEsophagus_1.rds','../downstream/input/scRNA/AdultEsophagus_2.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Esophagus_single']))
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultPancreas_1.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Pancreas_single']))
#../downstream/input/scRNA/AdultLiver_2.rds 
GR_merge=add_hyper_var(GR_merge,c('../downstream/input/scRNA/AdultLiver_1.rds','../downstream/input/scRNA/AdultLiver_2.rds',
                                  '../downstream/input/scRNA/AdultLiver_4.rds'),
                       unique(GR_merge$Sample[GR_merge$tissue=='Liver_single']))
