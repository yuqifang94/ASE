source('mainFunctions_sub.R')
# GWAS analysis -----------------------------------------------------------
#Check GWAS traits have hg19 coordinates
variant_trait=readRDS('../downstream/input/human_analysis/variant_traits.rds')
variant_trait_gr=GRanges(variant_trait)
# variant_trait_gr=do.call(c,variant_trait)
# variant_trait_gr=unique(variant_trait_gr)
# genome_freq=readRDS('../downstream/input/genome_1k_variant.rds')
# genome_freq$variant_id=gsub('dbSNP_153\\:','',genome_freq$Dbxref)
# olap=findOverlaps(variant_trait,genome_freq,type='equal')
#Check proportion of dNME SNP from GWAS
#get traits location
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
#variant_HetCpG_meta=readRDS(GR_merge_file)
mcols(variant_HetCpG_meta)=mcols(variant_HetCpG_meta)[,c('dNME','dMML','dNME_pval','dMML_pval')]
length(subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],variant_trait_gr,maxgap = 0))
olap=findOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],variant_trait_gr)
variant_HetCpG_meta$GWAS=NA
variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff][queryHits(olap)]$GWAS=variant_trait_gr[subjectHits(olap)]$`DISEASE/TRAIT`


#  23827/380526= 0.06261596
#dMML:1236/14006=0.088
length(subsetByOverlaps(variant_HetCpG_meta,variant_trait_gr,maxgap = 1000))# 601586/7967588=0.076
length(subsetByOverlaps(variant_HetCpG_meta[!is.na(variant_HetCpG_meta$genes_promoter)],variant_trait,maxgap = 1000))# 7264/230407=0.0315, #ntrait>1-> 
#735/230407=0.0032
length(subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff&!is.na(variant_HetCpG_meta$genes_promoter)],
                        variant_trait))#180/6351=0.028,ntrait>1->11/6351=0.0017
#check gwas near motif
motif_dir=readRDS('../downstream/output/motif_dirction_all_JASPAR_default.rds')
#motif_dir=readRDS('../downstream/output/motif_dir.rds')
motif_dir$qvalue=p.adjust(motif_dir$binom.pval,method='BH')
motif_dir_sig_ent=motif_dir[motif_dir$qval_binom<=0.1&motif_dir$prob>0.5,]
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene=readRDS(motif_gene_file)
motif_gene=motif_gene[motif_gene$geneSymbol%in%motif_dir_sig_ent$TF]
olap=findOverlaps(variant_HetCpG_meta,motif_gene)
variant_HetCpG_meta$Hign_ent_motif=FALSE
variant_HetCpG_meta$Hign_ent_motif[queryHits(olap)]=TRUE
variant_HetCpG_meta$dNME_pval[!variant_HetCpG_meta$Hign_ent_motif]=1
#GWAS enrichment at different ranges, mc.core settings
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=500)
saveRDS(dNME_traits,'../downstream/output/human_analysis/GWAS/dNME_traits_500_gr.rds')
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=1000)
saveRDS(dNME_traits,'../downstream/output/human_analysis/GWAS/dNME_traits_1000_motif_SNP.rds')
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=5000)
saveRDS(dNME_traits,'../downstream/output/human_analysis/GWAS/dNME_traits_5k_gr.rds')
dNME_traits=get_traits_GWAS(variant_HetCpG_meta,variant_trait,CMH=FALSE,maxgap=10000)
saveRDS(dNME_traits,'../downstream/output/human_analysis/GWAS/dNME_traits_10k_gr.rds')
dNME_traits=readRDS('../downstream/output/human_analysis/GWAS/dNME_traits_1000_gr2.rds')
#dNME_traits=readRDS('../downstream/output/dNME_traits_1000_motif_SNP.rds')
#downstream analysis
dNME_traits=readRDS('../downstream/output/human_analysis/GWAS/dNME_traits_500_gr.rds')
dNME_traits_df=do.call(rbind,lapply(dNME_traits,function(x) x[[1]]))
dNME_traits_gr=do.call(c,lapply(dNME_traits,function(x) x[[2]]))
dNME_traits_df=dNME_traits_df[dNME_traits_df$dNME_trait>=10&dNME_traits_df$OR>1,]
dNME_traits_df$qvalue=p.adjust(dNME_traits_df$p_value,method='BH')
dNME_traits_sig=dNME_traits_df[dNME_traits_df$qvalue<=0.1,]
dNME_traits_sig=dNME_traits_sig[order(dNME_traits_sig$OR,decreasing=F),]
dNME_traits_sig$trait=factor(dNME_traits_sig$trait,levels=dNME_traits_sig$trait)
write.csv(dNME_traits_sig,'../downstream/output/human_analysis/GWAS/GWAS_dNME.csv')
pdf('../downstream/output/human_analysis/GWAS/GWAS_traits_500_gr.pdf',width=18,height=5)
ggplot(dNME_traits_sig,aes(y=OR,x=trait))+
  geom_bar(stat="identity",color="black",position=position_dodge(0.9),fill='lightblue')+
  #geom_errorbar(aes(ymin=log(lower_CI),ymax=log(upper_CI)),width=0.2,position=position_dodge(0.9))+ 
  coord_flip()+ theme(axis.text.x = element_text(hjust = 1),legend.position = "none")+xlab("GWAS traits")+
  ggtitle("GWAS trait enrichment")+geom_text(aes(label=round(OR,digits = 2)), hjust=-0.5, color="black", size=3.5)
dev.off()