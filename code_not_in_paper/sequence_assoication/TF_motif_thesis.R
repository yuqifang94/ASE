#Plot example for ASCL1
source('mainFunctions_sub.R')
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                 axis.title.x=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.text.x=element_text(size=6),
                                 axis.text.y=element_text(size=6))
motif_gene <- readRDS(motif_gene_file)
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
stat="NME"
tf="ASCL1"
ASCL1_dNME=merge_SNP_motif(variant_HetCpG_meta[mcols(variant_HetCpG_meta)[[paste0('d',stat,'_pval')]]<=pval_cutoff],
                                motif_gene[motif_gene$geneSymbol==tf],stat=stat)
pdf("../downstream/output/human_analysis/motif_analysis/ASCL1_dNME.pdf")
ggplot(as.data.frame(mcols(ASCL1_dNME)),aes(x=alleleDiff,y=stat_diff))+geom_point()+xlab("Binding affinity difference")+ylab("dNME")+
    geom_vline(xintercept=0)+geom_hline(yintercept=0)
dev.off()
sum(ASCL1_dNME$alleleDiff>0&ASCL1_dNME$stat_diff>0)#44
sum(ASCL1_dNME$alleleDiff<0&ASCL1_dNME$stat_diff>0)#9
sum(ASCL1_dNME$alleleDiff<0&ASCL1_dNME$stat_diff<0)#19
sum(ASCL1_dNME$alleleDiff>0&ASCL1_dNME$stat_diff<0)#14

#P-value 

t.test(ASCL1_dNME$refNME,ASCL1_dNME$altNME)
#This can be in the GO analysis discussion
#GO analysis for dNME near TF binding sites
motif_dir_dNME=readRDS('../downstream/output/human_analysis/motif_analysis/dNME_all.rds')
variant_HetCpG_meta_motif=subsetByOverlaps(variant_HetCpG_meta,motif_gene[motif_gene$geneSymbol%in%  motif_dir_dNME[FDR<=0.1]$TF])
variant_HetCpG_meta_motif=convert_GR(variant_HetCpG_meta_motif,direction="DT")

gl=unique(unlist(variant_HetCpG_meta_motif[dNME_pval<=0.1]$genes_body))
back=unique(unlist(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=0.1]$genes_body))
GO_motif_only=GO_run(gl,back,cluster=1,ptcount=0,mapping="org.Hs.eg.db")
#GO for TF binding sites related SNP
gl=unique(unlist(variant_HetCpG_meta_motif[dNME_pval<=0.1]$genes_body))
back=unique(unlist(variant_HetCpG_meta_motif$genes_body))
GO_motif=GO_run(gl,back,cluster=1,ptcount=0,mapping="org.Hs.eg.db")44

#GO for all SNP
gl=unique(unlist(variant_HetCpG_meta_motif[dNME_pval<=0.1]$genes_body))
back=unique(unlist(variant_HetCpG_meta_motif$genes_body))
GO_all=GO_run(gl,back,cluster=1,ptcount=0,mapping="org.Hs.eg.db")

#Methyl-plus and Methyl-minus
methyl_sensitive_gene=as.data.table(read_csv("../downstream/input/human_analysis/motif_analysis/Methyl_call_SELEX.csv"),stringsASfactors=F)
methyl_plus=methyl_sensitive_gene[Call=="MethylPlus"]
methyl_minus=methyl_sensitive_gene[Call=="MethylMinus"]
motif_dir_dNME=readRDS('../downstream/output/human_analysis/motif_analysis/dNME_all.rds')
motif_dir_dMML=readRDS('../downstream/output/human_analysis/motif_analysis/dMML_all.rds')
#dNME methyl-plus
sig_plus=sum(motif_dir_dNME[FDR<=0.1 &proportion_low_NME<=0.5]$TF %in% methyl_plus$`TF name`)
non_sig_plus=sum(motif_dir_dNME[FDR>0.1|proportion_low_NME>0.5]$TF %in% methyl_plus$`TF name`)
sig_non_plus=sum(!motif_dir_dNME[FDR<=0.1&proportion_low_NME<=0.5]$TF %in% methyl_plus$`TF name`)
non_sig_non_plus=sum(!motif_dir_dNME[FDR>0.1|proportion_low_NME>0.5]$TF %in% methyl_plus$`TF name`)
fisher.test(matrix(c(sig_plus,non_sig_plus,sig_non_plus,non_sig_non_plus),nrow=2)) #0.5, P=0.04

sig_plus=sum(motif_dir_dNME[FDR<=0.1 &proportion_low_NME>0.5]$TF %in% methyl_plus$`TF name`)
non_sig_plus=sum(motif_dir_dNME[FDR>0.1|proportion_low_NME<=0.5]$TF %in% methyl_plus$`TF name`)
sig_non_plus=sum(!motif_dir_dNME[FDR<=0.1&proportion_low_NME>0.5]$TF %in% methyl_plus$`TF name`)
non_sig_non_plus=sum(!motif_dir_dNME[FDR>0.1|proportion_low_NME<=0.5]$TF %in% methyl_plus$`TF name`)
fisher.test(matrix(c(sig_plus,non_sig_plus,sig_non_plus,non_sig_non_plus),nrow=2)) #0.36, P=0.48

#dNME methyl-minus
sig_minus=sum(motif_dir_dNME[FDR<=0.1 &proportion_low_NME<=0.5]$TF %in% methyl_minus$`TF name`)
non_sig_minus=sum(motif_dir_dNME[FDR>0.1|proportion_low_NME>0.5]$TF %in% methyl_minus$`TF name`)
sig_non_minus=sum(!motif_dir_dNME[FDR<=0.1&proportion_low_NME<=0.5]$TF %in% methyl_minus$`TF name`)
non_sig_non_minus=sum(!motif_dir_dNME[FDR>0.1|proportion_low_NME>0.5]$TF %in% methyl_minus$`TF name`)
fisher.test(matrix(c(sig_minus,non_sig_minus,sig_non_minus,non_sig_non_minus),nrow=2)) #0.41, P=0.01

sig_minus=sum(motif_dir_dNME[FDR<=0.1 &proportion_low_NME>0.5]$TF %in% methyl_minus$`TF name`)
non_sig_minus=sum(motif_dir_dNME[FDR>0.1|proportion_low_NME<=0.5]$TF %in% methyl_minus$`TF name`)
sig_non_minus=sum(!motif_dir_dNME[FDR<=0.1&proportion_low_NME>0.5]$TF %in% methyl_minus$`TF name`)
non_sig_non_minus=sum(!motif_dir_dNME[FDR>0.1|proportion_low_NME<=0.5]$TF %in% methyl_minus$`TF name`)
fisher.test(matrix(c(sig_minus,non_sig_minus,sig_non_minus,non_sig_non_minus),nrow=2)) #6.51, P=0.0.001

#dMML methyl-plus
sig_plus=sum(motif_dir_dMML[FDR<=0.1 &Proportion_low_MML<=0.5]$TF %in% methyl_plus$`TF name`)
non_sig_plus=sum(motif_dir_dMML[FDR>0.1|Proportion_low_MML>0.5]$TF %in% methyl_plus$`TF name`)
sig_non_plus=sum(!motif_dir_dMML[FDR<=0.1&Proportion_low_MML<=0.5]$TF %in% methyl_plus$`TF name`)
non_sig_non_plus=sum(!motif_dir_dMML[FDR>0.1|Proportion_low_MML>0.5]$TF %in% methyl_plus$`TF name`)
fisher.test(matrix(c(sig_plus,non_sig_plus,sig_non_plus,non_sig_non_plus),nrow=2)) #1.93,0.34

sig_plus=sum(motif_dir_dMML[FDR<=0.1 &Proportion_low_MML>0.5]$TF %in% methyl_plus$`TF name`)
non_sig_plus=sum(motif_dir_dMML[FDR>0.1|Proportion_low_MML<=0.5]$TF %in% methyl_plus$`TF name`)
sig_non_plus=sum(!motif_dir_dMML[FDR<=0.1&Proportion_low_MML>0.5]$TF %in% methyl_plus$`TF name`)
non_sig_non_plus=sum(!motif_dir_dMML[FDR>0.1|Proportion_low_MML<=0.5]$TF %in% methyl_plus$`TF name`)
fisher.test(matrix(c(sig_plus,non_sig_plus,sig_non_plus,non_sig_non_plus),nrow=2)) #0.62, P=0.11

#dMML methyl-minus
sig_minus=sum(motif_dir_dMML[FDR<=0.1 &Proportion_low_MML<=0.5]$TF %in% methyl_minus$`TF name`)
non_sig_minus=sum(motif_dir_dMML[FDR>0.1|Proportion_low_MML>0.5]$TF %in% methyl_minus$`TF name`)
sig_non_minus=sum(!motif_dir_dMML[FDR<=0.1&Proportion_low_MML<=0.5]$TF %in% methyl_minus$`TF name`)
non_sig_non_minus=sum(!motif_dir_dMML[FDR>0.1|Proportion_low_MML>0.5]$TF %in% methyl_minus$`TF name`)
fisher.test(matrix(c(sig_minus,non_sig_minus,sig_non_minus,non_sig_non_minus),nrow=2)) #1.01, P=1

sig_minus=sum(motif_dir_dMML[FDR<=0.1 &Proportion_low_MML>0.5]$TF %in% methyl_minus$`TF name`)
non_sig_minus=sum(motif_dir_dMML[FDR>0.1|Proportion_low_MML<=0.5]$TF %in% methyl_minus$`TF name`)
sig_non_minus=sum(!motif_dir_dMML[FDR<=0.1&Proportion_low_MML>0.5]$TF %in% methyl_minus$`TF name`)
non_sig_non_minus=sum(!motif_dir_dMML[FDR>0.1|Proportion_low_MML<=0.5]$TF %in% methyl_minus$`TF name`)
fisher.test(matrix(c(sig_minus,non_sig_minus,sig_non_minus,non_sig_non_minus),nrow=2)) #0.753, P=0.3934

#ChromHMM motif

#find gnomic features enriched 
#OR: features vs non-feature, target vs non-targe, in dNME-ASM
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene <- readRDS(motif_gene_file)#See motif_break_array.R, default setting
genomic_features=readRDS(genomic_features_file)
selected_features=c("CpG island","CpG shore","CpG shelf","CpG open sea","gene body","exon","intron","intergenic","promoter","TSS")

motif_prefer_high_NME=fread('../downstream/output/graphs_tables/motif_preference_table/All_regions/motif_prefer_high_NME_only.csv')
motif_gene_high_NME=motif_gene[motif_gene$geneSymbol%in% motif_prefer_high_NME$TF]
variant_HetCpG_meta_dNME_ASM=variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff]
variant_HetCpG_meta_dNME_ASM$region=paste0(seqnames(variant_HetCpG_meta_dNME_ASM),":",start(variant_HetCpG_meta_dNME_ASM))
motif_gene_high_NME$region=paste0(seqnames(motif_gene_high_NME),':',start(motif_gene_high_NME))
# variant_HetCpG_meta_dNME_ASM_motif=variant_HetCpG_meta_dNME_ASM[variant_HetCpG_meta_dNME_ASM$region%in%motif_gene_high_NME_region]
# dNME_ASM_non_motif=variant_HetCpG_meta_dNME_ASM[!variant_HetCpG_meta_dNME_ASM$region%in%motif_gene_high_NME_region]

#Define colors 
feature_color=colorRampPalette(brewer.pal(9, "Set3"))(9)
#Extract motif binding regions
ah = AnnotationHub()
chromHMM_state=ENCODE_to_sample(unique(variant_HetCpG_meta$Sample))
chromHMM_list=list()
suppressMessages({for (sp in chromHMM_state$sample[!is.na(chromHMM_state$ENCODE)]){
  ah_num=names(AnnotationHub::query(ah, c("chromhmmSegmentations", chromHMM_state$ENCODE[chromHMM_state$sample==sp])))
  chromHMM_list[[sp]]=ah[[ah_num]]
  # }
}
})
cont_out_df=data.table()
chromHMM_motif_all_TF=data.table()
genomic_features_OR_out=data.table()
OR_cont_table=list()

for(TF in unique(motif_gene_high_NME$geneSymbol)){
  OR_cont_table[[TF]]=data.table()
  #ALT-REF
  motif_gene_high_NME_TF=motif_gene_high_NME[motif_gene_high_NME$geneSymbol==TF]
  variant_dNME_ASM_TF=variant_HetCpG_meta_dNME_ASM
  variant_dNME_ASM_TF$allele_diff=0
  olap=findOverlaps(variant_dNME_ASM_TF,motif_gene_high_NME_TF)
  variant_dNME_ASM_TF[queryHits(olap)]$allele_diff=motif_gene_high_NME_TF[subjectHits(olap)]$alleleDiff  
  
  variant_dNME_ASM_TF$allele_diff_NME=variant_dNME_ASM_TF$altNME-variant_dNME_ASM_TF$refNME
  #Fit the function
  variant_dNME_ASM_TF$ASM="No"
  variant_dNME_ASM_TF$ASM[sign(variant_dNME_ASM_TF$allele_diff_NME)==sign(variant_dNME_ASM_TF$allele_diff)]="Yes"
  genomic_features_OR=data.table()
  for(ft in selected_features){
    OR_feature=NULL
    OR_feature=testEnrichmentFeature_stat(variant_dNME_ASM_TF,genomic_features[[ft]],output_ft=1)
    if(!is.null(OR_feature)){
      OR_feature[[1]]$feature=ft
      OR_cont_table[[TF]]=rbind(OR_cont_table[[TF]],OR_feature[[1]])
      OR_feature=OR_feature[[2]]
      genomic_features_OR=rbind(genomic_features_OR,
                                data.table(feature=ft,OR_feature=OR_feature$estimate,p_value=OR_feature$p.value,
                                           lower_CI=OR_feature$conf.int[1],upper_CI=OR_feature$conf.int[2])
      )
    }
  }
  if(nrow(genomic_features_OR)>0){
    genomic_features_OR$TF=TF
    genomic_features_OR_out=rbind(genomic_features_OR_out,genomic_features_OR)
  }
  #Get samples for chromHMM analysis
  sample_chromHMM=names(table(variant_dNME_ASM_TF[variant_dNME_ASM_TF$ASM=="Yes"]$Sample))[table(variant_dNME_ASM_TF[variant_dNME_ASM_TF$ASM=="Yes"]$Sample)>=5]
  
  chromHMM_ls=list()
  for(sp in sample_chromHMM){
    chromHMM_in=chromHMM_list[[sp]]
    count_table=list()
    out_df=data.table()
    for(st in unique(chromHMM_in$name)){
      OR_chromHMM=NULL
      OR_chromHMM=testEnrichmentFeature_stat(variant_dNME_ASM_TF,chromHMM_in[chromHMM_in$name==st],output_ft=1)
      if(!is.null(OR_chromHMM)){
        count_table[[st]]=OR_chromHMM[[1]]
        OR_chromHMM=OR_chromHMM[[2]]
        out_df=rbind(out_df,
                     data.table(state=st,OR=OR_chromHMM$estimate,p_value=OR_chromHMM$p.value,
                                lower_CI=OR_chromHMM$conf.int[1],upper_CI=OR_chromHMM$conf.int[2]))
      }
      
      
    }
    out_df$Sample=sp
    
    chromHMM_ls[[sp]]=list(out_df,count_table)
  }
  if(length(chromHMM_ls)>0){chromHMM_motif_all=chromHMM_combine(chromHMM_ls)}
  if(nrow(chromHMM_motif_all)>0){
    chromHMM_motif_all$TF=TF
    out_df$TF=TF
    cont_out_df=rbind(cont_out_df,out_df)
    chromHMM_motif_all_TF=rbind(chromHMM_motif_all_TF,chromHMM_motif_all)
  }
  cat("Finish TF:",TF,'\n')
}

chromHMM_motif_all_TF$FDR=p.adjust(chromHMM_motif_all_TF$p_value,method="BH")
genomic_features_OR_out$FDR=p.adjust(genomic_features_OR_out$p_value,method='BH')
saveRDS(chromHMM_motif_all_TF,'../downstream/output/human_analysis/motif_analysis/chromHMM_motif_all_TF_NME.rds')
saveRDS(genomic_features_OR_out,'../downstream/output/human_analysis/motif_analysis/genomic_features_OR_motif_NME.rds')
saveRDS(OR_cont_table,'../downstream/output/human_analysis/motif_analysis/OR_cont_table_motif_NME.rds')
genomic_features_OR_out$feature =factor( genomic_features_OR_out$feature ,level=unique( genomic_features_OR_out$feature))
names(feature_color)=unique(genomic_features_OR_out$feature)
genomic_features_OR_out_sig=genomic_features_OR_out[FDR<=0.1]
pdf("../downstream/output/human_analysis/motif_analysis/OR_motif_preference_NME.pdf",width=4,height=4)
ggplot(genomic_features_OR_out_sig,aes(x=TF,y=OR_feature,group=feature,fill=feature))+
  geom_bar(stat ='identity',position=position_dodge())+ylab("OR")+
  scale_fill_manual(values=feature_color)+
  theme_glob+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom")
dev.off()
chromHMM_motif_all_TF_sig=chromHMM_motif_all_TF[FDR<=0.1&OR>=2]
chromHMM_motif_all_TF_sig_heatmap=dcast.data.table(chromHMM_motif_all_TF[TF%in%chromHMM_motif_all_TF_sig$TF],states~TF,value.var = "OR")
chromHMM_motif_all_TF_sig_FDR=as.matrix(dcast.data.table(chromHMM_motif_all_TF[TF%in%chromHMM_motif_all_TF_sig$TF],states~TF,value.var = "FDR")[,-1])
chromHMM_motif_all_TF_sig_FDR[chromHMM_motif_all_TF_sig_FDR <= 0.1]="*"
chromHMM_motif_all_TF_sig_FDR[chromHMM_motif_all_TF_sig_FDR > 0.1]=""
heatmap_mt=as.matrix(chromHMM_motif_all_TF_sig_heatmap[,-1])
rownames(heatmap_mt)=chromHMM_motif_all_TF_sig_heatmap$states
breaksList=seq(0,3,0.05)
heatmap_col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
pdf('../downstream/output/human_analysis/motif_analysis/chromHMM_SNP_feature_NME.pdf')
pheatmap(t(heatmap_mt),cluster_cols = F,cluster_rows = F,display_numbers=t(chromHMM_motif_all_TF_sig_FDR),
         color=heatmap_col,breaks = breaksList)
dev.off()

#MML

motif_prefer_low_MML=fread('../downstream/output/graphs_tables/motif_preference_table/All_regions/motif_prefer_low_MML_only.csv')
motif_gene_low_MML=motif_gene[motif_gene$geneSymbol%in% motif_prefer_low_MML$TF]
variant_HetCpG_meta_dMML_ASM=variant_HetCpG_meta[variant_HetCpG_meta$dMML_pval<=pval_cutoff]
variant_HetCpG_meta_dMML_ASM$region=paste0(seqnames(variant_HetCpG_meta_dMML_ASM),":",start(variant_HetCpG_meta_dMML_ASM))
motif_gene_low_MML$region=paste0(seqnames(motif_gene_low_MML),':',start(motif_gene_low_MML))
# variant_HetCpG_meta_dMML_ASM_motif=variant_HetCpG_meta_dMML_ASM[variant_HetCpG_meta_dMML_ASM$region%in%motif_gene_low_MML_region]
# dMML_ASM_non_motif=variant_HetCpG_meta_dMML_ASM[!variant_HetCpG_meta_dMML_ASM$region%in%motif_gene_low_MML_region]
#Extract motif binding regions
ah = AnnotationHub()
chromHMM_state=ENCODE_to_sample(unique(variant_HetCpG_meta$Sample))
chromHMM_list=list()
suppressMessages({for (sp in chromHMM_state$sample[!is.na(chromHMM_state$ENCODE)]){
  ah_num=names(AnnotationHub::query(ah, c("chromhmmSegmentations", chromHMM_state$ENCODE[chromHMM_state$sample==sp])))
  chromHMM_list[[sp]]=ah[[ah_num]]
  # }
}
})
cont_out_df=data.table()
chromHMM_motif_all_TF=data.table()
genomic_features_OR_out=data.table()
OR_cont_table=list()

for(TF in unique(motif_gene_low_MML$geneSymbol)){
  OR_cont_table[[TF]]=data.table()
  #ALT-REF
  motif_gene_low_MML_TF=motif_gene_low_MML[motif_gene_low_MML$geneSymbol==TF]
  variant_dMML_ASM_TF=variant_HetCpG_meta_dMML_ASM
  variant_dMML_ASM_TF$allele_diff=0
  olap=findOverlaps(variant_dMML_ASM_TF,motif_gene_low_MML_TF)
  variant_dMML_ASM_TF[queryHits(olap)]$allele_diff=motif_gene_low_MML_TF[subjectHits(olap)]$alleleDiff  
  
  variant_dMML_ASM_TF$allele_diff_MML=variant_dMML_ASM_TF$altMML-variant_dMML_ASM_TF$refMML
  #Fit the function
  variant_dMML_ASM_TF$ASM="No"
  variant_dMML_ASM_TF$ASM[sign(variant_dMML_ASM_TF$allele_diff_MML)==sign(variant_dMML_ASM_TF$allele_diff)]="Yes"
  genomic_features_OR=data.table()
  for(ft in selected_features){
    OR_feature=NULL
    OR_feature=testEnrichmentFeature_stat(variant_dMML_ASM_TF,genomic_features[[ft]],output_ft=1)
    if(!is.null(OR_feature)){
      OR_feature[[1]]$feature=ft
      OR_cont_table[[TF]]=rbind(OR_cont_table[[TF]],OR_feature[[1]])
      OR_feature=OR_feature[[2]]
      genomic_features_OR=rbind(genomic_features_OR,
                                data.table(feature=ft,OR_feature=OR_feature$estimate,p_value=OR_feature$p.value,
                                           lower_CI=OR_feature$conf.int[1],upper_CI=OR_feature$conf.int[2])
      )
    }
  }
  if(nrow(genomic_features_OR)>0){
    genomic_features_OR$TF=TF
    genomic_features_OR_out=rbind(genomic_features_OR_out,genomic_features_OR)
  }
  #Get samples for chromHMM analysis
  sample_chromHMM=names(table(variant_dMML_ASM_TF[variant_dMML_ASM_TF$ASM=="Yes"]$Sample))[table(variant_dMML_ASM_TF[variant_dMML_ASM_TF$ASM=="Yes"]$Sample)>=5]
  
  chromHMM_ls=list()
  for(sp in sample_chromHMM){
    chromHMM_in=chromHMM_list[[sp]]
    count_table=list()
    out_df=data.table()
    for(st in unique(chromHMM_in$name)){
      OR_chromHMM=NULL
      OR_chromHMM=testEnrichmentFeature_stat(variant_dMML_ASM_TF,chromHMM_in[chromHMM_in$name==st],output_ft=1)
      if(!is.null(OR_chromHMM)){
        count_table[[st]]=OR_chromHMM[[1]]
        OR_chromHMM=OR_chromHMM[[2]]
        out_df=rbind(out_df,
                     data.table(state=st,OR=OR_chromHMM$estimate,p_value=OR_chromHMM$p.value,
                                lower_CI=OR_chromHMM$conf.int[1],upper_CI=OR_chromHMM$conf.int[2]))
      }
      
      
    }
    out_df$Sample=sp
    
    chromHMM_ls[[sp]]=list(out_df,count_table)
  }
  if(length(chromHMM_ls)>0){chromHMM_motif_all=chromHMM_combine(chromHMM_ls)}
  if(nrow(chromHMM_motif_all)>0){
    chromHMM_motif_all$TF=TF
    out_df$TF=TF
    cont_out_df=rbind(cont_out_df,out_df)
    chromHMM_motif_all_TF=rbind(chromHMM_motif_all_TF,chromHMM_motif_all)
  }
  cat("Finish TF:",TF,'\n')
}

chromHMM_motif_all_TF$FDR=p.adjust(chromHMM_motif_all_TF$p_value,method="BH")
genomic_features_OR_out$FDR=p.adjust(genomic_features_OR_out$p_value,method='BH')
saveRDS(chromHMM_motif_all_TF,'../downstream/output/human_analysis/motif_analysis/chromHMM_motif_all_TF_MML.rds')
saveRDS(genomic_features_OR_out,'../downstream/output/human_analysis/motif_analysis/genomic_features_OR_motif_MML.rds')
saveRDS(OR_cont_table,'../downstream/output/human_analysis/motif_analysis/OR_cont_table_motif_MML.rds')
genomic_features_OR_out$feature =factor( genomic_features_OR_out$feature ,level=names(feature_color))
genomic_features_OR_out_sig=genomic_features_OR_out[FDR<=0.1]
pdf("../downstream/output/human_analysis/motif_analysis/OR_motif_preference_MML.pdf",width=5,height=4)
ggplot(genomic_features_OR_out_sig,aes(x=TF,y=OR_feature,group=feature,fill=feature))+
  geom_bar(stat ='identity',position=position_dodge())+ylab("OR")+
  scale_fill_manual(values=feature_color)+
  theme_glob+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position="bottom")
dev.off()
chromHMM_motif_all_TF_sig=chromHMM_motif_all_TF[FDR<=0.1&OR>=2]
chromHMM_motif_all_TF_sig_heatmap=dcast.data.table(chromHMM_motif_all_TF[TF%in%chromHMM_motif_all_TF_sig$TF],states~TF,value.var = "OR")
chromHMM_motif_all_TF_sig_FDR=as.matrix(dcast.data.table(chromHMM_motif_all_TF[TF%in%chromHMM_motif_all_TF_sig$TF],states~TF,value.var = "FDR")[,-1])
chromHMM_motif_all_TF_sig_FDR[chromHMM_motif_all_TF_sig_FDR <= 0.1]="*"
chromHMM_motif_all_TF_sig_FDR[chromHMM_motif_all_TF_sig_FDR > 0.1]=""
heatmap_mt=as.matrix(chromHMM_motif_all_TF_sig_heatmap[,-1])
rownames(heatmap_mt)=chromHMM_motif_all_TF_sig_heatmap$states
breaksList=seq(0,3,0.05)
heatmap_col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
pdf('../downstream/output/human_analysis/motif_analysis/chromHMM_SNP_feature_MML.pdf')
pheatmap(t(heatmap_mt),cluster_cols = F,cluster_rows = F,display_numbers=t(chromHMM_motif_all_TF_sig_FDR),
         color=heatmap_col,breaks = breaksList)
dev.off()