
# TBD  --------------------------------------------------------------------

#find gnomic features enriched 
#OR: features vs non-feature, target vs non-targe, in dNME-ASM
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene <- readRDS(motif_gene_file)#See motif_break_array.R, default setting
genomic_features=readRDS(genomic_features_file)
selected_features=c("CpG island","CpG shore","CpG shelf","CpG open sea","gene body","exon","intron","intergenic","promoter","TSS")
motif_prefer_high_NME=fread('../downstream/output/graphs/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
motif_gene_high_NME=motif_gene[motif_gene$geneSymbol%in% motif_prefer_high_NME$TF]
variant_HetCpG_meta_dNME_ASM=variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff]
variant_HetCpG_meta_dNME_ASM$region=paste0(seqnames(variant_HetCpG_meta_dNME_ASM),":",start(variant_HetCpG_meta_dNME_ASM))
motif_gene_high_NME$region=paste0(seqnames(motif_gene_high_NME),':',start(motif_gene_high_NME))
# variant_HetCpG_meta_dNME_ASM_motif=variant_HetCpG_meta_dNME_ASM[variant_HetCpG_meta_dNME_ASM$region%in%motif_gene_high_NME_region]
# dNME_ASM_non_motif=variant_HetCpG_meta_dNME_ASM[!variant_HetCpG_meta_dNME_ASM$region%in%motif_gene_high_NME_region]
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
saveRDS(chromHMM_motif_all_TF,'../downstream/output/chromHMM_motif_all_TF.rds')
saveRDS(genomic_features_OR_out,'../downstream/output/genomic_features_OR_motif.rds')
saveRDS(OR_cont_table,'../downstream/output/OR_cont_table_motif.rds')
genomic_features_OR_out_sig=genomic_features_OR_out[FDR<=0.1]
pdf("../downstream/output/OR_motif_preference.pdf",width=9,height=4)
ggplot(genomic_features_OR_out_sig,aes(x=TF,y=OR_feature,group=feature,fill=feature))+
  geom_bar(stat ='identity',position=position_dodge())+ylab("OR")+
  theme_glob+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
pdf('../downstream/output/chromHMM_SNP_feature.pdf')
pheatmap(t(heatmap_mt),cluster_cols = F,cluster_rows = F,display_numbers=t(chromHMM_motif_all_TF_sig_FDR),
         color=heatmap_col,breaks = breaksList)
dev.off()
#Find CTCF example

TF="CTCF"
motif_gene_high_NME_TF=motif_gene_high_NME[motif_gene_high_NME$geneSymbol==TF]
variant_dNME_ASM_TF=variant_HetCpG_meta_dNME_ASM
variant_dNME_ASM_TF$allele_diff=0
olap=findOverlaps(variant_dNME_ASM_TF,motif_gene_high_NME_TF)
variant_dNME_ASM_TF[queryHits(olap)]$allele_diff=motif_gene_high_NME_TF[subjectHits(olap)]$alleleDiff  
variant_dNME_ASM_TF$pctRef=NA
variant_dNME_ASM_TF[queryHits(olap)]$pctRef=motif_gene_high_NME_TF[subjectHits(olap)]$pctRef  
variant_dNME_ASM_TF$pctAlt=NA
variant_dNME_ASM_TF[queryHits(olap)]$pctAlt=motif_gene_high_NME_TF[subjectHits(olap)]$pctAlt    
variant_dNME_ASM_TF$allele_diff_NME=variant_dNME_ASM_TF$altNME-variant_dNME_ASM_TF$refNME
#Fit the function
variant_dNME_ASM_TF$ASM="No"
variant_dNME_ASM_TF$ASM[sign(variant_dNME_ASM_TF$allele_diff_NME)==sign(variant_dNME_ASM_TF$allele_diff)]="Yes"
cTCF_loc=variant_dNME_ASM_TF[variant_dNME_ASM_TF$ASM=="Yes"]
cTCF_loc=cTCF_loc[cTCF_loc$N>=4]
cTCF_loc=cTCF_loc[order(cTCF_loc$dNME,decreasing=T)]

source('plotMB.R')
pdf('../downstream/output/mouse_examples/motif_CTCF.pdf')
plotMB(subsetByOverlaps(motif_gene_high_NME_TF,cTCF_loc[2]),'STL003-874104')
dev.off()

#Mouse motif
# Motif analysis example in mouse -----------------------------------------
#read in significant motif
plot_motif_tissue<-function(tissue){
  FDR_cutoff=0.05
  motif_tissue=data.table()
  #target_regions_cluster_all=GRanges()
  target_regions=readRDS(paste0('../downstream/input/motif_target/',tissue,'_TF_motif_site.rds'))
  chromHMM_in=readRDS('../downstream/output/chromHMM_enhancer.rds')
  chromHMM_ts=chromHMM_in[chromHMM_in$tissue==tissue]
  region_in=fread(paste0('../downstream/input/mm10_cluster_all/',tissue,'.csv'))
  region_in=region_in[chromHMM_enhancer==TRUE]
  region_in_all=data.table()
  for(i in 1:10){
    motif_in=fread(paste0('../downstream/input/mouse_motif_cluster/',tissue,'/motif_',tissue,'_cluster_',i,'_enhancer.csv'))
    motif_in$cluster=i
    motif_in=motif_in[FDR<=FDR_cutoff]
    
    if(length(motif_in$motif)>0){
      #read in GO result
      GO_in=fread(paste0('../downstream/output/mm10_result/chromHMM_enhancer/cluster_GO/mm10_cluster_chromHMM/',tissue,'-',i,'_cluster_GO.csv'))
      GO_in=GO_in[FC>=1.5][1:5]
      GO_in_gene=unique(unlist(strsplit(GO_in$gene,";")))
      #read in all regions
      region_in_clu=region_in[cluster==i]
      target_regions_cluster=fastDoCall('c',lapply(motif_in$motif,function(x){
        target_regions_out=unique(target_regions[[x]][target_regions[[x]]$cluster==i])
        target_regions_out$motif=x
        return(target_regions_out)  }))
      target_regions_cluster=target_regions_cluster[target_regions_cluster$gene %in% GO_in_gene]
      olap=findOverlaps(convert_GR(region_in_clu$region),target_regions_cluster)
      region_in_clu=region_in_clu[queryHits(olap)]
      region_in_clu$gene_Jason=region_in_clu$gene
      region_in_clu$gene=NULL
      region_in_clu=cbind(region_in_clu,as.data.table(mcols(target_regions_cluster[subjectHits(olap)])))
      region_in_all=rbind(region_in_all,region_in_clu)
    }
    motif_tissue=rbind(motif_tissue,motif_in)
    
  }
  
  UC_raw=readRDS('../downstream/output/uc_matrix_DNase.rds')
  mml <- readRDS('../downstream/output/mml_matrix_DNase.rds')
  nme <- readRDS('../downstream/output/nme_matrix_DNase.rds')
  timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
  UC=UC_raw[[tissue]][,colnames(UC_raw[[tissue]])%in% timeorder]
  UC <-UC[unique(region_in_all$region),order(match(colnames(UC),timeorder))]
  stat_differential<-function(stat_in,UC,tissue){
    
    diff_out <- sapply(colnames(UC),function(i) {
      time <- strsplit(i,'-')
      sapply(time,function(time_in)
      {abs(stat_in[rownames(UC),paste0(tissue,'-',time_in[1],'-all')]-stat_in[rownames(UC),paste0(tissue,'-',time_in[2],'-all')])})
    })
    colnames(diff_out)=colnames(UC)
    rownames(diff_out)=rownames(UC)
    return(diff_out)
  }
  dmml=stat_differential(mml,UC,tissue)
  dnme=stat_differential(nme,UC,tissue)
  
  region_in_all_region=region_in_all[order(dNME_maxJSD,decreasing=T),list(cluster=paste(unique(cluster),collapse = ';'),
                                                                          motif=paste(unique(motif),collapse = ';'),
                                                                          dNME_maxJSD=dNME_maxJSD),
                                     by=list(region,gene,distance)]
  theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                   axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                   axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                   axis.text.x=element_text(size=16),
                   axis.text.y=element_text(size=16))+theme_classic()
  dmml_out=dmml[region_in_all$region,]
  colnames(dmml_out)=paste0('dMML-',colnames(dmml_out))
  dnme_out=dnme[region_in_all$region,]
  colnames(dnme_out)=paste0('dNME-',colnames(dnme_out))
  region_in_all=cbind(region_in_all,as.data.table(dmml_out),as.data.table(dnme_out))
  pdf(paste0('../downstream/output/graphs/Figure7/',tissue,'_motif.pdf'),width=5,height=5)
  for(rg in region_in_all_region$region){
    region_stat=data.table(stage=colnames(UC),
                           UC=UC[rg,],
                           dmml=dmml[rg,],
                           dnme=dnme[rg,]
    )
    region_stat=melt.data.table(region_stat,id.var='stage',variable.name = "stat")
    
    print(ggplot(region_stat,aes(x=stage,y=value,group=stat,color=stat))+geom_point(size=1)+geom_line(stat = "identity")+
            theme_glob+theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom",plot.title = element_text(hjust = 0.5))+
            scale_color_manual(values=c(UC="black", dmml="blue", dnme="red"))+
            ggtitle(paste(region_in_all_region[region==rg,],collapse='\n')))
    
    
  }
  dev.off()
  return(region_in_all)
}
heart_motif=plot_motif_tissue("heart")
limb_motif=plot_motif_tissue("limb")
forebrain_motif=plot_motif_tissue("forebrain")
#Find number of CG
CpG_mm10=getCpgSitesmm10()
heart_motif$N=countOverlaps(convert_GR(heart_motif$region),CpG_mm10)
heart_motif[order(N,decreasing=T)]


# Currently Not in use ----------------------------------------------------