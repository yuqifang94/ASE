source('mainFunctions_sub.R')
# GWAS analysis -----------------------------------------------------------
#Functions
get_traits_GWAS<-function(variant_in,trait_gr_in,pval_cutoff=0.1,count_cutoff=3,stat='dNME_pval',CMH=FALSE,maxgap=500,ncores=15){
  traits_ls=mclapply(trait_gr_in,get_traits_GWAS_all_trait,
                   variant_in=variant_in,pval_cutoff=pval_cutoff,count_cutoff=count_cutoff,stat=stat,CMH=CMH,maxgap=maxgap,
                   mc.cores=ncores)
  return(traits_ls)
  #list(do.call(rbind,lapply(traits_ls,function(x)x[[1]])),do.call(rbind,lapply(traits_ls,function(x)x[[2]])))
}
get_traits_GWAS_all_trait<-function(trait_gr,variant_in,pval_cutoff,count_cutoff,stat,CMH,maxgap=500){
  #trait_gr=trait_gr[trait_gr$trait%in%trait]
  trait=unique(trait_gr$`DISEASE/TRAIT`)
  OR_output=data.frame()
  CMH_df=data.frame()
  dNME_sig=variant_in[elementMetadata(variant_in)[,stat]<=pval_cutoff]
  dNME_non_sig=variant_in[elementMetadata(variant_in)[,stat]>pval_cutoff]
  dNME_traits_gr=findOverlaps(dNME_sig,trait_gr,maxgap = maxgap)
  # dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
  # non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
  # dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
  # non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
  dNME_trait=length(unique(queryHits(dNME_traits_gr)))
  dNME_non_trait=length(dNME_sig)-dNME_trait
  non_dNME_trait=length(subsetByOverlaps(dNME_non_sig,trait_gr,maxgap = maxgap))
  non_dNME_non_trait=length(dNME_non_sig)-non_dNME_trait
  trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
  journal=paste(unique(trait_gr$JOURNAL),collapse = ',')
  trait_gr_out=trait_gr[subjectHits(dNME_traits_gr)]
  mcols(trait_gr_out)=mcols(trait_gr_out)[c('MAPPED_GENE','STRONGEST SNP-RISK ALLELE','DISEASE/TRAIT')]
  names(mcols(trait_gr_out))=c('genes','risk allele','traits')
  rm(dNME_sig)
  rm(dNME_non_sig)
  rm(trait_gr)
  gc()
  if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
    
    CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
                      count=trait_count,subject='All')
    cat('processing trait',trait,'\n')
    
    if(length(CMH_df)>0){
      if (CMH){return(CMH_df)}
      else{
        
        OR=CMH_test(CMH_df)
        OR_output=cbind(t(CMH_df$count),data.frame(trait=trait,OR=OR$estimate,
                                                   p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
        colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
        OR_output$journal=journal
        
        return(list(OR_output,trait_gr_out))
      }
    }
  }
}
#Modify this to adapt trait analysis
get_traits_GWAS_all_trait_single<-function(variant_in,trait_gr,pval_cutoff,count_cutoff,stat,CMH,maxgap){
 
  OR_output=data.frame()
  CMH_df=data.frame()
  dNME_sig=variant_in[elementMetadata(variant_in)[,stat]<=pval_cutoff]
  dNME_non_sig=variant_in[elementMetadata(variant_in)[,stat]>pval_cutoff]
  
  # dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
  # non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
  # dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
  # non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
  dNME_trait=length(subsetByOverlaps(dNME_sig,trait_gr,maxgap = maxgap))
  dNME_non_trait=length(dNME_sig)-dNME_trait
  non_dNME_trait=length(subsetByOverlaps(dNME_non_sig,trait_gr,maxgap = maxgap))
  non_dNME_non_trait=length(dNME_non_sig)-non_dNME_trait
  trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
  rm(dNME_sig)
  rm(dNME_non_sig)
  rm(trait_gr)
  gc()
  if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
    
    CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
                      count=trait_count,subject='All')

    
    if(length(CMH_df)>0){
      if (CMH){return(CMH_df)}
      else{
        
        OR=CMH_test(CMH_df)
        OR_output=cbind(t(CMH_df$count),data.frame(OR=OR$estimate,
                                                   p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
        colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
        OR_output$journal=paste(unique(trait_gr$JOURNAL),collapse = ',')
        return(OR_output)
      }
    }
  }
}



get_traits_GWAS_trait<-function(trait,variant_in,pval_cutoff,count_cutoff,stat,CMH){
  traits_sp_ls=lapply(unique(variant_in$germlayer),get_traits_GWAS_sp_trait,trait=trait,
                      variant_in=variant_in,pval_cutoff=pval_cutoff,count_cutoff=count_cutoff,stat=stat,CMH=CMH)
  traits_sp=do.call(rbind,traits_sp_ls)
  if(CMH&length(traits_sp)>0){
    OR=CMH_test(traits_sp)
    OR_output=data.frame(trait=trait,OR=OR$estimate,p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2])
    return(OR_output)
  }else if(!CMH){
    return(traits_sp)
  }


}




get_traits_GWAS_sp_trait<-function(sp,trait,variant_in,pval_cutoff,count_cutoff,stat,CMH){
 
      OR_output=data.frame()
      CMH_df=data.frame()
      variant_in_sp=variant_in[variant_in$germlayer==sp&!is.na(variant_in$trait),]
 
      dNME_sig=variant_in_sp[,stat]<=pval_cutoff
      dNME_non_sig=variant_in_sp[,stat]>pval_cutoff
    
      dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
      non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
      dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
      non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
      trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
      if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
       
        CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
                                       count=trait_count,subject=sp)
        cat('processing',sp,'with trait',trait,'\n')

      if(length(CMH_df)>0){
        if (CMH){return(CMH_df)}
        else{
          OR=CMH_test(CMH_df)
          OR_output=cbind(t(CMH_df$count),data.frame(subject=sp,trait=trait,OR=OR$estimate,
                                                     p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
          colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
          return(OR_output)
        }
      }
     }
}

tissue_to_germlayer<-function(GR_input){
  GR_input$germlayer=NA
  tissue_ectoderm=c("foreskin_keratinocyte_paired",
                    "foreskin_melanocyte_paired",
                    "ectoderm_paired",
                    "brain_cerebellum_tissue_paired",
                    "brain_germinal_matrix_tissue_paired",
                    "Brain_substantia_nigra_paired",
                    "Brain_Hippocampus_middle_paired" )
  tissue_mesoderm=c("mesoderm_23_paired","Adipose_single",
                    "Left_Ventricle_single","Psoas_Muscle_single" ,
                    "Right_Ventricle_single","Right_Atrium_single","Spleen_single",
                    "Adrenal_Gland_single","Aorta_single","Ovary_single")
  tissue_endoderm=c("Small_Intestine_single","Lung_single","endoerm_27_paired",
                    "Bladder_single" ,"Gastric_single", "Sigmoid_Colon_single",
                    "Thymus_single","Esophagus_single", "Pancreas_single" ,"Liver_single")
  tissue_ESC=c("rep1","rep2","merged","42_embryonic_stem_cell_single" , "stem_27_undifferentiated_paired",'ESC')
  GR_input$germlayer[GR_input$tissue %in% tissue_ectoderm]='ectoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_mesoderm]='mesoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_endoderm]='endoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_ESC]='ESC'
  return(GR_input)
}

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

#Examine neurological traits
#Start with line 74
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
variant_HetCpG_meta_brain=variant_HetCpG_meta[grepl("brain|Brain",variant_HetCpG_meta$Sample)]

run_CMH_brain<-function(variant_HetCpG_meta_brain,SNP_trait,maxGap){
  SNP_trait$end=SNP_trait$start
  SNP_trait_gr=makeGRangesFromDataFrame(SNP_trait)
  seqlevels(SNP_trait_gr)=paste0("chr",seqlevels(SNP_trait_gr))

  cmh_out=c()
  for(sp in unique(variant_HetCpG_meta_brain$Sample)){
    variant_HetCpG_meta_sp=variant_HetCpG_meta_brain[variant_HetCpG_meta_brain$Sample==sp]
    dNME_trait_olap=findOverlaps(variant_HetCpG_meta_sp,SNP_trait_gr,maxgap=maxGap)
    variant_trait=variant_HetCpG_meta_sp[queryHits(dNME_trait_olap)]
    variant_non_trait=variant_HetCpG_meta_sp[-queryHits(dNME_trait_olap)]
    dNME_trait=sum(variant_trait$dNME_pval<=0.1)
    non_dNME_trait=sum(variant_trait$dNME_pval>0.1)
    dNME_non_trait=sum(variant_non_trait$dNME_pval<=0.1)
    non_dNME_non_trait=sum(variant_non_trait$dNME_pval>0.1)
    cmh_out=c(cmh_out,c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait))
  }
  cmh_out_run=array(cmh_out,dim=c(2,2,length(unique(variant_HetCpG_meta_brain$Sample))),
                  dimnames=list(dNME=c("dNME","non-dNME"),
                                trait=c("trait","non-trait"),
                                sample=unique(variant_HetCpG_meta_brain$Sample)))
  return(mantelhaen.test(cmh_out_run,exact=TRUE)) 
}                          

SNP_sch=fread('../downstream/input/human_analysis/PGC/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv',skip=73,header=T)
SNP_sch_trait=SNP_sch[PVAL<=5*10^-8]#22345
colnames(SNP_sch_trait)[c(1,3)]=c("chr","start")
sch=run_CMH_brain(variant_HetCpG_meta_brain,SNP_sch_trait,maxGap=1000)#2.946185

SNP_bp=fread('../downstream/input/human_analysis/PGC/pgc-bip2021-all.vcf.tsv',skip=72,header=T)
SNP_bp_trait=SNP_bp[PVAL<=5*10^-8]#22345
colnames(SNP_bp_trait)[c(1,2)]=c("chr","start")
bp=run_CMH_brain(variant_HetCpG_meta_brain,SNP_bp_trait,maxGap=1000)#3.215223

SNP_az=fread('../downstream/input/human_analysis/PGC/PGCALZ2sumstatsExcluding23andMe.txt',skip=0,header=T)
SNP_az_trait=SNP_az[p<=5*10^-8]#22345
colnames(SNP_az_trait)[c(1,2)]=c("chr","start")
az=run_CMH_brain(variant_HetCpG_meta_brain,SNP_az_trait,maxGap=1000)#0.22

SNP_adhd=fread('../downstream/input/human_analysis/PGC/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta',skip=0,header=T)
SNP_adhd_trait=SNP_adhd[P<=5*10^-8]#22345
colnames(SNP_adhd_trait)[c(1,3)]=c("chr","start")
adhd=run_CMH_brain(variant_HetCpG_meta_brain,SNP_adhd_trait,maxGap=1000)#0.8197396