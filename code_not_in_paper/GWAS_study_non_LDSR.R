rm(list=ls())
source('mainFunctions_sub.R')
traits_hg19 =makeCurrentGwascat(genome='GRCh37')
seqlevels(traits_hg19)=paste0('chr',seqlevels(traits_hg19))
saveRDS(traits_hg19,'../downstream/input/human_analysis/variant_traits.rds')

#Check GWAS traits have hg19 coordinates
variant_trait=readRDS('../downstream/input/human_analysis/variant_traits.rds')
# variant_trait_gr=do.call(c,variant_trait)
# variant_trait_gr=unique(variant_trait_gr)
# genome_freq=readRDS('../downstream/input/genome_1k_variant.rds')
# genome_freq$variant_id=gsub('dbSNP_153\\:','',genome_freq$Dbxref)
# olap=findOverlaps(variant_trait,genome_freq,type='equal')
#Check proportion of dNME SNP from GWAS
#get traits location

variant_trait_gr=do.call(c,variant_trait)

variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
#variant_HetCpG_meta=readRDS(GR_merge_file)
elementMetadata(variant_HetCpG_meta)=elementMetadata(variant_HetCpG_meta)[,c('dNME','dMML','dNME_pval','dMML_pval','tissue')]
length(subsetByOverlaps(variant_HetCpG_meta[variant_HetCpG_meta$dNME_pval<=pval_cutoff],variant_trait_gr,maxgap = 0))
olap=findOverlaps(variant_HetCpG_meta,variant_trait_gr)
olap_dt=data.table(qt=queryHits(olap),traits=variant_trait_gr[subjectHits(olap)]$`DISEASE/TRAIT`)
olap_dt=olap_dt[,list(traits=paste(traits,collapse = ';')),by=qt]
variant_HetCpG_meta$GWAS=NA
variant_HetCpG_meta$GWAS[olap_dt$qt]=olap_dt$traits
GWAS_enrich<-function(var_in,stat_in){
  diff_GWAS=sum(mcols(var_in)[[stat_in]]<=0.1&(!is.na(var_in$GWAS)))
  diff_non_GWAS=sum(mcols(var_in)[[stat_in]]<=0.1&(is.na(var_in$GWAS)))
  non_diff_GWAS=sum(mcols(var_in)[[stat_in]]>0.1&(!is.na(var_in$GWAS)))
  non_diff_non_GWAS=sum(mcols(var_in)[[stat_in]]>0.1&(is.na(var_in$GWAS)))
  return(fisher.test(matrix(c(diff_GWAS,diff_non_GWAS,non_diff_GWAS,non_diff_non_GWAS),nrow=2)))
  
}
#Not enriched in GWAS
GWAS_enrich(variant_HetCpG_meta,"dMML_pval")
GWAS_enrich(variant_HetCpG_meta,"dNME_pval")
#variant_HetCpG_meta=convert_GR(variant_HetCpG_meta[!is.na(variant_HetCpG_meta$GWAS)],direction="DT")
variant_HetCpG_meta$region=gsub('-.*','',variant_HetCpG_meta$region)#It's a SNP
variant_HetCpG_meta$GWAS=gsub(',',' ',variant_HetCpG_meta$GWAS)
write.csv(variant_HetCpG_meta,'../downstream/output/human_analysis/GWAS_traits_dNME.csv',row.names = F,quote=F)



# Use allele-agnostic data to get GWAS percent ----------------------------
NME_in=readRDS(NME_agnostic_DNase_file)
DNase_hg19=readRDS(DNase_hg19_file)
control_hg19=readRDS(control_hg19_file)
 variant_trait=readRDS('../downstream/input/human_analysis/variant_traits.rds')
 NME_in$DNase=NA
 
 NME_in_olap=findOverlaps(NME_in,DNase_hg19,type='equal')
 
 NME_in$DNase[queryHits(NME_in_olap)]="DNase"
 NME_in_olap=findOverlaps(NME_in,control_hg19,type='equal')
 NME_in$DNase[queryHits(NME_in_olap)]="control"

# Looking for enrichment for tissue-specific high NME regions -------------
 NME_in_DNase=NME_in[NME_in$DNase=="DNase"]
 high_NME=quantile(NME_in_DNase$NME,prob=0.9)
 NME_in_DNase_dt=convert_GR(NME_in_DNase,direction='DT')
 NME_in_DNase_dt$tissue=gsub(' - .*','',NME_in_DNase_dt$Sample)
 NME_in_DNase_dt$tissue=gsub('Adipose_Tissue_single','Adipose_single',NME_in_DNase_dt$tissue)
 NME_in_DNase_dt=tissue_to_germlayer(NME_in_DNase_dt)
 NME_in_DNase_dt_max_NME=NME_in_DNase_dt[,list(NME_max=max(NME)),by=list(region,tissue)]
 #NME_in_DNase_dt_max_NME=NME_in_DNase_dt[,list(NME_max=max(NME)),by=list(region,germlayer)]
 NME_in_DNase_dt_max_NME_unique=NME_in_DNase_dt_max_NME[,list(num_highNME=sum(NME_max>=high_NME),tissue_max=paste(tissue[NME_max>=high_NME],collapse=';')),by=list(region)]
 #NME_in_DNase_dt_max_NME_unique=NME_in_DNase_dt_max_NME[,list(num_highNME=sum(NME_max>=high_NME),tissue_max=paste(germlayer[NME_max>=high_NME],collapse=';')),by=list(region)]
 sum(NME_in_DNase_dt_max_NME_unique$num_highNME==1)/nrow(NME_in_DNase_dt_max_NME_unique)
 NME_in_DNase_dt=cbind(NME_in_DNase_dt,NME_in_DNase_dt_max_NME_unique[match(NME_in_DNase_dt$region,region),list(num_highNME,tissue_max)])
 NME_in_DNase_dt_max_NME_unique_gr=convert_GR(NME_in_DNase_dt_max_NME_unique$region)
 mcols(NME_in_DNase_dt_max_NME_unique_gr)=NME_in_DNase_dt_max_NME_unique
 variant_trait=GRanges(variant_trait)
 trait_NME_out_ts=mclapply(unique(variant_trait$`DISEASE/TRAIT`),function(trait){
   variant_trait_in=variant_trait[variant_trait$`DISEASE/TRAIT`==trait]
   #Overlap with high NME
   unique_high_NME_trait_gr=subsetByOverlaps(NME_in_DNase_dt_max_NME_unique_gr[NME_in_DNase_dt_max_NME_unique_gr$num_highNME==1],variant_trait_in)
   unique_high_NME_trait=length(unique_high_NME_trait_gr)
   unique_high_NME_non_trait=sum(NME_in_DNase_dt_max_NME_unique_gr$num_highNME==1)-unique_high_NME_trait
   all_low_NME_trait=length(subsetByOverlaps(NME_in_DNase_dt_max_NME_unique_gr[NME_in_DNase_dt_max_NME_unique_gr$num_highNME==0],variant_trait_in))
   all_low_NME_non_trait=sum(NME_in_DNase_dt_max_NME_unique_gr$num_highNME==0)-all_low_NME_trait
   non_unique_high_NME_trait=length(subsetByOverlaps(NME_in_DNase_dt_max_NME_unique_gr[NME_in_DNase_dt_max_NME_unique_gr$num_highNME>1],variant_trait_in))
   non_unique_high_NME_non_trait=sum(NME_in_DNase_dt_max_NME_unique_gr$num_highNME>1)-non_unique_high_NME_trait
   
   OR_lowNME=fisher.test(matrix(c(unique_high_NME_trait,unique_high_NME_non_trait,
                                  all_low_NME_trait,all_low_NME_non_trait),ncol=2))
   
   OR_highNME=fisher.test(matrix(c(unique_high_NME_trait,unique_high_NME_non_trait,
                                   non_unique_high_NME_trait,non_unique_high_NME_non_trait),ncol=2))
   return(data.table(unique_high_NME_trait=unique_high_NME_trait,
                     unique_high_NME_non_trait=unique_high_NME_non_trait,
                     all_low_NME_trait=all_low_NME_trait,
                     all_low_NME_non_trait=all_low_NME_non_trait,
                     non_unique_high_NME_trait=non_unique_high_NME_trait,
                     non_unique_high_NME_non_trait=non_unique_high_NME_non_trait,
                     OR_lowNME=OR_lowNME$estimate,
                     OR_lowNME_pval=OR_lowNME$p.value,
                     OR_lowNME_lowCI=OR_lowNME$conf.int[1],
                     OR_lowNME_highCI=OR_lowNME$conf.int[2],
                     OR_highNME=OR_highNME$estimate,
                     OR_highNME_pval=OR_highNME$p.value,
                     OR_highNME_lowCI=OR_highNME$conf.int[1],
                     OR_highNME_highCI=OR_highNME$conf.int[2],
                     trait=trait,
                     tissue_max=paste0(unique_high_NME_trait_gr$tissue_max,collapse=';')
   ))
  
   
},mc.cores=20)
 saveRDS(trait_NME_out_ts,'../downstream/output/human_analysis/GWAS/trait_NME_out_tissue.rds')
 trait_NME_out_ts=readRDS('../downstream/output/human_analysis/GWAS/trait_NME_out_tissue.rds')
 trait_NME_out_ts=do.call(rbind,trait_NME_out_ts)
 write.csv(trait_NME_out_ts[!is.infinite(OR_lowNME )&!is.infinite(OR_highNME)],'../downstream/output/human_analysis/GWAS/trait_enrich_tissue.csv')
 trait_NME_out_gm=readRDS('../downstream/output/human_analysis/GWAS/trait_NME_out_germlayer.rds')
 trait_NME_out_gm=do.call(rbind,trait_NME_out_gm)
 write.csv(trait_NME_out_gm[!is.infinite(OR_lowNME )&!is.infinite(OR_highNME)],'../downstream/output/human_analysis/GWAS/trait_enrich_germlayer.csv')
# SNP analysis ------------------------------------------------------------

 library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
 snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
NME_in_gr=unique(granges(NME_in))
 seqlevels(NME_in_gr)=gsub('chr','',seqlevels(NME_in_gr))
 SNP_NME_in=snpsByOverlaps(snps,NME_in_gr)
seqlevels(SNP_NME_in)=paste0('chr',seqlevels(SNP_NME_in))
SNP_NME_in=SNP_NME_in[seqnames(SNP_NME_in)!="chrMT"]
genome(SNP_NME_in)=gsub('.p13','',genome(SNP_NME_in))
SNP_NME_in=GRanges(SNP_NME_in)
variant_trait=GRanges(variant_trait)
SNP_NME_in$GWAS="Non GWAS"
variant_trait$cancer="Non-cancer"
variant_trait$cancer[grepl('lympho|cancer|leukemia|melanoma|myeloma|carcinoma|Meningioma|Adverse response to chemotherapy|Cutaneous nevi|Glioma|tumor',variant_trait$`DISEASE/TRAIT`,ignore.case = T)]=
  "cancer-related trait"
olap_non_cancer=findOverlaps(SNP_NME_in,variant_trait[variant_trait$cancer=="Non-cancer"])
SNP_NME_in$GWAS[queryHits(olap_non_cancer)]="Non-cancer"
olap_cancer=findOverlaps(SNP_NME_in,variant_trait[variant_trait$cancer=="cancer-related trait"])
SNP_NME_in$GWAS[queryHits(olap_cancer)]="cancer-related trait"

# variant_trait=do.call('c',variant_trait)
# 
# NME_ecdf=ecdf(NME_in$NME)
# high_NME=quantile(NME_in$NME,prob=0.75)


#Compare GWAS traits vs non GWAS SNP
NME_in_DNase=NME_in[NME_in$DNase=='DNase']
GWAS_NME_olap_cancer=findOverlaps(NME_in_DNase,SNP_NME_in[SNP_NME_in$GWAS== "cancer-related trait"],maxgap = 2000)
non_GWAS_NME_olap=findOverlaps(NME_in_DNase,granges(SNP_NME_in[SNP_NME_in$GWAS=="Non GWAS"]),maxgap = 2000)
GWAS_NME_olap_non_cancer=findOverlaps(NME_in_DNase,SNP_NME_in[SNP_NME_in$GWAS== "Non-cancer"],maxgap = 2000)

NME_in_DNase$GWAS_SNPs='Not SNP overlapping'

NME_in_DNase$GWAS_SNPs[queryHits(non_GWAS_NME_olap)]="non GWAS SNP"
NME_in_DNase$GWAS_SNPs[queryHits(GWAS_NME_olap_non_cancer)]="GWAS SNP non cancer"
NME_in_DNase$GWAS_SNPs[queryHits(GWAS_NME_olap_cancer)]="GWAS SNP cancer"

pdf('../downstream/output/human_analysis/GWAS_NME_density_2k.pdf')
ggplot(data=as.data.frame(mcols(NME_in_DNase)),aes(x=NME,group=GWAS_SNPs,color=GWAS_SNPs))+geom_density()+theme(legend.position = 'bottom')
dev.off()
table(NME_in_DNase$GWAS_SNPs)/length(NME_in_DNase)
as.data.table(mcols(NME_in_DNase))[,list(mean_NME=mean(NME),sd_NME=sd(NME)),by=list(GWAS_SNPs)]


NME_in_dt=convert_GR(NME_in,direction="DT")
#Do it for each type of germlayer
NME_in_dt$tissue=gsub(' - .*','',NME_in_dt$Sample)
NME_in_dt$tissue=gsub('Adipose_Tissue_single','Adipose_single',NME_in_dt$tissue)
NME_in_dt=tissue_to_germlayer(NME_in_dt)
trait_NME_out_ts=list()
#   
trait_NME_out_ts=lapply(unique(NME_in_dt$germlayer),function(ts){
  cat("Processing:",ts,'\n')
  trait_NME_out=data.table()
  NME_in_dt_ts=NME_in_dt[germlayer==ts]
  #NME_in_dt_ts_mean=NME_in_dt_ts[,list(NME_mean=mean(NME)),by=list(region,DNase)]
  NME_in_dt_ts_mean=NME_in_dt_ts
  NME_in_dt_ts_mean$NME_mean=NME_in_dt_ts_mean$NME
  #NME_ecdf=ecdf(NME_in_dt_ts_mean$NME_mean)
  high_NME=quantile(NME_in_dt_ts_mean$NME[NME_in_dt_ts_mean$DNase=="DNase"],prob=0.75)
  NME_in_gr_mean_ts=convert_GR(NME_in_dt_ts_mean$region)
  mcols(NME_in_gr_mean_ts)=NME_in_dt_ts_mean
  variant_trait=GRanges(variant_trait)
  print(length(unique(variant_trait$`DISEASE/TRAIT`)))
  trait_NME_out=do.call(rbind,mclapply(unique(variant_trait$`DISEASE/TRAIT`),trait_overlap,
                                       NME_in_gr_mean_ts=NME_in_gr_mean_ts,
                                       variant_trait=variant_trait,high_NME=high_NME,
                                       mc.cores=20))
  trait_NME_out$germlayer=ts
  print(trait_NME_out)
  return(trait_NME_out)
})
saveRDS(trait_NME_out_ts,'../downstream/output/human_analysis/trait_NME_out_ts_germlayer.rds')
trait_NME_out_ts_germlayer=readRDS('../downstream/output/human_analysis/GWAS/trait_NME_out_ts_germlayer.rds')
trait_NME_out_ts_germlayer=do.call(rbind,trait_NME_out_ts_germlayer)
# trait_NME_out_ts_merged=do.call(rbind,trait_NME_out_ts)
# trait_NME_out_ts_merged=trait_NME_out_ts_merged[DNase_total>0]
# p_exp=sum(trait_NME_out_ts_merged$DNase_high)/sum(trait_NME_out_ts_merged$DNase_total)
# trait_NME_out_ts_merged$p_value_binom=unlist(lapply(1:nrow(trait_NME_out_ts_merged),function(x) 
#   binom.test(trait_NME_out_ts_merged[x]$DNase_high,trait_NME_out_ts_merged[x]$DNase_total,p=p_exp,alternative = 'greater')$p.value))
# trait_NME_out_ts_merged$DNase_high_percent=trait_NME_out_ts_merged$DNase_high/trait_NME_out_ts_merged$DNase_total
# trait_NME_out_ts_merged$FDR=p.adjust(trait_NME_out_ts_merged$p_value_binom,method="BH")
# saveRDS(trait_NME_out_ts_merged,'../downstream/output/human_analysis/GWAS_agnostic_NME_germlayer.rds')
# trait_NME_out_ts_merged_sub=trait_NME_out_ts_merged[DNase_total>=20]
# NME_in_dt_mean=NME_in_dt[,list(NME_mean=mean(NME)),by=list(region,DNase)]


#Do it for each type of tissue
NME_in_dt=convert_GR(NME_in,direction="DT")
NME_in_dt$tissue=gsub(' - .*','',NME_in_dt$Sample)
NME_in_dt$tissue=gsub('Adipose_Tissue_single','Adipose_single',NME_in_dt$tissue)



trait_NME_out_ts=lapply(unique(NME_in_dt$tissue),function(ts){
  cat("Processing:",ts,'\n')
  trait_NME_out=data.table()
  NME_in_dt_ts=NME_in_dt[tissue==ts]
  #NME_in_dt_ts_mean=NME_in_dt_ts[,list(NME_mean=mean(NME)),by=list(region,DNase)]
  NME_in_dt_ts_mean=NME_in_dt_ts
  NME_in_dt_ts_mean$NME_mean=NME_in_dt_ts_mean$NME
  #NME_ecdf=ecdf(NME_in_dt_ts_mean$NME_mean)
  high_NME=quantile(NME_in_dt_ts_mean$NME[NME_in_dt_ts_mean$DNase=="DNase"],prob=0.75)
  NME_in_gr_mean_ts=convert_GR(NME_in_dt_ts_mean$region)
  mcols(NME_in_gr_mean_ts)=NME_in_dt_ts_mean
  variant_trait=GRanges(variant_trait)
  print(length(unique(variant_trait$`DISEASE/TRAIT`)))
  trait_NME_out=do.call(rbind,mclapply(unique(variant_trait$`DISEASE/TRAIT`),trait_overlap,
                                       NME_in_gr_mean_ts=NME_in_gr_mean_ts,
                                       variant_trait=variant_trait,high_NME=high_NME,
                                       mc.cores=20))
  trait_NME_out$tissue=ts
  print(trait_NME_out)
  return(trait_NME_out)
})
trait_NME_out_ts_merged=do.call(rbind,trait_NME_out_ts)
trait_NME_out_ts_merged=trait_NME_out_ts_merged[DNase_total>0]
trait_NME_out_ts_merged=readRDS('../downstream/output/human_analysis/GWAS/GWAS_agnostic_NME_DNase_mean.rds')


p_exp=sum(trait_NME_out_ts_merged$DNase_high)/sum(trait_NME_out_ts_merged$DNase_total)
trait_NME_out_ts_merged$p_value_binom=unlist(lapply(1:nrow(trait_NME_out_ts_merged),function(x) 
  binom.test(trait_NME_out_ts_merged[x]$DNase_high,trait_NME_out_ts_merged[x]$DNase_total,p=p_exp,alternative = 'greater')$p.value))
trait_NME_out_ts_merged$DNase_high_percent=trait_NME_out_ts_merged$DNase_high/trait_NME_out_ts_merged$DNase_total
trait_NME_out_ts_merged$FDR=p.adjust(trait_NME_out_ts_merged$p_value_binom,method="BH")
trait_NME_out_ts_merged_sub=trait_NME_out_ts_merged[DNase_total>=20]
NME_in_dt_mean=NME_in_dt[,list(NME_mean=mean(NME)),by=list(region,DNase)]
saveRDS(trait_NME_out_ts_merged,'../downstream/output/human_analysis/GWAS_agnostic_NME_tissue_mean_NME.rds')
trait_NME_out_ts_merged_sub=readRDS('../downstream/output/human_analysis/GWAS_agnostic_NME_tissue.rds')
NME_ecdf=ecdf(NME_in_dt_mean$NME_mean)
high_NME=quantile(NME_in_dt_mean$NME[NME_in_dt_mean$DNase=="DNase"],prob=0.75)
NME_in_gr_mean=convert_GR(NME_in_dt_mean$region)
mcols(NME_in_gr_mean)=NME_in_dt_mean
trait_NME_out=data.table()
variant_trait=GRanges(variant_trait)
#For mean 
for(trait in unique(variant_trait$`DISEASE/TRAIT`)){
  NME_trait=subsetByOverlaps(NME_in_gr_mean,variant_trait[variant_trait$`DISEASE/TRAIT`==trait])
  DNase_high=sum(NME_trait$NME_mean>=high_NME&NME_trait$DNase=="DNase")
  DNase_total=sum(NME_trait$DNase=="DNase")
  control_high=sum(NME_trait$NME_mean>=high_NME&NME_trait$DNase=="control")
  control_total=sum(NME_trait$DNase=="control")
  DNase_mean_quantile=mean(NME_ecdf(NME_trait[NME_trait$DNase=="DNase"]$NME_mean))
  control_mean_quantile=mean(NME_ecdf(NME_trait[NME_trait$DNase=="control"]$NME_mean))
  trait_NME_out=rbind(trait_NME_out,data.table(trait=trait,DNase_high=DNase_high,DNase_total=DNase_total,
                                               #control_high=control_high,control_total=control_total,
                                               DNase_mean_quantile=DNase_mean_quantile,
                                               control_mean_quantile=control_mean_quantile))
  
  
}

for(trait in unique(variant_trait$`DISEASE/TRAIT`)){
  NME_trait=subsetByOverlaps(NME_in,variant_trait[variant_trait$`DISEASE/TRAIT`==trait])
  DNase_high=sum(NME_trait$NME>=high_NME_mean&NME_trait$DNase=="DNase")
  DNase_total=sum(NME_trait$DNase=="DNase")
  control_high=sum(NME_trait$NME>=high_NME&NME_trait$DNase=="control")
  control_total=sum(NME_trait$DNase=="control")
  DNase_mean_quantile=mean(NME_ecdf(NME_trait[NME_trait$DNase=="DNase"]$NME))
  control_mean_quantile=mean(NME_ecdf(NME_trait[NME_trait$DNase=="control"]$NME))
  trait_NME_out=rbind(trait_NME_out,data.table(trait=trait,DNase_high=DNase_high,DNase_total=DNase_total,
                                               #control_high=control_high,control_total=control_total,
                                               DNase_mean_quantile=DNase_mean_quantile,
                                               control_mean_quantile=control_mean_quantile))
  
  
}
trait_NME_out$DNase_high_percent=trait_NME_out$DNase_high/trait_NME_out$DNase_total
trait_NME_out$control_high_percent=trait_NME_out$control_high/trait_NME_out$control_total 
saveRDS(trait_NME_out,'../downstream/output/human_analysis/GWAS_agnostic_NME_DNase_mean.rds')
#Do cancer only, do binomial test
trait_NME_out=readRDS('../downstream/output/human_analysis/GWAS_agnostic_NME_DNase.rds')
trait_NME_out=trait_NME_out[DNase_total>0]
p_exp=sum(trait_NME_out$DNase_high)/sum(trait_NME_out$DNase_total)
trait_NME_out$p_value_binom=unlist(lapply(1:nrow(trait_NME_out),function(x) binom.test(trait_NME_out[x]$DNase_high,trait_NME_out[x]$DNase_total,p=p_exp,alternative = 'greater')$p.value))
trait_NME_out$DNase_high_percent=trait_NME_out$DNase_high/trait_NME_out$DNase_total
trait_NME_out$FDR=p.adjust(trait_NME_out$p_value_binom,method="BH")
trait_NME_out$cancer=grepl('lympho|cancer|leukemia|melanoma|myeloma|carcinoma|Meningioma|Adverse response to chemotherapy|Cutaneous nevi|Glioma|tumor',trait_NME_out$trait,ignore.case = T)
write.csv(trait_NME_out[!is.na(DNase_high)],'../downstream/output/human_analysis/GWAS_agnostic_NME_DNase.csv')
ggplot(trait_NME_out,aes(x=p_value_binom))+geom_density()+xlab('P-values')+theme(legend.position = 'bottom')+
  geom_vline(xintercept =0.1)
t.test(trait_NME_out[cancer==TRUE]$DNase_high_percent,trait_NME_out[cancer==FALSE]$DNase_high_percent,alternative = "greater")
#Perfect control
#By chance
lapply(1:length(trait_NME_out))

#Color


#[X,Y,Z] = peaks(25);
#mesh(X,Y,Z)
# s.FaceColor = 'flat';
