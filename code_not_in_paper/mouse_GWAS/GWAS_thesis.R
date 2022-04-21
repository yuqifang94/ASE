source('mainFunctions_sub.R')
#This will be used in thesis---------------------------------------------------------------------------------------------------------------------------------------------------------
NME_in=readRDS(NME_agnostic_DNase_file)
DNase_hg19=readRDS(DNase_hg19_file)
control_hg19=readRDS(control_hg19_file)
variant_trait=readRDS('../downstream/input/human_analysis/variant_traits.rds')
NME_in$DNase=NA
NME_in_olap=findOverlaps(NME_in,DNase_hg19,type='equal')
NME_in$DNase[queryHits(NME_in_olap)]="DNase"
NME_in_olap=findOverlaps(NME_in,control_hg19,type='equal')
NME_in$DNase[queryHits(NME_in_olap)]="control"
variant_trait=GRanges(variant_trait)
SNP_NME_in$GWAS="Non GWAS"
variant_trait$cancer="Non-cancer"
variant_trait$cancer[grepl('lympho|cancer|leukemia|melanoma|myeloma|carcinoma|Meningioma|Adverse response to chemotherapy|Cutaneous nevi|Glioma|tumor',variant_trait$`DISEASE/TRAIT`,ignore.case = T)]=
  "cancer-related trait"
trait_enrich<-function(stat_in,variant_trait,ecdf=F,stat="NME"){
    if(ecdf){stat_ecdf=ecdf(mcols(stat_in)[[stat]]) }
    high_stat = quantile(mcols(stat_in)[[stat]],prob=0.75)
    trait_stat_out=data.table()
    for(trait in unique(variant_trait$`DISEASE/TRAIT`)){

        stat_trait=subsetByOverlaps(stat_in,variant_trait[variant_trait$`DISEASE/TRAIT`==trait])
        stat_high=sum(mcols(stat_trait)[[stat]]>=high_stat)
        stat_total=length(stat_trait)
        
        out = data.table(trait=trait,stat_high=stat_high,stat_total=stat_total)
        if(ecdf){out$mean_quantile=mean(stat_ecdf(stat_trait$stat))}
        trait_stat_out=rbind(trait_stat_out,out)
    
    
    }
    trait_stat_out$high_percent=trait_stat_out$stat_high/trait_stat_out$stat_total
    #Statistical testing
    trait_stat_out=trait_stat_out[stat_high>10]
    p_exp=sum(trait_stat_out$stat_high)/sum(trait_stat_out$stat_total)
  
    trait_stat_out$p_value_binom=unlist(lapply(1:nrow(trait_stat_out),function(x) binom.test(trait_stat_out[x]$stat_high,trait_stat_out[x]$stat_total,p=p_exp)$p.value))
    trait_stat_out$DNase_high_percent=trait_stat_out$stat_high/trait_stat_out$stat_total
    trait_stat_out$FDR=p.adjust(trait_stat_out$p_value_binom,method="BH")
    trait_stat_out$cancer=grepl('lympho|cancer|leukemia|melanoma|myeloma|carcinoma|Meningioma|Adverse response to chemotherapy|Cutaneous nevi|Glioma|tumor',trait_stat_out$trait,ignore.case = T)
return(trait_stat_out)
}
trait_NME_out=trait_enrich(NME_in[NME_in$DNase=="DNase"],variant_trait)
saveRDS(trait_NME_out,'../downstream/output/human_analysis/GWAS/GWAS_agnostic_NME_DNase2.rds')
write.csv(trait_NME_out[FDR<=0.2&DNase_high_percent>0.5][order(FDR,decreasing=F)],'../downstream/output/human_analysis/GWAS/GWAS_agnostic_NME_DNase.csv')
trait_NME_out_spAll=data.table()
for (sp in unique(NME_in$Sample)){
        trait_NME_out_sp=trait_enrich(NME_in[NME_in$DNase=="DNase"&NME_in$Sample==sp],variant_trait)
        trait_NME_out_sp$Sample=sp
        print(trait_NME_out_sp[FDR<=0.1])
        trait_NME_out_spAll=rbind(trait_NME_out_spAll,trait_NME_out_sp)
}
saveRDS(trait_NME_out_spAll,'../downstream/output/human_analysis/GWAS/GWAS_agnostic_NME_DNase_allsp.rds')
trait_NME_out_spAll_sig=trait_NME_out_spAll[FDR<=0.2&DNase_high_percent>0.25]
trait_NME_out_spAll_sig$pooled_analysis=trait_NME_out_spAll_sig$trait %in% trait_NME_out[FDR<=0.2]$trait
write.csv(trait_NME_out_spAll_sig,'../downstream/output/human_analysis/GWAS/GWAS_agnostic_NME_DNase_allsp.csv')
#Do cancer only, do binomial test
trait_NME_out=readRDS('../downstream/output/human_analysis/GWAS/GWAS_agnostic_NME_DNase2.rds')
trait_NME_out=trait_NME_out[stat_high>10]
p_exp=sum(trait_NME_out$DNase_high)/sum(trait_NME_out$DNase_total)
trait_NME_out$p_value_binom=unlist(lapply(1:nrow(trait_NME_out),function(x) binom.test(trait_NME_out[x]$stat_high,trait_NME_out[x]$stat_total,p=p_exp,alternative = 'greater')$p.value))
trait_NME_out$DNase_high_percent=trait_NME_out$stat_high/trait_NME_out$stat_total
trait_NME_out$FDR=p.adjust(trait_NME_out$p_value_binom,method="BH")
trait_NME_out$cancer=grepl('lympho|cancer|leukemia|melanoma|myeloma|carcinoma|Meningioma|Adverse response to chemotherapy|Cutaneous nevi|Glioma|tumor',trait_NME_out$trait,ignore.case = T)
write.csv(trait_NME_out[!is.na(stat_high)],'../downstream/output/human_analysis/GWAS/GWAS_agnostic_NME_DNase2.csv')