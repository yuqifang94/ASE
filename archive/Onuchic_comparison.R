##############################Onuchic paper###################################
if (!requireNamespace("jsonlite", quietly = TRUE))
{BiocManager::install("jsonlite")}
library(jsonlite)
Onuchic_SNP=fromJSON('../downstream/input/AllelicEpigenome-sigOnly-AllDocs.jsonld',flatten=FALSE)
#Extract data: in paper 241360 loci, here 240811 loci, there's pval missing
Onuchic_SNP_df=data.frame(REF=Onuchic_SNP$`Reference Allele`,
                          is_enc=Onuchic_SNP$`Is In Enhancer Region`,
                          GWAS=Onuchic_SNP$`Is Near GWAS Variant`,
                          chr=Onuchic_SNP$Chromosome,
                          start=Onuchic_SNP$Position,
                          ALT=Onuchic_SNP$`Alternative Allele`,
                          REF_met=Onuchic_SNP$`Combined Data Analysis`$`Ref Allele Methylated CpG Count`,
                          REF_unmet=Onuchic_SNP$`Combined Data Analysis`$`Ref Allele Unmethylated CpG Count`,
                          ALT_met=Onuchic_SNP$`Combined Data Analysis`$`Alt Allele Methylated CpG Count`,
                          ALT_unmet=Onuchic_SNP$`Combined Data Analysis`$`Alt Allele Unmethylated CpG Count`,
                          stringsAsFactors = F)
Onuchic_SNP_df$end=Onuchic_SNP_df$start
#After combining data, much more CpG in REF than in ALT, given more homozygous than heterozygous
Onuchic_SNP_df$REF_CpG=Onuchic_SNP_df$REF_met+Onuchic_SNP_df$REF_unmet
Onuchic_SNP_df$ALT_CpG=Onuchic_SNP_df$ALT_met+Onuchic_SNP_df$ALT_unmet
hist(log10(Onuchic_SNP_df$REF_CpG/Onuchic_SNP_df$ALT_CpG))
Onuchic_SNP_gr=makeGRangesFromDataFrame(Onuchic_SNP_df,keep.extra.columns = TRUE)
GR=readRDS('../downstream/output/GR.all.diff.H1.GM12878.rds')
GR=GR[!GR$Subject %in% c("H1","GM12878")]
GR_UC=GR[GR$Statistic=='UC' & GR$pvalue<=pval_cutoff]
GR_dMML=GR[GR$Statistic=='dMML' & GR$pvalue<=pval_cutoff]
GR_dNME=GR[GR$Statistic=='dNME']
#Find overlapped regions
Onuchic_SNP_gr=subsetByOverlaps(Onuchic_SNP_gr,GR)
Onuchic_SNP_gr$diff=abs(Onuchic_SNP_gr$REF_met/(Onuchic_SNP_gr$REF_met+Onuchic_SNP_gr$REF_unmet)-
                          Onuchic_SNP_gr$ALT_met/(Onuchic_SNP_gr$ALT_met+Onuchic_SNP_gr$ALT_unmet))

Onuchic_SNP_gr_dMML=subsetByOverlaps(Onuchic_SNP_gr,GR_merge_exclude_GM[GR_merge_exclude_GM$dMML_pval<=pval_cutoff])
Onuchic_SNP_gr_dMML=Onuchic_SNP_gr_dMML[Onuchic_SNP_gr_dMML$REF_CpG>100 &Onuchic_SNP_gr_dMML$ALT_CpG>100]
Onuchic_SNP_gr_dMML=Onuchic_SNP_gr_dMML[order(Onuchic_SNP_gr_dMML$diff,decreasing = T)]
GR_Onuchic_olap=findOverlaps(GR_merge_exclude_GM[GR_merge_exclude_GM$dMML_pval<=pval_cutoff],Onuchic_SNP_gr_dMML)
GR_dMML_Onuchic=GR_merge_exclude_GM[GR_merge_exclude_GM$dMML_pval<=pval_cutoff][queryHits(GR_Onuchic_olap)]
GR_dMML_Onuchic$Onuchic_diff=Onuchic_SNP_gr_dMML$diff[subjectHits(GR_Onuchic_olap)]
GR_dMML_Onuchic_top_dMML=GR_dMML_Onuchic[GR_dMML_Onuchic$Onuchic_diff>0.8 &GR_dMML_Onuchic$dMML>0.8]
unique(GR_dMML_Onuchic_top_dMML)[1:6]
#High dNME region can have little dMML


#Find overlaps 
length(unique(GR_dMML))#7699
length(unique(subsetByOverlaps(GR_UC,subsetByOverlaps(GR_dMML,GR_dNME,type='equal'),type='equal')))#2243
length(unique(subsetByOverlaps(GR_dMML,GR_dNME,type='equal')))#3254
length(unique(subsetByOverlaps(GR_dMML,GR_UC,type='equal')))#3749
length(unique(subsetByOverlaps(GR_dNME,GR_UC,type='equal')))#3063
length(unique(GR_UC))#7207
length(unique(GR_dNME))#207626
#Find overlap with Onuchic SNP
length(subsetByOverlaps(Onuchic_SNP_gr,GR_dMML))#1883
length(subsetByOverlaps(Onuchic_SNP_gr,GR_dNME))#5617
length(subsetByOverlaps(Onuchic_SNP_gr,GR_UC))#1677
length(subsetByOverlaps(Onuchic_SNP_gr,GR,maxgap = 100))
#Find dMML data
GR_dMML=GR[GR$Statistic=='dMML']
GR_dNME=GR[GR$Statistic=='dNME']
GR_dMML_SD_ASM=subsetByOverlaps(GR_dMML,Onuchic_SNP_gr)

#Perchr data summary
fp='../../../../../../allele_specific/AllelicEpigenome-AllChrs-AllDocs.jsonld/'
chr_df=list(Subject=c(),Chromosome=c(),Position=c(),`Reference Allele`=c(),`Alternative Allele`=c(),
            `Is In Enhancer Region`=c(), `Is Near GWAS Variant`=c(),`Is On Heterogenous CpG`=c(),
            `Is In Promoter Region`=c(),`Tissue Specific Analysis`=list(),
            `Transcription Factor`=c(),`1000 Genomes Allele Frequency`=c(),
            REF_met=c(), REF_unmet=c(),ALT_met=c(),ALT_unmet=c())
for (fin in dir(fp)){
  tt1=proc.time()[[1]]
  cat("Reading", fin,'\n')
  gc()
  chr_in=fromJSON(paste(fp,fin,sep=''),flatten=FALSE)
  cat("processing",fin,'\n')
  chr_in_list=chr_in[c("Subject","Chromosome","Position","Reference Allele","Alternative Allele",
                       "Is In Enhancer Region", "Is Near GWAS Variant","Is On Heterogenous CpG",
                       "Is In Promoter Region","Tissue Specific Analysis",
                       "Transcription Factor","1000 Genomes Allele Frequency")]
  chr_in_list$REF_met=chr_in$`Combined Data Analysis`$`Ref Allele Methylated CpG Count`
  chr_in_list$REF_unmet=chr_in$`Combined Data Analysis`$`Ref Allele Unmethylated CpG Count`
  chr_in_list$ALT_met=chr_in$`Combined Data Analysis`$`Alt Allele Methylated CpG Count`
  chr_in_list$ALT_unmet=chr_in$`Combined Data Analysis`$`Alt Allele Unmethylated CpG Count`
  chr_df=mapply(c,chr_df,chr_in_list,SIMPLIFY=F)
  cat(tail(chr_df$REF_met),'\n')
  cat("Finish processing",fin,"in",proc.time()[[1]]-tt1,'\n')
}
saveRDS(chr_df,"D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/Onuchic_SNP.rds")
Onuchic_SNP <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/Onuchic_SNP.rds")
#Make granges over Onuchic SNP
Onuchic_SNP_df=data.frame(REF=Onuchic_SNP$`Reference Allele`,
                          GWAS=Onuchic_SNP$`Is Near GWAS Variant`,
                          chr=Onuchic_SNP$Chromosome,
                          start=Onuchic_SNP$Position,
                          ALT=Onuchic_SNP$`Alternative Allele`,
                          REF_met=Onuchic_SNP$REF_met,
                          REF_unmet=Onuchic_SNP$REF_unmet,
                          ALT_met=Onuchic_SNP$ALT_met,
                          ALT_unmet=Onuchic_SNP$ALT_unmet,
                          stringsAsFactors = F)
Onuchic_SNP_df=Onuchic_SNP_df[!is.na(Onuchic_SNP_df$REF_met),]
Onuchic_SNP_df$end=Onuchic_SNP_df$start
#After combining data, much more CpG in REF than in ALT, given more homozygous than heterozygous
Onuchic_SNP_df$REF_CpG=Onuchic_SNP_df$REF_met+Onuchic_SNP_df$REF_unmet
Onuchic_SNP_df$ALT_CpG=Onuchic_SNP_df$ALT_met+Onuchic_SNP_df$ALT_unmet
hist(log10(Onuchic_SNP_df$REF_CpG/Onuchic_SNP_df$ALT_CpG))
Onuchic_SNP_gr=makeGRangesFromDataFrame(Onuchic_SNP_df,keep.extra.columns = TRUE)
Onuchic_SNP_gr$diff=abs(Onuchic_SNP_gr$REF_met/Onuchic_SNP_gr$REF_CpG-Onuchic_SNP_gr$ALT_met/Onuchic_SNP_gr$ALT_CpG)
#Find overlap with dNME region
dNME_Onuchic=subsetByOverlaps(GR_merge_exclude_GM[GR_merge_exclude_GM$dNME_pval<=pval_cutoff],Onuchic_SNP_gr)
dNME_Onuchic$N=NA
dNME_Onuchic$Subject=unlist(lapply(dNME_Onuchic$Sample,function(x) strsplit(x,' - ')[[1]][2]))
for (sp in unique(dNME_Onuchic$Subject)){
  olap=findOverlaps(dNME_Onuchic[dNME_Onuchic$Subject==sp],gff_in[gff_in$Subject==sp])
  dNME_Onuchic$N[dNME_Onuchic$Subject==sp][queryHits(olap)]=gff_in$N[gff_in$Subject==sp][subjectHits(olap)]
  
  
}
dNME_Onuchic=dNME_Onuchic[dNME_Onuchic$N>=10]
dNME_Onuchic_high_dNME=dNME_Onuchic[dNME_Onuchic$dNME>=0.8]
unique(dNME_Onuchic_high_dNME[order(dNME_Onuchic_high_dNME$dMML,decreasing = F)])#2nd one
subsetByOverlaps(Onuchic_SNP_gr, unique(dNME_Onuchic_high_dNME[order(dNME_Onuchic_high_dNME,decreasing = F)])[2])
dNME_Onuchic_low_dMML=dNME_Onuchic[dNME_Onuchic$dMML<=0.1 &dNME_Onuchic$dNME>=0.6]
dNME_Onuchic_low_dMML[order(dNME_Onuchic_low_dMML$dMML)]

#True is no data, false is have data
Onuchic_SNP_gr_raw$Tissue_data=is.na(Onuchic_SNP$`Tissue Specific Analysis`) | unlist(lapply(Onuchic_SNP$`Tissue Specific Analysis`,is.null))
SNP_sig=readRDS('../downstream/output/Onuchic_SNP_gr.rds')
#Number of SD-ASM overlapping Onuchic SNP
Onuchic_SNP_gr_sig=subsetByOverlaps(Onuchic_SNP_gr_raw,SNP_sig)
variant_HetCpG_new <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/input/variant_HetCpG_new.rds")
variant_HetCpG_new=variant_HetCpG_new[1:7]
names(variant_HetCpG_new)=NULL
variant_HetCpG_new=do.call(c,variant_HetCpG_new)
olap_SNP=findOverlaps(Onuchic_SNP_gr_raw,variant_HetCpG_new)
subsetByOverlaps(Onuchic_SNP_gr_raw,variant_HetCpG_new)
non_olap=which(!1:length(Onuchic_SNP_gr_raw) %in% queryHits(olap_SNP))
non_olap_not_na=non_olap[!non_olap %in% which(Onuchic_SNP_gr_raw$Tissue_data)]
#Which patient are those
non_olap_patient=do.call("c",lapply(Onuchic_SNP$`Tissue Specific Analysis`[non_olap_not_na],function(x) x$Patient))
#GR analyzed
GR <- readRDS("D:/Dropbox/JHU/Projects/Allele-spcific/code/downstream/output/GR.all.diff.H1.GM12878.rds")
GR=GR[!GR$Subject %in% c("H1","GM12878")]
Onuchic_SNP_gr_vcf=subsetByOverlaps(Onuchic_SNP_gr_raw,variant_HetCpG_new)
SNP_with_stat=subsetByOverlaps(Onuchic_SNP_gr_vcf,GR,maxgap = 100)#Might be some double counts, maxgap=100
#Check SNP for those we don't have stat
SNP_with_stat_olap=findOverlaps(Onuchic_SNP_gr_vcf,GR,maxgap = 100)
Onuchic_SNP_gr_vcf_noolap=Onuchic_SNP_gr_vcf[!(1:length(Onuchic_SNP_gr_vcf) %in% queryHits(SNP_with_stat_olap))]
Onuchic_SNP_gr_vcf_noolap$SNP=variants_collapase(Onuchic_SNP_gr_vcf_noolap)
SNP_with_stat$SNP=variants_collapase(SNP_with_stat)

SNP_are_sig=subsetByOverlaps(Onuchic_SNP_gr_sig,Onuchic_SNP_gr_vcf)
subsetByOverlaps(SNP_are_sig,SNP_with_stat)
sum(Onuchic_SNP_gr_sig$Tissue_data)
SNP_tissue_no_WGS=subsetByOverlaps(Onuchic_SNP_gr_sig[!Onuchic_SNP_gr_sig$Tissue_data],Onuchic_SNP_gr_raw[non_olap_not_na])
subsetByOverlaps(Onuchic_SNP_gr_vcf,GR.all.diff)
sum(1:length(Onuchic_SNP_gr_vcf) %in% queryHits(findOverlaps(Onuchic_SNP_gr_vcf,GR.all.diff)))
sum(!1:length(Onuchic_SNP_gr_vcf) %in% queryHits(findOverlaps(Onuchic_SNP_gr_vcf,GR.all.diff)))
#Finding SNPs that have Stats
Onuchic_SNP_gr_olap=findOverlaps(Onuchic_SNP_gr_raw,GR[GR$Statistic=='dMML'],maxgap = 100)
Onuchic_SNP_qt=unique(queryHits(Onuchic_SNP_gr_olap))#1163614, maxgap=200:2314666
Onuchic_SNP_gr_df=as.data.frame(Onuchic_SNP_gr_raw[Onuchic_SNP_qt])
#Getting ASM stats
library(data.table)
Onuchic_SNP_stats=lapply(Onuchic_SNP$`Tissue Specific Analysis`[Onuchic_SNP_qt],extract_stats)
Onuchic_SNP_gr_stats=lapply(1:length(Onuchic_SNP_stats), function(x) merge_stat(Onuchic_SNP_stats[[x]],Onuchic_SNP_gr_df[x,]))
Onuchic_SNP_gr_stats=rbindlist(Onuchic_SNP_gr_stats)
Onuchic_SNP_gr_stats=as.data.frame(Onuchic_SNP_gr_stats)
Onuchic_SNP_gr_stats$end=Onuchic_SNP_gr_stats$start
Onuchic_SNP_gr_stats_sig=Onuchic_SNP_gr_stats[Onuchic_SNP_gr_stats$pval<=pval_cutoff,]
dMML=GR[GR$Statistic=='dMML']
Onuchic_SNP_gr_stats=makeGRangesFromDataFrame(Onuchic_SNP_gr_stats,keep.extra.columns = T)

#Find number of overlap with each tissue-sample combination, both have data
olap_single_sample=data.frame()

for (subj in unique(dMML$Subject)){
  for (ts in unique(dMML$Tissue)){
    cat("processing", subj,ts,'\n')
    dMML_sp=dMML[dMML$Subject==subj & dMML$Tissue==ts]
    Onuchic_sp=Onuchic_SNP_gr_stats[Onuchic_SNP_gr_stats$subject == subj & Onuchic_SNP_gr_stats$sample==ts]
    if (length(dMML_sp)!=0){
      cat("Counting", subj,ts,'\n')
      olap_sig=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
      olap_sig_cpel=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
      olap_sig_Onuchic=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
      olap_nonsig=length(subsetByOverlaps(dMML_sp[dMML_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
      olap_single_sample=rbind(olap_single_sample,data.frame(sp=paste(subj,ts,sep='-'),olap_sig=olap_sig,
                                                             olap_sig_cpel=olap_sig_cpel,olap_sig_Onuchic=olap_sig_Onuchic,olap_nonsig=olap_nonsig))
    }
  }
}


saveRDS(Onuchic_SNP_gr_stats,'../downstream/output/Onuchic_SNP_df_within_GR.rds')
#reform all regions
Onuchic_SNP_gr_df=as.data.frame(Onuchic_SNP_gr_raw)
Onuchic_SNP_stats=lapply(Onuchic_SNP$`Tissue Specific Analysis`,extract_stats)
Onuchic_SNP_gr_stats=lapply(1:length(Onuchic_SNP_stats), function(x) merge_stat(Onuchic_SNP_stats[[x]],Onuchic_SNP_gr_df[x,]))
Onuchic_SNP_gr_stats=rbindlist(Onuchic_SNP_gr_stats)
Onuchic_SNP_gr_stats$end=Onuchic_SNP_gr_stats$start
Onuchic_SNP_gr_stats=makeGRangesFromDataFrame(as.data.frame(Onuchic_SNP_gr_stats),keep.extra.columns = T)
saveRDS(Onuchic_SNP_gr_stats,'../downstream/output/Onuchic_SNP_df_all.rds')
#Check if regions having stats
#Get rid of regions with unavailable dataset
Onuchic_SNP_gr_stats_ft=Onuchic_SNP_gr_stats[!(Onuchic_SNP_gr_stats$subject %in% c("112","149","150","HuFGM02") | 
                                                 Onuchic_SNP_gr_stats$sample == "Penis, Foreskin, Fibroblast Primary Cells")]
dMML_sp=single_sample(GR,Onuchic_SNP_gr_stats_ft,'dMML')
dNME_sp=single_sample(GR,Onuchic_SNP_gr_stats_ft,'dNME')
single_sample<-function(CPEL_in,Onuchic_in,stat_in){
  #Matching names
  GR_no_olap=GRanges()
  CPEL_in=CPEL_in[CPEL_in$Statistic==stat_in]
  CPEL_in$Tissue[CPEL_in$Tissue=='Adipose']="Adipose Tissue"
  CPEL_in$Tissue[CPEL_in$Tissue=='Embryonic Stem Cell']="H9 Cell Line"
  CPEL_in$Tissue[CPEL_in$Tissue=='Foreskin Melanocyte']="Penis, Foreskin, Melanocyte Primary Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Foreskin Keratinocyte']="Penis, Foreskin, Keratinocyte Primary Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Embyonic Stem Cell']="HUES64 Cell Line"
  CPEL_in$Tissue[CPEL_in$Tissue=='Ectoderm']="hESC Derived CD56+ Ectoderm Cultured Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Mesoderm']="hESC Derived CD56+ Mesoderm Cultured Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Endoderm']="hESC Derived CD184+ Endoderm Cultured Cells"
  CPEL_in$Tissue[CPEL_in$Tissue=='Liver']="Adult Liver"
  olap_single_sample=data.frame()
  for (subj in unique(CPEL_in$Subject)){
    for (ts in unique(CPEL_in$Tissue)){
      cat("processing", subj,ts,'\n')
      CPEL_in_sp=CPEL_in[CPEL_in$Subject==subj & CPEL_in$Tissue==ts]
      Onuchic_sp=Onuchic_in[Onuchic_in$subject == subj & Onuchic_in$sample==ts]
      if (length(CPEL_in_sp)!=0){
        cat("Counting", subj,ts,'\n')
        olap=findOverlaps(Onuchic_sp,CPEL_in_sp,maxgap = 100)
        GR_no_olap=c(GR_no_olap,Onuchic_sp[!(1:length(Onuchic_sp) %in% queryHits(olap))])
        olap_GR=length(subsetByOverlaps(Onuchic_sp,CPEL_in_sp,maxgap = 100))
        olap_sig=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
        olap_sig_cpel=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue<=pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
        olap_sig_Onuchic=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval<=pval_cutoff],maxgap = 100))
        olap_nonsig=length(subsetByOverlaps(CPEL_in_sp[CPEL_in_sp$pvalue>pval_cutoff],Onuchic_sp[Onuchic_sp$pval>pval_cutoff],maxgap = 100))
        olap_single_sample=rbind(olap_single_sample,data.frame(olap=olap_GR,non_olap=length(Onuchic_sp)-olap_GR,sp=paste(subj,ts,sep='-'),olap_sig=olap_sig,
                                                               olap_sig_cpel=olap_sig_cpel,olap_sig_Onuchic=olap_sig_Onuchic,olap_nonsig=olap_nonsig))
      }
    }
  }
  return(list(olap_single_sample,GR_no_olap))
}
Onuchic_SNP_gr_stats_ft_with_stat=subsetByOverlaps(Onuchic_SNP_gr_stats_ft,GR,maxgap = 100)
sum(Onuchic_SNP_gr_stats$subject %in% c("112","149","150","HuFGM02") | Onuchic_SNP_gr_stats$sample == "Penis, Foreskin, Fibroblast Primary Cells" )
#Extract stats from Onuchic paper
extract_stats<-function(dat_in){
  dat_in=dat_in[dat_in$Experiment$`foaf:name`=="Bisulfite-Seq",]
  return(data.frame(subject=dat_in$Patient,sample=dat_in$Tissue$`rdfs:label`,diff=dat_in$`Methylation Difference`,pval=dat_in$`FDR P-value`,
                    REF_unmet=dat_in$`Ref Allele Unmethylated CpG Count`,
                    REF_met=dat_in$`Ref Allele Methylated CpG Count`,
                    ALT_unmet=dat_in$`Alt Allele Unmethylated CpG Count`,
                    ALT_met=dat_in$`Alt Allele Methylated CpG Count`,
                    stringsAsFactors = F))
}
#Extract location for each stat
merge_stat<-function(stat_df,gr_df){
  if(nrow(stat_df)>0){
    stat_df$chr=gr_df$seqnames
    stat_df$start=gr_df$start
    return(stat_df)
  }
  
}