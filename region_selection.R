#Region selection
#Looking for NME separately
NME_allele_calc=readRDS('../downstream/output/NME_allele_calc.rds')
NME_allele_all_calc=NME_allele_calc[[1]]
genomic_features=readRDS('../downstream/input/genomic_features.rds')
#HUES64 stem
HUES64_stem=NME_allele_all_calc[NME_allele_all_calc$Sample=='stem_27_undifferentiated_paired - HUES64']
HUES64_stem_sub=HUES64_stem[which(width(HUES64_stem)>10 & HUES64_stem$pval<=0.1)]
#HUES64 ectoderm
HUES64_ectoderm=NME_allele_all_calc[NME_allele_all_calc$Sample=='ectoderm_paired - HUES64']
HUES64_ectoderm_sub=HUES64_ectoderm[which(width(HUES64_ectoderm)>10 & HUES64_ectoderm$pval<=0.1)]
#HUES64 mesoderm
HUES64_mesoderm=NME_allele_all_calc[NME_allele_all_calc$Sample=='mesoderm_23_paired - HUES64']
HUES64_mesoderm_sub=HUES64_mesoderm[which(width(HUES64_mesoderm)>10 & HUES64_mesoderm$pval<=0.1)]
#HUES64 endoerm
HUES64_endoerm=NME_allele_all_calc[NME_allele_all_calc$Sample=='endoerm_27_paired - HUES64']
HUES64_endoerm_sub=HUES64_endoerm[which(width(HUES64_endoerm)>10 & HUES64_endoerm$pval<=0.1)]

#Find regions shared by all
region_all=subsetByOverlaps(HUES64_stem,HUES64_ectoderm)
region_all=subsetByOverlaps(region_all,HUES64_endoerm)
region_all=subsetByOverlaps(region_all,HUES64_mesoderm)
#Add pval and diff of each
region_all_sum=granges(region_all)
region_all_sum$N=region_all$N
olap_stem=findOverlaps(HUES64_stem,region_all)
olap_mesoderm=findOverlaps(HUES64_mesoderm,region_all)
olap_endoerm=findOverlaps(HUES64_endoerm,region_all)
olap_ectoderm=findOverlaps(HUES64_ectoderm,region_all)
#Add pval
region_all_sum$stem_pval[subjectHits(olap_stem)]=HUES64_stem$pval[queryHits(olap_stem)]
region_all_sum$mesoderm_pval[subjectHits(olap_mesoderm)]=HUES64_mesoderm$pval[queryHits(olap_mesoderm)]
region_all_sum$endoerm_pval[subjectHits(olap_endoerm)]=HUES64_endoerm$pval[queryHits(olap_endoerm)]
region_all_sum$ectoderm_pval[subjectHits(olap_ectoderm)]=HUES64_ectoderm$pval[queryHits(olap_ectoderm)]
#Add diff
region_all_sum$stem_diff[subjectHits(olap_stem)]=HUES64_stem$diff[queryHits(olap_stem)]
region_all_sum$mesoderm_diff[subjectHits(olap_mesoderm)]=HUES64_mesoderm$diff[queryHits(olap_mesoderm)]
region_all_sum$endoerm_diff[subjectHits(olap_endoerm)]=HUES64_endoerm$diff[queryHits(olap_endoerm)]
region_all_sum$ectoderm_diff[subjectHits(olap_ectoderm)]=HUES64_ectoderm$diff[queryHits(olap_ectoderm)]
#find most different regions
region_all_sum$diff_sum=region_all_sum$stem_diff+region_all_sum$mesoderm_diff+
  region_all_sum$endoerm_diff+region_all_sum$ectoderm_diff
regions_all_sum_diff_sorted=region_all_sum[order(abs(region_all_sum$diff_sum),decreasing=T)]
#Shared by all 
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.1 & region_all_sum$mesoderm_pval<=0.1 &
                                                       region_all_sum$endoerm_pval<=0.1 & region_all_sum$ectoderm_pval<=0.1)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$stem_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=5]
#Only in mesoderm
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval>0.1 & region_all_sum$mesoderm_pval<=0.05 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$mesoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in mesoderm and stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval<=0.05 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$mesoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]

#Only in ectoderm
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval>0.1 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval<=0.05)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$ectoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in ectoderm and stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval<=0.05)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$ectoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in endoderm
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval>0.1 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval<=0.05 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$endoerm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in endoderm & stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval<=0.05 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$endoerm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]


#H1
H1=NME_allele_all_calc[NME_allele_all_calc$Sample=='merged - H1']
H1_sub=H1[which(H1$N>4 & abs(H1$diff)>0.65 & H1$pval<=0.1)]
H1_sub[order(abs(H1_sub$diff),decreasing = T)]
#High entropy region
H1_sub_diff=H1_sub[order(abs(H1_sub$diff),decreasing = T)]
H1_sub_diff_out=GRanges(H1_sub_diff)
H1_sub_diff_out$score=abs(H1_sub_diff$diff)
#Within 5k of promoter region
H1_sub_diff_TSS=H1_sub_diff[precede(H1_sub_diff,genomic_features$TSS)<=10000| follow(H1_sub_diff,genomic_features$TSS)<=10000]
olap=findOverlaps(H1_sub_diff_TSS,genomic_features$TSS,maxgap = 10000)
H1_sub_diff_TSS$gene=NA
H1_sub_diff_TSS$gene[queryHits(olap)]=genomic_features$TSS$gene_name[subjectHits(olap)]
H1_sub_diff_TSS_gene=H1_sub_diff_TSS[!is.na(H1_sub_diff_TSS$gene)]
H1_sub_diff_TSS_gene[order(H1_sub_diff_TSS_gene$N,decreasing = T)]
#Measure average distance, find nearst ones
H1_sub_diff_st=as.data.frame(sort(H1_sub_diff))
distance_gr=abs(H1_sub_diff_st$end[1:nrow(H1_sub_diff_st)-1]-H1_sub_diff_st$start[2:nrow(H1_sub_diff_st)])
sort(H1_sub_diff)[which(distance_gr<1000)]
#List of targets
#chr6 126072755-126072802 Hey2
# SNX18, AKAP11
#A gene expression profile of stem cell pluripotentiality and differentiation is conserved across diverse solid and hematopoietic cancers
GSA=read.table('../downstream/bamfiles/H1/GSA.txt',stringsAsFactors = F)
#With CpG shores
H1_sub_diff_shores=subsetByOverlaps(H1_sub_diff,genomic_features$`CpG shore`)
#H9
H9=NME_allele_all_calc[NME_allele_all_calc$Sample=='42_embryonic_stem_cell_single - H9']
H9_sub=H9[which(width(H9)>10 & H9$pval<=0.1)]
#GM12878
GM12878=NME_allele_all_calc[NME_allele_all_calc$Sample=='merged - GM12878']
GM12878_sub=GM12878[which(width(GM12878)>10 & GM12878$pval<=0.1)]

#Region selection MML
#Looking for NME separately
MML_allele_all_calc=readRDS('../downstream/output/MML_allele_calc.rds')
MML_allele_all_calc=MML_allele_all_calc[[1]]
genomic_features=readRDS('../downstream/input/genomic_features.rds')
#HUES64 stem
HUES64_stem=MML_allele_all_calc[MML_allele_all_calc$Sample=='stem_27_undifferentiated_paired - HUES64']
HUES64_stem_sub=HUES64_stem[which(width(HUES64_stem)>10 & HUES64_stem$pval<=0.1)]
#HUES64 ectoderm
HUES64_ectoderm=MML_allele_all_calc[MML_allele_all_calc$Sample=='ectoderm_paired - HUES64']
HUES64_ectoderm_sub=HUES64_ectoderm[which(width(HUES64_ectoderm)>10 & HUES64_ectoderm$pval<=0.1)]
#HUES64 mesoderm
HUES64_mesoderm=MML_allele_all_calc[MML_allele_all_calc$Sample=='mesoderm_23_paired - HUES64']
HUES64_mesoderm_sub=HUES64_mesoderm[which(width(HUES64_mesoderm)>10 & HUES64_mesoderm$pval<=0.1)]
#HUES64 endoerm
HUES64_endoerm=MML_allele_all_calc[MML_allele_all_calc$Sample=='endoerm_27_paired - HUES64']
HUES64_endoerm_sub=HUES64_endoerm[which(width(HUES64_endoerm)>10 & HUES64_endoerm$pval<=0.1)]

#Find regions shared by all
region_all=subsetByOverlaps(HUES64_stem,HUES64_ectoderm)
region_all=subsetByOverlaps(region_all,HUES64_endoerm)
region_all=subsetByOverlaps(region_all,HUES64_mesoderm)
#Add pval and diff of each
region_all_sum=granges(region_all)
region_all_sum$N=region_all$N
olap_stem=findOverlaps(HUES64_stem,region_all)
olap_mesoderm=findOverlaps(HUES64_mesoderm,region_all)
olap_endoerm=findOverlaps(HUES64_endoerm,region_all)
olap_ectoderm=findOverlaps(HUES64_ectoderm,region_all)
#Add pval
region_all_sum$stem_pval[subjectHits(olap_stem)]=HUES64_stem$pval[queryHits(olap_stem)]
region_all_sum$mesoderm_pval[subjectHits(olap_mesoderm)]=HUES64_mesoderm$pval[queryHits(olap_mesoderm)]
region_all_sum$endoerm_pval[subjectHits(olap_endoerm)]=HUES64_endoerm$pval[queryHits(olap_endoerm)]
region_all_sum$ectoderm_pval[subjectHits(olap_ectoderm)]=HUES64_ectoderm$pval[queryHits(olap_ectoderm)]
#Add diff
region_all_sum$stem_diff[subjectHits(olap_stem)]=HUES64_stem$diff[queryHits(olap_stem)]
region_all_sum$mesoderm_diff[subjectHits(olap_mesoderm)]=HUES64_mesoderm$diff[queryHits(olap_mesoderm)]
region_all_sum$endoerm_diff[subjectHits(olap_endoerm)]=HUES64_endoerm$diff[queryHits(olap_endoerm)]
region_all_sum$ectoderm_diff[subjectHits(olap_ectoderm)]=HUES64_ectoderm$diff[queryHits(olap_ectoderm)]
#find most different regions
region_all_sum$diff_sum=region_all_sum$stem_diff+region_all_sum$mesoderm_diff+
  region_all_sum$endoerm_diff+region_all_sum$ectoderm_diff
regions_all_sum_diff_sorted=region_all_sum[order(abs(region_all_sum$diff_sum),decreasing=T)]
#Shared by all 
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.1 & region_all_sum$mesoderm_pval<=0.1 &
                                                       region_all_sum$endoerm_pval<=0.1 & region_all_sum$ectoderm_pval<=0.1)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$stem_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=5]
#Only in mesoderm
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval>0.1 & region_all_sum$mesoderm_pval<=0.05 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$mesoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in mesoderm and stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval<=0.05 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$mesoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]

#Only in ectoderm
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval>0.1 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval<=0.05)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$ectoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in ectoderm and stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval>0.1 & region_all_sum$ectoderm_pval<=0.05)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$ectoderm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in endoderm
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval>0.1 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval<=0.05 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$endoerm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#Only in endoderm & stem
regions_all_sum_diff_sorted_sig=region_all_sum[which(region_all_sum$stem_pval<=0.05 & region_all_sum$mesoderm_pval>0.1 &
                                                       region_all_sum$endoerm_pval<=0.05 & region_all_sum$ectoderm_pval>0.1)]
regions_all_sum_diff_sorted_sig=regions_all_sum_diff_sorted_sig[order(abs(regions_all_sum_diff_sorted_sig$endoerm_diff),decreasing=T)]
regions_all_sum_diff_sorted_sig[regions_all_sum_diff_sorted_sig$N>=10]
#H1
H1_MML=MML_allele_all_calc[MML_allele_all_calc$Sample=='merged - H1']
H1_MML_sub=H1_MML[which(width(H1_MML)>10 & H1_MML$pval<=0.1)]
#High dMML region
H1_MML_sub_diff=H1_MML_sub[order(abs(H1_MML_sub$diff),decreasing = T)]
H1_MML_sub_diff_out=GRanges(H1_MML_sub_diff)
H1_MML_sub_diff_out$score=abs(H1_MML_sub_diff$diff)
export.bedGraph(H1_MML_sub_diff_out,'../downstream/output/dMML/H1_MML.bedGraph')
#Within 5k of promoter region
H1_MML_sub_diff_TSS=H1_MML_sub_diff[precede(H1_MML_sub_diff,genomic_features$TSS)<=10000| follow(H1_MML_sub_diff,genomic_features$TSS)<=10000]
olap=findOverlaps(H1_MML_sub_diff_TSS,genomic_features$TSS,maxgap = 10000)
H1_MML_sub_diff_TSS$gene=NA
H1_MML_sub_diff_TSS$gene[queryHits(olap)]=genomic_features$TSS$gene_name[subjectHits(olap)]
H1_MML_sub_diff_TSS_gene=H1_MML_sub_diff_TSS[!is.na(H1_MML_sub_diff_TSS$gene)]
H1_MML_sub_diff_TSS_gene[order(H1_MML_sub_diff_TSS_gene$N,decreasing = T)]
#Measure average distance, find nearst ones
H1_MML_sub_diff_st=as.data.frame(sort(H1_MML_sub_diff))
distance_gr=abs(H1_MML_sub_diff_st$end[1:nrow(H1_MML_sub_diff_st)-1]-H1_MML_sub_diff_st$start[2:nrow(H1_MML_sub_diff_st)])
sort(H1_MML_sub_diff)[which(distance_gr<1000)]
#List of targets
#chr6 126072755-126072802 Hey2
# SNX18, AKAP11
#A gene expression profile of stem cell pluripotentiality and differentiation is conserved across diverse solid and hematopoietic cancers
GSA=read.table('../downstream/bamfiles/H1_MML/GSA.txt',stringsAsFactors = F)
#With CpG shores
H1_MML_sub_diff_shores=subsetByOverlaps(H1_MML_sub_diff,genomic_features$`CpG shore`)
#H9
H9=MML_allele_all_calc[MML_allele_all_calc$Sample=='42_embryonic_stem_cell_single - H9']
H9_sub=H9[which(width(H9)>10 & H9$pval<=0.1)]
#GM12878
GM12878=MML_allele_all_calc[MML_allele_all_calc$Sample=='merged - GM12878']
GM12878_sub=GM12878[which(width(GM12878)>10 & GM12878$pval<=0.1)]
