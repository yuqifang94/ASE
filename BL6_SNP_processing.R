
#separate BL6 comparison and the mm10 comparison  
#read in vcf
BL6DBA_SNP=read.table('../vcf/DBA_2J.mgp.v5.snps.dbSNP142.chr.vcf')
BL6DBA_SNP=BL6DBA_SNP[,c(1,2,4,5)]
colnames(BL6DBA_SNP)[c(1,2)]=c('chr','start')
BL6DBA_SNP$SNP=paste(BL6DBA_SNP$V4,BL6DBA_SNP$V5,sep='/')
BL6DBA_SNP$V4=BL6DBA_SNP$V5=NULL
BL6DBA_SNP=makeGRangesFromDataFrame(BL6DBA_SNP,end.field='start',keep.extra.columns = T)
mm10_region=import.gff('mm10_allele_agnostic_analysis.gff')
#doing it for each chromosome
mm10_region_het=mclapply(as.character( unique(seqnames(mm10_region))),function(x,gff_in,vcf_in){
  vcf_in=vcf_in[seqnames(vcf_in)==x]
  gff_in=gff_in[seqnames(gff_in)==x]
  gff_in_sub=subsetByOverlaps(gff_in,vcf_in)
  #form CpG Granges
  CpG=as.numeric(unlist(lapply(gff_in_sub$CpGs,function(x){
    x=gsub('\\]','',x)
    x=gsub('\\[','',x)
    gsub(' ','',x)
  })))
  CpG=GRanges(seqnames = x,ranges=IRanges(start=CpG,end=CpG+1))
  hetCpG=subsetByOverlaps(CpG,vcf_in)
  hetCpG=resize(hetCpG,width=1,fix='start')
  olap=findOverlaps(gff_in,hetCpG)
  gff_in$hetCpG=FALSE
  gff_in$hetCpG[queryHits(olap)]=TRUE
  return(gff_in)
},gff_in=mm10_region,vcf_in=BL6DBA_SNP,mc.cores=10)
mm10_region_het=do.call('c',mm10_region_het)
olap=findOverlaps(UC_in,mm10_region_het[mm10_region_het$hetCpG],type='equal')
UC_in=UC_in[-queryHits(olap)]
#sbatch cpelasm_allele_agnostic_uc.slrm mm10_liver_day13_5_merged2 mm10_liver_day14_5_merged2


UC_in_Epiblast=UC_in[UC_in$tissue1=='Epiblast']
UC_in_ICM=UC_in[UC_in$tissue1=='ICM']
saveRDS(UC_in,'UC_agnostic_mouse_all_merged_BL6DBA.rds')
saveRDS(UC_in_Epiblast,'UC_agnostic_mouse_all_merged_BL6DBA_Epiblast.rds')
saveRDS(UC_in_ICM,'UC_agnostic_mouse_all_merged_BL6DBA_ICM.rds')

UC_in_matrix_ls_Epiblast=agnostic_matrix_conversion(UC_in_Epiblast,'UC')
UC_in_Epiblast= UC_in_Epiblast[UC_in_Epiblast$N>=2]
UC_in_matrix_ls_Epi=mclapply(unique(UC_in_Epiblast$tissue2),
                             function(x) agnostic_matrix_conversion(UC_in_Epiblast[UC_in_Epiblast$tissue2==x],'UC'),mc.cores=6)
UC_in_ICM=UC_in_ICM[UC_in_ICM$N>=2]
UC_in_matrix_ls_ICM=mclapply(unique(UC_in_ICM$tissue2),
                             function(x) agnostic_matrix_conversion(UC_in_ICM[UC_in_ICM$tissue2==x],'UC'),mc.cores=6)
saveRDS(UC_in_matrix_ls_Epiblast,'UC_agnostic_mouse_all_matrix_dedup_N2_all_merged_ls_Epiblast.rds')#74% regiOn have all data
saveRDS(UC_in_matrix_ls_ICM,'UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_ICM.rds')#74% regiOn have all data

saveRDS(UC_in_DNase[UC_in_DNase$tissue1=='Epiblast'],'UC_agnostic_mouse_all_merged_BL6DBA_Epiblast_DNase.rds')
saveRDS(UC_in_Epiblast_PRC,'UC_agnostic_mouse_all_merged_BL6DBA_Epiblast_PRC.rds')
saveRDS(UC_in_ICM_DNase,'UC_agnostic_mouse_all_merged_BL6DBA_ICM_DNase.rds')
saveRDS(UC_in_ICM_PRC,'UC_agnostic_mouse_all_merged_BL6DBA_ICM_PRC.rds')