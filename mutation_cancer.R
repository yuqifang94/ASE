source('mainFunctions_sub.R')
pan_mutation=fread('../downstream/input/human_analysis/cancer_SNP/ICGC_pan_cancer.tsv.gz')
nc=c("A","T","C","G")
pan_mutation_single=pan_mutation[(ref%in%nc)&(alt %in% nc),list(chr,pos,ref,alt,gene,driver,driver_statement,category)]#27568469
pan_mutation_single=makeGRangesFromDataFrame(pan_mutation_single,seqnames.field = "chr",start.field = "pos",end.field = "pos")
seqlevels(pan_mutation_single)=paste0("chr",seqlevels(pan_mutation_single))
variant_human=readRDS(variant_HetCpG_meta_file)
subsetByOverlaps(variant_human,pan_mutation_single)#No overlap

#From https://cancer.sanger.ac.uk/cosmic/download
#See scripted download starting with
#echo "email@example.com:mycosmicpassword" | base64


pan_mutation_coding=fread(paste0(cancer_input_dir,'CosmicMutantExport.tsv.gz'))
pan_mutation_coding=pan_mutation_coding[,list(`HGNC ID`,HGVSG,`Primary site`,`Primary histology`,`Genome-wide screen`,
                                              GENOMIC_MUTATION_ID,GRCh,`Mutation genome position`,`Mutation Description`,
                                              `FATHMM score`,`FATHMM prediction`)]
pan_mutation_noncoding=fread(paste0(cancer_input_dir,'CosmicNCV.tsv.gz'))
pan_mutation_noncoding=pan_mutation_noncoding[,list(HGVSG,`Primary site`,`Primary histology`,
                          GENOMIC_MUTATION_ID,GRCh,`genome position`,`Mutation somatic status`,
                          `FATHMM_MKL_NON_CODING_SCORE`)]
pan_mutation_noncoding$`HGNC ID`="Non-coding"
pan_mutation_noncoding$`Genome-wide screen`="Non-coding"
pan_mutation_noncoding$`Mutation genome position`=pan_mutation_noncoding$`genome position`
pan_mutation_noncoding$`genome position`=NULL
pan_mutation_noncoding$`Mutation Description`=pan_mutation_noncoding$`Mutation somatic status`
pan_mutation_noncoding$`Mutation somatic status`=NULL
pan_mutation_noncoding$`FATHMM score`=pan_mutation_noncoding$`FATHMM_MKL_NON_CODING_SCORE`
pan_mutation_noncoding$`FATHMM_MKL_NON_CODING_SCORE`=NULL
#FATHMM-MKL non-coding score. A p-value ranging from 0 to 1 where >= 0.7 is functionally significant.
pan_mutation_noncoding$`FATHMM prediction`="Not Functional"
pan_mutation_noncoding[`FATHMM score`>=0.7]$`FATHMM prediction`="Functional"
pan_mutation=rbind(pan_mutation_coding,pan_mutation_noncoding)
pan_mutation=readRDS('../downstream/input/pan_cancer_mutation_coding_non_coding.rds')
pan_mutation=pan_mutation[!is.na(`Mutation genome position`)& `Mutation genome position`!=""& !is.na(`FATHMM score`)]#60763076/65719644
pan_mutation_gr=GRanges(seqnames = sub(":.*", "", pan_mutation$`Mutation genome position`), 
                        IRanges(start = as.numeric(sub("-.*","", sub(".*:", "", pan_mutation$`Mutation genome position`))), 
                                end = as.numeric(sub(".*-", "",pan_mutation$`Mutation genome position`))), strand = "*")
mcols(pan_mutation_gr)=pan_mutation
saveRDS(pan_mutation_gr,'../downstream/input/pan_cancer_mutation_coding_non_coding_gr.rds')
pan_mutation_gr=readRDS('../downstream/input/human_analysis/cancer_SNP/pan_cancer_mutation_coding_non_coding_gr.rds')
rm(pan_mutation_coding)
rm(pan_mutation_noncoding)
ch = import.chain('../downstream/input/human_analysis/liftOver/hg38ToHg19.over.chain')
seqlevelsStyle(pan_mutation_gr) = "UCSC"  # necessary
pan_mutation_gr_19 = liftOver(pan_mutation_gr, ch)
pan_mutation_gr_19=unlist(pan_mutation_gr_19)
saveRDS(pan_mutation_gr_19,'../downstream/input/pan_cancer_mutation_coding_non_coding_gr19.rds')
human_variant=readRDS(variant_HetCpG_meta_file)
human_pan=subsetByOverlaps(human_variant,pan_mutation_gr_19)#442710/5357609
saveRDS(human_pan,'../downstream/output/human_pan.rds')
human_pan=readRDS('../downstream/output/human_analysis/cancer_analysis/human_pan.rds')
human_pan_dNME=human_pan[human_pan$dNME_pval<=pval_cutoff]
pan_mutation_gr_19=readRDS('../downstream/input/human_analysis/cancer_SNP/pan_cancer_mutation_coding_non_coding_gr19.rds')
motif_gene <- readRDS(motif_gene_file)
high_NME_variant=read.csv('../downstream/output/graphs/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
motif_gene_NME_sig=motif_gene[motif_gene$geneSymbol %in%high_NME_variant$TF]
human_pan_dNME_olap=findOverlaps(human_pan_dNME,motif_gene_NME_sig)
human_pan_dNME$alleleDiff=0 #Alt - ref
human_pan_dNME$alleleDiff[queryHits(human_pan_dNME_olap)]=motif_gene_NME_sig$alleleDiff[subjectHits(human_pan_dNME_olap)]
human_pan_dNME$NME_diff=human_pan_dNME$altNME-human_pan_dNME$refNME
human_pan_dNME_same_sign=human_pan_dNME[sign(human_pan_dNME$NME_diff)==sign(human_pan_dNME$alleleDiff)]#454/784
olap_pan_info=findOverlaps(human_pan_dNME_same_sign,pan_mutation_gr_19)
human_pan_dNME_same_sign$FATHMM_pred="NA"
human_pan_dNME_same_sign$FATHMM_pred[queryHits(olap_pan_info)]=pan_mutation_gr_19$`FATHMM prediction`[subjectHits(olap_pan_info)]
human_pan_dNME_same_sign$tissue_cancer="NA"
human_pan_dNME_same_sign$tissue_cancer[queryHits(olap_pan_info)]=pan_mutation_gr_19$`Primary site`[subjectHits(olap_pan_info)]
human_pan_dNME_same_sign=human_pan_dNME_same_sign[order(round(abs(human_pan_dNME_same_sign$NME_diff),digits = 2),decreasing=T)]
human_pan_dNME_same_sign_st=
  human_pan_dNME_same_sign[human_pan_dNME_same_sign$FATHMM_pred%in%c("Functional","PATHOGENIC")&
                             (human_pan_dNME_same_sign$refCG!=0|human_pan_dNME_same_sign$altCG!=0)]
write.csv(human_pan_dNME_same_sign_st,'../downstream/output/human_analysis/cancer_analysis/pan_cancer_dNME_SNP.csv')
subsetByOverlaps(pan_mutation_gr_19, human_pan_dNME_same_sign_st[10])
source('plotMB.R')
plotMB(subsetByOverlaps(motif_gene_NME_sig,human_pan_dNME_same_sign_st[10]),rsid="H9-311767")

