rm(list=ls())
source("mainFunctions_sub.R")

stat_merge<-function(gr_in,allele_in,vcf_in,CpG){
  #Check merge behavior, 
  dMML=gr_in[gr_in$Statistic=="dMML"]
  dNME=gr_in[gr_in$Statistic=="dNME"]
  UC=gr_in[gr_in$Statistic=="UC"]
  gr=unique(granges(gr_in))
  olap_dMML=findOverlaps(gr,dMML,type="equal")
  gr$dMML[queryHits(olap_dMML)]=dMML$Value[subjectHits(olap_dMML)]
  gr$dMML_pval[queryHits(olap_dMML)]=dMML$pvalue[subjectHits(olap_dMML)]
  gr$Sample=unique(gr_in$Sample)
  
  olap_dNME=findOverlaps(gr,dNME,type="equal")
  gr$dNME[queryHits(olap_dNME)]=dNME$Value[subjectHits(olap_dNME)]
  gr$dNME_pval[queryHits(olap_dNME)]=dNME$pvalue[subjectHits(olap_dNME)]
  
  olap_UC=findOverlaps(gr,UC,type="equal")
  gr$UC[queryHits(olap_UC)]=UC$Value[subjectHits(olap_UC)]
  gr$UC_pval[queryHits(olap_UC)]=UC$pvalue[subjectHits(olap_UC)]
  
  gr$NME1=gr$MML1=gr$NME2=gr$MML2=gr$N=NA
  #genome 1 NME
  olap_NME1=findOverlaps(gr,allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME'],type="equal")
  gr$NME1[queryHits(olap_NME1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME']$Value[subjectHits(olap_NME1)]
  gr$N[queryHits(olap_NME1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME']$N[subjectHits(olap_NME1)]
  #genome 2 NME
  olap_NME2=findOverlaps(gr,allele_in[allele_in$Genome==2 & allele_in$Statistic=='NME'],type="equal")
  gr$NME2[queryHits(olap_NME2)]=allele_in[allele_in$Genome==2 & allele_in$Statistic=='NME']$Value[subjectHits(olap_NME2)]
  #genome 1 MML
  olap_MML1=findOverlaps(gr,allele_in[allele_in$Genome==1 & allele_in$Statistic=='MML'],type="equal")
  gr$MML1[queryHits(olap_MML1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='MML']$Value[subjectHits(olap_MML1)]
  #genome 1 MML
  olap_MML2=findOverlaps(gr,allele_in[allele_in$Genome==2 & allele_in$Statistic=='MML'],type="equal")
  gr$MML2[queryHits(olap_MML2)]=allele_in[allele_in$Genome==2 & allele_in$Statistic=='MML']$Value[subjectHits(olap_MML2)]
  gr$Subject=unique(gr_in$Subject)
  gr$tissue=unique(gr_in$tissue)
  #Add g1cg etc information to the gr
  gr=add_hetCPG(gr,vcf_in,CpG)
  return(gr)
}
#Put hetCpG count into each sample in gr_allele

##redo instead of adding allele information to allele CpG, add it to GR_merge
#This is mainly modifed to fit the output of CPEL ASM

add_hetCPG<-function(gr,vcf_in,CpG){
  olap=findOverlaps(vcf_in,gr,type='within',select='all')
  #Count number of CG, here g1CG=genome1 and g2CG=genome2, however, we need to split based on the GT, also ref
  vcf_in$g1CG=as.numeric(grepl("CG",vcf_in$genome1_tri))
  vcf_in$g2CG=as.numeric(grepl("CG",vcf_in$genome2_tri))
  vcf_in$refCG=as.numeric(grepl("CG",vcf_in$REF_tri))
  vcf_in$altCG=as.numeric(grepl("CG",vcf_in$ALT_tri))
  #Count number of Het CpG here *Check 
  df_sub=data.table(subjHits=subjectHits(olap),
                    g1CG=vcf_in$g1CG[queryHits(olap)],g2CG=vcf_in$g2CG[queryHits(olap)],
                    refCG=vcf_in$refCG[queryHits(olap)],altCG=vcf_in$altCG[queryHits(olap)])
  
  agg_sub=df_sub[,.(g1CG=sum(g1CG),g2CG=sum(g2CG),refCG=sum(refCG),altCG=sum(altCG)),by=subjHits]
  gr$g1CG[agg_sub$subjHits]=agg_sub$g1CG#agg_sub$subjHits is the unique subject hits
  gr$g2CG[agg_sub$subjHits]=agg_sub$g2CG
  gr$refCG[agg_sub$subjHits]=agg_sub$refCG#agg_sub$subjHits is the unique subject hits
  gr$altCG[agg_sub$subjHits]=agg_sub$altCG
  gr$N_hg19=countOverlaps(gr,CpG)
  #count number of CG lost from reference
  gr$N_nonhet=gr$N_hg19-countOverlaps(gr,vcf_in[(vcf_in$refCG-vcf_in$altCG)>0])
  
  return(gr)
}


#add genes to GR
add_gene_GR<-function(GR_merge_in,feature_in,feature_names){
  mcols(GR_merge_in)[[feature_names]]=NA
  olap_promoter=findOverlaps(GR_merge_in,feature_in)
  df_idx=data.table(qt=queryHits(olap_promoter),
                    genes=as.character(feature_in$gene_name[subjectHits(olap_promoter)]),
                    stringsAsFactors = F)
  df_idx=df_idx[,list(genes=list(genes)),by=list(df_idx$qt)]
  colnames(df_idx)=c('qt','genes')
  df_idx$genes=lapply(df_idx$genes,as.character)
  mcols(GR_merge_in)[[feature_names]][df_idx$qt]=df_idx$genes
  return(GR_merge_in)
  
}
# get all hg19 CpG site ---------------------------------------------------

CpG_hg19=readRDS('../downstream/input/human_analysis/CpG_hg19.rds')
GR=readRDS(GR_file)
GR_allele=readRDS(GR_allele_file)
variant_HetCpG=readRDS(variant_HetCpG_file)
# creating merged object --------------------------------------------------
GR_merge=GRanges()
for(sp in unique(GR$Sample)){
  cat("Processing",sp,'\n')
  ts=gsub('.*- ','',sp)
  GR_merge=c(GR_merge,stat_merge(GR[GR$Sample==sp],
                                 GR_allele[GR_allele$Sample==sp],
                                 variant_HetCpG[[ts]],CpG_hg19))
}
#Count Number of Het CpG at extened regions, each extend 500 bp
hetCGallele_merged<-function(sub,gr_merge,CpG,vcf_in,gene_size=500){
  cat('Analyzing',sub,'\n')
  t1=proc.time()[[3]]
  #Import vcf file
  sub_vcf=vcf_in[[sub]]
  sub_allele=gr_merge[gr_merge$Subject==sub]
  #gr_allele got resized with start and end +1, use +2 to include equal start & end, no longer needed for new output?
  #sub_het=resize(sub_het, width(sub_het) + 4, fix="center")
  sub_allele$gff_size=width(sub_allele)
  gr_out=GR_resize_merged(sub_allele,CpG,sub_vcf,gene_size=gene_size,sub)#Check for calculating density
  cat('Finish analyzing',sub,proc.time()[[3]]-t1,'\n')
  return(gr_out)
}
#Generate density from GR allele, at least need 200 bp for density
#Count number of hetCG at each allele GR from result, add gff size, and N
GR_resize_merged<-function(GR.in,CpG_sites,hetCpG,sub,gene_size=500){
  ##From definitino of CpG island, use 200 bp regions
  GR.extend=resize(GR.in,width=width(GR.in)+gene_size,fix='center')
  GR.in$CG_hg19_extend=countOverlaps(GR.extend,CpG_sites)
  #any SNP contain CG in ref genome
  GR.in$CG_nonhet_extend=GR.in$CG_hg19_extend-countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$REF_tri)])
  gr_seq=getSeq(Hsapiens,GR.extend,as.character=T)
  GR.in$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
  #Count CpG in genome 1
  GR.in$CG_allele_extend_g1=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$genome1_tri)])
  GR.in$CG_allele_extend_g2=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$genome2_tri)])
  GR.in$gff_size_extend=width(GR.extend)
  return(GR.in)
}
# Counting splitting events -----------------------------------------------
#Note Jordi splitted regions with large CpG numbers into 2 regions. 
#There's possibility that one region have SNP the other region does not have SNP but we can still separate allele 
#Because allele-separation happened before splitting regions
#The event is rare and only happened in small portion 0.07% of the regions
#They're treated as allele without CG difference in CpG analysis and without binding site difference in motif analysis 
#Per-sample data, highest is 0.56%:   Brain_substantia_nigra_paired - 112
as.data.table(mcols(GR_merge))[,sum(is.na(g1CG))/length(g1CG),by=list(Sample)]
GR_merge$g1CG[is.na(GR_merge$g1CG)]<-GR_merge$g2CG[is.na(GR_merge$g2CG)]<-GR_merge$refCG[is.na(GR_merge$refCG)]<-GR_merge$altCG[is.na(GR_merge$altCG)]<-0
GR_merge=GR_merge[GR_merge$N>=2&!(GR_merge$Sample%in%(c('rep1 - H1','rep2 - H1')))]
genomic_features=readRDS(genomic_features_file)
# Adding gene information to object ---------------------------------------
GR_merge=add_gene_GR(GR_merge,genomic_features$promoter,'genes_promoter')
GR_merge=add_gene_GR(GR_merge,genomic_features$`gene body`,'genes_body')
GR_merge=add_gene_GR(GR_merge,genomic_features$TSS,'TSS')
GR_merge$Sample[GR_merge$Sample=="merged - H1"] = "ESC - H1"
# loading the varibility file for each sample -----------------------------

agnostic_dir="../downstream/input/human_analysis/NME_expression_var/scRNA/"
GR_merge$hyper_var_fn=NA
GR_merge$tissue[GR_merge$tissue=='Adipose_Tissue_single']='Adipose_single'
GR_merge$hyper_var_fn[GR_merge$Sample %in% c('42_embryonic_stem_cell_single - H9','stem_27_undifferentiated_paired - HUES64',"ESC - H1")]=
  paste0(agnostic_dir,'HESC_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Adipose_single']=paste0(agnostic_dir,'AdultAdipose_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Bladder_single']=
  paste0(agnostic_dir,'AdultBladder_1.rds',';',agnostic_dir,'AdultBladder_2.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Small_Intestine_single']=paste0(agnostic_dir,'AdultIleum_2.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Gastric_single']=
  paste0(agnostic_dir,'AdultStomach_2.rds',';',agnostic_dir,'AdultStomach_1.rds',';',
                                  agnostic_dir,'AdultStomach_3.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Left_Ventricle_single']=
  paste0(agnostic_dir,'AdultHeart_1.rds',';',agnostic_dir,'AdultHeart_2.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Lung_single']=
  paste0(agnostic_dir,'AdultLung_1.rds',';',agnostic_dir,'AdultLung_2.rds',';',
                                  agnostic_dir,'AdultLung_3.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Psoas_Muscle_single']=paste0(agnostic_dir,'AdultMuscle_1.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Sigmoid_Colon_single']=paste0(agnostic_dir,'AdultColon_1.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Spleen_single']=paste0(agnostic_dir,'AdultSpleen_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Adrenal_Gland_single']=
  paste0(agnostic_dir,'AdultAdrenalGland_2.rds',';',agnostic_dir,'AdultAdrenalGland_3.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Aorta_single']=paste0(agnostic_dir,'AdultArtery_1.rds')


GR_merge$hyper_var_fn[GR_merge$tissue=='Esophagus_single']=
  paste0(agnostic_dir,'AdultEsophagus_1.rds',';',agnostic_dir,'AdultEsophagus_2.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Pancreas_single']=paste0(agnostic_dir,'AdultPancreas_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Liver_single']=
  paste0(agnostic_dir,'AdultLiver_1.rds',';',agnostic_dir,'AdultLiver_2.rds',';',
                                  agnostic_dir,'AdultLiver_4.rds')

saveRDS(GR_merge,GR_merge_file)

# add CpG information how many CG in each allele etc ----------------------

GR_merge_CpG=GRanges()
for (subj in subjects){GR_merge_CpG=c(GR_merge_CpG,hetCGallele_merged(subj,GR_merge,CpG_hg19,variant_HetCpG,gene_size=500))}
GR_merge_CpG$density=GR_merge_CpG$CG_hg19_extend/GR_merge_CpG$gff_size_extend
GR_merge_CpG$density_diff=(GR_merge_CpG$CG_allele_extend_g1-GR_merge_CpG$CG_allele_extend_g2)/
  GR_merge_CpG$gff_size_extend
saveRDS(GR_merge_CpG,GR_merge_file)