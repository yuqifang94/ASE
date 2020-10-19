rm(list=ls())
source("mainFunctions_sub.R")
##############################Updated CPEL pipeline modify it accordingly###############
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003","STL011","H1","HuFGM02","112","149","150")
#Reading in regions to define the region analyzed, check if all x,y,MT are excluded
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
saveRDS(gff_in,gff_in_file)
#read in vcf files
#Fix the issue with genome file order does not equal to the ref/alt order, calculate the CG or het CG in genome 1 based on actual ref/alt order
variant_HetCpG=lapply(subjects,function(x) extractHetCpG('../downstream/data/vcfFiles/',x)) 
names(variant_HetCpG)=subjects
saveRDS(variant_HetCpG,variant_HetCpG_file)
#Note the first line cannot be NA
GR=import.subject('../downstream/data/bedGraph_diff/')
saveRDS(GR,GR_file)
GR_allele=import.subject('../downstream/data/bedGraph_allele/',calc='allele')
saveRDS(GR_allele,GR_allele_file)
###Checking if GR_allele match gff matched###########################
for(subj in subjects){
  cat(paste(subj,':\n'))
  cat(length(subsetByOverlaps(GR_allele[GR_allele$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-
        length(GR_allele[GR_allele$Subject==subj]),'\n')
  cat(length(subsetByOverlaps(GR[GR$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-length(GR[GR$Subject==subj]),'\n')
}
#Extracting genomic features
genomic_features=getGeneralFeats_CpG("../downstream/input/")
saveRDS(genomic_features,genomic_features_file)
#gff with HetCpG information
hetCpG_gff=lapply(subjects,gff_hetCpG_count,gff_in=gff_in,vcf_in=variant_HetCpG,CpG=CpG_hg19)
#Checking if agree with gff file
names(hetCpG_gff)=subjects
saveRDS(hetCpG_gff,hetCpG_gff_file)
#create merged object
GR_merge=GRanges()
for(sp in unique(GR$Sample)){
  cat("Processing",sp,'\n')
  GR_merge=c(GR_merge,stat_merge(GR[GR$Sample==sp],GR_allele[GR_allele$Sample==sp]))
}

GR_merge=GR_merge[GR_merge$N>=2]
#get rid of not unique direction, 10%, 2059 genes, closest to middle point, hypervaribility use mean.
GR_merge=add_gene_GR(GR_merge,genomic_features$promoter,'genes_promoter')
GR_merge=add_gene_GR(GR_merge,genomic_features$`gene body`,'genes_body')
GR_merge=add_gene_GR(GR_merge,genomic_features$TSS,'TSS')
GR_merge$Sample[GR_merge$Sample=="merged - H1"] = "ESC - H1"
# loading the varibility data for each sample -----------------------------
# regions with high NME at one allele is more likely to have hyper vari or with preferential motif binding
GR_merge$hyper_var_fn=NA
GR_merge$tissue[GR_merge$tissue=='Adipose_Tissue_single']='Adipose_single'
GR_merge$hyper_var_fn[GR_merge$tissue %in% c('42_embryonic_stem_cell_single - H9','stem_27_undifferentiated_paired - HUES64',"ESC - H1")]=
  paste('../downstream/input/scRNA/HESC_1.rds',sep = ';')

GR_merge$hyper_var_fn[GR_merge$tissue=='Adipose_single']=paste('../downstream/input/scRNA/AdultAdipose_1.rds',sep = ';')
#../downstream/input/scRNA/AdultBladder_2.rds 
GR_merge$hyper_var_fn[GR_merge$tissue=='Bladder_single']=
  paste('../downstream/input/scRNA/AdultBladder_1.rds','../downstream/input/scRNA/AdultBladder_2.rds',sep = ';')
GR_merge$hyper_var_fn[GR_merge$tissue=='Small_Intestine_single']=paste('../downstream/input/scRNA/AdultIleum_2.rds',sep = ';')
#../downstream/input/scRNA/AdultStomach_3.rds 
GR_merge$hyper_var_fn[GR_merge$tissue=='Gastric_single']=
  paste('../downstream/input/scRNA/AdultStomach_2.rds','../downstream/input/scRNA/AdultStomach_1.rds',
                                  '../downstream/input/scRNA/AdultStomach_3.rds',sep = ';')
#../downstream/input/scRNA/AdultHeart_2.rds
GR_merge$hyper_var_fn[GR_merge$tissue=='Left_Ventricle_single']=
  paste('../downstream/input/scRNA/AdultHeart_1.rds','../downstream/input/scRNA/AdultHeart_2.rds',sep = ';')
#../downstream/input/scRNA/AdultLung_3.rds 
GR_merge$hyper_var_fn[GR_merge$tissue=='Lung_single']=
  paste('../downstream/input/scRNA/AdultLung_1.rds','../downstream/input/scRNA/AdultLung_2.rds',
                                  '../downstream/input/scRNA/AdultLung_3.rds',sep = ';')

GR_merge$hyper_var_fn[GR_merge$tissue=='Psoas_Muscle_single']=paste('../downstream/input/scRNA/AdultMuscle_1.rds',sep = ';')
GR_merge$hyper_var_fn[GR_merge$tissue=='Sigmoid_Colon_single']=paste('../downstream/input/scRNA/AdultColon_1.rds',sep = ';')
GR_merge$hyper_var_fn[GR_merge$tissue=='Spleen_single']=paste('../downstream/input/scRNA/AdultSpleen_1.rds',sep = ';')

GR_merge$hyper_var_fn[GR_merge$tissue=='Adrenal_Gland_single']=
  paste('../downstream/input/scRNA/AdultAdrenalGland_2.rds','../downstream/input/scRNA/AdultAdrenalGland_3.rds',sep = ';')

GR_merge$hyper_var_fn[GR_merge$tissue=='Aorta_single']=paste('../downstream/input/scRNA/AdultArtery_1.rds',sep = ';')


GR_merge$hyper_var_fn[GR_merge$tissue=='Esophagus_single']=
  paste('../downstream/input/scRNA/AdultEsophagus_1.rds','../downstream/input/scRNA/AdultEsophagus_2.rds',sep = ';')
GR_merge$hyper_var_fn[GR_merge$tissue=='Pancreas_single']=paste('../downstream/input/scRNA/AdultPancreas_1.rds',sep = ';')

GR_merge$hyper_var_fn[GR_merge$tissue=='Liver_single']=
  paste('../downstream/input/scRNA/AdultLiver_1.rds','../downstream/input/scRNA/AdultLiver_2.rds',
                                  '../downstream/input/scRNA/AdultLiver_4.rds',sep = ';')

saveRDS(GR_merge,GR_merge_file)
##Merged object with CpG information############
#GR_merge=readRDS(GR_merge_file)
GR_merge_CpG=GRanges()
for (subj in subjects){GR_merge_CpG=c(GR_merge_CpG,hetCGallele_merged(subj,GR_merge,hetCpG_gff,CpG_hg19,variant_HetCpG,gene_size=500))}
GR_merge_CpG$density=GR_merge_CpG$CG_hg19_extend/GR_merge_CpG$gff_size_extend
GR_merge_CpG$density_diff=(GR_merge_CpG$CG_allele_extend_g1-GR_merge_CpG$CG_allele_extend_g2)/
  GR_merge_CpG$gff_size_extend
saveRDS(GR_merge_CpG,GR_merge_file)
##############################Variant analysis: Counting enrichment of each variant from all variant ##############################################
variant_HetCpG_meta=do.call('c',lapply(names(variant_HetCpG),variant_meta,variant_in=variant_HetCpG,GR_in=GR_merge))
#Trinucleotide analysis
#variant_HetCpG_meta$mask_tri=unlist(lapply(variant_HetCpG_meta$REF_tri,mask_tri))
saveRDS(variant_HetCpG_meta,variant_HetCpG_meta_file)

in_dir='../downstream/data/allele_agnostic_2kb/'
NME_hypervar_calc=GRanges()
MML_hypervar_calc=GRanges()
NME_meanvar_calc=GRanges()
MML_meanvar_calc=GRanges()
GR_calc=data.frame()
GR_merge$fn=paste(GR_merge$Subject,GR_merge$tissue,sep='_')
scRNA_result=data.frame()
#For the mean expression: log2 or not?
for (fn in unique(GR_merge$fn)){
  
  fn_mml=paste(in_dir,fn,'_phased_allele_agnostic_mml.bedGraph',sep='')
  fn_nme=paste(in_dir,fn,'_phased_allele_agnostic_nme.bedGraph',sep='')
  hyper_var_file=unlist(strsplit(unique(GR_merge$hyper_var_fn[GR_merge$fn==fn]),';'))
  cat('Processing',fn,'\n')
  if(file.exists(fn_mml)&file.exists(fn_nme)&all(file.exists(hyper_var_file))){
    cat('reading in:',fn_mml,'\n')
    cat('reading in:',fn_nme,'\n')
    sp_hyper_var=read_hypervar(hyper_var_file)
    #scRNA_result=rbind(scRNA_result,sp_hyper_var)
    sp_hyper_var$log2mean=log2(sp_hyper_var$mean)
    #Add hypervaribility inforamtion
    #Add hypervaribility inforamtion
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(fn_nme,sp_hyper_var,
                                                         genomic_features,GR_merge[GR_merge$fn==fn],exp_stat=quote(hypervar_logvar)))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(fn_mml,sp_hyper_var,
                                                         genomic_features,GR_merge[GR_merge$fn==fn],exp_stat=quote(hypervar_logvar)))
    MML_meanvar_calc=c(MML_meanvar_calc,dist_plot_calc(fn_mml,sp_hyper_var,
                                                       genomic_features,GR_merge[GR_merge$fn==fn],exp_stat=quote(log2mean)))
    NME_meanvar_calc=c(NME_meanvar_calc,dist_plot_calc(fn_nme,sp_hyper_var,
                                                       genomic_features,GR_merge[GR_merge$fn==fn],exp_stat=quote(log2mean)))
    
  }else{cat("file not exist for:",fn,'\n')}
}


saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc,
             MML_meanvar_calc=MML_meanvar_calc,
             NME_meanvar_calc=NME_meanvar_calc),'../downstream/output/allele_agnostic_var.rds')
#human allele-agnostic analysis is from +/- 20k of TSS
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir)){
  cat('Reading in',fn,'\n')
  stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
  sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
  subject_in=sub('_.*','',sample_in)
  tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
  sample_in=paste0(tissue_in,' - ',subject_in)

  if(sample_in=="merged - H1"){sample_in="ESC - H1"}
  if(stat_in=="NME"){
  NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$Sample==sample_in]))}
  else if(stat_in=="MML"){
  MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$Sample==sample_in]))}else
  {cat("Error stat_in:", stat_in,'\n')}
}

NME_in$NME=NME_in$score
MML_in$MML=MML_in$score
NME_in=NME_in[!(NME_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
MML_in=MML_in[!(MML_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
NME_in=NME_in[NME_in$N>=2]
MML_in=MML_in[MML_in$N>=2]
saveRDS(NME_in,NME_agnostic_file)
saveRDS(MML_in,MML_agnostic_file)


