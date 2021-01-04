rm(list=ls())
source("mainFunctions_sub.R")
##############################Updated CPEL pipeline modify it accordingly###############
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003","STL011","H1","HuFGM02","112","149","150")
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
#Reading in regions to define the region analyzed, check if all x,y,MT are excluded
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
saveRDS(gff_in,gff_in_file)
 for(subj in subjects){
  cat(paste(subj,':\n'))
  cat(length(subsetByOverlaps(GR_allele[GR_allele$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-
        length(GR_allele[GR_allele$Subject==subj]),'\n')
  cat(length(subsetByOverlaps(GR[GR$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-length(GR[GR$Subject==subj]),'\n')
}
#Extracting genomic features
genomic_features=getGeneralFeats_CpG("../downstream/input/")
saveRDS(genomic_features,genomic_features_file)

#create merged object
GR_merge=GRanges()
for(sp in unique(GR$Sample)){
  cat("Processing",sp,'\n')
  ts=gsub('.*- ','',sp)
  GR_merge=c(GR_merge,stat_merge(GR[GR$Sample==sp],
                                 GR_allele[GR_allele$Sample==sp],
                                 variant_HetCpG[[ts]],CpG_hg19))
}
#Note Jordi splitted regions with large CpG numbers into 2 regions. 
#There's possibility that one region have SNP the other region does not have SNP but we can still separate allele 
#Because allele-separation happened before splitting regions
#The event is rare and only happened in small portion 0.07% of the regions
#They're treated as allele without CG difference in CpG analysis and without binding site difference in motif analysis 
#Per-sample data, highest is 0.56%:   Brain_substantia_nigra_paired - 112
as.data.table(mcols(GR_merge))[,sum(is.na(g1CG))/length(g1CG),by=list(Sample)]
GR_merge$g1CG[is.na(GR_merge$g1CG)]<-GR_merge$g2CG[is.na(GR_merge$g2CG)]<-GR_merge$ref1CG[is.na(GR_merge$refCG)]<-GR_merge$altCG[is.na(GR_merge$altCG)]<-0
GR_merge=GR_merge[GR_merge$N>=2]
GR_merge=GR_merge[!GR_merge$Sample%in%c("rep1 - H1","rep2 - H1")]
#get rid of not unique direction, 10%, 2059 genes, closest to middle point, hypervaribility use mean.
GR_merge=add_gene_GR(GR_merge,genomic_features$promoter,'genes_promoter')
GR_merge=add_gene_GR(GR_merge,genomic_features$`gene body`,'genes_body')
GR_merge=add_gene_GR(GR_merge,genomic_features$TSS,'TSS')
GR_merge$Sample[GR_merge$Sample=="merged - H1"] = "ESC - H1"
# loading the varibility data for each sample -----------------------------
# regions with high NME at one allele is more likely to have hyper vari or with preferential motif binding
agnostic_dir="../downstream/input/scRNA/"
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
##Merged object with CpG information############
#GR_merge=readRDS(GR_merge_file)
GR_merge_CpG=GRanges()
for (subj in subjects){GR_merge_CpG=c(GR_merge_CpG,hetCGallele_merged(subj,GR_merge,CpG_hg19,variant_HetCpG,gene_size=500))}
GR_merge_CpG$density=GR_merge_CpG$CG_hg19_extend/GR_merge_CpG$gff_size_extend
GR_merge_CpG$density_diff=(GR_merge_CpG$CG_allele_extend_g1-GR_merge_CpG$CG_allele_extend_g2)/
  GR_merge_CpG$gff_size_extend
saveRDS(GR_merge_CpG,GR_merge_file)
#processing variant analysis
variant_HetCpG_meta=fastDoCall('c',lapply(names(variant_HetCpG),variant_meta,variant_in=variant_HetCpG,GR_in=GR_merge))
#Trinucleotide analysis
#variant_HetCpG_meta$mask_tri=unlist(lapply(variant_HetCpG_meta$REF_tri,mask_tri))
saveRDS(variant_HetCpG_meta,variant_HetCpG_meta_file)

#human allele-agnostic analysis is from +/- 20k of TSS
GR_merge=readRDS(GR_merge_file)
in_dir='../downstream/data/allele_agnostic_20kb/'
all_regions=import.gff3('../downstream/output/human_20kb_allele_agnostic_250bp.gff')
mcols(all_regions)=mcols(all_regions)[,c("N","score")]
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir)){
  cat('Reading in',fn,'\n')
  stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
  sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
  subject_in=sub('_.*','',sample_in)
  tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
  sample_in=paste0(tissue_in,' - ',subject_in)
  if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
  fn_in=unique(GR_merge[GR_merge$Sample==sample_in]$hyper_var_fn)
if(length(fn_in)==0){fn_in=NA}

  if(stat_in=="NME"){
  #Remove overlapped regions
  NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}
  else if(stat_in=="MML"){
  MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,allele_include=F,
                                sample_in=sample_in,hyper_var_file=fn_in))}else
   {cat("Error stat_in:", stat_in,'\n')}
}

NME_in$NME=NME_in$score
NME_in=NME_in[NME_in$N>=2]
NME_in=NME_in[!(NME_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
#number of regions check:77457572
saveRDS(NME_in,NME_agnostic_file)
MML_in$MML=MML_in$score
MML_in=MML_in[!(MML_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
MML_in=MML_in[MML_in$N>=2]
#number check: 77460740
saveRDS(MML_in,MML_agnostic_file)


# add hypervar to allele-agnostic data ------------------------------------

GR_merge=readRDS(GR_merge_file)
NME_hypervar_calc=GRanges()
NME_in=readRDS(NME_agnostic_file)
MML_hypervar_calc=GRanges()
MML_in=readRDS(MML_agnostic_file)

GR_calc=data.frame()
scRNA_result=data.frame()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){

    sp_hyper_var=read_hypervar(hyper_var_file)
    #scRNA_result=rbind(scRNA_result,sp_hyper_var)
    #Add hypervaribility inforamtion
    #Add hypervaribility inforamtion
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))

    
  }else{cat("file not exist for:",sp,'\n')}
}

#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),'../downstream/output/allele_agnostic_var_homogeneous2.rds')
# # human_FANTOM ------------------------------------------------------------
# GR_merge=readRDS('GR_merge_final12_ls.rds')
# in_dir='../allele_agnostic_hg19_cov5_3kb_FANTOM/'
# NME_in=GRanges()
# MML_in=GRanges()
# for(fn in  dir(in_dir)){
#   cat('Reading in',fn,'\n')
#   stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
#   sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
#   subject_in=sub('_.*','',sample_in)
#   tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
#   sample_in=paste0(tissue_in,' - ',subject_in)
#   
#   if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
#   if(stat_in=="NME"){
#     NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$Sample==sample_in]))}
#   else if(stat_in=="MML"){
#     MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$Sample==sample_in]))}else
#     {cat("Error stat_in:", stat_in,'\n')}
# }
# 
# NME_in$NME=NME_in$score
# MML_in$MML=MML_in$score
# NME_in=NME_in[!(NME_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
# MML_in=MML_in[!(MML_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
# NME_in=NME_in[NME_in$N>=2]
# MML_in=MML_in[MML_in$N>=2]
# saveRDS(NME_in,"allele_agnostic_hg19_cov10_3kb_FANTOM_NME.rds")
# saveRDS(MML_in,"allele_agnostic_hg19_cov10_3kb_FANTOM_MML.rds")
# NME_in$hyper_var_fn[NME_in$Sample %in% c('42_embryonic_stem_cell_single - H9','stem_27_undifferentiated_paired - HUES64',"ESC - H1")]=
#   paste('../downstream/input/scRNA/HESC_1.rds',sep = ';')

# DNase and control ------------------------------------------------------------
GR_merge=readRDS(GR_merge_file)
in_dir='../downstream/data/agnostic_DNase/'
all_regions=granges(import.gff3('../downstream/output/human_DNase_allele_agnostic_250bp.gff'))
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir,pattern="bedGraph")){
  cat('Reading in',fn,'\n')
  stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
  sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
  subject_in=sub('_.*','',sample_in)
  tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
  sample_in=paste0(tissue_in,' - ',subject_in)
  if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
  
  fn_in=unique(GR_merge[GR_merge$Sample==sample_in]$hyper_var_fn)
  if(length(fn_in)==0){fn_in=NA}
  if(stat_in=="NME"&sum(GR_merge$Sample==sample_in)>0){
    NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                  allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}
  else if(stat_in=="MML"&sum(GR_merge$Sample==sample_in)>0){
    MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}else
    {cat("Error stat_in:", stat_in,'in sample:',sample_in,'\n')}
}

NME_in$NME=NME_in$score
MML_in$MML=MML_in$score
NME_in=NME_in[!(NME_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
MML_in=MML_in[!(MML_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
NME_in=NME_in[NME_in$N>=2]
MML_in=MML_in[MML_in$N>=2]
#check: NME 57678420
#check: MML 57681712
saveRDS(NME_in,"../downstream/output/allele_agnostic_hg19_DNase_NME_homogeneous_excluding_dMML2.rds")
saveRDS(MML_in,"../downstream/output/allele_agnostic_hg19_DNase_MML_homogeneous2.rds")

# allele-specific regions with agnostic --------------------------------------------------------
in_dir='../downstream/data/allele_specific_region_agnostic/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir,pattern="bedGraph")){
  cat('Reading in',fn,'\n')
  stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
  sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
  subject_in=sub('_.*','',sample_in)
  tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
  sample_in=paste0(tissue_in,' - ',subject_in)

  if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
  

  if(stat_in=="NME"){
    NME_in_sp=read.bedGraph.informME(paste0(in_dir,fn))
    NME_in_sp$Sample=sample_in
    NME_in_sp$statistics=stat_in
    NME_in=c(NME_in,NME_in_sp)}
  else if(stat_in=="MML"){
    MML_in_sp=read.bedGraph.informME(paste0(in_dir,fn))
    MML_in_sp$Sample=sample_in
    MML_in_sp$statistics=stat_in
    MML_in=c(MML_in,MML_in_sp)}else
    {cat("Error stat_in:", stat_in,'\n')}
}
rm(NME_in_sp)
rm(MML_in_sp)
#NME in check: 51459969
#MML in check: 51459969
saveRDS(NME_in,'../downstream/output/NME_agnostic_ASM.rds')
saveRDS(MML_in,'../downstream/output/MML_agnostic_ASM.rds')
NME_in=readRDS('../downstream/output/NME_agnostic_ASM.rds')
MML_in=readRDS('../downstream/output/MML_agnostic_ASM.rds')
NME_in_tb=as.data.table(mcols(NME_in))
NME_in_tb$region=paste0(seqnames(NME_in),':',start(NME_in),'-',end(NME_in))
MML_in_tb=as.data.table(mcols(MML_in))
MML_in_tb$region=paste0(seqnames(MML_in),':',start(MML_in),'-',end(MML_in))
GR_merge_tb=rbind(NME_in_tb,MML_in_tb)
rm(NME_in_tb)
rm(MML_in_tb)
rm(NME_in)
rm(MML_in)
GR_merge_tb$K=NULL
GR_merge_tb$N=NULL
GR_merge=readRDS(GR_merge_file)
GR_merge_tb_asm=rbind(
  data.table(score=GR_merge$dNME,Sample=GR_merge$Sample,statistics='dNME',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dNME_pval,Sample=GR_merge$Sample,statistics='dNME_pval',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dMML,Sample=GR_merge$Sample,statistics='dMML',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=GR_merge$dMML_pval,Sample=GR_merge$Sample,statistics='dMML_pval',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=(GR_merge$MML1+GR_merge$MML2)/2,Sample=GR_merge$Sample,statistics='MML_ASM',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge))),
  data.table(score=(GR_merge$NME1+GR_merge$NME2)/2,Sample=GR_merge$Sample,statistics='NME_ASM',region=paste0(seqnames(GR_merge),':',start(GR_merge),'-',end(GR_merge)))
)
rm(GR_merge)

GR_merge_tb=rbind(GR_merge_tb,GR_merge_tb_asm)
GR_merge_tb=dcast.data.table(GR_merge_tb,Sample+region~statistics,value.var = "score")
saveRDS(GR_merge_tb,'../downstream/output/GR_merge_ASM_comp.rds')
GR_merge_tb=readRDS('../downstream/output/GR_merge_ASM_comp.rds')

#Read in mouse NME
NME_in=readRDS('../downstream/input/NME_agnostic_mouse_all_merged.rds')
dir='../downstream/data/Mouse_C1/'
NME_in=NME_in[NME_in$tissue=="limb"&NME_in$N>=2]
gtf <- fread('../downstream/input/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
genes <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
genes$gene_name <- gn
NME_in=dist_calc(NME_in,genes)
#Percent gene covered?
length(unique(NME_in[abs(NME_in$dist)<=3000]$gene))/length(genes[seqnames(genes)!="chrM"])#85%
NME_in_dt=as.data.table(mcols(NME_in))
NME_in_dt$region=paste0(seqnames(NME_in),':',start(NME_in),'-',end(NME_in))
NME_in_dt$hyper_var=-100
NME_in_dt$var=-100
NME_in_dt$mean=-100
for(st in unique(NME_in$stage)){
  tt1=proc.time()[[3]]
  if(file.exists(paste0(dir,sub('E','',st),'.rds'))){
    scRNA_in=readRDS(paste0(dir,sub('E','',st),'.rds'))
    scRNA_in=scRNA_in[rownames(scRNA_in)%in% unique(c(NME_in_dt[(stage==st)]$gene)),]
    if(nrow(scRNA_in)>0){
      #Add hypervar to TSS 
      NME_in_dt[(stage==st)]$hyper_var=scRNA_in[NME_in_dt[(stage==st)]$gene,"hypervar_logvar"]
      NME_in_dt[(stage==st)]$var=scRNA_in[NME_in_dt[(stage==st)]$gene,"var"]
      NME_in_dt[(stage==st)]$mean=scRNA_in[NME_in_dt[(stage==st)]$gene,"mean"]
      
    }
  }else{cat("File not exist for ",st,'\n')}
  cat('Finish processing ',sub('E','',st),'in: ',proc.time()[[3]]-tt1,'\n')
  
}
saveRDS(NME_in_dt,'../downstream/output/NME_in_limb_ENOCD3_imputed.rds')

# Read in mouse enhancer --------------------------------------------------

rm(list=ls())
library(data.table)
library(rtracklayer)
library(Gmisc)
#Read in FeDMR 
dir_in='../downstream/input/FeDMR/'
FeDMR_fn=dir(dir_in,pattern='tsv')
FeDMR_in=fastDoCall('c',lapply(dir(dir_in,pattern='tsv'),function(fn){
  sp=gsub('.tsv','',fn)
  sp=gsub("feDMR_","",sp)
  tissue=gsub(".*_","",sp)
  stage=gsub(paste0("_",tissue),"",sp)
  stage=gsub("_5",".5",stage)
  tissue=gsub("craniofacial","EFP",tissue)
  tissue=gsub("tube","NT",tissue)
  fn_in=fread(paste0(dir_in,fn))
  fn_in=makeGRangesFromDataFrame(fn_in,seqnames.field = "chrom",keep.extra.columns = T)
  fn_in$tissue=tissue
  fn_in$stage=stage
  fn_in$sample=paste(tissue,stage,sep='-')
  return(fn_in)
  
})
)
FeDMR_in_mcols=as.data.table(mcols(FeDMR_in))
FeDMR_in_mcols$tissue_exist=TRUE
FeDMR_in_gr=unique(FeDMR_in)
FeDMR_in_gr=FeDMR_in_gr[order(FeDMR_in_gr$dmr_id,decreasing=T)]
mcols(FeDMR_in_gr)=mcols(FeDMR_in_gr)[,c("dmr_id","score","tissue_specificity","stage_specificity")]
FeDMR_in_gr$dmr_id_original=FeDMR_in_gr$dmr_id
FeDMR_in_gr$dmr_id=NULL
FeDMR_in_mcols_dc=dcast.data.table(FeDMR_in_mcols,dmr_id~sample,value.var = "tissue_exist",fill=FALSE)
FeDMR_in_mcols_dc=FeDMR_in_mcols_dc[order(FeDMR_in_mcols_dc$dmr_id,decreasing = T)]
mcols(FeDMR_in_gr)=cbind(mcols(FeDMR_in_gr),FeDMR_in_mcols_dc)
print(which(FeDMR_in_gr$dmr_id_original!=FeDMR_in_gr$dmr_id))
FeDMR_in_gr$dmr_id_original=NULL
saveRDS(FeDMR_in_gr,'../downstream/output/FeDMR.rds')
# percent_cov=NULL
# UC_in=readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')
# for(tissue in unique(FeDMR_in_mcols$tissue)){
# percent_cov=rbind(percent_cov,data.frame(tissue=tissue,
#   percent=length(subsetByOverlaps(FeDMR_in_gr[apply(mcols(FeDMR_in_gr)[,gsub('-.*','',colnames(mcols(FeDMR_in_gr)))==tissue],1,any)],UC_in[[tissue]]))/
#     length(FeDMR_in_gr[apply(mcols(FeDMR_in_gr)[,gsub('-.*','',colnames(mcols(FeDMR_in_gr)))==tissue],1,any)])
#   ))
# }

#Read in chromHMM

#ChromHMM

chromHMM_pooled= read_chromHMM_bed('../downstream/input/chromHMM_mm10/pooled/','pooled')
saveRDS(chromHMM_pooled,'../downstream/output/chromHMM.rds')
chromHMM_enhancer=chromHMM_pooled[sub('-.*','',as.character(chromHMM_pooled$chrom_state))=="En"]
saveRDS(chromHMM_enhancer,'../downstream/output/chromHMM_enhancer.rds')

