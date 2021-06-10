rm(list=ls())
source("mainFunctions_sub.R")

# get all hg19 CpG site ---------------------------------------------------
CpG_hg19=getCpgSitesH19()
saveRDS(CpG_hg19,'../downstream/input/human_analysis/CpG_hg19.rds')
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003",
           "STL011","H1","HuFGM02","112","149","150")
# reading in vcf files ----------------------------------------------------
#Linux for converting vcf to vcf.gz in ../downstream/data/vcfFiles/:
#for fn in *.vcf; do bgzip -c $fn > $fn.gz; done
variant_HetCpG=lapply(subjects,function(x) extractHetCpG('../downstream/data/vcfFiles/',x)) 
names(variant_HetCpG)=subjects
saveRDS(variant_HetCpG,variant_HetCpG_file)

# reading in stat for each sample -------------------------------
#Note here we're using coverage cutoff=5 and boundary check == true
GR=import.subject('../downstream/data/bedGraph_diff/')
saveRDS(GR,GR_file)
GR_allele=import.subject('../downstream/data/bedGraph_allele/',calc='allele')
saveRDS(GR_allele,GR_allele_file)

# reading in analyzed regions for each sample -----------------------------
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
saveRDS(gff_in,gff_in_file)

# Sanity check output regions not in original gff file --------------------
for(subj in subjects){
  cat(paste(subj,':\n'))
  cat(length(subsetByOverlaps(GR_allele[GR_allele$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-
        length(GR_allele[GR_allele$Subject==subj]),'\n')
  cat(length(subsetByOverlaps(GR[GR$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-length(GR[GR$Subject==subj]),'\n')
}
#Result all 0

# Extracting genomic features ---------------------------------------------

genomic_features=getGeneralFeats_CpG("../downstream/input/human_analysis/")
saveRDS(genomic_features,genomic_features_file)


# creating merged object --------------------------------------------------
GR_merge=GRanges()
for(sp in unique(GR$Sample)){
  cat("Processing",sp,'\n')
  ts=gsub('.*- ','',sp)
  GR_merge=c(GR_merge,stat_merge(GR[GR$Sample==sp],
                                 GR_allele[GR_allele$Sample==sp],
                                 variant_HetCpG[[ts]],CpG_hg19))
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

# Processing variant based result -----------------------------------------

variant_HetCpG_meta=fastDoCall('c',lapply(names(variant_HetCpG),variant_meta,variant_in=variant_HetCpG,GR_in=GR_merge))
#Trinucleotide analysis
#variant_HetCpG_meta$mask_tri=unlist(lapply(variant_HetCpG_meta$REF_tri,mask_tri))
saveRDS(variant_HetCpG_meta,variant_HetCpG_meta_file)


# reading in allele-agnostic analysis -------------------------------------
GR_merge=readRDS(GR_merge_file)
in_dir='../downstream/data/allele_agnostic_20kb/'
#all_regions=import.gff3('../downstream/output/human_20kb_allele_agnostic_250bp.gff')
#mcols(all_regions)=mcols(all_regions)[,c("N","score")]
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

# DNase and control ------------------------------------------------------------
GR_merge=readRDS(GR_merge_file)
in_dir='../downstream/data/agnostic_DNase/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir,pattern="[mn]m[le].bedGraph")){
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
saveRDS(NME_in,NME_agnostic_DNase_file)
saveRDS(MML_in,MML_agnostic_DNase_file)

# allele-specific regions with agnostic --------------------------------------------------------
in_dir='../downstream/data/allele_specific_region_agnostic/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir,pattern="[mn]m[le].bedGraph")){
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
saveRDS(NME_in[NME_in$N>=2],NME_agnostic_ASM_file)
saveRDS(MML_in[MML_in$N>=2],MML_agnostic_ASM_file)

# Allele-agnostic analysis for rest of regions --------------------------------------------------------
in_dir='../downstream/data/compliment_MML_NME_human/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir,pattern="[mn]m[le].bedGraph")){
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
#NME in check: 202721894
#MML in check: 202721894
saveRDS(NME_in[NME_in$N>=2],NME_agnostic_comp_file)
saveRDS(MML_in[MML_in$N>=2],MML_agnostic_comp_file)


saveRDS(c(readRDS(NME_agnostic_file),
          readRDS(NME_agnostic_DNase_file),
          readRDS(NME_agnostic_ASM_file),
          readRDS(NME_agnostic_comp_file)
                  ),NME_agnostic_all_file)
saveRDS(c(readRDS(MML_agnostic_file),
          readRDS(MML_agnostic_DNase_file),
          readRDS(MML_agnostic_ASM_file),
          readRDS(MML_agnostic_comp_file)
),MML_agnostic_all_file)
# reading in mouse MML and NME --------------------------------------------
#Complimentary regions
dir_comp='../downstream/data/compliment_MML_NME_model_mouse/'
MML_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*mml"),
                               read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
MML_in$MML=MML_in$score
MML_in$score=NULL
NME_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*nme"),read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
NME_in$NME=NME_in$score
NME_in$score=NULL
#Analyzed PRC, DNase and control
dir_analyzed='../downstream/data/DNase_control_PRC_MML_NME_model_mouse/'
MML_in_analyzed=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*mml"),
                               read.agnostic.mouse,in_dir=dir_analyzed,mc.cores=20))
MML_in_analyzed$MML=MML_in_analyzed$score
MML_in_analyzed$score=NULL
NME_in_analyzed=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*nme"),read.agnostic.mouse,in_dir=dir_analyzed,mc.cores=20))
NME_in_analyzed$NME=NME_in_analyzed$score
NME_in_analyzed$score=NULL
#Combine two datasets
NME_in=c(NME_in,NME_in_analyzed)
MML_in=c(MML_in,MML_in_analyzed)
#Convert to matrix: note here I didn't filter out the N<=17 since all NME and MML are intersect with UC regions in later analysis, the filtering in only done in UC
NME_in_matrix=agnostic_matrix_conversion(NME_in[NME_in$N>=2])
MML_in_matrix=agnostic_matrix_conversion(MML_in[MML_in$N>=2],'MML')

saveRDS(NME_in_matrix,NME_matrix_file)
saveRDS(MML_in_matrix,MML_matrix_file)
rm(MML_in)
rm(NME_in)
rm(MML_in_matrix)
rm(NME_in_matrix)
gc()
# reading in mouse UC --------------------------------------------
#Complimentary regions
UC_in_dir='../downstream/data/compliment_UC_non_MDS_mouse/'
UC_in=fastDoCall('c',mclapply(dir(UC_in_dir,pattern = '.*uc.bedGraph'),function(x){UC_in=read.agnostic.mouse.uc(paste(UC_in_dir,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=20))
UC_in$tissue=sub('-.*','',UC_in$Sample)
UC_in$Sample=sub('.5-.*-E1','.5-E1',UC_in$Sample)
#DNase,control,PRC
UC_in_dir_analyzed='../downstream/data/DNase_control_PRC_non_MDS_mouse/'
UC_in_analyzed=fastDoCall('c',mclapply(dir(UC_in_dir,pattern = paste0(paste(unique(UC_in$tissue),collapse='|'),'.*uc.bedGraph')),
                                       function(x){
  UC_in=read.agnostic.mouse.uc(paste(UC_in_dir_analyzed,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=20))
UC_in_analyzed$Sample=sub('.5-.*-E1','.5-E1',UC_in_analyzed$Sample)
UC_in_analyzed=UC_in_analyzed[UC_in_analyzed$Sample %in% unique(UC_in$Sample)]
#Merging data
UC_in=c(UC_in_analyzed[UC_in_analyzed$N>=2&UC_in_analyzed$N<=17],UC_in[UC_in$N>=2&UC_in$N<=17])
#Convert to matrix
UC_in_matrix_ls=mclapply(unique(UC_in$tissue),function(x) agnostic_matrix_conversion(UC_in[UC_in$tissue==x],'UC'),mc.cores=20)
names(UC_in_matrix_ls)=unique(UC_in$tissue)
saveRDS(UC_in_matrix_ls,UC_in_matrix_ls_file)


# UC for mouse MDS comparison ---------------------------------------------
gff_in_compliment=import.gff3(mouse_compliment_gff_file)
gff_in_compliment=paste0(seqnames(gff_in_compliment),':',start(gff_in_compliment),'-',end(gff_in_compliment))
UC_in_MDS_comp=data.table(region=gff_in_compliment)
compliment_MDS_dir='../downstream/data/compliment_UC_MDS_mouse/'
UC_in_MDS_comp_UC=fastDoCall('cbind',
                             mclapply(dir(compliment_MDS_dir,pattern = '.*uc.bedGraph'),function(x){
                               read.agnostic.mouse.uc(paste(compliment_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_compliment)},mc.cores=10))

UC_in_MDS_comp=cbind(UC_in_MDS_comp,UC_in_MDS_comp_UC)
saveRDS(UC_in_MDS_comp,'../downstream/output/UC_in_MDS_comp.rds')#This is temporary, delete afterwards
#Filter based on N first to save space
DNase_conrol_MDS_dir='../downstream/data/DNase_control_PRC_MDS_mouse/'
gff_in_DNase=import.gff3(mouse_DNase_control_gff_file)
gff_in_DNase=paste0(seqnames(gff_in_DNase),':',start(gff_in_DNase),'-',end(gff_in_DNase))
UC_in_analyzed_MDS=data.table(region=gff_in_DNase)
UC_in_analyzed_MDS_UC=fastDoCall('cbind',
                                 mclapply(dir(compliment_MDS_dir,pattern = '.*uc.bedGraph'),function(x){
                                   read.agnostic.mouse.uc(paste(DNase_conrol_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_DNase)},mc.cores=10))
UC_in_analyzed_MDS=cbind(UC_in_analyzed_MDS,UC_in_analyzed_MDS_UC)
saveRDS(UC_in_analyzed_MDS,'../downstream/output/UC_in_analyzed_MDS.rds')#This is temporary, delete afterwards
UC_in_MDS_all=rbind(UC_in_MDS_comp,UC_in_analyzed_MDS)
saveRDS(UC_in_MDS_all,UC_in_MDS_all_file)


#Read in mouse NME and scRNA
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

