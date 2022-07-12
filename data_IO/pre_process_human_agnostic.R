source('mainFunctions_sub.R')
# reading in allele-agnostic analysis -------------------------------------
GR_merge=readRDS(GR_merge_file)
in_dir_20kb='../downstream/data/allele_agnostic_20kb/'
#all_regions=import.gff3('../downstream/output/human_20kb_allele_agnostic_250bp.gff')
#mcols(all_regions)=mcols(all_regions)[,c("N","score")]
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir_20kb)){
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
  NME_in=c(NME_in,read.agnostic(paste0(in_dir_20kb,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}
  else if(stat_in=="MML"){
  MML_in=c(MML_in,read.agnostic(paste0(in_dir_20kb,fn),GR_merge_in=NULL,allele_include=F,
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
in_dir_DNase='../downstream/data/agnostic_DNase/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir_DNase,pattern="[mn]m[le].bedGraph")){
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
    NME_in=c(NME_in,read.agnostic(paste0(in_dir_DNase,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                  allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}
  else if(stat_in=="MML"&sum(GR_merge$Sample==sample_in)>0){
    MML_in=c(MML_in,read.agnostic(paste0(in_dir_DNase,fn),GR_merge_in=NULL,allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}else
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
in_dir_ASM='../downstream/data/allele_specific_region_agnostic/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir_ASM,pattern="[mn]m[le].bedGraph")){
  cat('Reading in',fn,'\n')
  stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
  sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
  subject_in=sub('_.*','',sample_in)
  tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
  sample_in=paste0(tissue_in,' - ',subject_in)

  if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
  

  if(stat_in=="NME"){
    NME_in_sp=read.bedGraph.informME(paste0(in_dir_ASM,fn))
    NME_in_sp$Sample=sample_in
    NME_in_sp$statistics=stat_in
    NME_in=c(NME_in,NME_in_sp)}
  else if(stat_in=="MML"){
    MML_in_sp=read.bedGraph.informME(paste0(in_dir_ASM,fn))
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
in_dir_compliment='../downstream/data/compliment_MML_NME_human/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir_compliment,pattern="[mn]m[le].bedGraph")){
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
  NME_in=c(NME_in,read.agnostic(paste0(in_dir_compliment,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))
 }
  else if(stat_in=="MML"){
  MML_in=c(MML_in,read.agnostic(paste0(in_dir_compliment,fn),GR_merge_in=NULL,allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}else
    {cat("Error stat_in:", stat_in,'\n')}
}

#NME in check: 202721894
#MML in check: 202721894
saveRDS(NME_in[NME_in$N>=2],NME_agnostic_comp_file)
saveRDS(MML_in[MML_in$N>=2],MML_agnostic_comp_file)


saveRDS(c(readRDS(NME_agnostic_file),
          readRDS(NME_agnostic_DNase_file),
          #readRDS(NME_agnostic_ASM_file),
          readRDS(NME_agnostic_comp_file)
                  ),NME_agnostic_all_file)
saveRDS(c(readRDS(MML_agnostic_file),
          readRDS(MML_agnostic_DNase_file),
          #readRDS(MML_agnostic_ASM_file),
          readRDS(MML_agnostic_comp_file)
),MML_agnostic_all_file)

#Check each file for size: difference due to removing ASM regions
NME_agnostic=readRDS(NME_agnostic_file)#77457572
MML_agnostic=readRDS(MML_agnostic_file)#77460740
NME_DNase=readRDS(NME_agnostic_DNase_file)#57678420
MML_DNase=readRDS(MML_agnostic_DNase_file)#57681712
NME_compliment=readRDS(NME_agnostic_comp_file)#122068847
MML_compliment=readRDS(MML_agnostic_comp_file)#122068847