rm(list=ls())
source("mainFunctions_sub.R")
#Read in allele-agnositc model
read.agnostic<-function(file_in,GR_merge_in=NULL,allele_include=T,olap_type="any",all_regions=NA,sample_in=NA,hyper_var_file=NA){
  stat=toupper(sub('.*_','',sub('.bedGraph','',file_in)))
  informME_in=read.bedGraph.informME(file_in)
  if(length(GR_merge_in)>0){
    GR_merge_in=GR_merge_in[GR_merge_in$Sample==sample_in]
  #Find overlapped region
    olap=findOverlaps(informME_in,GR_merge_in,type=olap_type)
    if(length(olap)>0) {informME_in=informME_in[-queryHits(olap)]}
  }
  #add GR_merge data
  if(allele_include){
  
    cat('Percent overlap with dNME region:',length(unique(queryHits(olap)))/length(informME_in)*100,'%\n')
    informME_in$score_original=informME_in$score
 
    #replace value instead of remove regions
    olap=findOverlaps(all_regions,GR_merge_in)
    asm_replacement=data.table(idx=queryHits(olap), 
                               score= rowMeans(as.matrix(elementMetadata(GR_merge_in)[paste(stat,c('1','2'),sep='')]))[subjectHits(olap)],
                               dMML=GR_merge_in$dMML[subjectHits(olap)],
                               dMML_pval=GR_merge_in$dMML_pval[subjectHits(olap)],
                               N=GR_merge_in$N[subjectHits(olap)])
    asm_replacement=asm_replacement[,list(score=mean(score),dMML=mean(dMML),dMML_pval=mean(dMML_pval)),by=list(idx)]
    all_regions=all_regions[asm_replacement$idx]
    all_regions$score=asm_replacement$score
    all_regions$K=NA
    all_regions$dMML=asm_replacement$dMML
    all_regions$dMML_pval=asm_replacement$dMML_pval
    informME_in=c(informME_in,all_regions)
  }
 informME_in=informME_in[!is.infinite(informME_in$score)]
  informME_in$Sample=sample_in
  informME_in$hyper_var_fn=hyper_var_file
  return(informME_in)

}
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
  NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))
 }
  else if(stat_in=="MML"){
  MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}else
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
          #readRDS(NME_agnostic_ASM_file),
          readRDS(NME_agnostic_comp_file)
                  ),NME_agnostic_all_file)
saveRDS(c(readRDS(MML_agnostic_file),
          readRDS(MML_agnostic_DNase_file),
          #readRDS(MML_agnostic_ASM_file),
          readRDS(MML_agnostic_comp_file)
),MML_agnostic_all_file)
# #Unique analyzed region
# unique_gr=unique(granges(MML_all))
# CG_hg19=getCpgSitesHg19()
# CG_hg19_autosome=CG_hg19[seqnames(CG_hg19) %in% c(paste0('chr',c(1:22,'X','Y')))]
# length(subsetByOverlaps(CG_hg19_autosome,unique_gr,minoverlap=2))/length(CG_hg19_autosome)#0.7910136
# CG_hg19_covered=countOverlaps(CG_hg19_autosome,unique_gr,minoverlap=2)
# CG_hg19_covered_tb=table(CG_hg19_covered)
# CG_hg19_covered_tb_cov=CG_hg19_covered_tb[names(CG_hg19_covered_tb)!="0"]
# CG_hg19_covered_tb_cov[1]/sum(CG_hg19_covered_tb_cov)#0.7925
