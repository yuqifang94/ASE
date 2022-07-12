source('mainFunctions_sub.R')

in_dir='../downstream/input/human_analysis/pseComb/'
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
  NME_in_sp=read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,
                                allele_include = F,sample_in=sample_in,hyper_var_file=NA)
  NME_in_sp$tissue=tissue_in
  NME_in=c(NME_in,NME_in_sp)
 }
  else if(stat_in=="MML"){
  MML_in_sp=read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,allele_include = F,sample_in=sample_in,hyper_var_file=NA)
  MML_in_sp$tissue=tissue_in
  MML_in=c(MML_in,MML_in_sp)}else
    {cat("Error stat_in:", stat_in,'\n')}
}
rm(NME_in_sp)
rm(MML_in_sp)
#NME in check: 202721894
#MML in check: 202721894
saveRDS(convert_GR(NME_in[NME_in$N>=2],direction="DT"),"../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_pse.rds")




mix_dt=data.table(tissueName=
                c("endoerm_01_mesoderm_09_paired","endoerm_03_mesoderm_07_paired",
                    "endoerm_05_mesoderm_05_paired", "endoerm_07_mesoderm_03_paired", "endoerm_09_mesoderm_01_paired"),
                  endoermProp=c(0.1,0.3,0.5,0.7,0.9),
                  mesodermProp=c(0.9,0.7,0.5,0.3,0.1)

)
NME_in=readRDS("../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_pse.rds")
NME_all=readRDS(NME_agnostic_all_file)
NME_all=convert_GR(NME_all[NME_all$Sample %in% c("endoerm_27_paired - HUES64","mesoderm_23_paired - HUES64")],direction="DT")
NME_in=cbind(NME_in,mix_dt[match(NME_in$tissue,tissueName)])
NME_in$pureMesoderm=NME_all[Sample=="mesoderm_23_paired - HUES64"][match(NME_in$region,region)]$score
NME_in$pureEndoerm=NME_all[Sample=="endoerm_27_paired - HUES64"][match(NME_in$region,region)]$score
NME_in$predScore=NME_in$pureMesoderm*NME_in$mesodermProp+NME_in$pureEndoerm*NME_in$endoermProp
NME_in=NME_in[!is.na(score)&!is.na(predScore)]
cor(NME_in$score,NME_in$predScore)##Correlation: 0.9257687
saveRDS(NME_in,"../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_pse_calc.rds")
pdf("../downstream/output/human_analysis/cell_deconvolution/NME_pred_score.pdf")
ggplot(NME_in,aes(x=predScore,y=score))+geom_bin2d(bins=100)+xlab("predicted NME")+ylab("actual NME")
dev.off()
MML_in=readRDS("../downstream/output/human_analysis/CPEL_outputs/MML_agnostic_pse.rds")
MML_all=readRDS(MML_agnostic_all_file)
MML_all=convert_GR(MML_all[MML_all$Sample %in% c("endoerm_27_paired - HUES64","mesoderm_23_paired - HUES64")],direction="DT")
MML_in=cbind(MML_in,mix_dt[match(MML_in$tissue,tissueName)])
MML_in$pureMesoderm=MML_all[Sample=="mesoderm_23_paired - HUES64"][match(MML_in$region,region)]$score
MML_in$pureEndoerm=MML_all[Sample=="endoerm_27_paired - HUES64"][match(MML_in$region,region)]$score
MML_in$predScore=MML_in$pureMesoderm*MML_in$mesodermProp+MML_in$pureEndoerm*MML_in$endoermProp
MML_in=MML_in[!is.na(score)&!is.na(predScore)]
cor(MML_in$score,MML_in$predScore)##Correlation: 0.993693
saveRDS(MML_in,"../downstream/output/human_analysis/CPEL_outputs/MML_agnostic_pse_calc.rds")#Correlation: 0.9936582

pdf("../downstream/output/human_analysis/cell_deconvolution/MML_pred_score.pdf")
ggplot(MML_in,aes(x=predScore,y=score))+geom_bin2d(bins=100)+xlab("predicted MML")+ylab("actual MML")
dev.off()
