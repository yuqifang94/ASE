source('mainFunctions_sub.R')
CG_exp_agnostic_hg19_file='../downstream/output/human_analysis/CpG_density/analyzed_region_CG_hg19.rds'
#NME
NME_in=readRDS(NME_agnostic_file)
analyzed_region=unique(granges(NME_in))
gr_seq=getSeq(Hsapiens,analyzed_region,as.character=T)
analyzed_region$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
saveRDS(analyzed_region,CG_exp_agnostic_hg19_file)

analyzed_region=readRDS(CG_exp_agnostic_hg19_file)
NME_in=readRDS(NME_agnostic_file)
MML_in=readRDS(MML_agnostic_file)
NME_in_dt=convert_GR(NME_in,direction="DT")
NME_in_dt$region_Sample=paste0(NME_in_dt$region,'-',NME_in_dt$Sample)
MML_in_dt=convert_GR(MML_in,direction="DT")
MML_in_dt$region_Sample=paste0(MML_in_dt$region,'-',MML_in_dt$Sample)
NME_in_dt$MML=MML_in_dt[match(NME_in_dt$region_Sample,region_Sample)]$MML
analyzed_region=convert_GR(analyzed_region,direction="DT")
NME_in_dt[,NME_res:=resid(loess(NME~MML)),by=list(Sample)]
saveRDS(NME_in_dt,"../downstream/output/human_analysis/CPEL_outputs/NME_in_dt_mean_corrected.rds")
NME_in_dt$CGcont_exp=analyzed_region_sub[match(NME_in_dt$region,region)]$CGcont_exp
NME_in_dt_gr=convert_GR(NME_in_dt$region)
mcols(NME_in_dt_gr)=NME_in_dt
saveRDS(NME_in_dt_gr,"../downstream/output/human_analysis/CPEL_outputs/NME_in_dt_mean_corrected_gr.rds")
CpG_hg19=readRDS(CpG_hg19_file)
NME_in_dt_gr$CG_hg19=countOverlaps(NME_in_dt_gr,CpG_hg19)
NME_in_dt_gr$density=NME_in_dt_gr$CG_hg19/NME_in_dt_gr$CGcont_exp
cor.test(NME_in_dt_gr$density,NME_in_dt_gr$NME,method='pearson')
cor.test(NME_in_dt_gr$density,NME_in_dt_gr$NME_res,method='pearson')

# Mean methylation level is not predictive of normalized methylation entropy.
# Genome-wide scatter plots of pairs of mean methylation level (MML) and normalized methylation entropy
# (NME) values within genomic units in colon normal/cancer and lung normal/cancer, as well as scatter
# plots of pairs of differential mean methylation level (dMML) and differential normalized methylation
# entropy (dNME) values. Although some NME values can be predicted from MML values using the red
# curve determined by the relationship NME = −MML × log2
# (MML) − (1 − MML) × log2
# (1 − MML), in
# general, the NME cannot be inferred from the MML using this formula, which is only valid within genomic
# units containing one CpG site. Notably, genomic units with the same MML values may be characterized by
# different NME values. Moreover, genomic units with zero dMML values are not necessarily characterized
# by zero dNME values. In addition, genomic units with dMML values of ±0.5 are characterized by the
# widest range of dNME values, whereas genomic units with dMNL values of ±1 are characterized by zero
# dNME values.

#linear fit between NME and MML
summary(lm(NME_in_dt$NME~NME_in_dt$MML))

#NME  MAV

NME_in=readRDS("../downstream/output/human_analysis/CPEL_outputs/NME_in_dt_gastric_STL002_gr.rds")
genomic_features=readRDS(genomic_features_file)
GR_calc=data.frame()
scRNA_result=data.frame()
NME_hypervar_calc=GRanges()
NME_in$NME=NME_in_dt_mean_corrected_gr$NME_res
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  hyper_var_file=gsub(scRNA_dir,scRNA_dir,hyper_var_file)
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){
    
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))

    
    
  }else{cat("file not exist for:",sp,'\n')}
}
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir='../downstream/output/human_analysis/NME_MAV/mean_adjust_NME/')
#Small sample STL002 Gastric
NME_in=readRDS(NME_agnostic_file)
MML_in=readRDS(MML_agnostic_file)
NME_in_dt=convert_GR(NME_in,direction="DT")
NME_in_dt=NME_in_dt[Sample=="Gastric_single - STL002"]
NME_in_dt$region_Sample=paste0(NME_in_dt$region,'-',NME_in_dt$Sample)
MML_in_dt=convert_GR(MML_in,direction="DT")
MML_in_dt$region_Sample=paste0(MML_in_dt$region,'-',MML_in_dt$Sample)
NME_in_dt$MML=MML_in_dt[match(NME_in_dt$region_Sample,region_Sample)]$MML
NME_in_dt[,NME_res:=resid(loess(NME~MML)),by=list(Sample)]
saveRDS(NME_in_dt,"../downstream/output/human_analysis/CPEL_outputs/NME_in_dt_gastric_STL002.rds")
NME_in_dt=readRDS("../downstream/output/human_analysis/CPEL_outputs/NME_in_dt_gastric_STL002.rds")
NME_in_dt_gr=convert_GR(NME_in_dt$region)
mcols(NME_in_dt_gr)=NME_in_dt
saveRDS(NME_in_dt_gr,"../downstream/output/human_analysis/CPEL_outputs/NME_in_dt_gastric_STL002_gr.rds")
#Plot MML and NME
NME_in=readRDS(NME_agnostic_file)
MML_in=readRDS(MML_agnostic_file)
NME_in_dt=convert_GR(NME_in,direction="DT")
MML_in_dt=convert_GR(MML_in,direction="DT")
NME_in_dt$region_Sample=paste0(NME_in_dt$region,'-',NME_in_dt$Sample)
MML_in_dt$region_Sample=paste0(MML_in_dt$region,'-',MML_in_dt$Sample)
NME_in_dt$MML=MML_in_dt[match(NME_in_dt$region_Sample,region_Sample)]$MML
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.text.x=element_text(size=10),
                                 axis.text.y=element_text(size=10))
pdf("../downstream/output/human_analysis/QC/NME_vs_MML_bin2d.pdf")
for (sp in unique(NME_in_dt$Sample)){
    print(ggplot(NME_in_dt[Sample==sp],
                  aes(x=MML, y=NME))+
                  xlim(c(0,1))+ylim(c(0,1))+ggtitle(paste0("MML and NME relationship\n",sp))+
                  #geom_smooth(method="loess",se=FALSE)+
                  xlab("MML")+ ylab("NME")+geom_bin2d(bins=100)+theme_glob+
                  #geom_smooth(method="lm",color="red")+
                  scale_fill_gradient(name = "count", trans = "log10",high="#132B43",low="#56B1F7"))
}
dev.off()
#add line geom_function(fun = function(x) 0.5*exp(-abs(x)))

for (sp in unique(NME_in_dt$Sample)){
    png(paste0("../downstream/output/human_analysis/QC/NME_vs_MML_dot/",sp,".png"))
    print(ggplot(NME_in_dt[Sample==sp],
                  aes(x=MML, y=NME))+
                  xlim(c(0,1))+ylim(c(0,1))+ggtitle(paste0("MML and NME relationship\n",sp))+
                  #geom_smooth(method="loess",se=FALSE)+
                  xlab("MML")+ ylab("NME")+geom_point(size=0.01)+theme_glob)
    dev.off()
}


