calc_OR<-function(feature,GR_all,GR_hit,CpG){
  olap_gr=findOverlaps(GR_all,GR_hit,type='equal')
  GR_not_hit=GR_all[-queryHits(olap_gr)]
  GR_hit=subsetByOverlaps(CpG,GR_hit,type='within')
  GR_not_hit=subsetByOverlaps(CpG,GR_not_hit,type='within')
  GR_hit_feature=length(subsetByOverlaps(GR_hit,feature))
  
  GR_hit_notfeature=length(GR_hit)-GR_hit_feature
  GR_not_hit_feature=length(subsetByOverlaps(GR_not_hit,feature))
  GR_not_hit_not_feature=length(GR_not_hit)-GR_not_hit_feature
  cont_table=matrix(c(GR_hit_feature,GR_hit_notfeature,GR_not_hit_feature,GR_not_hit_not_feature),nrow=2)
  colnames(cont_table)=c('Hit','Not Hit')
  rownames(cont_table)=c('Feature','Not feature')
  print(fisher.test(cont_table))
  return(cont_table)
}

MML=GR[GR$Statistic=='dMML'&GR$ASM=='Yes']
#MML_SNP=dat_check_MML[dat_check_MML$diffMML>0]
MML_GR=subsetByOverlaps(MML,dat_check_MML)
calc_OR(outGR[['CpG open sea']],MML,GR_sub_ASM_MML,CpG_all)
calc_OR(outGR[['CpG island']],MML,GR_sub_ASM_MML,CpG_all)
calc_OR(outGR[['CpG shore']],MML,GR_sub_ASM_MML,CpG_all)
calc_OR(outGR[['CpG shelf']],MML,GR_sub_ASM_MML,CpG_all)
calc_OR(outGR[['promoter']],MML,GR_sub_ASM_MML,CpG_all)
calc_OR(outGR[['gene body']],MML,GR_sub_ASM_MML,CpG_all)



NME=GR[GR$Statistic=='dNME'&GR$ASM=='Yes']
#NME_SNP=dat_check_NME[dat_check_NME$diffNME>0]
NME_GR=subsetByOverlaps(NME,dat_check_NME)
calc_OR(outGR[['CpG open sea']],NME,GR_sub_ASM_NME,CpG_all)
calc_OR(outGR[['CpG island']],NME,GR_sub_ASM_NME,CpG_all)
calc_OR(outGR[['CpG shore']],NME,GR_sub_ASM_NME,CpG_all)
calc_OR(outGR[['CpG shelf']],NME,GR_sub_ASM_NME,CpG_all)
calc_OR(outGR[['promoter']],NME,GR_sub_ASM_NME,CpG_all)


####Density Distribution: balanced CG vs imbalanced CG####
#Form data frame, using ASM regions 
NME_ASM_CG_df=data.frame(density=NME_ASM$density,CG_type=NME_ASM$CG_type)
NME_ASM_CG_df_sub=NME_ASM_CG_df[NME_ASM_CG_df$density<=0.02,]
aggregate(NME_ASM_CG_df_sub,by=list(NME_ASM_CG_df_sub$CG_type),FUN=mean)
ggplot(NME_ASM_CG_df[NME_ASM_CG_df$density<=0.025,],aes(x=density,fill=CG_type))+#scale_fill_manual(values = c("blue","red"))+
  geom_density(stat = "density",alpha=0.6)+xlab('CpG density')+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
  ggtitle('dNME ASM CpG density at balanced CG vs unbalanced CG')

#####Density plot ASM vs non-ASM####
NME_CG_df=data.frame(density=NME_het$density,ASM=NME_het$ASM)
NME_CG_df=NME_CG_df[!is.na(NME_CG_df$ASM),]
ggplot(NME_CG_df[NME_CG_df$density<0.05,],aes(x=density,fill=ASM))+scale_fill_manual(values = c("blue","red"))+
  geom_density(alpha=0.6)+xlab('CpG density')+ylim(c(0,250))+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
  ggtitle('dNME CpG density at ASM vs non ASM (density <0.05)')

#####Plot log10 density vs difference####
density_df=data.frame(density=as.factor(round(round(log10(NME_ASM_het$density),digits = 1),digits = 4)),diff=NME_ASM_het$diff)
density_df=data.frame(density=as.factor(NME_ASM_het$density),diff=NME_ASM_het$diff)
#Violin plot is showing something if properly binned
ggplot(density_df,aes(x=density,y=diff))+geom_violin()+ stat_summary(fun.y=median, geom="point", size=2, color="red")+
  xlab('local CpG density')+ylab('dNME (More -less)')+ggtitle('dNME distribution at different CpG density')

####Plot aggregated graph using density levels and average dNME in that level####
#For ASM with CpG difference only
density_df=data.frame(density=round(NME_ASM_het$CGcount_diff,digits = 2),diff=NME_ASM_het$diff)
density_agg=aggregate(density_df,by=list(density_df$density),FUN=mean)
plot(density_agg$Group.1,density_agg$diff,pch=1,xlab='log10(CpG density)',ylab='dNME',main='dNME distribution with CpG density at ASM')
#For all regions at ASM
density_df=data.frame(density=round(log10(NME_ASM$density),digits = 1),diff=NME_ASM_het$diff)
density_agg=aggregate(density_df,by=list(density_df$density),FUN=mean)
plot(density_agg$Group.1,density_agg$diff,pch=1,xlab='log10(CpG density)',ylab='dNME',main='dNME distribution with CpG density at ASM')