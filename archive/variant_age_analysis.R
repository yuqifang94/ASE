source('mainFunctions_sub.R')
# checking Age ------------------------------------------------------------
age_out=readRDS('../downstream/input/age_atlas_sub.rds')
qual_score=0.5
#10055361, 3575815 have SNP,867652 SNP
age_out=age_out[age_out$QualScore_Jnt>=qual_score&age_out$QualScore_Mut>=qual_score&age_out$QualScore_Rec>=qual_score]
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
genome_freq=readRDS('../downstream/input/genome_1k_variant.rds')
olap=findOverlaps(variant_HetCpG_meta,age_out)
variant_HetCpG_meta$variant_age=NA
variant_HetCpG_meta$variant_age[queryHits(olap)]=age_out$AgeMedian_Jnt[subjectHits(olap)]
olap=findOverlaps(variant_HetCpG_meta,genome_freq)
variant_HetCpG_meta$DAF=NA
variant_HetCpG_meta$DAF[queryHits(olap)]=genome_freq$MAF[subjectHits(olap)]
genomic_features=readRDS(genomic_features_file)
olap=findOverlaps(variant_HetCpG_meta,genomic_features$TSS,maxgap = 500)
variant_HetCpG_meta$TSS=NA
variant_HetCpG_meta$TSS[queryHits(olap)]=genomic_features$TSS$gene_name[subjectHits(olap)]
rm(age_out)
rm(genome_freq)
# olap=findOverlaps(variant_HetCpG_meta,motif_gene_sig_dir)
# variant_HetCpG_meta$motif="No"
# variant_HetCpG_meta$motif[queryHits(olap)]='Yes'
NME_quant=quantile(c(variant_HetCpG_meta$NME1,variant_HetCpG_meta$NME2),prob=0.75)
olap=findOverlaps(variant_HetCpG_meta,genomic_features$promoter,maxgap = 500)
variant_HetCpG_meta$promoter=NA
variant_HetCpG_meta$promoter[queryHits(olap)]=genomic_features$promoter$gene_name[subjectHits(olap)]
GR_merge=readRDS(GR_merge_file)
olap=findOverlaps(variant_HetCpG_meta,GR_merge)
GR_merge$density=GR_merge$CG_hg19_extend/GR_merge$gff_size_extend
variant_HetCpG_meta$density=NA
variant_HetCpG_meta$density[queryHits(olap)]=GR_merge$density[subjectHits(olap)]
#summary analysis
median(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff],na.rm = T)
median(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dMML_pval<=pval_cutoff],na.rm = T)
median(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dMML_pval>pval_cutoff&variant_HetCpG_meta$dMML_pval>pval_cutoff],na.rm = T)
median(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$motif=='Yes'],na.rm = T)
median(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$motif=='No'],na.rm = T)
#HetCpG is older than nonHet CpG, but still the dNME-Hap is younger than non-dNME-Hap
median(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$HetCpG],na.rm = T)
median(variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&!variant_HetCpG_meta$HetCpG],na.rm = T)
#HetCpG
all_df_HetCpG=rbind(data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&!variant_HetCpG_meta$HetCpG],
                               type='Non-HetCpG-dNME'),
                    data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$HetCpG],
                               type='HetCpG-dNME'),
                    data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$HetCpG],
                               type='HetCpG-nondNME'), 
                    data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&!variant_HetCpG_meta$HetCpG],
                               type='non-HetCpG-non-dNME'))

aggregate(all_df_HetCpG$age[!is.na(all_df_HetCpG$age)],by=list(all_df_HetCpG$type[!is.na(all_df_HetCpG$age)]),FUN=median)
ggplot(all_df_HetCpG,aes(x=age,color=type))+geom_density(size=1)+xlab('variant age')+theme(legend.position = 'bottom')
ggplot(all_df_HetCpG,aes(x=age,color=type))+stat_ecdf(size=1)+xlab('variant age')+theme(legend.position = 'bottom')
#Motif binding
all_df_motif=rbind(data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$motif=='No'],
                              type='Non-Motif-dNME'),
                   data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$motif=='Yes'],
                              type='Motif-dNME'),
                   data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$motif=='Yes'],
                              type='Motif-nondNME'), 
                   data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$motif=='No'],
                              type='non-Motif-non-dNME'))
aggregate(all_df_motif$age[!is.na(all_df_motif$age)],by=list(all_df_motif$type[!is.na(all_df_motif$age)]),FUN=median)
ggplot(all_df_motif[!is.na(all_df_motif$age),],aes(x=age,color=type))+geom_density(size=1)+xlab('variant age')+theme(legend.position = 'bottom')

#Density
all_df_density=rbind(data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&log10(variant_HetCpG_meta$density)>=-2],
                                type='high-density-dNME'),
                     data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval<=pval_cutoff&log10(variant_HetCpG_meta$density)< -2],
                                type='low-density-dNME'),
                     data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&log10(variant_HetCpG_meta$density)>=-2],
                                type='high-density-nondNME'), 
                     data.frame(age=variant_HetCpG_meta$variant_age[variant_HetCpG_meta$dNME_pval>pval_cutoff&log10(variant_HetCpG_meta$density)< -2],
                                type='low-density-non-dNME'))
aggregate(all_df_density$age[!is.na(all_df_density$age)],by=list(all_df_density$type[!is.na(all_df_density$age)]),FUN=median)
ggplot(all_df_density[!is.na(all_df_density$age),],aes(x=age,color=type))+geom_density(size=1)+xlab('variant age')+theme(legend.position = 'bottom')

#HetCpG location is older

#Large haplotype ->larger density?
variant_HetCpG_meta=variant_HetCpG_meta[variant_HetCpG_meta$N>=2]
dNME=which(variant_HetCpG_meta$dNME_pval<=pval_cutoff)
dMML=which(variant_HetCpG_meta$dMML_pval<=pval_cutoff)
non_dNME_dMML=variant_HetCpG_meta$dMML_pval>pval_cutoff& variant_HetCpG_meta$dNME_pval

dNME_in=variant_HetCpG_meta$variant_age[dNME]
dNME_in=dNME_in[!is.na(dNME_in)]
dMML_in=variant_HetCpG_meta$variant_age[dMML]
dMML_in=dMML_in[!is.na(dMML_in)]
non_dNME_dMML_in=variant_HetCpG_meta$variant_age[non_dNME_dMML]
non_dNME_dMML_in=non_dNME_dMML_in[!is.na(non_dNME_dMML_in)]


#ecdf takes too long run it later
fit1.dNME <- ecdf(dNME_in)
fit1.dMML<-ecdf(dMML_in)
fit1.non_dNME_dMML<-ecdf(non_dNME_dMML_in)
tt=proc.time()[[3]]
# fit2.dNME <- replicate(10000, { x <- sample(dNME_in, replace=TRUE);
# ecdf(x)(unique(dNME_in))})
# saveRDS(fit2.dNME,'../downstream/output/variant_age/edcf.dNME2.rds')
# fit2.dMML <- replicate(10000, { x <- sample(dMML_in, replace=TRUE);
# ecdf(x)(unique(dMML_in))})
# saveRDS(fit2.dMML,'../downstream/output/variant_age/edcf.dMML2.rds')
# fit2.non_dNME_dMML<- replicate(10000, { x <- sample(non_dNME_dMML_in, replace=TRUE);
# ecdf(x)(unique(non_dNME_dMML_in))})
# saveRDS(fit2.non_dNME_dMML,'../downstream/output/variant_age/edcf.non_dNME_dMML2.rds')
# cat('Finishing bootstrap',proc.time()[[3]]-tt,'\n')
fit2.dNME=readRDS('../downstream/output/variant_age/edcf.dNME.rds')
fit2.dMML=readRDS('../downstream/output/variant_age/edcf.dMML.rds')
fit2.non_dNME_dMML=readRDS('../downstream/output/variant_age/edcf.non_dNME_dMML.rds')

fit3.dNME <- apply(fit2.dNME, 1, quantile, c(0.025,0.975) )
fit3.dMML <- apply(fit2.dMML, 1, quantile, c(0.025,0.975) )
fit3.non_dNME_dMML <- apply(fit2.non_dNME_dMML, 1, quantile, c(0.025,0.975) )
jpeg('../downstream/output/variant_age/ecdf_N1.jpg')
plot(fit1.dNME, ylim=range(0,1))

polygon(c(fit1.dMML(unique(dMML_in)), rev(fit1.dMML(unique(dMML_in)))), c(fit3.dMML[1,], rev(fit3.dMML[2,])), col='lightblue', border=F)
polygon(c(fit1.non_dNME_dMML(unique(non_dNME_dMML_in)), rev(fit1.non_dNME_dMML(unique(non_dNME_dMML_in)))), c(fit3.non_dNME_dMML[1,], rev(fit3.non_dNME_dMML[2,])), col='grey', border=F)
polygon(c(fit1.dNME(unique(dNME_in)), rev(fit1.dNME(unique(dNME_in)))), c(fit3.dNME[1,], rev(fit3.dNME[2,])), col='pink', border=F)
lines(fit1.dMML,col='blue')
lines(fit1.dNME,col='red')
lines(fit1.non_dNME_dMML,col='black')
dev.off()



fit1.dNME <- density(dNME_in)
fit1.dMML<-density(dMML_in)
fit1.non_dNME_dMML<-density(non_dNME_dMML_in)
# tt=proc.time()[[3]]
# fit2.dNME <- mcreplicate(10000, {density(sample(dNME_in, replace=TRUE), 
#                                          from=min(fit1.dNME$x), to=max(fit1.dNME$x))$y },mc.cores=8)
# saveRDS(fit2.dNME,'../downstream/output/variant_age/density.dNME.rds')
# fit2.dMML <- mcreplicate(10000, {density(sample(dMML_in, replace=TRUE),
#                                          from=min(fit1.dMML$x), to=max(fit1.dMML$x))$y},mc.cores=8)
# saveRDS(fit2.dMML,'../downstream/output/variant_age/density.dMML.rds')
# fit2.non_dNME_dMML <- mcreplicate(10000, {density(sample(non_dNME_dMML_in, replace=TRUE), 
#                                                   from=min(fit1.non_dNME_dMML$x), to=max(fit1.non_dNME_dMML$x))$y},mc.cores=8)
# saveRDS(fit2.non_dNME_dMML,'../downstream/output/variant_age/density.non_dNME_dMML.rds')
# cat('Finishing bootstrap',proc.time()[[3]]-tt,'\n')

fit2.dNME=readRDS('../downstream/output/variant_age/density.dNME.rds')
fit2.dMML=readRDS('../downstream/output/variant_age/density.dMML.rds')
fit2.non_dNME_dMML=readRDS('../downstream/output/variant_age/density.non_dNME_dMML.rds')

fit3.dNME <- apply(fit2.dNME, 1, quantile, c(0.025,0.975) )
fit3.dMML <- apply(fit2.dMML, 1, quantile, c(0.025,0.975) )
fit3.non_dNME_dMML <- apply(fit2.non_dNME_dMML, 1, quantile, c(0.025,0.975) )

plot(fit1.dNME, ylim=range(c(fit3.dMML,fit3.dNME,fit3.non_dNME_dMML)))

polygon(c(fit1.dMML$x, rev(fit1.dMML$x)), c(fit3.dMML[1,], rev(fit3.dMML[2,])), col='lightblue', border=F)
polygon(c(fit1.non_dNME_dMML$x, rev(fit1.non_dNME_dMML$x)), c(fit3.non_dNME_dMML[1,], rev(fit3.non_dNME_dMML[2,])), col='grey', border=F)
polygon(c(fit1.dNME$x, rev(fit1.dNME$x)), c(fit3.dNME[1,], rev(fit3.dNME[2,])), col='pink', border=F)
lines(fit1.dMML,col='blue')
lines(fit1.dNME,col='red')
lines(fit1.non_dNME_dMML,col='black')

all_df=rbind(data.frame(age=variant_HetCpG_meta$variant_age[dNME],
                        DAF=variant_HetCpG_meta$DAF[dNME],
                        type='dNME-Hap'),
             data.frame(age=variant_HetCpG_meta$variant_age[dMML],
                        DAF=variant_HetCpG_meta$DAF[dMML],
                        type='dMML-Hap'),
             data.frame(age=variant_HetCpG_meta$variant_age[non_dNME_dMML], 
                        DAF=variant_HetCpG_meta$DAF[non_dNME_dMML],
                        type='non-dMML-dNME-Hap')
)
all_df=all_df[!is.na(all_df$DAF),]
ggplot(all_df,aes(x=age,color=type))+geom_density(size=1)+xlab('variant age')+theme(legend.position = 'bottom')
ggplot(all_df,aes(x=DAF,color=type))+stat_ecdf(size=1)+xlab('Derived Allele frequency')+theme(legend.position = 'bottom')
#Age per sample 
pdf('../downstream/output/variant_age/variant_age_sp_density_ecdf.pdf')
variant_HetCpG_meta=variant_HetCpG_meta[!is.na(variant_HetCpG_meta$variant_age)]
variant_age_out=data.frame()
all_df_out=data.frame()
for (sp in unique(variant_HetCpG_meta$Sample)){
  cat('processing sample:',sp,'\n')
  dNME_sp=which(variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$Sample==sp)
  dMML_sp=which(variant_HetCpG_meta$dMML_pval<=pval_cutoff&variant_HetCpG_meta$Sample==sp)
  non_dNME_dMML_sp=variant_HetCpG_meta$dMML_pval>pval_cutoff& variant_HetCpG_meta$dNME_pval>pval_cutoff&variant_HetCpG_meta$Sample==sp
  all_df_sp=rbind(data.frame(age=variant_HetCpG_meta$variant_age[dNME_sp],
                             DAF=variant_HetCpG_meta$DAF[dNME_sp],
                             type='dNME-Hap',stringsAsFactors = F),
                  data.frame(age=variant_HetCpG_meta$variant_age[dMML_sp],
                             DAF=variant_HetCpG_meta$DAF[dMML_sp],
                             type='dMML-Hap',stringsAsFactors = F),
                  data.frame(age=variant_HetCpG_meta$variant_age[non_dNME_dMML_sp], 
                             DAF=variant_HetCpG_meta$DAF[non_dNME_dMML_sp],
                             type='non-dMML-dNME-Hap',stringsAsFactors = F)
                 
  )

  all_df_sp=all_df_sp[!is.na(all_df_sp$DAF)&!is.na(all_df_sp$age),]
  all_df_sp$sample=sp
  all_df_sp$tissue=unique(variant_HetCpG_meta$tissue[dNME_sp])
  all_df_sp$subject=strsplit(sp,' - ')[[1]][2]
  all_df_out=rbind(all_df_out,all_df_sp)
  
  cat('plotting\n')
# 
#  print(ggplot(all_df_sp,aes(x=DAF,y=age,color=type))+geom_smooth(size=1)+xlab('DAF')+ylab('variant age')+theme(legend.position = 'bottom')+
#  ggtitle(sp))
#   # ggplot(all_df)+geom_histogram(aes(x=round(DAF,digits = 1),y = ..density..,fill=type,color=type),binwidth=0.1,position="dodge")+xlab('DAF')+
#   #   theme(legend.position = 'bottom')
#   print(ggplot(all_df_sp,aes(x=age,color=type))+geom_density(size=1)+xlab('variant age')+theme(legend.position = 'bottom')+ggtitle(sp))
#   print(ggplot(all_df_sp,aes(x=age,color=type))+stat_ecdf(size=1)+xlab('variant age')+theme(legend.position = 'bottom')+ggtitle(sp))
  # if(sum(all_df_sp$type=='dMML-Hap')>2){dMML_age=t.test(all_df_sp$age[all_df_sp$type=='dMML-Hap'])$estimate}else(dMML_age=NA)
  # dNME_age=t.test(all_df_sp$age[all_df_sp$type=='dNME-Hap'])$estimate
  # non_dMML_dNME_age=t.test(all_df_sp$age[all_df_sp$type=='non-dMML-dNME-Hap'])$estimate
  # if(sum(variant_HetCpG_meta$Sample==sp& variant_HetCpG_meta$dNME_pval<=pval_cutoff)>2&
  #    sum(variant_HetCpG_meta$Sample==sp& variant_HetCpG_meta$dMML_pval<=pval_cutoff)>2){
  # variant_age_out=rbind(variant_age_out,
  #                       data.frame(dMML_age=dMML_age,
  #                                  dNME_age=dNME_age,
  #                                  non_dMML_dNME_age=non_dMML_dNME_age,
  #                                  cor_dNME_sig=cor.test(variant_HetCpG_meta$dNME[variant_HetCpG_meta$Sample==sp&
  #                                                                                   variant_HetCpG_meta$dNME_pval<=pval_cutoff],
  #                                                    variant_HetCpG_meta$variant_age[variant_HetCpG_meta$Sample==sp&
  #                                                                                      variant_HetCpG_meta$dNME_pval<=pval_cutoff])$estimate,
  #                                  cor_dMML_sig=cor.test(variant_HetCpG_meta$dMML[variant_HetCpG_meta$Sample==sp&
  #                                                                                   variant_HetCpG_meta$dMML_pval<=pval_cutoff],
  #                                                    variant_HetCpG_meta$variant_age[variant_HetCpG_meta$Sample==sp&
  #                                                                                      variant_HetCpG_meta$dMML_pval<=pval_cutoff])$estimate,
  #                                  cor_dNME=cor.test(variant_HetCpG_meta$dNME[variant_HetCpG_meta$Sample==sp],
  #                                                        variant_HetCpG_meta$variant_age[variant_HetCpG_meta$Sample==sp])$estimate,
  #                                  cor_dMML=cor.test(variant_HetCpG_meta$dMML[variant_HetCpG_meta$Sample==sp],
  #                                                        variant_HetCpG_meta$variant_age[variant_HetCpG_meta$Sample==sp])$estimate,
  #                                  cor_dNME_pval=cor.test(variant_HetCpG_meta$dNME[variant_HetCpG_meta$Sample==sp],
  #                                                    variant_HetCpG_meta$variant_age[variant_HetCpG_meta$Sample==sp])$p.value,
  #                                  cor_dMML_pval=cor.test(variant_HetCpG_meta$dMML[variant_HetCpG_meta$Sample==sp],
  #                                                    variant_HetCpG_meta$variant_age[variant_HetCpG_meta$Sample==sp])$p.value,
  #                                  sample_sp=sp,stringsAsFactors = F))
  # }
}   
dev.off()
#get the TSS
#try stat_ecdf maybe
all_df_out$sample_type=paste(all_df_out$sample,all_df_out$type,sep='_')
pdf('../downstream/output/variant_age/ecdf_by_sample_subject.pdf')
ggplot(all_df_out,aes(x=age,group=sample_type))+stat_ecdf(size=0.1,aes(color=subject),alpha=0.5)+xlab('Age')+theme(legend.position = 'bottom')
dev.off()
pdf('../downstream/output/variant_age/ecdf_daf_by_sample_type.pdf')
ggplot(all_df_out,aes(x=age,group=sample_type))+stat_ecdf(size=0.1,aes(color=type),alpha=0.5)+xlab('Age')+theme(legend.position = 'bottom')
dev.off()
ggplot(all_df_out,aes(x=age,group=type))+geom_density(size=1,aes(color=type))+xlab('Age')+theme(legend.position = 'bottom')
ggplot(all_df_out,aes(x=age,group=type))+stat_ecdf(size=1,aes(color=type))+xlab('Age')+theme(legend.position = 'bottom')
#Rank genes
variant_HetCpG_meta_uq=unique(variant_HetCpG_meta[!is.na(variant_HetCpG_meta$genes_promoter)])
variant_HetCpG_meta_uq$dNME_sum=NA
variant_HetCpG_meta_uq$dMML_sum=NA
olap=findOverlaps(variant_HetCpG_meta,variant_HetCpG_meta_uq,type='equal')
variant_HetCpG_meta_promoter=variant_HetCpG_meta[!is.na(variant_HetCpG_meta$genes_promoter)&!is.na(variant_HetCpG_meta$variant_age)]
#variant_HetCpG_meta_promoter=variant_HetCpG_meta_promoter[variant_HetCpG_meta_promoter$dNME_pval<=pval_cutoff]
summary_stat_df=data.frame(dNME=variant_HetCpG_meta_promoter$dNME,
                           dMML=variant_HetCpG_meta_promoter$dMML,genes=variant_HetCpG_meta_promoter$genes_promoter,
                           age=variant_HetCpG_meta_promoter$variant_age)
summary_stat_agg=aggregate(summary_stat_df,by=list(summary_stat_df$genes),FUN=mean)
summary_stat_agg=summary_stat_agg[,c(1,2,3,5)]
names(summary_stat_agg)=c('gene_name','dNME_sum','dMML_sum','variant_age')
age_rank=rank(-summary_stat_agg$variant_age)
dNME_rank=rank(-summary_stat_agg$dNME_sum)
#N_rank=rank(-summary_stat_agg$N)
genes_rank_gsea=rank(age_rank+dNME_rank)
names(genes_rank_gsea)=summary_stat_agg$gene_name

#Merge to a single rank

gsea_anno=read.gmt('../downstream/input/gsea/c2.all.v7.1.symbols.gmt')
gsea_anno_background=lapply(gsea_anno,function(x) x[x%in%unique(variant_HetCpG_meta$genes_promoter)])
fgseaRes <- fgsea(pathways = gsea_anno_background, 
                  stats    = genes_rank_gsea,
                  nperm     = 1000,
                  minSize  = 10,
                  maxSize  = 500,
                  nproc=20,BPPARAM=SnowParam())

genes_rank=variant_HetCpG_meta_uq$genes_promoter[order(age_rank+dNME_rank+N_rank)]

write(genes_rank,'../downstream/output/variant_age/genes_rank_dNME_age_promoter.txt')
#TSS
variant_HetCpG_meta_uq=unique(variant_HetCpG_meta[!is.na(variant_HetCpG_meta$TSS)])
variant_HetCpG_meta_uq$dNME_sum=NA
variant_HetCpG_meta_uq$dMML_sum=NA
olap=findOverlaps(variant_HetCpG_meta,variant_HetCpG_meta_uq,type='equal')
summary_stat_df=data.frame(dNME=variant_HetCpG_meta$dNME[queryHits(olap)],dMML=variant_HetCpG_meta$dMML[queryHits(olap)],sht=subjectHits(olap))
summary_stat_agg=aggregate(summary_stat_df,by=list(summary_stat_df$sht),FUN=mean)
variant_HetCpG_meta_uq$dNME_sum[summary_stat_agg$sht]=summary_stat_agg$dNME
variant_HetCpG_meta_uq$dMML_sum[summary_stat_agg$sht]=summary_stat_agg$dMML
age_rank=rank(-variant_HetCpG_meta_uq$variant_age)
dNME_rank=rank(-variant_HetCpG_meta_uq$dNME_sum)
N_rank=rank(-variant_HetCpG_meta_uq$N)

genes_rank=variant_HetCpG_meta_uq$genes_promoter[order(age_rank+dNME_rank+N_rank)]


write(genes_rank,'../downstream/output/variant_age/genes_rank_dNME_age_TSS.txt')


genes_sig=unique(variant_HetCpG_meta$genes_promoter[variant_HetCpG_meta$variant_age<=20000&variant_HetCpG_meta$dNME_pval<=pval_cutoff&variant_HetCpG_meta$N>=1])
genes_all=unique(variant_HetCpG_meta$genes_promoter[variant_HetCpG_meta$variant_age<=20000&variant_HetCpG_meta$N>=1])
GO_young=GO_anno(genes_sig,genes_all)
GO_young_genes=GO_young[[1]]
GO_young_genes_sig=GO_young_genes[GO_young_genes$Significant>=10&GO_young_genes$FC>=2,]
GO_young_genes_sig$qval=p.adjust(GO_young_genes_sig$classicFisher)
GO_young_genes_sig=GO_young_genes_sig[GO_young_genes_sig$qval<=0.1,]
GO_young_obj=GO_young[[2]]
gene_out=genesInTerm(GO_young_obj, GO_young_genes_sig$GO.ID)
names(gene_out)=GO_young_genes_sig$Term
gene_out=lapply(gene_out,function(x) x[x%in%genes_sig])

cat(unique(unlist(gene_out)),sep=',')

