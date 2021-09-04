
# compare nme with other tissues
nme=readRDS('../downstream/output/nme_matrix_DNase.rds')
nme_prob=nme[,c("EFP-E15.5-all","forebrain-E15.5-all","heart-E15.5-all","liver-E15.5-all","midbrain-E15.5-all")]
nme_no_prob=nme[,c("EFP-E14.5-all","forebrain-E14.5-all","heart-E14.5-all","liver-E14.5-all","midbrain-E14.5-all")]
nme_no_prob=nme[,c("forebrain-E16.5-all","heart-E16.5-all","liver-E16.5-all","midbrain-E16.5-all")]
nme_comb=rbind(data.table(nme=as.numeric(nme_prob[!is.na(nme_prob)]),region_type="bias problem"),
               data.table(nme=as.numeric(nme_no_prob[!is.na(nme_no_prob)]),region_type="no bias problem"))
ggplot(nme_comb,aes(x=nme,color=region_type,group=region_type))+geom_density()

nme_prob_difference=rowMeans(nme_prob,na.rm = T)-rowMeans(nme_no_prob,na.rm=T)
GO_out_all_region=unique(do.call('c',lapply(GO_out_all$all,function(x) {
  do.call('c',lapply(x,function(clu){
    return(clu$csv_in_ts_clu$region)
    
    
    
  }))})))  
GO_out_all_region_forebrain=unique(do.call('c',lapply(GO_out_all$all$forebrain,function(clu) {
  
  return(clu$csv_in_ts_clu$region)
  
  
  
}))) 
ggplot(data.frame(difference=nme_prob_difference[GO_out_all_region]),aes(difference))+geom_density(size=1)
nme_prob=nme[GO_out_all_region,c("EFP-E15.5-all","forebrain-E15.5-all","heart-E15.5-all","liver-E15.5-all","midbrain-E15.5-all")]
nme_no_prob=nme[GO_out_all_region,c("EFP-E14.5-all","forebrain-E14.5-all","heart-E14.5-all","liver-E14.5-all","midbrain-E14.5-all")]
nme_comb=rbind(data.table(nme=as.numeric(nme_prob[!is.na(nme_prob)]),region_type="bias problem"),
               data.table(nme=as.numeric(nme_no_prob[!is.na(nme_no_prob)]),region_type="no bias problem"))
ggplot(nme_comb,aes(x=nme,color=region_type,group=region_type))+geom_density(size=1)

UC_all=readRDS('../downstream/output/UC_merge_max_loc.rds')
hist(nme[,"forebrain-E15.5-all"]-nme[,"forebrain-E16.5-all"])
forebrain_dnme_15_16=nme[,"forebrain-E15.5-all"]-nme[,"forebrain-E16.5-all"]
forebrain_dnme_14_16=nme[,"forebrain-E14.5-all"]-nme[,"forebrain-E16.5-all"]
plot(density(forebrain_dnme[!is.na(forebrain_dnme)]))
plot(density(forebrain_dnme[GO_out_all_region_forebrain],na.rm=T))
forebrain_dnme=rbind(data.table(dnme=forebrain_dnme_15_16,region_type="dnme_15_16"),
                     data.table(dnme=forebrain_dnme_14_16,region_type="dnme_14_16"))
ggplot(forebrain_dnme,aes(x=dnme,color=region_type))+geom_density(size=1)

#Compare trim 30 vs trim 5
#Forebrain
nme=readRDS('../downstream/output/nme_matrix_DNase.rds')
forebrain_30_15=import.bedGraph('../downstream/input/trim15/mm10_forebrain_day15_5_all_allele_agnostic_nme.bedGraph')
start(forebrain_30_15)=start(forebrain_30_15)-1

forebrain_comp=data.table(region=paste0(seqnames(forebrain_30_15),":",start(forebrain_30_15),'-',end(forebrain_30_15)),
                          nme_30=forebrain_30_15$score)
forebrain_comp=forebrain_comp[region%in%rownames(nme)]
forebrain_comp$nme_5=nme[forebrain_comp$region,grepl("forebrain-E15.5-all",colnames(nme))]
forebrain_comp$nme_E16=nme[forebrain_comp$region,grepl("forebrain-E16.5-all",colnames(nme))]
forebrain_comp$nme_E14=nme[forebrain_comp$region,grepl("forebrain-E14.5-all",colnames(nme))]
forebrain_comp$relative_diff=forebrain_comp$nme_30-forebrain_comp$nme_5
ggplot(forebrain_comp,aes(x=relative_diff))+geom_density()+xlim(c(-0.1,0.1))
ggplot(forebrain_comp,aes(x=relative_diff))+geom_density()
forebrain_comp_melt=melt.data.table(forebrain_comp[,list(region,nme_30,nme_5,nme_E14)],id.vars = "region")
ggplot(forebrain_comp_melt,aes(x=value,color=variable))+geom_density(size=1)+ggtitle("forebrain E15.5")
GO_out_all=readRDS('../downstream/output/GO_out_all_dMML_dNME_0rm_FC_enhancer.rds')
GO_out_forebrain_region=do.call('c',lapply(GO_out_all$all$forebrain,function(x) x$csv_in_ts_clu$region))
ggplot(forebrain_comp[region%in% GO_out_forebrain_region],aes(x=relative_diff))+geom_density()+xlim(c(-0.5,0.5))
ggplot(forebrain_comp[region%in% GO_out_forebrain_region],aes(x=relative_diff))+geom_density()+xlim(c(-0.1,0.1))
#EFP
EFP_30_15=import.bedGraph('../downstream/input/trim15/mm10_EFP_day15_5_all_allele_agnostic_nme.bedGraph')
start(EFP_30_15)=start(EFP_30_15)-1

EFP_comp=data.table(region=paste0(seqnames(EFP_30_15),":",start(EFP_30_15),'-',end(EFP_30_15)),
                    nme_30=EFP_30_15$score)
EFP_comp=EFP_comp[region%in%rownames(nme)]
EFP_comp$nme_5=nme[EFP_comp$region,grepl("EFP-E15.5-all",colnames(nme))]
EFP_comp$nme_E14=nme[EFP_comp$region,grepl("EFP-E14.5-all",colnames(nme))]
EFP_comp$relative_diff=EFP_comp$nme_30-EFP_comp$nme_5
ggplot(EFP_comp,aes(x=relative_diff))+geom_density()
EFP_comp_melt=melt.data.table(EFP_comp[,list(region,nme_30,nme_5,nme_E14)],id.vars = "region")
ggplot(EFP_comp_melt,aes(x=value,color=variable))+geom_density(size=1)+ggtitle("EFP E15.5")

#heart
heart_30_15=import.bedGraph('../downstream/input/trim15/mm10_heart_day15_5_all_allele_agnostic_nme.bedGraph')
start(heart_30_15)=start(heart_30_15)-1

heart_comp=data.table(region=paste0(seqnames(heart_30_15),":",start(heart_30_15),'-',end(heart_30_15)),
                      nme_30=heart_30_15$score)
heart_comp=heart_comp[region%in%rownames(nme)]
heart_comp$nme_5=nme[heart_comp$region,grepl("heart-E15.5-all",colnames(nme))]
heart_comp$nme_E14=nme[heart_comp$region,grepl("heart-E14.5-all",colnames(nme))]
heart_comp$relative_diff=heart_comp$nme_30-heart_comp$nme_5
ggplot(heart_comp,aes(x=relative_diff))+geom_density()
heart_comp_melt=melt.data.table(heart_comp[,list(region,nme_30,nme_5,nme_E14)],id.vars = "region")
ggplot(heart_comp_melt,aes(x=value,color=variable))+geom_density(size=1)+ggtitle("heart E15.5")

#midbrain
midbrain_30_15=import.bedGraph('../downstream/input/trim15/mm10_midbrain_day15_5_all_allele_agnostic_nme.bedGraph')
start(midbrain_30_15)=start(midbrain_30_15)-1

midbrain_comp=data.table(region=paste0(seqnames(midbrain_30_15),":",start(midbrain_30_15),'-',end(midbrain_30_15)),
                         nme_30=midbrain_30_15$score)
midbrain_comp=midbrain_comp[region%in%rownames(nme)]
midbrain_comp$nme_5=nme[midbrain_comp$region,grepl("midbrain-E15.5-all",colnames(nme))]
midbrain_comp$nme_E14=nme[midbrain_comp$region,grepl("midbrain-E14.5-all",colnames(nme))]
midbrain_comp$relative_diff=midbrain_comp$nme_30-midbrain_comp$nme_5
ggplot(midbrain_comp,aes(x=relative_diff))+geom_density()
midbrain_comp_melt=melt.data.table(midbrain_comp[,list(region,nme_30,nme_5,nme_E14)],id.vars = "region")
ggplot(midbrain_comp_melt,aes(x=value,color=variable))+geom_density(size=1)+ggtitle("midbrain E15.5")

#hindbrain

hindbrain_comp=data.table(region=rownames(nme),nme_E16=nme[,"hindbrain-E16.5-all"],nme_E15=nme[,"hindbrain-E15.5-all"])
hindbrain_comp_melt=melt.data.table(hindbrain_comp[,list(region,nme_E15,nme_E16)],id.vars = "region")
ggplot(hindbrain_comp_melt,aes(x=value,color=variable))+geom_density(size=1)+ggtitle("hindbrain")


#heart E14.5
heart_30_15=import.bedGraph('../downstream/input/trim15/mm10_heart_day14_5_all_allele_agnostic_nme.bedGraph')
start(heart_30_15)=start(heart_30_15)-1

heart_comp=data.table(region=paste0(seqnames(heart_30_15),":",start(heart_30_15),'-',end(heart_30_15)),
                      nme_30=heart_30_15$score)
heart_comp=heart_comp[region%in%rownames(nme)]
heart_comp$nme_5=nme[heart_comp$region,grepl("heart-E14.5-all",colnames(nme))]
heart_comp$nme_E15=nme[heart_comp$region,grepl("heart-E15.5-all",colnames(nme))]
heart_comp$relative_diff=heart_comp$nme_30-heart_comp$nme_5
ggplot(heart_comp,aes(x=relative_diff))+geom_density()
heart_comp_melt=melt.data.table(heart_comp[,list(region,nme_30,nme_5)],id.vars = "region")
ggplot(heart_comp_melt,aes(x=value,color=variable))+geom_density(size=1)+ggtitle("heart E14.5")
