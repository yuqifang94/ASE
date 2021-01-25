rm(list=ls())
source("mainFunctions_sub.R")
theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()

# Repeat NME vs CpG density analysis --------------------------------------
#Get CpG location for mm10
CpG_mm10=getCpgSitesmm10()
NME_in=readRDS('../downstream/input/NME_agnostic_mouse_all_merged.rds')
NME_in=NME_in[!grepl('-P0-',NME_in$Sample)]
NME_in=NME_in[NME_in$N>=2]
dnase<- readRDS('../downstream/input/mm10_DNase.rds')
NME_in=subsetByOverlaps(NME_in,dnase,type='equal')
NME_in_gr=unique(granges(NME_in))
gr_seq=getSeq(Mmusculus,NME_in_gr,as.character=T)
NME_in_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
olap=findOverlaps(NME_in,NME_in_gr,type='equal')
NME_in$CG_exp=NME_in_gr$CGcont_exp[subjectHits(olap)]               
NME_in$CpG_number=countOverlaps(NME_in,CpG_mm10)
saveRDS(NME_in,'../downstream/input/NME_agnostic_mouse_all_merged_DNase_CpG.rds')
NME_in=readRDS('../downstream/input/NME_agnostic_mouse_all_merged_DNase_CpG.rds')
NME_in=NME_in[NME_in$CpG_number>0]
NME_in$density=NME_in$CpG_number/NME_in$CG_exp
NME_in$density_quant=findInterval(NME_in$density,seq(0,1,0.1))
#NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
NME_in$density_quant=factor(quant_conv[NME_in$density_quant],levels=quant_conv)

pdf('../downstream/output/graphs/FigureS12/CpG_density_NME_boxplot_CG_exp_mouse.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in)),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off() 
cor.test(NME_in$density,NME_in$NME)


# MAV vs NME -------------------------------------------------
NME_in_dt=readRDS('../downstream/output/NME_in_limb_ENOCD3_imputed.rds')
NME_in_dt=NME_in_dt[(!is.na(hyper_var)&hyper_var!=-100)]

# matrix and quantile normalization ---------------------------------------
#No need quantile normalization since we calculate correlation for each sample it for each sample
# hyper_var_dc=matrix_conv(NME_in_dt,"hyper_var")
# NME_dc=matrix_conv(NME_in_dt,"NME")
# NME_dc=NME_dc[rowSums(is.na(NME_dc))==0,]
# hyper_var_dc=hyper_var_dc[rowSums(is.na(hyper_var_dc))==0,]
# rn=intersect(rownames(NME_dc),rownames(hyper_var_dc))
# NME_dc=NME_dc[rn,]
# hyper_var_dc=hyper_var_dc[rn,]
# #Test quantile normalization
# hyper_var_dc_nm=normalize.quantiles(hyper_var_dc)
# rownames(hyper_var_dc_nm)=rownames(hyper_var_dc)
# colnames(hyper_var_dc_nm)=colnames(hyper_var_dc)
# #After quantile normalization, check plot
# hyper_var_dc_nm_dt=data.table(region=rownames(hyper_var_dc_nm))
# hyper_var_dc_nm_dt=cbind(hyper_var_dc_nm_dt,as.data.table(hyper_var_dc_nm))
# hyper_var_dc_nm_dt=melt.data.table(hyper_var_dc_nm_dt,id.vars="region",variable.name = "stage",value.name = "hyper_var")
# #After
# ggplot(hyper_var_dc_nm_dt,aes(x=hyper_var,color=stage))+geom_density(size=1)+theme(legend.position = "bottom")
# #Before
# ggplot(NME_in_dt,aes(x=hyper_var,color=stage))+geom_density(size=1)+theme(legend.position = "bottom")
# ggplot(NME_in_dt,aes(x=NME,color=stage))+geom_density(size=1)+theme(legend.position = "bottom")
# #assign to orignal values
# NME_in_dt$hyper_var=NULL
# NME_in_dt$hyper_var=hyper_var_dc_nm_dt[match(paste0(NME_in_dt$region,NME_in_dt$stage),paste0(region,stage))]$hyper_var
# NME_in_dt=NME_in_dt[!is.na(hyper_var)]
# saveRDS(NME_in_dt,'../downstream/output/NME_in_dt_limb_ENCODE_C1_nrom.rds')
# NME_in_dt=readRDS('../downstream/output/NME_in_dt_limb_ENCODE_C1_nrom.rds')

dist_plot_run(NME_in_dt,theme_glob,ylab="NME",stat_in="hyper_var",dir='../downstream/output/graphs/FigureS13/')
dist_plot_run(NME_in_dt,theme_glob,ylab="NME",stat_in="var",dir='../downstream/output/graphs/FigureS13/')
dist_plot_run(NME_in_dt,theme_glob,ylab="NME",stat_in="mean",dir='../downstream/output/graphs/FigureS13/')
# motif preprocessing for Ken ----------------------------------------------------------
#See mouse_motif_processing.R
# add information to Ken's list -------------------------------------------
#Motif prefer high ent
motif_human=fread('../downstream/output/graphs/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
shared_motif=fread('../downstream/input/mouse_motif_enriched_Ken/perfer_high_NME_overlap_motif_dNME.csv')$V2[-1]
shared_motif=data.table(shared_motif=shared_motif,prob_human= motif_human[match(shared_motif,motif_human$TF)]$Proportion)
ken_dir='../downstream/input/mouse_motif_enriched_Ken/'
file_in=dir(ken_dir,pattern="OR_residual.csv")
for(fn in file_in){
  Ken_in=fread(paste0(ken_dir,fn))
  tissue=gsub("_OR_residual.csv",'',fn)
  shared_motif[[tissue]]=Ken_in[match(shared_motif$shared_motif,gsub('.*_','',Ken_in$V1))]$residual     

}
plot_motif_binding<-function(motif,tissue,UC_raw,mml,nme,CpG_mm10){
  motif_tissue=data.table()
  #target_regions_cluster_all=GRanges()
  target_regions=readRDS(paste0('../downstream/input/motif_target/',tissue,'_motif_site_dNME.rds'))
  names(target_regions)=gsub('.*_','',names(target_regions))
  #chromHMM_in=readRDS('../downstream/output/chromHMM_enhancer.rds')
  #chromHMM_ts=chromHMM_in[chromHMM_in$tissue==tissue]
  region_in=fread(paste0('../downstream/input/mm10_cluster_all/',tissue,'.csv'))
  region_in=region_in[chromHMM_enhancer==TRUE]
  region_in_all=data.table()
  for(i in 1:10){
    
    if(length(motif)>0){
      #read in GO result
      GO_in=fread(paste0('../downstream/output/mm10_result/chromHMM_enhancer/cluster_GO/mm10_cluster_chromHMM/',tissue,'-',i,'_cluster_GO.csv'))
      GO_in=GO_in[FC>=1.5&FDR<=0.1][1:5]
      GO_in_gene=unique(unlist(strsplit(GO_in$gene,";")))
      #read in all regions
      region_in_clu=region_in[cluster==i]
      target_regions_cluster=do.call(c,unlist(lapply(motif,function(x){
        target_regions_out=unique(target_regions[[x]][target_regions[[x]]$cluster==i])
        if(length(target_regions_out)>0){
        target_regions_out$motif=x
        
        return(target_regions_out)}
        })))
      #target_regions_cluster=target_regions_cluster[target_regions_cluster$gene %in% GO_in_gene]
      olap=findOverlaps(convert_GR(region_in_clu$region),target_regions_cluster)
      region_in_clu=region_in_clu[queryHits(olap)]
      #region_in_clu$gene_Jason=region_in_clu$gene
      region_in_clu$gene_top_GO=region_in_clu$gene %in% GO_in_gene
      #region_in_clu$gene=NULL
      region_in_clu=cbind(region_in_clu,as.data.table(mcols(target_regions_cluster[subjectHits(olap)])))
      region_in_all=rbind(region_in_all,region_in_clu)
    }
    
    
  }
  
  timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
  UC=UC_raw[[tissue]][,colnames(UC_raw[[tissue]])%in% timeorder]
  UC <-UC[unique(region_in_all$region),order(match(colnames(UC),timeorder))]
  stat_differential<-function(stat_in,UC,tissue){
    
    diff_out <- sapply(colnames(UC),function(i) {
      time <- strsplit(i,'-')
      sapply(time,function(time_in)
      {abs(stat_in[rownames(UC),paste0(tissue,'-',time_in[1],'-all')]-stat_in[rownames(UC),paste0(tissue,'-',time_in[2],'-all')])})
    })
    colnames(diff_out)=colnames(UC)
    rownames(diff_out)=rownames(UC)
    return(diff_out)
  }
  dmml=stat_differential(mml,UC,tissue)
  dnme=stat_differential(nme,UC,tissue)
  
  region_in_all_region=region_in_all[order(dNME_maxJSD,decreasing=T),list(motif=paste(unique(motif),collapse = ';')),
                                     by=list(region,gene,distance,dNME_maxJSD,gene_top_GO)]
  region_in_all_region$N=countOverlaps(convert_GR(region_in_all_region$region),CpG_mm10)
  theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                   axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                   axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                   axis.text.x=element_text(size=16),
                   axis.text.y=element_text(size=16))+theme_classic()
  dmml_out=dmml[region_in_all$region,]
  colnames(dmml_out)=paste0('dMML-',colnames(dmml_out))
  dnme_out=dnme[region_in_all$region,]
  colnames(dnme_out)=paste0('dNME-',colnames(dnme_out))
  region_in_all=cbind(region_in_all,as.data.table(dmml_out),as.data.table(dnme_out))
  pdf(paste0('../downstream/output/graphs/Figure7/',tissue,'_motif_dNME_only.pdf'),width=5,height=5)
  for(rg in region_in_all_region$region){
    region_stat=data.table(stage=colnames(UC),
                           #UC=scale(UC[rg,])[,1],
                           #dmml=scale(dmml[rg,])[,1],
                           dnme=dnme[rg,]
    )
    #nme_cor=cor(region_stat$UC,region_stat$dnme,method='spearman')
    #mml_cor=cor(region_stat$UC,region_stat$dmml,method='spearman')
    # region_stat=melt.data.table(region_stat,id.var='stage',variable.name = "stat")
    # 
    # print(ggplot(region_stat,aes(x=stage,y=value,group=stat,color=stat))+geom_point(size=1)+geom_line(stat = "identity")+
    #         theme_glob+theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom",plot.title = element_text(hjust = 0.5))+
    #         scale_color_manual(values=c(UC="black", dmml="blue", dnme="red"))+ylab('Scaled value')+
    #         ggtitle(paste0('dnme cor:',round(nme_cor,digits = 3),'\n',
    #                        'dmml cor:',round(mml_cor,digits = 3),'\n',
    #                        paste(region_in_all_region[region==rg,],collapse='\n'))))
    
    print(ggplot(region_stat,aes(x=stage,y=dnme))+geom_point(size=1)+geom_line(stat = "identity",aes(group=1))+
            theme_glob+theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom",plot.title = element_text(hjust = 0.5))+
            ylab('dNME')+
            ggtitle(paste(region_in_all_region[region==rg,],collapse='\n')))
  }
  dev.off()
  return(list(region_in_all_region,region_stat))
}
UC_raw=readRDS('../downstream/output/uc_matrix_DNase.rds')
UC_raw_control=readRDS('../downstream/output/uc_matrix_.rds')
mean(unlist(lapply(UC_raw,function(x) quantile(x,prob=0.95))))
mml <- readRDS('../downstream/output/mml_matrix_DNase.rds')
nme <- readRDS('../downstream/output/nme_matrix_DNase.rds')
CpG_mm10=getCpgSitesmm10()

heart_regions=plot_motif_binding(shared_motif$shared_motif,"heart",UC_raw,mml,nme,CpG_mm10)

#saveRDS(heart_regions,'../downstream/output/heart_regions.rds')

limb_regions=plot_motif_binding(shared_motif$shared_motif,"limb",UC_raw,mml,nme,CpG_mm10)
#saveRDS(limb_regions,'../downstream/output/limb_regions.rds')
forebrain_regions=plot_motif_binding(shared_motif$shared_motif,"forebrain",UC_raw,mml,nme,CpG_mm10)
#saveRDS(forebrain_regions,'../downstream/output/forebrain_regions.rds')
# Getting dNME in each targeted motif -------------------------------------
dnme=readRDS('../downstream/output/dnme_matrix_DNase.rds')
motif_target_dir='../downstream/input/motif_target/'
nme <- readRDS('../downstream/output/nme_matrix_DNase.rds')
nme=nme[,!grepl("P0",colnames(nme))]
motif_dir_human_high_ent=fread('../downstream/input/Ken_agnostic_motif_accessible_ent.txt',header = F)
motif_dir_human_high_ent$TF=gsub('.*_','',motif_dir_human_high_ent$V1)
motif_dir_human_low_ent=fread('../downstream/input/Ken_agnostic_motif_accessible_no_ent.txt',header = F)
motif_dir_human_low_ent$TF=gsub('.*_','',motif_dir_human_low_ent$V1)
# motif_dir_human_high_ent=fread('../downstream/output/graphs/tableS1_motif_prefer_ent_OMIM.csv')
# motif_dir_human_high_ent=motif_dir_human_high_ent[qval_binom<=0.05]
# motif_dir_human_low_ent=fread('../downstream/output/graphs/table3_motif_not_prefer_ent_OMIM.csv')
# motif_dir_human_low_ent=motif_dir_human_low_ent[qval_binom<=0.05]
for (fn in dir(motif_target_dir)){
motif_target_in=readRDS(paste0(motif_target_dir,fn))
names(motif_target_in)=sub('.*_','',names(motif_target_in))
#For each motif_check dNME
tissue=sub('_.*','',fn)
cat('processing:',tissue,'\n')
#for each motif, it's enriched in some cluster, 
#we want to check in that cluster dNME_maxJSD is higher than regions not overlapping that motif
#Or if that's higher than motif showing no dNME preference
# motif_target_in_subset=motif_target_in
# names(motif_target_in_subset)=NULL
# motif_target_in_subset=unique(do.call(c,motif_target_in_subset))
# nme_motif=nme[queryHits(findOverlaps(convert_GR(rownames(nme)),motif_target_in_subset)),]
nme_tissue=nme[,grepl(tissue,colnames(nme))]
clu_sig=do.call(rbind,lapply(1:10,function(x){
  csv_in=fread(paste0('../downstream/input/mouse_motif_cluster/',tissue,'/motif_',tissue,'_cluster_',x,'_enhancer.csv'))
  csv_in$cluster=x
  return(csv_in)
}))
#clu_sig=clu_sig[FDR<=0.05&odds_ratio>=1.2]
clu_sig$motif=gsub('.*_','',clu_sig$motif)

#use all regions, not only enhancer maybe
NME_motif_target=mclapply(names(motif_target_in),function(x){
  #make a data.table, showing cluster, mean dNME, mean dMML at UC
  if(x %in% clu_sig$motif){
    motif_target=motif_target_in[[x]][motif_target_in[[x]]$cluster %in% clu_sig[motif ==x]$cluster]
    nme_motif=nme_tissue[queryHits(findOverlaps(convert_GR(rownames(nme_tissue)),motif_target)),]
    motif_out=as.data.table(nme_motif)
    motif_out$regions=rownames(nme_motif)
    motif_out$motif=x
    motif_out=melt.data.table(motif_out,id.vars=c('regions','motif'),variable.name = 'tissue_stage',value.name = "NME")
   return(motif_out)
  }
},mc.cores=20)
NME_motif_target=fastDoCall('rbind',NME_motif_target)

NME_motif_target$motif_type="Not in human analysis"
NME_motif_target[motif %in% motif_dir_human_low_ent$TF]$motif_type="Decreased NME in accessible region"
NME_motif_target[motif %in% motif_dir_human_high_ent$TF]$motif_type="Increased NME in accessible region"
NME_motif_target_mean=NME_motif_target[,list(mean_nme=median(NME,na.rm=T)),by=list(motif_type,motif)]
pdf(paste0('../downstream/output/graphs/Figure6/motif_accessible/mean_NME_',tissue,'.pdf'),width=3.5,height = 7)
print(ggplot(NME_motif_target_mean,aes(x=motif_type,y=mean_nme))+geom_boxplot()+geom_text(aes(label=motif),size=1)+xlab('')+ylab("NME")+
  theme_glob+theme(axis.text.x = element_text(angle = 90)))
dev.off()

pdf(paste0('../downstream/output/graphs/Figure6/motif_accessible/all_NME_',tissue,'.pdf'),width=3.5,height = 7)
print(ggplot(NME_motif_target,aes(x=motif_type,y=NME))+geom_boxplot()+xlab('')+ylab("NME")+
  theme_glob+theme(axis.text.x = element_text(angle = 90)))
dev.off()
#dNME at UC
dnme_tissue=fread(paste0('../downstream/input/mm10_cluster/',tissue,'.csv'))
dNME_motif_target=mclapply(names(motif_target_in),function(x){
  #make a data.table, showing cluster, mean dNME, mean dMML at UC
  if(x %in% clu_sig$motif){
    clusters=clu_sig[motif ==x]$cluster
    dnme_motif_out=data.table()
    for (clu in clusters){
      motif_target=motif_target_in[[x]][motif_target_in[[x]]$cluster ==clu]
      dnme_tissue_clu=dnme_tissue[dnme_tissue$cluster==clu]
      dnme_motif_clu=dnme_tissue_clu[queryHits(findOverlaps(convert_GR(dnme_tissue_clu$region ),motif_target))]
      dnme_motif_out=rbind(dnme_motif_out,data.table(dNME=dnme_motif_clu$dNME_maxJSD,region=dnme_motif_clu$region,cluster=clu))
    }

    dnme_motif_out$motif=x
    return(dnme_motif_out)
  }
},mc.cores=20)
dNME_motif_target=fastDoCall('rbind',dNME_motif_target)
dNME_motif_target$motif_type="Not in human analysis"
dNME_motif_target[motif %in% motif_dir_human_low_ent$TF]$motif_type="Decreased NME in accessible region"
dNME_motif_target[motif %in% motif_dir_human_high_ent$TF]$motif_type="Increased NME in accessible region"
dNME_motif_target_mean=dNME_motif_target[,list(mean_dnme=mean(dNME,na.rm=T)),by=list(motif_type,motif)]
pdf(paste0('../downstream/output/graphs/Figure6/motif_accessible/motif_dNME_',tissue,'.pdf'),width=3.5,height = 7)
print(ggplot(dNME_motif_target_mean,aes(x=motif_type,y=mean_dnme))+geom_boxplot()+xlab('')+ylab("dNME")+
        geom_text(aes(label=motif),size=1)+theme_glob+theme(axis.text.x = element_text(angle = 90)))
dev.off()

pdf(paste0('../downstream/output/graphs/Figure6/motif_accessible/all_dNME_',tissue,'.pdf'),width=3.5,height = 7)
print(ggplot(dNME_motif_target,aes(x=motif_type,y=dNME))+geom_boxplot()+xlab('')+ylab("dNME")+
  theme_glob+theme(axis.text.x = element_text(angle = 90)))
dev.off()

}

# prepare regions for Ken -------------------------------------------------
#Prepare for Jason
UC_raw=readRDS('../downstream/output/uc_matrix_DNase.rds')
cutoffs=do.call(rbind,lapply(UC_raw,function(x) quantile(x,prob=c(0.7,0.8,0.9,0.95,0.975,0.99,0.999,0.9999))))
colMeans(cutoffs)

#prediction
cor_dt_classify=list()
lv=0.50
FDR_cor_cutoff=0.25
for(cor_in in cor_dt){

  #Use selected regions to predict
  model <- lm(dNME_cor~dMML_cor+0, data = cor_in)
  predict_out=predict(model, newdata =cor_in, interval = "prediction",level=lv)
  colnames(predict_out)=paste0('dNME_',colnames(predict_out))
  cor_in=cbind(cor_in,predict_out)
  cor_in_fitting=data.table(dMML_cor=unique(round(cor_in$dMML_cor,digits = 3)))
  cor_in_fitting=cbind(cor_in_fitting,predict(model, newdata =cor_in_fitting, interval = "prediction",level=lv))
  cor_in[dNME_upr>1]$dNME_upr=1
  cor_in[dNME_lwr<(-1)]$dNME_lwr=-1
  tissue=unique(cor_in$tissue)
  dNME_cutoff=min(cor_in[dNME_FDR<=FDR_cor_cutoff]$dNME_cor)
  dMML_cutoff=min(cor_in[dMML_FDR<=FDR_cor_cutoff]$dMML_cor)
  cor_in$significant="NA"
  cor_in[(dNME_FDR>FDR_cor_cutoff&dMML_FDR>FDR_cor_cutoff)]$significant="Not Significant"
  cor_in[(dNME_FDR<=FDR_cor_cutoff|dMML_FDR<=FDR_cor_cutoff)&(dNME_cor>dNME_lwr&dNME_cor<dNME_upr)]$significant="Correlated One significant"
  cor_in[(dNME_FDR<=FDR_cor_cutoff&dMML_FDR<=FDR_cor_cutoff)&(dNME_cor>dNME_lwr&dNME_cor<dNME_upr)]$significant="Correlated both significant"
  cor_in[(dNME_FDR<=FDR_cor_cutoff)&(dNME_cor>=dNME_upr)]$significant="dNME correlated"
  cor_in[(dMML_FDR<=FDR_cor_cutoff)&(dNME_cor<=dNME_lwr)]$significant="dMML correlated"
  sig_levels=c("Not Significant",
               "Correlated One significant",
               "Correlated both significant",
               "dNME correlated",
               "dMML correlated")
  cor_in$significant=factor(cor_in$significant,levels=sig_levels)

  col_theme=c("grey","purple","green","red","blue")
  png(paste0('../downstream/output/dmml_dnme_cor/',tissue,'_dNME_dMML_dots_0_inter_',lv,'.png'),width=7,height=7,units = "in",res=96)
  print(ggplot(cor_in,aes(x=dMML_cor,y=dNME_cor))+
          geom_point(aes(color=significant),alpha=0.1)+
          theme(legend.position = "bottom")+ 
          scale_color_manual(breaks= sig_levels, values = col_theme)+
          geom_line(data=cor_in_fitting,aes(x=dMML_cor,y=fit),size=1)+
          geom_line(data=cor_in_fitting,aes(x=dMML_cor,y=lwr),size=1,linetype='dashed')+
          geom_line(data=cor_in_fitting,aes(x=dMML_cor,y=upr),size=1,linetype='dashed')+
         geom_hline(yintercept = dNME_cutoff,size=1,color="black")+geom_vline(xintercept = dMML_cutoff,size=1,color="black"))
 
  
  dev.off()
  cor_dt_classify[[tissue]]=cor_in
}
#correlation cutoff
cor_dt_classify_plot=fastDoCall('rbind',cor_dt_classify)
pdf('../downstream/output/dmml_dnme_cor/bar_plot_0_inter_50.pdf')
ggplot(cor_dt_classify_plot,aes(x=tissue))+geom_bar(stat="count",aes(group=significant,fill=significant))+theme(legend.position = 'bottom')
dev.off()
mml <- readRDS('../downstream/output/mml_matrix_DNase.rds')
nme <- readRDS('../downstream/output/nme_matrix_DNase.rds')
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
for(tissue in names(UC_raw)){
  UC=UC_raw[[tissue]][,colnames(UC_raw[[tissue]])%in% timeorder]
  
  stat_differential<-function(stat_in,UC,tissue){
    
    diff_out <- sapply(colnames(UC),function(i) {
      time <- strsplit(i,'-')
      sapply(time,function(time_in)
      {abs(stat_in[rownames(UC),paste0(tissue,'-',time_in[1],'-all')]-stat_in[rownames(UC),paste0(tissue,'-',time_in[2],'-all')])})
    })
    colnames(diff_out)=colnames(UC)
    rownames(diff_out)=rownames(UC)
    return(diff_out)
  }
  dmml=stat_differential(mml,UC,tissue)
  dnme=stat_differential(nme,UC,tissue)
  #Merge them into single data table
  dmml_dt=as.data.table(scalematrix(dmml))
  dmml_dt$region=rownames(dmml)
  dnme_dt=as.data.table(scalematrix(dnme))
  dnme_dt$region=rownames(dnme)
  UC_dt=as.data.table(UC)
  UC_dt$region=rownames(UC)
  dmml_dt=melt.data.table(dmml_dt,id='region')
  dnme_dt=melt.data.table(dnme_dt,id='region')
  UC_dt=melt.data.table(UC_dt,id='region')
  dmml_dt=dmml_dt[order(region,variable)]
  dnme_dt=dnme_dt[order(region,variable)]
  UC_dt=UC_dt[order(region,variable)]
  plot_dt=data.table(region=dmml_dt$region,
                     stage=dmml_dt$stage,
                     dmml=dmml_dt$value,
                     dnme=dnme_dt$value,
                     UC=UC_dt$value)
  plot_dt$dmml_cor=dMML_cor[[tissue]][plot_dt$region]
  plot_dt$dnme_cor=dNME_cor[[tissue]][plot_dt$region]
  plot_dt=plot_dt[!is.na(dmml_cor)]
  # jpeg('../downstream/output/dmml_UC_scale2.jpeg')
  # ggplot(plot_dt[dmml_cor>=0.7&region=="chr10:100015944-100016193"],aes(x=dmml,y=UC))+geom_point(alpha=1)
  # dev.off()
  png(paste0('../downstream/output/dmml_dnme_cor/dnme_UC_scale_cor',tissue,'.png'))
  print(ggplot(plot_dt[dnme_cor>=0.7],aes(x=dnme,y=UC))+geom_point(alpha=0.01)+geom_smooth())
  dev.off()
  png(paste0('../downstream/output/dnme_UC_scale_all',tissue,'.png'))
  print(ggplot(plot_dt,aes(x=dnme,y=UC))+geom_point(alpha=0.01)+geom_smooth())
  dev.off()
  png(paste0('../downstream/output/dmml_UC_scale_cor',tissue,'.png'))
  print(ggplot(plot_dt[dmml_cor>=0.7],aes(x=dmml,y=UC))+geom_point(alpha=0.01)+geom_smooth())
  dev.off()
  png(paste0('../downstream/output/dmml_UC_scale_all',tissue,'.png'))
  print(ggplot(plot_dt,aes(x=dmml,y=UC))+geom_point(alpha=0.01)+geom_smooth())
  dev.off()
}

dir_in='../downstream/input/mm10_cluster_all/'
region_out=list()
for(fn in dir(dir_in,pattern='.csv')){
  csv_in=fread(paste0(dir_in,fn))
  csv_in=csv_in[,list(cluster,region ,chromHMM_enhancer,dNME_maxJSD,dMML_maxJSD)]
  tissue=gsub('.csv','',fn)
  csv_in$mml_cor=dMML_cor[[tissue]][csv_in$region]
  csv_in$nme_cor=dNME_cor[[tissue]][csv_in$region]
  csv_in$cor_type="Not assigned"
  csv_in[nme_cor>=0.7&mml_cor<0.7]$cor_type="dNME_only"
  csv_in[mml_cor>=0.7&nme_cor<0.7]$cor_type="dMML_only"
  csv_in[mml_cor<0.7&nme_cor<0.7]$cor_type="Not_dNME_dMML"
  csv_in[nme_cor>=0.7&mml_cor>=0.7]$cor_type="dMML_dNME"
  region_out[[tissue]]=csv_in
}
saveRDS(region_out,'../downstream/output/region_annotation.rds')

#Correlation setting:
#Try truncated mltivariate:tmvtnorm: A Package for the Truncated Multivariate Normal Distribution
#Getting fitted data and correlation outside the fitted data
dMML_cor=readRDS('../downstream/input/dmmlcor.rds')
dNME_cor=readRDS('../downstream/input/dNMEcor.rds')
dmml_perm=readRDS('../downstream/input/permudmmlcor.rds')
dnme_perm=readRDS('../downstream/input/permudnmecor.rds')
FDR_cor<-function(real_value,perm){
  cor_null=ecdf(perm)
  cor_pval=1-cor_null(real_value)
  print(length(perm))
  return(p.adjust(cor_pval,method='BH'))
}
convert_matrix<-function(perm_in,tissue,stat){
  perm_in_dt=as.data.table(perm_in)
  perm_in_dt$region=rownames(perm_in)
  perm_in_dt$tissue="tissue"
  perm_in_dt=melt.data.table(perm_in_dt,id.vars=c('region','tissue'))
  perm_in_dt$variable=paste0(stat,"_perm")
  return(perm_in_dt)
  
}
#Merge into data table
cor_dt=lapply(names(dNME_cor),function(x) {
  regions=names(dNME_cor[[x]])
  out_dt=data.table(region=regions,
                    dMML_cor=dMML_cor[[x]][regions],
                    dNME_cor=dNME_cor[[x]][regions],
                    tissue=x)
  dMML_perm_in=convert_matrix(dmml_perm[[x]],x,'dmml')
  dNME_perm_in=convert_matrix(dnme_perm[[x]],x,'dnme')
  out_dt$dMML_FDR=FDR_cor(out_dt$dMML_cor,dMML_perm_in$value)
  out_dt$dNME_FDR=FDR_cor(out_dt$dNME_cor,dNME_perm_in$value)
  dNME_cutoff=min(out_dt[dNME_FDR<=0.25]$dNME_cor)
  dMML_cutoff=min(out_dt[dMML_FDR<=0.25]$dMML_cor)
  out_dt_mt=melt.data.table(out_dt,id.vars = c('region','tissue'))
  out_dt_mt=rbind(out_dt_mt,dMML_perm_in,dNME_perm_in)
  # pdf(paste0('../downstream/output/dmml_dnme_cor/',x,'_dMML_permute_real.pdf'))
  # print(ggplot(out_dt_mt[variable%in%c("dMML_cor","dmml_perm")],aes(x=value,group=variable,color=variable))+geom_density(size=1)+
  #         theme(legend.position = "bottom")+geom_vline(xintercept = dMML_cutoff,size=1,color="purple"))
  # 
  # dev.off()
  # 
  # pdf(paste0('../downstream/output/dmml_dnme_cor/',x,'_dNME_permute_real.pdf'))
  # print(ggplot(out_dt_mt[variable%in%c("dNME_cor","dnme_perm")],aes(x=value,group=variable,color=variable))+geom_density(size=1)+
  #         theme(legend.position = "bottom")+geom_vline(xintercept = dNME_cutoff,size=1,color="purple"))
  # 
  # dev.off()
  return(out_dt)
  
  
})

library(tmvtnorm)
library(metR)
cor_dt=lapply(cor_dt,function(x) x[!grepl("X",region)&!grepl("Y",region)])

# density estimation ------------------------------------------------------

#Estimate bandwidth using ucv:
library(MASS)
tissue_in=cor_dt[[4]]
tissue_in_mt=as.matrix(tissue_in[,list(dMML_cor,dNME_cor)])
bw_opt=ucv(tissue_in_mt)

#Using mle on all data below FDR=0.25
set.seed(12345)
FDR_cor_cutoff=0.25
contour_data_out=list()
estimage_out=list()
x_seq <- seq(-1, 1, by=0.001)
y_seq <- seq(-1, 1, by=0.001)
for(cor_mv_norm in cor_dt){
  tissue=unique(cor_mv_norm$tissue)
  cor_mv_norm_trunc=cor_mv_norm[(dMML_FDR>FDR_cor_cutoff)&(dNME_FDR>FDR_cor_cutoff)]
  #cor_mv_norm_trunc=cor_mv_norm[(dMML_cor<=0)&(dNME_cor<=0)]
  #Use range instead of max due to decimal problem
  upper=c(range(cor_mv_norm_trunc$dMML_cor)[2],range(cor_mv_norm_trunc$dNME_cor)[2])
  #There's some digit problem?
  cor_mv_norm_trunc[dMML_cor<(-1)]$dMML_cor=-1
  cor_mv_norm_trunc[dNME_cor<(-1)]$dNME_cor=-1
  cor_mv_norm_trunc[dMML_cor>(1)]$dMML_cor=1
  cor_mv_norm_trunc[dNME_cor>(1)]$dNME_cor=1
  lower=c(-1,-1)
  trunc_matrix=as.matrix(cor_mv_norm_trunc[,list(dMML_cor,dNME_cor)])
  mle.fit1 <- mle.tmvnorm(trunc_matrix, lower=lower, upper=upper,method = "L-BFGS-B",lower.bounds=c(-1,-1),upper.bounds = c(1,1))
  mle.fit1_sum=summary(mle.fit1)
  
  mu= mle.fit1_sum@coef[c('mu_1','mu_2'),"Estimate"]
  sigma= matrix(mle.fit1_sum@coef[c('sigma_1.1','sigma_1.2',"sigma_1.2","sigma_2.2"),"Estimate"],nrow=2,byrow=T)

  contour_data=as.data.table(expand.grid(x_seq,y_seq))
  colnames(contour_data)=c("x_seq","y_seq")
  contour_data$MD=mahalanobis(as.matrix(contour_data[,.(x_seq,y_seq)]),center=mu,cov=sigma)
  contour_data$density=NA
  # for (i in 1:nrow(contour_data)){
  # 
  #     #prepare for contour
  #   contour_data$density[i] = dmvnorm(contour_data[i,.(x_seq,y_seq)], mean=mu,sigma=sigma)
  # 
  # }
  # jpeg(paste0('../downstream/output/dmml_dnme_cor/truncated_nrom_small_0_',tissue,'.jpg'))
  # print(ggplot(cor_mv_norm[(dMML_FDR>FDR_cor_cutoff)&(dNME_FDR>FDR_cor_cutoff)],aes(x=dMML_cor,y=dNME_cor))+ geom_point(alpha=0.01)+
  #   geom_vline(xintercept = upper[1],size=1)+geom_hline(yintercept = upper[2],size=1)+
  #   geom_contour(mapping=aes(x=x_seq,y=y_seq,z=MD),data=contour_data,breaks=seq(1,5,1),color='red')+ 
  #   geom_label_contour(mapping=aes(x=x_seq,y=y_seq,z=MD),data=contour_data,breaks=seq(1,5,1))+
  #   geom_point(mapping=aes(x=mu1,y=mu2),data=data.frame(mu1=mu[1],mu2=mu[2]),size=2,color='blue'))
  #  dev.off()
   contour_data_out[[tissue]]=contour_data
   estimage_out[[tissue]]=list(mu=mu,sigma=sigma)
   gc()
  #geom_contour_fill()+ scale_fill_gradient2(high = "blue",low="white")+
}
#Use distance cutoff: 4
MD_cutoff=4
FDR_cutoff=0.25
cor_mv_norm_MD_out=list()
col_theme=c("red","blue","green","purple")
sig_levels=c("dMML_only","Both","dNME_only","None")
for(cor_mv_norm_MD in cor_dt){
  tissue=unique(cor_mv_norm_MD$tissue)
  estimage_tissue=estimage_out[[tissue]]
  cor_mv_norm_MD$MD=mahalanobis(as.matrix(cor_mv_norm_MD[,.(dMML_cor,dNME_cor)]),center=estimage_tissue$mu,cov=estimage_tissue$sigma)
  cor_mv_norm_MD$cor_type="NA"
  cor_mv_norm_MD[MD<=MD_cutoff&(dMML_FDR<=FDR_cutoff|dNME_FDR<=FDR_cutoff)]$cor_type="Both"
  cor_mv_norm_MD[(dMML_FDR>FDR_cutoff&dNME_FDR>FDR_cutoff)]$cor_type="None"
  cor_mv_norm_MD[MD>MD_cutoff&dNME_FDR<=FDR_cutoff]$cor_type="dNME_only"
  cor_mv_norm_MD[MD>MD_cutoff&dMML_FDR<=FDR_cutoff]$cor_type="dMML_only"
  dNME_cutoff=min(cor_mv_norm_MD[dNME_FDR<=FDR_cutoff]$dNME_cor)
  dMML_cutoff=min(cor_mv_norm_MD[dMML_FDR<=FDR_cutoff]$dMML_cor)

  cor_mv_norm_MD$cor_type=factor(cor_mv_norm_MD$cor_type,levels=sig_levels)
  

  png(paste0('../downstream/output/dmml_dnme_cor/cor_type_',tissue,'_025_4_md.png'),width=7,height=7,units = "in",res=96)
  print(
    ggplot(cor_mv_norm_MD,aes(x=dMML_cor,y=dNME_cor))+
          geom_point(aes(color=cor_type),alpha=0.1)+
          theme(legend.position = "bottom")+ 
          scale_color_manual(breaks= sig_levels, values = col_theme)+
          geom_hline(yintercept = dNME_cutoff,size=1,color="black")+geom_vline(xintercept = dMML_cutoff,size=1,color="black")+
          geom_contour(mapping=aes(x=x_seq,y=y_seq,z=MD),data=contour_data_out[[tissue]],breaks=MD_cutoff,color='black')
    )
  
dev.off()
cor_mv_norm_MD_out[[tissue]]=cor_mv_norm_MD
}
 
cor_mv_norm_MD_out=fastDoCall('rbind',cor_mv_norm_MD_out)
cor_mv_norm_MD_out_plot=cor_mv_norm_MD_out[,.N,by=c("tissue","cor_type")][,list(prop=N/sum(N),cor_type=cor_type),by=tissue]
pdf('../downstream/output/dmml_dnme_cor/bar_plot_025_4_trunc.pdf')
ggplot(cor_mv_norm_MD_out_plot,aes(x=tissue,y=prop))+geom_bar(stat="identity",aes(group=cor_type,fill=cor_type))+theme(legend.position = 'bottom')+
  scale_fill_manual(breaks= sig_levels, values = col_theme)
dev.off()
library(fpc)
for(i in 1:length(cor_dt)){
  df=cor_dt[[i]][,.(dMML_cor,dNME_cor)]
  tissue=unique(cor_dt[[i]]$tissue)
  for(eps in seq(0.001,0.01,0.001)){
    png(paste0('../downstream/output/dmml_dnme_cor/DBscan_',eps,'_',tissue,'.png'),width=7,height=7,units = "in",res=96)
    print(plot(fpc::dbscan(df, eps = eps, MinPts = 10),df, main = "DBSCAN", frame = FALSE))
    dev.off()
    
  }
}
library(mclust)
cor_dt_mc_out_all=list()
set.seed(123)
for(i in 1:length(cor_dt)){
  cor_dt_ts=cor_dt[[i]]
  cor_dt_mt=as.matrix(cor_dt_ts[,.(dMML_cor,dNME_cor)])
  tissue=unique(cor_dt_ts$tissue)
  cor_dt_mc_out=list()
  for(g in 3:12){
    cor_dt_mc_out[[as.character(g)]]=Mclust(cor_dt_mt,G=g)
    png(paste0('../downstream/output/dmml_dnme_cor/mcluster_',tissue,'_',g,'.png'),width=7,height=7,units = "in",res=96)
    plot( cor_dt_mc_out[[as.character(g)]], what = 'classification')
    dev.off()
    
  }
  cor_dt_mc_out_all[[tissue]]=cor_dt_mc_out
}

