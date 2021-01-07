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
dMML_cor=readRDS('../downstream/input/dmmlcor.rds')
dNME_cor=readRDS('../downstream/input/dNMEcor.rds')
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
