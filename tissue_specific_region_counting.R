 source('mainFunctions_sub.R')
 UC_in=readRDS(UC_in_matrix_cluster_file)
 cut=0.1
 #Non tissue-specific high UC
aid <- sapply(names(UC_in),function(i) {
 names(which(rowSums(UC_in[[i]] > cut) > 0))
  })  
 ts_aid <- sapply(names(UC_in),function(i) {
    #This is for one and only one
    return(setdiff(aid[[i]],unlist(aid[names(aid)!=i])))
    
})
ts_aid_dt_fn='../downstream/output/mouse_analysis/tissue_specific_enhancer/ts_aid.rds'
saveRDS(ts_aid,ts_aid_dt_fn)
length(unique(unlist(aid)))
length(unique(unlist(ts_aid)))
aid_tb_count=table(unlist(aid))
aid_tb_count_dt=data.table(region=aid_tb_count)
colnames(aid_tb_count_dt)=c('region','n_tissue')
pdf(paste0(figure_path,'n_tissue_pass_filter_ts.pdf'))
for(ts in names(aid)){
    aid_tb_count_dt_ts=data.table(count=table(aid_tb_count_dt[region %in% aid[[ts]]]$n_tissue))
    colnames(aid_tb_count_dt_ts)=c('n_tissue','count')
    aid_tb_count_dt_ts$prop=aid_tb_count_dt_ts$count/sum(aid_tb_count_dt_ts$count)
    print(ggplot(aid_tb_count_dt_ts,aes(x=n_tissue,y=prop))+geom_bar(stat='identity', fill="steelblue")+
    xlab('Number of tissue with regions having UC>=0.1')+ylab('Proportion of regions')+
    geom_text(aes(label=round(prop,digits=2)),vjust=-1.6)+ggtitle(ts))



}

dev.off()
aid_tb=table(table(unlist(aid)))
#    1      2      3      4      5      6      7
# 588189 281624 192486 202206 278883 324387 218109
aid_tb_dt=data.table(n_region=aid_tb)
colnames(aid_tb_dt)=c('n_tissue','n_region')
aid_tb_dt$n_region_pro=aid_tb_dt$n_region/sum(aid_tb_dt$n_region)
aid_tb_dt$n_tissue=factor(aid_tb_dt$n_tissue,levels=as.character(1:10))
pdf(paste0(figure_path,'number_tissue_pass_filter.pdf'))
ggplot(aid_tb_dt,aes(x=n_tissue,y=n_region_pro))+geom_bar(stat='identity', fill="steelblue")+
xlab('Number of tissue with regions having UC>=0.1')+ylab('Proportion of regions')+
geom_text(aes(label=round(n_region_pro,digits=2)),vjust=-1.6)
dev.off()
#How do they overlap enhancers
aid_gr=lapply(aid,convert_GR,direction='GR')
enhancer_bin=readRDS(bin_enhancer_rds)
aid_gr_enhancer=lapply(aid_gr,function(x) subsetByOverlaps(x,enhancer_bin))
aid_gr_enhancer_dt=lapply(aid_gr_enhancer,function(x) convert_GR(x,direction="DT")$region)

aid_tb_enhancer=table(table(unlist(aid_gr_enhancer_dt)))
#  1     2     3     4     5     6     7
#26436 13077  8453  6645  7896  9589  7659
aid_tb_dt_enhancer=data.table(n_region=aid_tb_enhancer)
colnames(aid_tb_dt_enhancer)=c('n_tissue','n_region')
aid_tb_dt_enhancer$n_region_pro=aid_tb_dt_enhancer$n_region/sum(aid_tb_dt_enhancer$n_region)
aid_tb_dt_enhancer$n_tissue=factor(aid_tb_dt_enhancer$n_tissue,levels=as.character(1:10))
pdf(paste0(figure_path,'number_tissue_pass_filter_enhancer.pdf'))
ggplot(aid_tb_dt_enhancer,aes(x=n_tissue,y=n_region_pro))+geom_bar(stat='identity', fill="steelblue")+
xlab('Number of tissue with regions having UC>=0.1')+ylab('Proportion of regions')+
geom_text(aes(label=round(n_region_pro,digits=2)),vjust=-1.6)
dev.off()


#Count different dMML regions
#read in dMML
source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
dMML_in=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=20)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=20)
# dMML_vec=unlist(mclapply(dMML_in,function(x) as.vector(x),mc.cores=20))
# quantile(,prob=0.95,na.rm=T)#0.90 quant: 0.109, 0.95 quant:0.149
# 
# UC_in_vec=unlist(mclapply(UC_in_only,function(x) as.vector(x),mc.cores=20))
# quantile(UC_in_vec,prob=0.95,na.rm=T)#0.95quant: 0.1284, 0.90 quant: 0.066
# quantile(unlist(mclapply(UC_in,function(x) as.vector(x[,grepl('dNME',colnames(x))]),mc.cores=20)),prob=0.9)

 #Non tissue-specific high UC
 cut=0.1
aid_UC <- sapply(names(UC_in_only),function(i) {
 names(which(rowSums(UC_in_only[[i]] > cut) > 0))
  })  
mean(unlist(lapply(aid_UC,length)))#1074546
#calculate dMML between any possible pairs
 #Non tissue-specific high dMML
 #Find a cut that generate same number of regions as UC

 cut_dMML=0.1
 low_dMML=0.1
 high_dMML=0.5
 diff_log=data.table()
 diff=1
 while(abs(diff)>0.01){
     if(diff>0){low_dMML=cut_dMML}else(high_dMML=cut_dMML)
    cut_dMML=(high_dMML+low_dMML)/2
    aid_dMML <- sapply(names(dMML_in),function(i) {

        names(which(rowSums(dMML_in[[i]] > cut_dMML) > 0))
    })  
    diff=(mean(unlist(lapply(aid_dMML[names(aid_dMML)!='liver'],length)))-1074546)/1074546
    #diff >0: increase cut, replace low dMML to current cut
    
    log_tb=data.table(cut_dMML=cut_dMML,diff=diff)
    print(log_tb)
    diff_log=rbind(diff_log,log_tb)
}
#Without liver cut_dMML=0.134375
aid_tb=table(table(unlist(aid_dMML[names(aid_dMML)!='liver'])))
#      1       2       3       4       5       6       7
#2051100  864257  495128  313922  189435  109431   44978
aid_tb_dt=data.table(n_region=aid_tb)
colnames(aid_tb_dt)=c('n_tissue','n_region')
aid_tb_dt$n_region_pro=aid_tb_dt$n_region/sum(aid_tb_dt$n_region)
aid_tb_dt$n_tissue=factor(aid_tb_dt$n_tissue,levels=as.character(1:10))
pdf(paste0(figure_path,'number_tissue_pass_filter_dMML_1p_with_liver_cutoff_prop.pdf'))
ggplot(aid_tb_dt,aes(x=n_tissue,y=n_region_pro))+geom_bar(stat='identity', fill="steelblue")+
xlab('Number of tissue (except liver)')+ylab('Proportion of regions')+
geom_text(aes(label=round(n_region_pro,digits=2)),vjust=-1.6)
dev.off()
plot_dt=data.table(value=c(dMML_vec,UC_in_vec),
                    dat_type=c(rep('dMML',length(dMML_vec)),rep('UC',length(UC_in_vec))))
pdf(paste0(figure_path,'dMML_UC_dist.pdf'))

ggplot(plot_dt, aes(x=value,group=dat_type,color=dat_type))+geom_density()+xlim(c(0,0.15))
dev.off()

#generate Venn diagrame showing overlap between dMML region and UC in each tissue
#Filter no include liver
venn_out=data.table()
for(ts in names(aid_UC)){
    UC_specific=setdiff(aid_UC[[ts]],aid_dMML[[ts]])
    dMML_specific=setdiff(aid_dMML[[ts]],aid_UC[[ts]])
    shared=intersect(aid_UC[[ts]],aid_dMML[[ts]])
    venn_out=rbind(venn_out,data.table(tissue=ts,
                                        UC_specific=length(UC_specific),
                                        dMML_specific=length(dMML_specific),
                                        shared=length(shared)))
    print(venn_out)
}
venn_out_prop=cbind(venn_out[,1],venn_out[,-1]/rowSums(venn_out[,-1]))

#Count dNME regions
source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
dNME_in=mclapply(UC_in,function(x) x[,grepl('dNME',colnames(x))],mc.cores=20)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=20)
 cut_dNME=0.1
 low_dNME=0.1
 high_dNME=0.5
 diff_log=data.table()
 diff=1
 while(abs(diff)>0.01){
     if(diff>0){low_dNME=cut_dNME}else(high_dNME=cut_dNME)
    cut_dNME=(high_dNME+low_dNME)/2
    aid_dNME <- sapply(names(dNME_in),function(i) {

        names(which(rowSums(dNME_in[[i]] > cut_dNME) > 0))
    })  
    diff=(mean(unlist(lapply(aid_dNME,length)))-1074546)/1074546
    #diff >0: increase cut, replace low dNME to current cut
    
    log_tb=data.table(cut_dNME=cut_dNME,diff=diff)
    print(log_tb)
    diff_log=rbind(diff_log,log_tb)
}
#Without liver cut_dNME=0.31875
aid_tb=table(table(unlist(aid_dNME)))
#  1       2       3       4       5       6       7
# 1164341  750930  527652  355675  209875   96894   26770
aid_tb_dt=data.table(n_region=aid_tb)
colnames(aid_tb_dt)=c('n_tissue','n_region')
aid_tb_dt$n_region_pro=aid_tb_dt$n_region/sum(aid_tb_dt$n_region)
aid_tb_dt$n_tissue=factor(aid_tb_dt$n_tissue,levels=as.character(1:10))
pdf(paste0(figure_path,'number_tissue_pass_filter_dNME_1p_cutoff_prop.pdf'))
ggplot(aid_tb_dt,aes(x=n_tissue,y=n_region_pro))+geom_bar(stat='identity', fill="steelblue")+
xlab('Number of tissue')+ylab('Proportion of regions')+
geom_text(aes(label=round(n_region_pro,digits=2)),vjust=-1.6)
dev.off()
#generate Venn diagrame showing overlap between dNME region and UC in each tissue
#Filter no include liver
venn_out=data.table()
for(ts in names(aid_UC)){
    UC_specific=setdiff(aid_UC[[ts]],aid_dNME[[ts]])
    dNME_specific=setdiff(aid_dNME[[ts]],aid_UC[[ts]])
    shared=intersect(aid_UC[[ts]],aid_dNME[[ts]])
    venn_out=rbind(venn_out,data.table(tissue=ts,
                                        UC_specific=length(UC_specific),
                                        dNME_specific=length(dNME_specific),
                                        shared=length(shared)))
    print(venn_out)
}
venn_out_prop=cbind(venn_out[,1],venn_out[,-1]/rowSums(venn_out[,-1]))
cut=0.1
aid_UC <- sapply(names(UC_in_only),function(i) {
 names(which(rowSums(UC_in_only[[i]] > cut) > 0))
  })  
venn_out_all=list()

#Try different cut of dNME and see overlap
for(cut_dNME in seq(0.001,1,0.001)){
    cat("Processing:",cut_dNME,'\n')
      aid_dNME <- sapply(names(dNME_in),function(i) {

        names(which(rowSums(dNME_in[[i]] > cut_dNME) > 0))
    })  
    venn_out=data.table()
    for(ts in names(aid_UC)){
    UC_specific=setdiff(aid_UC[[ts]],aid_dNME[[ts]])
    dNME_specific=setdiff(aid_dNME[[ts]],aid_UC[[ts]])
    shared=intersect(aid_UC[[ts]],aid_dNME[[ts]])
    venn_out=rbind(venn_out,data.table(tissue=ts,
                                        UC_specific=length(UC_specific),
                                        dNME_specific=length(dNME_specific),
                                        shared=length(shared)))
    
}
venn_out_all[[as.character(cut_dNME)]]=venn_out

}
saveRDS(venn_out_all,'../downstream/output/mouse_analysis/dNME_UC_venn_out_all.rds')
#Find region overlap example
install_version("rJava", version = "0.9.12", repos = "http://cran.us.r-project.org")