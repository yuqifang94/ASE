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