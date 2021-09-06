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
dNME_in=mclapply(UC_in,function(x) x[,grepl('dNME',colnames(x))],mc.cores=7)
dMML_in=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=7)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=7)
rm(UC_in)
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
#cut_dNME=0.31875
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

#dMML

venn_out_all_dMML=list()

#Try different cut of dNME and see overlap
for(cut_dMML in seq(0.001,1,0.001)){
    cat("Processing:",cut_dMML,'\n')
      aid_dMML <- sapply(names(dMML_in),function(i) {

        names(which(rowSums(dMML_in[[i]] > cut_dMML) > 0))
    })  
    venn_out=data.table()
    for(ts in names(aid_UC)){
    UC_specific=setdiff(aid_UC[[ts]],aid_dMML[[ts]])
    dMML_specific=setdiff(aid_dMML[[ts]],aid_UC[[ts]])
    shared=intersect(aid_UC[[ts]],aid_dMML[[ts]])
    venn_out=rbind(venn_out,data.table(tissue=ts,
                                        UC_specific=length(UC_specific),
                                        dMML_specific=length(dMML_specific),
                                        shared=length(shared)))
    
}
venn_out_all_dMML[[as.character(cut_dMML)]]=venn_out

}
saveRDS(venn_out_all_dMML,'../downstream/output/mouse_analysis/UC_dNME_olap/dMML_UC_venn_out_all.rds')

#Fixing dNME and dMML cutoff, vary UC cutoff
#Try different cut of dNME and see overlap
venn_out_all_UC_all=data.table()
for(cut_dMML in seq(0,1,0.01)){
    for(cut_dNME in seq(0,1,0.01)){
        aid_dMML <- sapply(names(dMML_in),function(i) {

                names(which(rowSums(dMML_in[[i]] > cut_dMML) > 0))
        })  
        aid_dNME <- sapply(names(dNME_in),function(i) {

                names(which(rowSums(dNME_in[[i]] > cut_dNME) > 0))
        })  
        set.seed(12345)
        venn_out_all_UC=list()
        tt1=proc.time()
        venn_out_all_UC=mclapply(seq(0,1,0.01),function(cut_UC){
            cat("Processing:",cut_UC,'\n')
            aid_UC <- sapply(names(UC_in_only),function(i) {

                names(which(rowSums(UC_in_only[[i]] > cut_UC) > 0))
            })  
            venn_out=data.table()
            

            for(ts in names(aid_UC)){

            UC_specific_dMML=setdiff(aid_UC[[ts]],aid_dMML[[ts]])
            shared_dMML=intersect(aid_UC[[ts]],aid_dMML[[ts]])
            dMML_specific=setdiff(aid_dMML[[ts]],aid_UC[[ts]])

            shared_dNME=intersect(aid_UC[[ts]],aid_dNME[[ts]])
            dNME_specific=setdiff(aid_dNME[[ts]],aid_UC[[ts]])
            UC_specific_dNME=setdiff(aid_UC[[ts]],aid_dNME[[ts]])

            #Get random control for UC
            aid_UC_random=sample(rownames(UC_in_only[[ts]]),length(aid_UC[[ts]]))
            UC_specific_dMML_rand=setdiff(aid_UC_random,aid_dMML[[ts]])
            shared_dMML_rand=intersect(aid_UC_random,aid_dMML[[ts]])
            dMML_specific_rand=setdiff(aid_dMML[[ts]],aid_UC_random)

            shared_dNME_rand=intersect(aid_UC_random,aid_dNME[[ts]])
            dNME_specific_rand=setdiff(aid_dNME[[ts]],aid_UC_random)
            UC_specific_dNME_rand=setdiff(aid_UC_random,aid_dNME[[ts]])


            venn_out=rbind(venn_out,data.table(tissue=ts,
                                                UC_specific_dMML=length(UC_specific_dMML),
                                                shared_dMML=length(shared_dMML),
                                                dMML_specific=length(dMML_specific),
                                                shared_dNME=length(shared_dNME),
                                                dNME_specific=length(dNME_specific),
                                                UC_specific_dNME=length(UC_specific_dNME),
                                                cutoff=cut_UC,
                                                UC_specific_dMML_rand=length(UC_specific_dMML_rand),
                                                shared_dMML_rand=length(shared_dMML_rand),
                                                dMML_specific_rand=length(dMML_specific_rand),
                                                shared_dNME_rand=length(shared_dNME_rand),
                                                dNME_specific_rand=length(dNME_specific_rand),
                                                UC_specific_dNME_rand=length(UC_specific_dNME_rand)
                                            ))

            
        }
        print(head(venn_out))
        return(venn_out)
        },mc.cores=20
        )
        venn_out_all_UC=do.call(rbind,venn_out_all_UC)
        venn_out_all_UC$dNME_cutoff=dNME_cutoff
        venn_out_all_UC$dMML_cutoff=dMML_cutoff
       venn_out_all_UC_all=rbind(venn_out_all_UC_all,))
         print(head(venn_out_all_UC_all))
        cat("Finishing processing: dMML=",dMML_cutoff," dNME=",dNME_cutoff,"in:",proc.time()[[3]]-tt1,'\n')
    }
}

saveRDS(venn_out_all_UC_all,'../downstream/output/mouse_analysis/UC_dNME_olap/UC_var_venn_out_all_dNME_dMML_cutoff.rds')

#Combine dMML and dNME
venn_out_all_dMML=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/dMML_UC_venn_out_all.rds')
venn_out_all_dNME=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_UC_venn_out_all.rds')
venn_out_all_dMML=lapply(venn_out_all_dMML,function(x) {
                                            x$stat_type='dMML'
                                            colnames(x)[colnames(x)=='dMML_specific']='stat_specific'
                                            return(x)
                                            })
venn_out_all_dNME=lapply(venn_out_all_dNME,function(x) {
                                            x$stat_type='dNME'
                                            colnames(x)[colnames(x)=='dNME_specific']='stat_specific'
                                            return(x)
                                            })
venn_out_all=list()
for(cutoff in names(venn_out_all_dNME)){
    venn_out_all[[cutoff]]=rbind(venn_out_all_dMML[[cutoff]],venn_out_all_dNME[[cutoff]])
    venn_out_all[[cutoff]]$cutoff=as.numeric(cutoff)


}
venn_out_all=do.call(rbind,venn_out_all)

venn_out_all_prop1=data.table(stat_specific_prop=venn_out_all$stat_specific/(venn_out_all$stat_specific+venn_out_all$shared),
                                shared_prop=venn_out_all$shared/(venn_out_all$stat_specific+venn_out_all$shared))
venn_out_all_prop2=venn_out_all[,list(UC_specific,stat_specific,shared)]/rowSums(venn_out_all[,list(UC_specific,stat_specific,shared)])
colnames(venn_out_all_prop2)=paste0(colnames(venn_out_all_prop2),'_prop_all')
venn_out_all=cbind(venn_out_all,venn_out_all_prop2)
venn_out_all=cbind(venn_out_all,venn_out_all_prop1)
# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_olap_prop.pdf')
# ggplot(venn_out_all[!is.na(shared_prop)],aes(x=cutoff,y=shared_prop,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion regions in UC")+xlim(c(1,0))
# dev.off()
#Remove the filter that gives 0 region in dNME or dMML
venn_out_all=venn_out_all[stat_specific+shared!=0]
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_specific_prop2.pdf')
ggplot(venn_out_all[!is.na(stat_specific_prop)],aes(x=cutoff,y=stat_specific_prop,color=stat_type))+geom_smooth(se=TRUE)+
        xlab("cutoff")+ylab("Proportion stat specific regions")+xlim(c(0,1))
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_shared_prop2.pdf')
ggplot(venn_out_all[!is.na(shared_prop)],aes(x=cutoff,y=shared_prop,color=stat_type))+geom_smooth(se=TRUE)+
        xlab("cutoff")+ylab("Proportion shared regions")+xlim(c(0,1))
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_specific_prop_all.pdf')
ggplot(venn_out_all[!is.na(stat_specific_prop_all)],aes(x=cutoff,y=stat_specific_prop_all,color=stat_type))+geom_smooth(se=TRUE)+
        xlab("cutoff")+ylab("Proportion stat specific regions")+xlim(c(0,1))
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_shared_prop_all.pdf')
ggplot(venn_out_all[!is.na(shared_prop_all)],aes(x=cutoff,y=shared_prop_all,color=stat_type))+geom_smooth(se=TRUE)+
        xlab("cutoff")+ylab("Proportion shared regions")+xlim(c(0,1))
dev.off()
#Generating random control
set.seed(12345)
venn_out_all_rand=data.table()
pct_all=seq(0,100,5)
tt1=proc.time()[[3]]
for(i in 1:nrow(venn_out_all)){
    pct=round(i/nrow(venn_out_all)*100)
    if(pct %in% pct_all){
        #remove that percent
        pct_all=pct_all[pct_all!=pct]
        cat("Finishing processing:",pct,'% in',proc.time()[[3]]-tt1,'\n')
    }
    venn_in= venn_out_all[i]
    if(venn_in$stat_type=='dMML'){
        stat_mat=dMML_in[[venn_in$tissue]]
    }else 
    if(venn_in$stat_type=='dNME'){
        stat_mat=dNME_in[[venn_in$tissue]]
    }
    rand_region=sample(rownames(stat_mat),venn_in$stat_specific+venn_in$shared)
    venn_out=data.table(
                            UC_specific_rand=length(setdiff(aid_UC[[venn_in$tissue]],rand_region)),
                            shared_rand=length(intersect(aid_UC[[venn_in$tissue]],rand_region)),
                            stat_specific_rand=length(setdiff(rand_region,aid_UC[[venn_in$tissue]]))
    )
    venn_out_prop_all=venn_out/rowSums(venn_out)
    colnames(venn_out_prop_all)=paste0(colnames(venn_out_prop_all),'_prop_all')
    venn_out_prop=venn_out/(venn_out$shared_rand+venn_out$stat_specific_rand)
    colnames(venn_out_prop)=paste0(colnames(venn_out_prop),'_prop')
    venn_out_all_rand=rbind(venn_out_all_rand,
                        cbind(venn_in,venn_out,venn_out_prop_all,venn_out_prop))
}
saveRDS(venn_out_all_rand,'../downstream/output/mouse_analysis/UC_dNME_olap/venn_out_all_rand.rds')
venn_out_all_rand=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/venn_out_all_rand.rds')
venn_out_all_rand$total_stat=venn_out_all_rand$stat_specific+venn_out_all_rand$shared
#Plot actual numbers
venn_out_all_rand_dNME=venn_out_all_rand[stat_type=='dNME']

venn_out_all_rand_dNME_mt=melt.data.table(venn_out_all_rand_dNME[,list(stat_specific,
                                                                        shared,total_stat,
                                                                        stat_specific_rand,
                                                                        shared_rand,
                                                                        cutoff,
                                                                        tissue,
                                                                        stat_specific_prop,
                                                                        shared_prop,
                                                                        stat_specific_rand_prop,
                                                                        shared_rand_prop)],
                                                        id.vars=c('cutoff','tissue'))
venn_out_all_rand_dNME_mt$value_type=gsub('_rand','',venn_out_all_rand_dNME_mt$variable)
venn_out_all_rand_dNME_mt$random_control=grepl("_rand",venn_out_all_rand_dNME_mt$variable)
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_UC_count_tissue_shared_dNME.pdf')
for(ts in unique(venn_out_all_rand_dNME_mt$tissue)){
    print(ggplot(venn_out_all_rand_dNME_mt[tissue==ts&variable%in%c("shared","shared_rand")],aes(x=cutoff,y=value,color=random_control))+
    #geom_smooth(aes(linetype=value_type),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Number of regions")+xlim(c(0,1))+ggtitle(ts))
}
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_UC_count_tissue_specific_dNME.pdf')
for(ts in unique(venn_out_all_rand_dNME_mt$tissue)){
    print(ggplot(venn_out_all_rand_dNME_mt[tissue==ts&variable%in%c("stat_specific","stat_specific_rand")],aes(x=cutoff,y=value,color=random_control))+
    #geom_smooth(aes(linetype=value_type),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Number of regions")+xlim(c(0,1))+ggtitle(ts))
}
dev.off()

pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_UC_prop_tissue_shared_dNME.pdf')
for(ts in unique(venn_out_all_rand_dNME_mt$tissue)){
    print(ggplot(venn_out_all_rand_dNME_mt[tissue==ts&variable%in%c("shared_prop","shared_rand_prop")],aes(x=cutoff,y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Porportion overlapped of regions")+xlim(c(0,1))+ggtitle(ts))+ylim(c(0,1))
}
dev.off()

pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_UC_prop_tissue_shared_dNME_quant.pdf')
for(ts in unique(venn_out_all_rand_dNME_mt$tissue)){
    dNME_ecdf=ecdf(UC_in_max_loc[[ts]]$dNME_max_pair)
    print(ggplot(venn_out_all_rand_dNME_mt[tissue==ts&variable%in%c("shared_prop","shared_rand_prop")],aes(x=dNME_ecdf(cutoff),y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Porportion overlapped of regions")+xlim(c(0,1))+ggtitle(ts))+ylim(c(0,1))
}
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_UC_prop_tissue_specific_dNME.pdf')
for(ts in unique(venn_out_all_rand_dNME_mt$tissue)){
    print(ggplot(venn_out_all_rand_dNME_mt[tissue==ts&variable%in%c("stat_specific_prop","stat_specific_rand_prop")],aes(x=cutoff,y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
            xlab("cutoff")+ylab("Porportion of regions")+xlim(c(0,1))+ggtitle(ts))+ylim(c(0,1))
}
dev.off()
#dMML
venn_out_all_rand_dMML_mt=melt.data.table(venn_out_all_rand[stat_type=='dMML',list(stat_specific,
                                                                        shared,
                                                                        #total_stat,
                                                                        stat_specific_rand,
                                                                        shared_rand,
                                                                        cutoff,
                                                                        tissue,
                                                                        stat_specific_prop,
                                                                        shared_prop,
                                                                        stat_specific_rand_prop,
                                                                        shared_rand_prop)],
                                                        id.vars=c('cutoff','tissue'))
venn_out_all_rand_dMML_mt$value_type=gsub('_rand','',venn_out_all_rand_dMML_mt$variable)
venn_out_all_rand_dMML_mt$random_control=grepl("_rand",venn_out_all_rand_dMML_mt$variable)
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dMML_UC_count_tissue_shared_dMML.pdf')
for(ts in unique(venn_out_all_rand_dMML_mt$tissue)){
    print(ggplot(venn_out_all_rand_dMML_mt[tissue==ts&variable%in%c("shared","shared_rand")],aes(x=cutoff,y=value,color=random_control))+
    #geom_smooth(aes(linetype=value_type),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Number of regions")+xlim(c(0,1))+ggtitle(ts))
}
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dMML_UC_count_tissue_specific_dMML.pdf')
for(ts in unique(venn_out_all_rand_dMML_mt$tissue)){
    print(ggplot(venn_out_all_rand_dMML_mt[tissue==ts&variable%in%c("stat_specific","stat_specific_rand")],aes(x=cutoff,y=value,color=random_control))+
    #geom_smooth(aes(linetype=value_type),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Number of regions")+xlim(c(0,1))+ggtitle(ts))
}
dev.off()

pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dMML_UC_prop_tissue_shared_dMML.pdf')
for(ts in unique(venn_out_all_rand_dMML_mt$tissue)){
    print(ggplot(venn_out_all_rand_dMML_mt[tissue==ts&variable%in%c("shared_prop","shared_rand_prop")],aes(x=cutoff,y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Porportion of overlapped regions")+xlim(c(0,1))+ggtitle(ts))+ylim(c(0,1))
}
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dMML_UC_prop_tissue_specific_dMML.pdf')
for(ts in unique(venn_out_all_rand_dMML_mt$tissue)){
    print(ggplot(venn_out_all_rand_dMML_mt[tissue==ts&variable%in%c("stat_specific_prop","stat_specific_rand_prop")],aes(x=cutoff,y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Porportion of regions")+xlim(c(0,1))+ggtitle(ts))+ylim(c(0,1))
}
dev.off()


pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dMML_UC_prop_tissue_shared_dMML_quant.pdf')
for(ts in unique(venn_out_all_rand_dMML_mt$tissue)){
    dMML_ecdf=ecdf(UC_in_max_loc[[ts]]$dMML_max_pair)
    print(ggplot(venn_out_all_rand_dMML_mt[tissue==ts&variable%in%c("shared_prop","shared_rand_prop")],aes(x=dMML_ecdf(cutoff),y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=random_control),size=0.05)+
    xlab("cutoff")+ylab("Porportion of overlapped regions")+xlim(c(0,1))+ggtitle(ts))+ylim(c(0,1))
}
dev.off()

#Compare dMML and dNME
venn_out_all_rand_mt=melt.data.table(venn_out_all_rand[,list(stat_specific,
                                                                        shared,
                                                                        stat_type,
                                                                        total_stat,
                                                                        stat_specific_rand,
                                                                        shared_rand,
                                                                        cutoff,
                                                                        tissue,
                                                                        stat_specific_prop,
                                                                        shared_prop,
                                                                        stat_specific_rand_prop,
                                                                        shared_rand_prop)],
                                                        id.vars=c('cutoff','tissue','total_stat','stat_type'))
venn_out_all_rand_mt$value_type=gsub('_rand','',venn_out_all_rand_mt$variable)
venn_out_all_rand_mt$random_control=""
venn_out_all_rand_mt$random_control[grepl("_rand",venn_out_all_rand_mt$variable)]=" random control"
venn_out_all_rand_mt$stat_type_random=paste0(venn_out_all_rand_mt$stat_type,venn_out_all_rand_mt$random_control)
library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dMML_dNME_comp_same_region_new_cutoff_no_vline.pdf')
for(ts in unique(venn_out_all_rand_mt$tissue)){
 
    print(ggplot(venn_out_all_rand_mt[tissue==ts&variable%in%c("shared_prop","shared_rand_prop")],aes(x=total_stat,y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=stat_type_random),size=0.05)+
    xlab("Total selected regions")+ylab("Porportion of overlapped regions")+ggtitle(ts)+ylim(c(0,1))+
    #geom_vline(xintercept=sum(tissue_out_filtered[[ts]]$region_type %in% c("NME only","Both")),color='blue',alpha=0.1,size=0.5)+
    #geom_vline(xintercept=sum(tissue_out_filtered[[ts]]$region_type %in% c("MML only","Both")),color='red',alpha=0.1,size=0.5)+
    guides(color=guide_legend(title="",override.aes = list(size=10)))+ scale_x_continuous(trans =reverselog_trans(10)))
}
dev.off()
#Choose a cutoff that gives same number of regions with dNME only and both
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
venn_out_cat=data.table()
for(ts in names(tissue_out_filtered)){
    venn_out_cat=rbind(venn_out_cat,
    venn_out_all_rand_mt[stat_type=='dNME'&variable=="shared_prop"&tissue==ts][which.min(abs(total_stat-sum(tissue_out_filtered[[ts]]$region_type %in% c("NME only","Both"))))],
    venn_out_all_rand_mt[stat_type=='dMML'&variable=="shared_prop"&tissue==ts][which.min(abs(total_stat-sum(tissue_out_filtered[[ts]]$region_type %in% c("MML only","Both"))))]
    )

}
#Fix UC compare dNME and dMML

UC_var_venn_out_all=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/UC_var_venn_out_all.rds')
UC_var_venn_out_all$total_UC=UC_var_venn_out_all$UC_specific_dMML+UC_var_venn_out_all$shared_dMML
UC_var_venn_out_all$total_dMML=UC_var_venn_out_all$dMML_specific+UC_var_venn_out_all$shared_dMML
UC_var_venn_out_all$total_dNME=UC_var_venn_out_all$dNME_specific+UC_var_venn_out_all$shared_dNME
UC_var_venn_out_all$total_dMML_rand=UC_var_venn_out_all$dMML_specific_rand+UC_var_venn_out_all$shared_dMML_rand
UC_var_venn_out_all$total_dNME_rand=UC_var_venn_out_all$dNME_specific_rand+UC_var_venn_out_all$shared_dNME_rand
UC_var_venn_out_all$dMML_shared_prop=UC_var_venn_out_all$shared_dMML/UC_var_venn_out_all$total_dMML
UC_var_venn_out_all$dNME_shared_prop=UC_var_venn_out_all$shared_dNME/UC_var_venn_out_all$total_dNME
UC_var_venn_out_all$dMML_shared_prop_rand=UC_var_venn_out_all$shared_dMML_rand/UC_var_venn_out_all$total_dMML_rand
UC_var_venn_out_all$dNME_shared_prop_rand=UC_var_venn_out_all$shared_dNME_rand/UC_var_venn_out_all$total_dNME_rand

UC_var_venn_out_all_mt=melt.data.table(UC_var_venn_out_all,id.vars=c('tissue','cutoff','total_UC','total_dNME','total_dMML','total_dMML_rand','total_dNME_rand'))
UC_var_venn_out_all_mt$stat_var=gsub('_.*','',UC_var_venn_out_all_mt$variable)
UC_var_venn_out_all_mt_prop_shared=UC_var_venn_out_all_mt[grepl('shared_prop',UC_var_venn_out_all_mt$variable)]
UC_var_venn_out_all_mt_prop_shared$random=""
UC_var_venn_out_all_mt_prop_shared[grepl('rand',UC_var_venn_out_all_mt_prop_shared$variable)]$random="random control"
UC_var_venn_out_all_mt_prop_shared$data_type=paste(UC_var_venn_out_all_mt_prop_shared$stat_var,UC_var_venn_out_all_mt_prop_shared$random)
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/fix_UC_var_dNME.pdf')
for(ts in unique(UC_var_venn_out_all_mt_prop_shared$tissue)){
 
    print(ggplot(UC_var_venn_out_all_mt_prop_shared[tissue==ts],aes(x=cutoff,y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=data_type),size=0.05)+
    xlab("UC cutoff")+ylab("Porportion of overlapped regions")+ggtitle(ts)+ylim(c(0,1))+
    guides(color=guide_legend(title="",override.aes = list(size=10)))+
    scale_x_reverse()
    #scale_x_continuous(trans =reverselog_trans(10))
    )
}
dev.off()
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/fix_UC_var_dNME_region.pdf')
for(ts in unique(UC_var_venn_out_all_mt_prop_shared$tissue)){
 
    print(ggplot(UC_var_venn_out_all_mt_prop_shared[tissue==ts],aes(x=total_UC,y=value))+
    #geom_smooth(aes(linetype=random_control),se=TRUE)+
    geom_point(aes(color=data_type),size=0.05)+
    xlab("Total selected UC regions")+ylab("Porportion of overlapped regions")+ggtitle(ts)+ylim(c(0,1))+
    guides(color=guide_legend(title="",override.aes = list(size=10)))
    #scale_x_continuous(trans ='reverse')
    )
}
dev.off()

# #Plot figure for randon control

# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_specific_prop_rand.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=stat_specific_rand_prop,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion stat specific regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()

# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_specific_prop.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=stat_specific_prop,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion stat specific regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()

# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_shared_prop_rand.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=shared_rand_prop,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion shared regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()

# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_shared_prop.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=shared_prop,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion shared regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()

# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_specific_prop_all_rand.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=stat_specific_rand_prop_all,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion stat specific regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()
# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_specific_prop_all.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=stat_specific_prop_all,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion stat specific regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()


# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_shared_prop_all_rand.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=shared_rand_prop_all,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion shared regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()
# pdf('../downstream/output/mouse_analysis/UC_dNME_olap/dNME_dMML_UC_shared_prop_all.pdf')
# ggplot(venn_out_all_rand,aes(x=cutoff,y=shared_prop_all,color=stat_type))+geom_smooth(se=TRUE)+
#         xlab("cutoff")+ylab("Proportion shared regions")+xlim(c(0,1))+ylim(c(0,1))
# dev.off()

#Find dNME value at overlapping and non-overlapping regions
cut_UC=0.1
aid_UC <- sapply(names(UC_in_only),function(i) {
 names(which(rowSums(UC_in_only[[i]] > cut_UC) > 0))
  })  
cut_dNME=0.31875
aid_dNME <- sapply(names(dNME_in),function(i) {

    names(which(rowSums(dNME_in[[i]] > cut_dNME) > 0))
})  
max_dat_out=data.table()
for(ts in names(aid_UC)){
    UC_specific=setdiff(aid_UC[[ts]],aid_dNME[[ts]])
    dNME_specific=setdiff(aid_dNME[[ts]],aid_UC[[ts]])
    shared=intersect(aid_UC[[ts]],aid_dNME[[ts]])
    max_dat_out=rbind(max_dat_out,
        data.table(
            tissue=ts,
            region=UC_specific,
            region_type='UC_specific',
            maxdNME=rowMaxs(dNME_in[[ts]][UC_specific,],na.rm=T),
            maxUC=rowMaxs(UC_in_only[[ts]][UC_specific,],na.rm=T)
        ),
        data.table(
            tissue=ts,
            region=shared,
            region_type='shared',
            maxdNME=rowMaxs(dNME_in[[ts]][shared,],na.rm=T),
            maxUC=rowMaxs(UC_in_only[[ts]][shared,],na.rm=T)
        ),
        data.table(
            tissue=ts,
            region=dNME_specific,
            region_type='dNME_specific',
            maxdNME=rowMaxs(dNME_in[[ts]][dNME_specific,],na.rm=T),
            maxUC=rowMaxs(UC_in_only[[ts]][dNME_specific,],na.rm=T)
        )
    )
}
pdf('../downstream/output/mouse_analysis/UC_dNME_olap/maxdNME_dist.pdf')
ggplot(max_dat_out,aes(x=maxdNME,color=region_type))+geom_density()+xlab("dNME")
dev.off()
#Find region overlap example
set.seed(123)
max_dat_rand=max_dat_out[region_type=='dNME_specific'&tissue=='heart'][sample(1:nrow(max_dat_out[region_type=='dNME_specific'&tissue=='heart']),15)]
dNME_rand=dNME_in[['heart']][max_dat_rand$region,]

colnames(dNME_rand)[apply(dNME_rand,1,which.max)]
rownames(dNME_rand)
UC_in[['heart']]["chr8:47268997-47269246",grepl('E13.5-E14.5',colnames(UC_in[['heart']]))]
UC_in[['heart']]["chr12:79071734-79071983",grepl('E12.5-E14.5',colnames(UC_in[['heart']]))]
UC_in[['heart']]["chr4:131255401-131255650",grepl('E14.5-E15.5',colnames(UC_in[['heart']]))]
UC_in[['heart']]["chr12:83702876-83703125",grepl('E12.5-E13.5',colnames(UC_in[['heart']]))]
UC_in[['heart']]["chr7:98180174-98180426",grepl('E15.5-E16.5',colnames(UC_in[['heart']]))]

#QQ plot

source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
UC_in_max_loc=readRDS(UC_merge_max_loc_file)
dNME_in=mclapply(UC_in,function(x) x[,grepl('dNME',colnames(x))],mc.cores=7)
dMML_in=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=7)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=7)
rm(UC_in)
for (ts in names(UC_in_max_loc)){
    max_ts=UC_in_max_loc[[ts]][,c('UC_max_pair','dNME_max_UC_pair','dMML_max_UC_pair')]
    max_ts=cbind(data.table(regions=rownames(max_ts)),as.data.table(max_ts))
    #type='cairo' is needed for JHPCE
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dNME_max_UC_pair,x=UC_max_pair))+geom_point(alpha=0.01,color='blue')+geom_vline(xintercept=0.1,color='red')+
    xlab('UC')+ylab('dNME'))
    dev.off()
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dMML_max_UC_pair,x=UC_max_pair))+geom_point(alpha=0.01,color='blue')+geom_vline(xintercept=0.1,color='red')+
    xlab('UC')+ylab('dMML'))
    dev.off()
    max_ts$dNME_quant=ecdf(max_ts$dNME_max_UC_pair)(max_ts$dNME_max_UC_pair)
    max_ts$dMML_quant=ecdf(max_ts$dMML_max_UC_pair)(max_ts$dMML_max_UC_pair)
    max_ts$UC_quant=ecdf(max_ts$UC_max_pair)(max_ts$UC_max_pair)
    #ecdf(max_ts$UC_max_pair)(0.1)=0.8060806
    UC_cutoff_quant=ecdf(max_ts$UC_max_pair)(0.1)
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_quant.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dNME_quant*100,x=UC_quant*100))+geom_point(alpha=0.001,color='blue')+
    geom_vline(xintercept=UC_cutoff_quant,color='red')+xlab('UC quantile')+ylab('dNME quantile'))
    dev.off()
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_quant.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dMML_quant*100,x=UC_quant*100))+geom_point(alpha=0.001,color='blue')+
    geom_vline(xintercept=UC_cutoff_quant,color='red')+xlab('UC quantile')+ylab('dMML quantile'))
    dev.off()
    max_ts$high_UC="smaller than 0.1"
    max_ts[UC_max_pair>0.1]$high_UC="greater than 0.1"
      png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_density_04.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts[dMML_max_UC_pair>0.4])+geom_density(aes(x=dMML_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0.4,1)
    dev.off()
        png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_density_all.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts)+geom_density(aes(x=dMML_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0,1)
    dev.off()
         png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_density_04.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts[dNME_max_UC_pair>0.4])+geom_density(aes(x=dNME_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0.4,1)
    dev.off()
        png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_density_all.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts)+geom_density(aes(x=dNME_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0,1)
    dev.off()
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_QQ.png'),type='cairo',width=1000,height=1000)
    qqplot(x=max_ts$UC_max_pair, y=max_ts$dNME_max_UC_pair,xlab='UC',ylab='dNME')
    abline(0,1,col='black')
    abline(v=0.1,col='red')
    dev.off()
      png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_QQ.png'),type='cairo',width=1000,height=1000)
    qqplot(x=max_ts$UC_max_pair, y=max_ts$dMML_max_UC_pair,xlab='UC',ylab='dMML')
    abline(0,1,col='black')
    abline(v=0.1,col='red')
    dev.off()
}
#Do the Q-Q plot for the UC>0.1 and catogrization
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
for(ts in names(UC_in_max_loc)){
    max_ts=UC_in_max_loc[[ts]][,c('UC_max_pair','dNME_max_UC_pair','dMML_max_UC_pair')]
    max_ts=cbind(data.table(regions=rownames(max_ts)),as.data.table(max_ts))
    cat_ts=tissue_out_filtered[[ts]]
    max_ts=max_ts[match(cat_ts$region,regions)]
    max_ts$region_type=cat_ts$region_type
    col_region_type=list(`NME only` ='red',`MML only`='blue',`Neither`='green',`Both`='orange')
    for (stat_type in c('dMML','dNME')){
        png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_',stat_type,'_QQ_all_cat.png'),type='cairo',width=1000,height=1000)
        qqplot(x=max_ts$UC_max_pair, y=max_ts[[paste0(stat_type,'_max_UC_pair')]],xlab='UC',ylab=stat_type,xlim=c(0,1),ylim=c(0,1))
        abline(0,1,col='black')
        dev.off()
        png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_',stat_type,'_QQ_cat.png'),type='cairo',width=1000,height=1000)
        par(new=F)
        for(rt in unique(max_ts$region_type)){
            qqplot(x=max_ts[region_type==rt]$UC_max_pair, y=max_ts[region_type==rt][[paste0(stat_type,'_max_UC_pair')]],
                    xlab='UC',ylab=stat_type,xlim=c(0,1),ylim=c(0,1),col=col_region_type[[rt]])
            par(new=T)
        
            }
        abline(0,1,col='black')
        legend("topleft", legend=c("NME only", "MML only","Both","Neither"),
        col=c("red", "blue","green","orange"),pch=1)
        dev.off()

    }
}
#Heart ecdf
dNME_in_heart_ecdf=ecdf(dNME_in$heart)
UC_in_heart_ecdf=ecdf(UC_in_only$heart)
dMML_in_heart_ecdf=ecdf(dMML_in$heart)

random_sample_low_UC=UC_in_max_loc[['heart']][sample(which(UC_in_max_loc[['heart']]$UC_max_pair<=0.1),15),c("UC_max_pair","dNME_max_UC_pair")]
colnames(random_sample_low_UC)=c("UC","dNME")
random_sample_low_UC$UC_quantile=UC_in_heart_ecdf(random_sample_low_UC$UC)
random_sample_low_UC$dNME_quantile=UC_in_heart_ecdf(random_sample_low_UC$dNME)
random_sample_low_UC$dNME_quantile_log10=-log10(1-UC_in_heart_ecdf(random_sample_low_UC$dNME))
rownames(random_sample_low_UC)=NULL
random_sample_high_UC=UC_in_max_loc[['heart']][sample(which(UC_in_max_loc[['heart']]$UC_max_pair>0.1),15),c("UC_max_pair","dNME_max_UC_pair")]
colnames(random_sample_high_UC)=c("UC","dNME")
random_sample_high_UC$UC_quantile=UC_in_heart_ecdf(random_sample_high_UC$UC)
random_sample_high_UC$dNME_quantile=UC_in_heart_ecdf(random_sample_high_UC$dNME)
random_sample_high_UC$dNME_quantile_log10=-log10(1-UC_in_heart_ecdf(random_sample_high_UC$dNME))
rownames(random_sample_high_UC)=NULL
#Get boxplot, etc for tissue
source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
UC_in_max_loc=readRDS(UC_merge_max_loc_file)
GO_out_ts_low_UC_high_dNME_out=list()
GO_out_ts_low_UC_high_dMML_out=list()
cut_dNME=0.31875
cut_dMML=0.134375
enhancer=readRDS(bin_enhancer_rds)#21441
enhancer_bg=subsetByOverlaps(enhancer,convert_GR(unique(unlist(lapply(UC_in_max_loc,rownames)))))
for(ts in names(UC_in_max_loc)){
    cat("Processing:",ts,"\n")
    tissue_UC=UC_in_max_loc[[ts]]
    tissue_UC$filter_UC="UC<=0.1"
    tissue_UC$filter_UC[tissue_UC$UC_max_pair>0.1]="UC>0.1"
    tissue_UC$dNME_quant=ecdf(tissue_UC$dNME_max_UC_pair)(tissue_UC$dNME_max_UC_pair)

    pdf(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_dNME_at_highUC.pdf'))
    print(ggplot(tissue_UC,aes(y=dNME_max_UC_pair,x=filter_UC))+geom_boxplot(outlier.shape=NA)+
        xlab('')+ylab('dNME')+ylim(c(0,0.75)))
    dev.off()
    pdf(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_dNME_quant_at_highUC.pdf'))
    print(ggplot(tissue_UC,aes(y=dNME_quant,x=filter_UC))+geom_boxplot(outlier.shape=NA)+
        xlab('')+ylab('dNME quantile'))
    dev.off()

    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_dNME_QQ_all.png'),type='cairo',width=1000,height=1000)
    qqplot(x=tissue_UC$dNME_max_UC_pair, y=tissue_UC$dNME_max_UC_pair[tissue_UC$filter_UC=="UC>0.1"],xlab='Expected dNME',ylab='Observed dNME at UC>0.1',xlim=c(0,1),ylim=c(0,1))
    abline(0,1,col='black')
    dev.off()
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_dMML_QQ_all.png'),type='cairo',width=1000,height=1000)
    qqplot(x=tissue_UC$dMML_max_UC_pair, y=tissue_UC$dMML_max_UC_pair[tissue_UC$filter_UC=="UC>0.1"],xlab='Expected dMML',ylab='Observed dMML at UC>0.1',xlim=c(0,1),ylim=c(0,1))
    abline(0,1,col='black')
    dev.off()
    #GO analysis with UC<=0.1 but not dNME

    tissue_UC$dNME_filter=FALSE
    tissue_UC$dNME_filter[tissue_UC$dNME_max_pair>cut_dNME]=TRUE
    tissue_UC$dMML_filter=FALSE
    tissue_UC$dMML_filter[tissue_UC$dMML_max_pair>cut_dMML]=TRUE
    tissue_enhancer_lowUC_high_dNME=subsetByOverlaps(enhancer,convert_GR(rownames(tissue_UC)[tissue_UC$filter_UC=="UC<=0.1"&tissue_UC$dNME_filter]))#9468
    tissue_enhancer_lowUC_high_dMML=subsetByOverlaps(enhancer,convert_GR(rownames(tissue_UC)[tissue_UC$filter_UC=="UC<=0.1"&tissue_UC$dMML_filter]))#18400
    tissue_enhancer_high_UC=subsetByOverlaps(enhancer,convert_GR(rownames(tissue_UC)[tissue_UC$filter_UC=="UC>0.1"]))#15654
    GO_out_ts_low_UC_high_dNME_out[[ts]]=GO_run(unique(tissue_enhancer_lowUC_high_dNME$`Target Gene`),enhancer_bg$`Target Gene`,cluster=NA)
    GO_out_ts_low_UC_high_dMML_out[[ts]]=GO_run(unique(tissue_enhancer_lowUC_high_dMML$`Target Gene`),enhancer_bg$`Target Gene`,cluster=NA)
    GO_out_ts_high_UC[[ts]]=GO_run(unique(tissue_enhancer_high_UC$`Target Gene`),enhancer_bg$`Target Gene`,cluster=NA)
    print(GO_out_ts_low_UC_high_dNME_out[[ts]][FDR<=0.2&FC>=1.5])
    print(GO_out_ts_low_UC_high_dMML_out[[ts]][FDR<=0.2&FC>=1.5])
    print(GO_out_ts_high_UC[[ts]][FDR<=0.2&FC>=1.5])
}
saveRDS(GO_out_ts_high_UC,'../downstream/output/mouse_analysis/UC_dNME_olap/GO_out_ts_high_UC.rds')
saveRDS(GO_out_ts_low_UC_high_dNME_out,'../downstream/output/mouse_analysis/UC_dNME_olap/GO_out_ts_low_UC_high_dNME_out.rds')
saveRDS(GO_out_ts_low_UC_high_dMML_out,'../downstream/output/mouse_analysis/UC_dNME_olap/GO_out_ts_low_UC_high_dMML_out.rds')

do.call(rbind, lapply(GO_out_ts_low_UC_high_dNME_out,function(x) x[FDR<=0.2&FC>=1.5]))
do.call(rbind, lapply(GO_out_ts_low_UC_high_dMML_out,function(x) x[FDR<=0.2&FC>=1.5]))

#For each tissue, find a cutoff of dNME or dMML that gives similar number of region and check overlap, do it for non-tissue-specific now

source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
dNME_in=mclapply(UC_in,function(x) x[,grepl('dNME',colnames(x))],mc.cores=7)
dMML_in=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=7)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=7)
rm(UC_in)
cutoff_selection<-function(dat_in,n_regions,diff_threshold=0,cut_start=0,cut_end=1){
#Initialize key parameters
    diff=10000000
    cutoff=cut_start
    diff_log=data.table()
    while(abs(diff)>diff_threshold){
        if(diff>0){cut_start=cutoff}else(cut_end=cutoff)
        cutoff=(cut_start+cut_end)/2
        dat_region=names(which(rowSums(dat_in> cutoff) > 0))
        diff=(length(dat_region)-n_regions)
        #diff >0: increase cut, replace low dNME to current cut
        
        log_tb=data.table(cutoff=cutoff,diff=diff)
        print(log_tb)
        diff_log=rbind(diff_log,log_tb)
    }
    return(list(cutoff=cutoff,region_sel=dat_region,diff=diff,diff_log=diff_log))

}

tissue_out_filtered=readRDS(tissue_out_filtered_fn)

 cut=0.1
aid_UC <- sapply(names(UC_in_only),function(i) {
 names(which(rowSums(UC_in_only[[i]] > cut) > 0))
  })  

cutoff_out=list()
summary_output=data.table()
for(ts in names(tissue_out_filtered)){
    all_region=tissue_out_filtered[[ts]]$region
    dNME_regions=tissue_out_filtered[[ts]][region_type%in%c("NME only","Both")]$region
    dMML_regions=tissue_out_filtered[[ts]][region_type%in%c("MML only","Both")]$region
    dNME_prop=length(dNME_regions)/length(all_region)
    dMML_prop=length(dMML_regions)/length(all_region)
    all_UC=aid_UC[[ts]]
    dNME_cutoff=cutoff_selection(dNME_in[[ts]],round(length(all_UC)*dNME_prop),diff_threshold=10)
    dMML_cutoff=cutoff_selection(dMML_in[[ts]],round(length(all_UC)*dMML_prop),diff_threshold=10)
    summary_output=rbind(summary_output,
                        data.table(stat_type="dNME",cutoff=dNME_cutoff$cutoff,tissue=ts,
                                    UC_only=length(setdiff(all_UC,dNME_cutoff$region_sel)),
                                    shared=length(intersect(all_UC,dNME_cutoff$region_sel)),
                                    stat_specific=length(setdiff(dNME_cutoff$region_sel,all_UC))),
                         data.table(stat_type="dMML",cutoff=dMML_cutoff$cutoff,tissue=ts,
                                    UC_only=length(setdiff(all_UC,dMML_cutoff$region_sel)),
                                    shared=length(intersect(all_UC,dMML_cutoff$region_sel)),
                                    stat_specific=length(setdiff(dMML_cutoff$region_sel,all_UC)))
                    )            
    cutoff_out[[ts]]=list(dNME_cutoff=dNME_cutoff,dMML_cutoff=dMML_cutoff)

}
saveRDS(summary_output,'../downstream/output/mouse_analysis/UC_dNME_olap/summary_output_non_ts.rds')
saveRDS(cutoff_out,'../downstream/output/mouse_analysis/UC_dNME_olap/cutoff_out_non_ts.rds')
summary_output$overlap_prop=summary_output$shared/(summary_output$shared+summary_output$stat_specific)

#NME only & MML only 
cutoff_out=list()
summary_output=data.table()
for(ts in names(tissue_out_filtered)){
    all_region=tissue_out_filtered[[ts]]$region
    dNME_regions=tissue_out_filtered[[ts]][region_type%in%c("NME only")]$region
    dMML_regions=tissue_out_filtered[[ts]][region_type%in%c("MML only")]$region
    dNME_prop=length(dNME_regions)/length(all_region)
    dMML_prop=length(dMML_regions)/length(all_region)
    all_UC=aid_UC[[ts]]
    dNME_cutoff=cutoff_selection(dNME_in[[ts]],round(length(all_UC)*dNME_prop),diff_threshold=10)
    dMML_cutoff=cutoff_selection(dMML_in[[ts]],round(length(all_UC)*dMML_prop),diff_threshold=10)
    summary_output=rbind(summary_output,
                        data.table(stat_type="dNME",cutoff=dNME_cutoff$cutoff,tissue=ts,
                                    UC_only=length(setdiff(all_UC,dNME_cutoff$region_sel)),
                                    shared=length(intersect(all_UC,dNME_cutoff$region_sel)),
                                    stat_specific=length(setdiff(dNME_cutoff$region_sel,all_UC))),
                         data.table(stat_type="dMML",cutoff=dMML_cutoff$cutoff,tissue=ts,
                                    UC_only=length(setdiff(all_UC,dMML_cutoff$region_sel)),
                                    shared=length(intersect(all_UC,dMML_cutoff$region_sel)),
                                    stat_specific=length(setdiff(dMML_cutoff$region_sel,all_UC)))
                    )            
    cutoff_out[[ts]]=list(dNME_cutoff=dNME_cutoff,dMML_cutoff=dMML_cutoff)

}
summary_output$overlap_prop=summary_output$shared/(summary_output$shared+summary_output$stat_specific)
summary_output[,list(mean(overlap_prop)),by=list(stat_type)]
saveRDS(summary_output,'../downstream/output/mouse_analysis/UC_dNME_olap/summary_output_non_ts_only.rds')
saveRDS(cutoff_out,'../downstream/output/mouse_analysis/UC_dNME_olap/cutoff_out_non_ts_only.rds')


#One and only one overlap analysis
cutoff_selection_only<-function(dat_in,ts#,
    #n_regions,
    #diff_threshold=0.1,
    #step=0.1
    ){
    #Initialize key parameters
    #cutoff=0.6
    #diff=1
    #diff_log=data.table()
    #This is no longer linear try different threshold
    #while(abs(diff)>diff_threshold){
    
    return(mclapply(seq(0,1,0.001),function(cutoff,ts) {
        #cutoff=cutoff+step*diff
        cat("Processing:",cutoff,'\n')
        aid <- sapply(names(dat_in),function(i) {
                    names(which(rowSums(dat_in[[i]] > cutoff) > 0))
                    })          
        out=list()
        return(out[[as.character(cutoff)]]=setdiff(aid[[ts]],unlist(aid[names(aid)!=ts])))
        
        #diff=(length(dat_region)-n_regions)/n_regions
        #diff >0: increase cut, replace low dNME to current cut
        
        #log_tb=data.table(cutoff=cutoff,diff=diff)
        #print(log_tb)
        #diff_log=rbind(diff_log,log_tb)
    },ts=ts,mc.cores=20
    ))

    #return(list(cutoff=diff_log[which.min(abs(diff))]$cutoff,region_sel=dat_region[which.min(abs(diff_log$diff))],diff=diff_log[which.min(abs(diff))]$diff,diff_log=diff_log))
    #return(list(cutoff=diff_log,region_sel=dat_region,diff=diff,diff_log=diff_log))

}
#Relax UC cutoff for one and only one
cutoff_out=list()
summary_output=data.table()
dNME_cutoff=list()
dMML_cutoff=list()
UC_cutoff=list()
for(ts in names(dNME_in)){
    cat("Processing:",ts,"\n")
    #all_region=tissue_out_filtered[[ts]]$region
    #dNME_regions=tissue_out_filtered[[ts]][region_type%in%c("NME only","Both")]$region
    #dMML_regions=tissue_out_filtered[[ts]][region_type%in%c("MML only","Both")]$region
    dNME_cutoff[[ts]]=cutoff_selection_only(dNME_in,ts)
    
    dMML_cutoff[[ts]]=cutoff_selection_only(dMML_in,ts)
    UC_cutoff[[ts]]=cutoff_selection_only(UC_in_only,ts)
    # all_UC=aid_UC[[ts]]
    # summary_output=rbind(summary_output,
    #                     data.table(stat_type="dNME",olap_type='stat only',cutoff=dNME_cutoff$cutoff,tissue=ts,
    #                                 UC_only=length(setdiff(dNME_regions,dNME_cutoff$region_sel)),
    #                                 shared=length(intersect(dNME_regions,dNME_cutoff$region_sel)),
    #                                 stat_specific=length(setdiff(dNME_cutoff$region_sel,dNME_regions))),
    #                     data.table(stat_type="dNME",olap_type='all',cutoff=dNME_cutoff$cutoff,tissue=ts,
    #                                 UC_only=length(setdiff(all_region,dNME_cutoff$region_sel)),
    #                                 shared=length(intersect(all_region,dNME_cutoff$region_sel)),
    #                                 stat_specific=length(setdiff(dNME_cutoff$region_sel,all_region))),
    #                     data.table(stat_type="dMML",olap_type='stat only',cutoff=dMML_cutoff$cutoff,tissue=ts,
    #                                 UC_only=length(setdiff(dMML_regions,dMML_cutoff$region_sel)),
    #                                 shared=length(intersect(dMML_regions,dMML_cutoff$region_sel)),
    #                                 stat_specific=length(setdiff(dMML_cutoff$region_sel,dMML_regions))),
    #                     data.table(stat_type="dMML",olap_type='all',cutoff=dMML_cutoff$cutoff,tissue=ts,
    #                                 UC_only=length(setdiff(all_region,dMML_cutoff$region_sel)),
    #                                 shared=length(intersect(all_region,dMML_cutoff$region_sel)),
    #                                 stat_specific=length(setdiff(dMML_cutoff$region_sel,all_region))),
    #                     data.table(stat_type="dMML",olap_type='all UC',cutoff=dMML_cutoff$cutoff,tissue=ts,
    #                                 UC_only=length(setdiff(all_UC,dMML_cutoff$region_sel)),
    #                                 shared=length(intersect(all_UC,dMML_cutoff$region_sel)),
    #                                 stat_specific=length(setdiff(dMML_cutoff$region_sel,all_UC))),
    #                     data.table(stat_type="dNME",olap_type='all UC',cutoff=dNME_cutoff$cutoff,tissue=ts,
    #                                 UC_only=length(setdiff(all_UC,dNME_cutoff$region_sel)),
    #                                 shared=length(intersect(all_UC,dNME_cutoff$region_sel)),
    #                                 stat_specific=length(setdiff(dNME_cutoff$region_sel,all_UC)))
                                    
    #                     )   
    # cutoff_out[[ts]]=list(dNME_cutoff=dNME_cutoff,dMML_cutoff=dMML_cutoff)

}
#summary_output$overlap_prop=summary_output$shared/(summary_output$shared+summary_output$stat_specific)
assign_name<-function(x) {names(x)=as.character(seq(0,1,0.001)); return(x)}
saveRDS(list(UC_cutoff=lapply(UC_cutoff,assign_name),
            dMML_cutoff=lapply(dMML_cutoff,assign_name),
            dNME_cutoff=lapply(dNME_cutoff,assign_name)
            ),
UC_in=readRDS(UC_merge_file)
#For each possible cutoff, find the overlap
region_rand_pool=lapply(UC_in,rownames)
rm(UC_in)
cutoff_ts_only_out=readRDS('../downstream/output/mouse_analysis/UC_dNME_olap/regions_ts_only.rds')
cutoff_ts_num_output_all=list()
set.seed(12345)
for(ts in names(region_rand_pool)){
    cat("Processing:",ts,'\n')
    tt1=proc.time()[[3]]
    cutoff_ts_num_output=as.data.table(expand.grid(seq(0,1,0.002),seq(0,1,0.01),seq(0,1,0.01)))
    colnames(cutoff_ts_num_output)=c("UC","dNME","dMML")

    cutoff_ts_only_out_ts=lapply(cutoff_ts_only_out,function(x) x[[ts]])
    cutoff_ts_num_output_all[[ts]]=cbind(cutoff_ts_num_output,
    fastDoCall("rbind",
    mclapply(1:nrow(cutoff_ts_num_output),
    function(x) {
         dNME_region=cutoff_ts_only_out_ts$dNME_cutoff[[as.character(cutoff_ts_num_output$dNME[[x]])]]
         dMML_region=cutoff_ts_only_out_ts$dMML_cutoff[[as.character(cutoff_ts_num_output$dMML[[x]])]]
         UC_region=cutoff_ts_only_out_ts$UC_cutoff[[as.character(cutoff_ts_num_output$UC[[x]])]]
        #if(length(UC_region)>0){
         UC_region_rand=sample(region_rand_pool[[ts]],length(UC_region))
         p5_num=round(nrow(cutoff_ts_num_output)*seq(0,1,0.01))
        if (x %in% p5_num){cat("Finishing:",seq(0,1,0.01)*100[p5_num==x],'\\% finished\n')}
         return(data.table(        
         dNME_specific=length(setdiff(dNME_region,UC_region)),
         dMML_specific=length(setdiff(dMML_region,UC_region)),
         UC_specific_dNME=length(setdiff(UC_region,dNME_region)),
         UC_specific_dMML=length(setdiff(UC_region,dMML_region)),
         shared_dNME=length(intersect(dNME_region,UC_region)),
         shared_dMML=length(intersect(dMML_region,UC_region)),
         dNME_specific_rand=length(setdiff(dNME_region,UC_region_rand)),
         dMML_specific_rand=length(setdiff(dMML_region,UC_region_rand)),
         UC_specific_dNME_rand=length(setdiff(UC_region_rand,dNME_region)),
         UC_specific_dMML_rand=length(setdiff(UC_region_rand,dMML_region)),
         shared_dNME_rand=length(intersect(dNME_region,UC_region_rand)),
         shared_dMML_rand=length(intersect(dMML_region,UC_region_rand)),
         dNME_cutoff=cutoff_ts_num_output$dNME[[x]],
         dMML_cutoff=cutoff_ts_num_output$dMML[[x]],
         UC_cutoff=cutoff_ts_num_output$UC[[x]]
         )
         )

         #}
    },mc.cores=20)))
  
    cat("Finish processing:",ts,'in',proc.time()[[3]]-tt1,'\n')

}