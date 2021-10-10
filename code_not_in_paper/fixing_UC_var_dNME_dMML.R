# nolint start
setwd('../')
source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
dNME_in=mclapply(UC_in,function(x) x[,grepl('dNME',colnames(x))],mc.cores=7)
dMML_in=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=7)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=7)
#Data generation

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


#Fixing dNME and dMML cutoff, vary UC cutoff
#Try different cut of dNME and see overlap
venn_out_all_UC_all=data.table()
cut_dMML=0.134375
cut_dNME=0.31875
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


saveRDS(venn_out_all_UC_all,'../downstream/output/mouse_analysis/UC_dNME_olap/UC_var_venn_out_all_dNME_dMML_cutoff.rds')
# nolint end