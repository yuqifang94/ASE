source('mainFunctions_sub.R')
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
UC_merge=readRDS(UC_merge_file)
UC_merge_cat=fastDoCall("rbind",lapply(names(UC_merge),function(x) {
            tissue_out_filtered_ts=tissue_out_filtered[[x]]
            UC_merge_ts =as.data.table(UC_merge[[x]][tissue_out_filtered_ts$region,grepl("UC-|dNME-|dMML-",colnames( UC_merge[[x]]))])
            UC_merge_ts=cbind(UC_merge_ts,tissue_out_filtered_ts[,list(region,dMML_cor,dNME_cor,tissue,cluster,region_type)])
            id_col=c("region","dMML_cor","dNME_cor","tissue","cluster","region_type")
            UC_merge_ts=melt.data.table(UC_merge_ts,id.vars=id_col,variable.name="stage")
            UC_merge_ts$stat=gsub('-.*','',UC_merge_ts$stage)
            UC_merge_ts$stage=gsub(paste0("UC-|dNME-|dMML-|-all|",x,'-'),"", UC_merge_ts$stage)
            return(dcast.data.table(UC_merge_ts,
            region+dMML_cor+dNME_cor+tissue+cluster+region_type+stage~stat,
            value.var="value")
            )


}))
saveRDS(UC_merge_cat,'../downstream/output/mouse_analysis/QC/UC_merge_cat.rds')
region_type_cor=do.call(rbind,lapply(unique(UC_merge_cat$region_type),function(rt) {
    return(rbind(
        cbind(corToDt(cor.test(UC_merge_cat[region_type==rt]$UC,UC_merge_cat[region_type==rt]$dMML)), 
            data.table(region_type=rt,correlation_type="UC-dMML")),
         cbind(corToDt(cor.test(UC_merge_cat[region_type==rt]$UC,UC_merge_cat[region_type==rt]$dNME)), 
            data.table(region_type=rt,correlation_type="UC-dNME")))
)

}))
plot_cor<-function(plot_dat,title){
    return(plot_dat+
    geom_bin2d(bins=100)+
    ggtitle(title)+ 
    scale_fill_gradient(name = "count", trans = "log10",high="#132B43",low="#56B1F7")+ylim(c(0,1))+
    geom_smooth())


}
pdf("../downstream/output/mouse_analysis/QC/correlation.pdf",height=5,width=7)
dMML_UC=plot_cor(ggplot(UC_merge_cat,aes(x=dMML,y=UC)),"dMML-UC")
dNME_UC=plot_cor(ggplot(UC_merge_cat,aes(x=dNME,y=UC)),"dNME-UC")
dNME_dMML=plot_cor(ggplot(UC_merge_cat,aes(x=dMML,y=dNME)),"dMML-dNME")
print(grid.arrange(dMML_UC,dNME_UC,dNME_dMML,nrow=2,ncol=2, top = textGrob("All regions")))
dev.off()
pdf("../downstream/output/mouse_analysis/QC/correlation_sep.pdf",height=5,width=7)
lapply(unique(UC_merge_cat$region_type),function(rt){
    dMML_UC=plot_cor(ggplot(UC_merge_cat[region_type==rt],aes(x=dMML,y=UC)),"dMML-UC")
    dNME_UC=plot_cor(ggplot(UC_merge_cat[region_type==rt],aes(x=dNME,y=UC)),"dNME-UC")
    dNME_dMML=plot_cor(ggplot(UC_merge_cat[region_type==rt],aes(x=dMML,y=dNME)),"dMML-dNME")
    print(grid.arrange(dMML_UC,dNME_UC,dNME_dMML,nrow=2,ncol=2, top = textGrob(rt)))


})
dev.off()

plot_cor_dist_single<-function(datIn){
  datIn=melt.data.table(datIn,id.vars=c("region","tissue","cluster","region_type"),measure.vars=c("dMML_cor","dNME_cor"))
  pdf(paste0("../downstream/output/mouse_analysis/QC/correlation_dist_",unique(datIn$tissue),".pdf"))
  lapply(unique(datIn$region_type),function(x) {
    print(ggplot(datIn[region_type==x],aes(x=value,color=variable))+geom_density()+xlab("correlation")+ggtitle(x)+xlim(c(0,1)))
  })
  dev.off()
}
lapply(tissue_out_filtered,plot_cor_dist_single)

#Example
pdf(paste0("../downstream/output/mouse_analysis/QC/correlation_example.pdf"))
    for(rt in unique(tissue_out_filtered$limb$region_type)){
        limb_region=tissue_out_filtered$limb[region_type==rt][order(cor_diff)]
        for(i in 1:20){
            datPlot=limb_region[i]
            UC_stat=UC_merge_cat[region==datPlot$region]
            UC_stat=melt.data.table(UC_stat,id.vars="UC",measure.vars=c("dMML","dNME"))
            print(
                ggplot(UC_stat,aes(x=value,y=UC,color=variable,group=variable))+geom_point()+geom_smooth(se=F)+ggtitle(rt)
            )

        }
    }

dev.off()