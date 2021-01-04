
# Compare old UC and fixed UC ---------------------------------------------
UC_matrix_ls_old=readRDS('../UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls.rds')
UC_matrix_ls=readRDS('../UC_run_before_MDS/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
DNase=readRDS('../mm10_DNase.rds')
UC_matrix_ls_old=lapply(UC_matrix_ls_old,function(x) subsetByOverlaps(x,DNase[DNase$region_type=="DNase"],type='equal'))
UC_in_matrix_ls=lapply(UC_in_matrix_ls,function(x) subsetByOverlaps(x,DNase[DNase$region_type=="DNase"],type='equal'))
comp_UC_fix=lapply(names(UC_matrix_ls_old),function(x){
  total_regions=length(UC_matrix_ls_old[[x]])
  UC_matrix_old=as.matrix(mcols(UC_matrix_ls_old[[x]]))
  rownames(UC_matrix_old)=paste0(seqnames(UC_matrix_ls_old[[x]]),':',start(UC_matrix_ls_old[[x]]),'-',end(UC_matrix_ls_old[[x]]))
  cat(x,'\n')
  total_regions_01=sum(rowSums(UC_matrix_old>0.1,na.rm=T)>0)
  UC_matrix_old=UC_matrix_old[rowSums(UC_matrix_old==1,na.rm = T)>0,]
  UC_matrix_fix=as.matrix(mcols(UC_in_matrix_ls[[x]]))
  rownames(UC_matrix_fix)=paste0(seqnames(UC_in_matrix_ls[[x]]),':',start(UC_in_matrix_ls[[x]]),'-',end(UC_in_matrix_ls[[x]]))
  UC_matrix_fix=UC_matrix_fix[rownames(UC_matrix_old),]
  #Number of regions
  diff=rowSums(UC_matrix_old==1,na.rm = T)-rowSums(UC_matrix_fix==1,na.rm = T)
  N_region=sum(diff!=0)
  #Number of sample x regions
  N_sample_region=sum(diff)
  return(list(data.table(tissue=x,N_region=N_region,N_sample_region=N_sample_region,total_regions=total_regions,
                         total_regions_01=total_regions_01),
              regions=rownames(UC_matrix_old)))
})
saveRDS(comp_UC_fix,'../downstream/output/UC_jsd_fix.rds')
#Bug Fix
UC_fix=readRDS('../downstream/output/UC_jsd_fix.rds')
tissue=unlist(lapply(UC_fix,function(x) x[[1]]$tissue))
names(UC_fix)=tissue
tissue=sub('_.*','',dir('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/'))
csv_error=lapply(tissue,function(x) {
  csv_in=fread(paste0('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/',x,'_all.csv'))
  UC_region=UC_fix[[x]][[2]]
  return(csv_in[csv_in$region%in%UC_region])
})
names(csv_error)=tissue
csv_error_dNME=do.call(rbind,lapply(names(csv_error),function(x) return(data.table(csv_error[[x]]$dNME_maxJSD,tissue=x))))
ggplot(csv_error_dNME,aes(x=tissue,y=V1))+geom_boxplot()+theme_glob+ylab('dNME')
lapply(UC_fix,function(x) x[[1]]$N_region/x[[1]]$total_regions_01)




# density plot ------------------------------------------------------------
UC_gr=lapply(UC,function(x) {
  gr=GRanges(seqnames=sub(':.*','',rownames(x)),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',rownames(x)))),
                     end=as.numeric(sub('.*-','',rownames(x)))))
  gr$maxUC=rowMax(x)
  gr$maxUC_loc=apply(x,1,which.max)
  return(gr)
})

dnme_gr=lapply(dnme,function(x) {
  gr=GRanges(seqnames=sub(':.*','',rownames(x)),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',rownames(x)))),
                     end=as.numeric(sub('.*-','',rownames(x)))))
  mcols(gr)=x
  return(gr)
})

dmml_gr=lapply(dmml,function(x) {
  gr=GRanges(seqnames=sub(':.*','',rownames(x)),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',rownames(x)))),
                     end=as.numeric(sub('.*-','',rownames(x)))))
  mcols(gr)=x
  return(gr)
})
for(sp in names(dnme_gr)){
  dnme_gr[[sp]]$maxUC_loc=UC_gr[[sp]]$maxUC_loc
  dnme_gr[[sp]]$dNME_max_UC= apply(mcols(dnme_gr[[sp]]),1,function(x) x[x["maxUC_loc"]])
  dmml_gr[[sp]]$maxUC_loc=UC_gr[[sp]]$maxUC_loc
  dmml_gr[[sp]]$dMML_max_UC= apply(mcols(dmml_gr[[sp]]),1,function(x) x[x["maxUC_loc"]])
}

UC_enhancer=unlist(lapply(UC_gr,function(x) subsetByOverlaps(x,chromHMM)$maxUC))
UC_promoter=unlist(lapply(UC_gr,function(x) subsetByOverlaps(x,promoters)$maxUC))

dnme_enhancer=unlist(lapply(dnme_gr,function(x) subsetByOverlaps(x,chromHMM)$dNME_max_UC))
dnme_promoter=unlist(lapply(dnme_gr,function(x) subsetByOverlaps(x,promoters)$dNME_max_UC))


dmml_enhancer=unlist(lapply(dmml_gr,function(x) subsetByOverlaps(x,chromHMM)$dMML_max_UC))
dmml_promoter=unlist(lapply(dmml_gr,function(x) subsetByOverlaps(x,promoters)$dMML_max_UC))

# plot_density_dt=rbind(data.table(UC=UC_enhancer,region="enhancer"),
#                        data.table(UC=UC_promoter,region="promoter"))
# pdf('../downstream/output/graphs/Figure5/Figure5B_UC_enhancer.pdf',width=3.5,height=3.5)
# print(ggplot(plot_density_dt,aes(x=UC,fill=region))+geom_density(alpha=0.5)+theme_glob+theme(legend.position = "bottom"))
# dev.off()


pdf('../downstream/output/graphs/Figure5/Figure5B_dNME_dMML_density_enhancer.pdf',width=3.5,height=3.5)
plot_density_dt=rbind(data.table(dnme=dnme_enhancer,region="enhancer"),
                      data.table(dnme=dnme_promoter,region="promoter"))
dNME_density=ggplot(plot_density_dt,aes(x=dnme,fill=region))+geom_density(alpha=0.5)+theme_glob+theme(legend.position = "bottom")

plot_density_dt=rbind(data.table(dmml=dmml_enhancer,region="enhancer"),
                      data.table(dmml=dmml_promoter,region="promoter"))

dMML_density=ggplot(plot_density_dt,aes(x=dmml,fill=region))+geom_density(alpha=0.5)+theme_glob+theme(legend.position = "bottom")
ggarrange(dNME_density,dMML_density,nrow=2,ncol=1,common.legend=T)
dev.off()
