library(GenomicRanges)
library(matrixStats)
library(data.table)
library(Gmisc)
library(ggplot2)
d <- readRDS('../downstream/output/UC_agnostic_mouse_dnase.rds')
mml <- readRDS('../downstream/output/mml.rds')
nme <- readRDS('../downstream/output/nme.rds')
#subset DNase region
chromHMM_enhancer=readRDS('../downstream/output/mm10_chromHMM_enhancer.rds')
#Max dNME, dMML
stat_in=mml
#Find dNME between any pair for each tissue
dNME_max=lapply(names(d),function(ts,mat_in){
  mat_in=mat_in[,gsub('-.*','',colnames(mat_in))==ts]
  mat_diff_all=matrix(nrow=nrow(mat_in),ncol = 0)
  rownames(mat_diff_all)=rownames(mat_in)
  nsp=ncol(mat_in)
  mat_diff_all_cn=c()
  for (i in 1:(nsp-1)){
    for (j in (i+1):nsp){
      mat_diff_all=cbind(mat_diff_all,matrix(abs(mat_in[,i]-mat_in[,j])))
      stage1=sub('-all','',sub(ts,'',colnames(mat_in)[i]))
      stage2=sub('-all','',sub(ts,'',colnames(mat_in)[j]))
      mat_diff_all_cn=c(mat_diff_all_cn,paste0(ts,stage1,stage2))
    }
  }
  colnames(mat_diff_all)=mat_diff_all_cn
  mat_out_gr=GRanges(seqnames=sub(':.*','',rownames(mat_diff_all)),
                     IRanges(start= as.numeric(sub('.*:','',sub('-.*','',rownames(mat_diff_all)))),
                                                                           end=as.numeric(sub('.*-','',rownames(mat_diff_all)))))
  mat_out_gr$maxstat=rowMaxs(mat_diff_all)
  mat_out_gr$tissue=ts
  return(mat_out_gr)
},mat_in=stat_in)
#saveRDS(dNME_max,'../downstream/output/dMML_max_mm10.rds')

#Max UC
UC_max=lapply(names(d),function(ts,mat_in){
  mat_in=mat_in[[ts]]
  mat_out_gr=GRanges(seqnames=sub(':.*','',rownames(mat_in)),IRanges(start= as.numeric(sub('.*:','',sub('-.*','',rownames(mat_in)))),
                                                                           end=as.numeric(sub('.*-','',rownames(mat_in)))))
  mat_out_gr$maxstat=rowMaxs(mat_in)
  mat_out_gr$tissue=ts
  return(mat_out_gr)
},mat_in=d)

gtf <- fread('../downstream/input/grcm38.gtf.gz',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
names(gr) <- gn
pro <- promoters(gr,upstream=2000,downstream=1000)
seqlevels(pro)=paste0("chr",seqlevels(pro))
dNME_max_enc=lapply(dNME_max,function(gr_in,enhancer,promoter){
  olap=findOverlaps(gr_in,enhancer[enhancer$tissue==unique(gr_in$tissue)])
  gr_in$enhancer=""
  gr_in$enhancer[queryHits(olap)]="enhancer"
  olap=findOverlaps(gr_in,pro)
  gr_in$promoter=""
  gr_in$promoter[queryHits(olap)]="promoter"
  return(gr_in)
},enhancer=chromHMM_enhancer,promoter=pro)

dNME_max_enc=fastDoCall('c',dNME_max_enc)
dNME_max_enc=dNME_max_enc[!is.na(dNME_max_enc$maxstat)]
dNME_max_enc$anno=paste0(dNME_max_enc$enhancer,dNME_max_enc$promoter)
dNME_max_enc$anno[dNME_max_enc$anno=="enhancerpromoter"]="enhancer"
dNME_max_enc$anno[dNME_max_enc$anno==""]="Not enhancer/promoter"
ecdf_df_all=data.frame()
for (anno in unique(dNME_max_enc$anno)){
  ecdf_df_all=rbind(ecdf_df_all,
                    data.frame(max_stat=seq(0,1,0.001),
                               estimate=ecdf(dNME_max_enc$maxstat[dNME_max_enc$anno==anno])(seq(0,1,0.001)),
                               anno=anno))
  
}
pdf('../downstream/output/max_UC_enhancer.pdf')
ggplot(ecdf_df_all,aes(x=max_stat,y=estimate,color=anno))+geom_line(size=1)+xlab('Max dMML')+ylab('Cumulative prob')+
  theme(legend.position = "bottom")
dev.off()
