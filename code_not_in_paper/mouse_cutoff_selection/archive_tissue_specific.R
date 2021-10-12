#Convert aid_tis region into data.table
region_tissue_assign=do.call(rbind,lapply(names(ts_aid),function(x) data.table(region=ts_aid[[x]],tissue=x)))
olap=findOverlaps(convert_GR(region_tissue_assign$region,direction='GR'),convert_GR(unique(H3K27ac_output_dt$region),direction='GR'))
region_tissue_assign_enc=region_tissue_assign[queryHits(olap)]
region_tissue_assign_enc$enhancer_region=unique(H3K27ac_output_dt$region)[subjectHits(olap)]
#prepare matrix for histone marker info and RNA info
#RNA#Refine here
H3K27ac_output_dt_analyzed=H3K27ac_output_dt[region %in% region_tissue_assign_enc$enhancer_region]

H3K27ac_output_dt_analyzed$analyzed_region=region_tissue_assign_enc[match(H3K27ac_output_dt$region,enhancer_region)]$region
heatmap_dt=dcast.data.table(H3K27ac_output_dt_analyzed,analyzed_region~sample,value.var='log2FPKM',fun.aggregate=mean)
heatmap_dt_mt=as.matrix(heatmap_dt[,-1])
rownames(heatmap_dt_mt)=heatmap_dt$region

colann= data.frame(tissue=gsub('-.*','',colnames(heatmap_dt_mt)))
rownames(colann)=colnames(heatmap_dt_mt)
rowann <- data.frame(tissue_r=H3K27ac_output_dt)
rownames(rowann)=rownames(heatmap_dt_mt)
  c1 <- mouse_color()
 tiff(paste0(figure_path,'expression_tissue_specific_raw.tiff'),width=5000,height=5000,res=500)
  pheatmap(heatmap_dt_mt,cluster_rows = F,annotation_row = rowann,cluster_cols = F,
           annotation_col = colann,show_colnames = F,show_rownames = F,
           gaps_row = cumsum(rle(as.character(rowann$tissue_r))$lengths),
           annotation_colors = list(tissue=c1,tissue_r=c1)
           )
  dev.off()
#Histone

pdf(paste0(figure_path,'correlation_rank_tissue_enhancer_stacked_bar.pdf'))
all_tissue=names(ts_aid)
#Also prepare the heatmap
for(ts in all_tissue){
#Zero filtered H3K27ac_output_dt_cor
  olap=findOverlaps(convert_GR(H3K27ac_output_dt_cor$region,dir="GR"),convert_GR(ts_aid[[ts]]))
  H3K27ac_output_dt_cor_ts=H3K27ac_output_dt_cor[queryHits(olap)]
  H3K27ac_output_dt_cor_ts=H3K27ac_output_dt_cor_ts[,list(rank=rank(-cor),tissue=tissue,cor=cor),by=list(region)]
  H3K27ac_output_dt_cor_ts$tissue=factor(H3K27ac_output_dt_cor_ts$tissue,levels=c(ts,all_tissue[all_tissue!=ts]))
  plot_out=ggplot(H3K27ac_output_dt_cor_ts,aes(x=tissue,y=as.numeric(rank)))+geom_bar()+
          ggtitle(ts)+ylim(c(0,8))+coord_flip()
  print(plot_out)

}

dev.off()
#Show denisty of the background and each tissue
ts_correaltion_all=H3K27ac_output_dt_cor[region%in%ts_aid_out$enhancer_region]
pdf(paste0(figure_path,'correlation_all_enhancer.pdf'))
ggplot(ts_correaltion_all,aes(x=cor))+geom_density(color="darkblue",fill="lightblue")
dev.off()
pdf(paste0(figure_path,'correlation_tissue_all.pdf'))
    plot_dt=rbind(data.table(correlation=ts_aid_out$correlation,cor_type='tissue_specific_UC'),
                  data.table(correlation=ts_correaltion_all$cor,cor_type='all_correlation'))

 print(ggplot(plot_dt,aes(x=correlation,group=cor_type,color=cor_type))+geom_density()+ggtitle(ts))


dev.off()
t.test(plot_dt[cor_type=="tissue_specific_UC"]$correlation,plot_dt[cor_type=="all_correlation"]$correlation,alternative="greater")
  pdf(paste0(figure_path,'correlation_tissue_enhancer.pdf'))

for(ts in unique(ts_aid_out$tissue)){
    plot_dt=rbind(data.table(correlation=ts_aid_out[tissue==ts]$correlation,cor_type='tissue_specific_UC'),
                  data.table(correlation=ts_correaltion_all$cor,cor_type='all_correlation'))

 print(ggplot(plot_dt,aes(x=correlation,group=cor_type,color=cor_type))+geom_density()+ggtitle(ts))


}
dev.off()
pdf(paste0(figure_path,'correlation_tissue_rank.pdf'))
for(ts in unique(ts_aid_out$tissue)){
    plot_dt=rbind(data.table(correlation=ts_aid_out[tissue==ts]$rank_correlation,cor_type='tissue_specific_UC'),
                  data.table(correlation=as.numeric(H3K27ac_output_dt_cor_dc_rank),cor_type='all_correlation'))
    plot_dt$correlation=factor(round(plot_dt$correlation),levels=as.character(1:7))#Round tiles
    
 print(ggplot(plot_dt,aes(x=correlation,group=cor_type,fill=cor_type))+geom_bar(aes(y = ..prop..),)+
      ggtitle(ts)+xlab("correlation rank")+ylim(c(0,0.25)))


}
dev.off()

#Sanity check with Bin's paper
 H3K27ac_output_dt=readRDS(H3K27ac_output_dt_fn)
 Bin_enhancer=readRDS(bin_enhancer_rds)
 Bin_enhancer=convert_GR(Bin_enhancer,direction='DT')
 H3K27ac_output_dt_cor_check=H3K27ac_output_dt[,list(cor_pearson=cor(log2RPKM,log2FPKM,method='pearson'),
   cor_spearman=cor(log2RPKM,log2FPKM,method='spearman'),
    std_log2RPKM=sd(log2RPKM),
    std_log2FPKM=sd(log2FPKM)),by=list(region,replicates)]
H3K27ac_output_dt_cor_check=H3K27ac_output_dt_cor_check[std_log2RPKM!=0&std_log2FPKM!=0]
H3K27ac_output_dt_cor_check_rep1=H3K27ac_output_dt_cor_check[replicates==1]
H3K27ac_output_dt_cor_check_rep1$SCC=Bin_enhancer[match(H3K27ac_output_dt_cor_check_rep1$region,region)]$SCC
pdf(paste0(figure_path,'spearman_SCC_difference.pdf'))
method_diff=abs(H3K27ac_output_dt_cor_check_rep1$cor_spearman-H3K27ac_output_dt_cor_check_rep1$SCC)
hist(method_diff,xlab='difference',main=paste0("mean diff spearman:",round(mean(method_diff),digits=3)))

dev.off()
pdf(paste0(figure_path,'spearman_SCC_difference_percent.pdf'))
method_diff=H3K27ac_output_dt_cor_check_rep1[,list(percent_diff=abs((cor_spearman-SCC)/SCC))]$percent_diff
hist(method_diff,xlab='difference',main=paste0("mean diff spearman:",round(mean(method_diff),digits=3)))
dev.off()
pdf(paste0(figure_path,'spearman_SCC_difference.pdf'))
method_diff=abs(H3K27ac_output_dt_cor_check_rep1$cor_spearman-H3K27ac_output_dt_cor_check_rep1$SCC)
hist(method_diff,xlab='difference',main=paste0("mean diff spearman:",round(mean(method_diff),digits=3)))
dev.off()
pdf(paste0(figure_path,'spearman_SCC_difference_percent.pdf'))
method_diff=H3K27ac_output_dt_cor_check_rep1[,list(percent_diff=abs((cor_spearman-SCC)/SCC))]$percent_diff
hist(method_diff,xlab='difference',main=paste0("mean diff spearman:",round(mean(method_diff),digits=3)))
dev.off()
pdf(paste0(figure_path,'example_correlation.pdf'))
example_dt=H3K27ac_output_dt[region=="chr1:6733000-6735000"&tissue=="EFP"]
plot(example_dt$log2RPKM,example_dt$log2FPKM)
dev.off()
pdf(paste0(figure_path,'log2FPKM_dis.pdf'))
hist(H3K27ac_output_dt$log2FPKM)
dev.off()