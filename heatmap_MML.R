library(RColorBrewer)
library(pheatmap)
library(gplots)
source('mainFunctions_sub.R')
MML_in=readRDS('../downstream/input/mouse_analysis/MML_matrix_mouse_all_dedup_N2_all_regions.rds')
MML_in_dt=convert_GR(MML_in,direction='matrix')
colnames(MML_in_dt)=gsub('-all','',colnames(MML_in_dt))
cluster_region_out_fn='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_1/cluster_all_region_assignment_filtered_0_1.rds'
clu=readRDS(cluster_region_out_fn)
clu_all=do.call(rbind,lapply(clu,function(x) x[order(as.numeric(cluster))]))
MML_in_dt=MML_in_dt[clu_all$regions,]
rowann_out=data.frame(tissue_r=clu_all$tissue,cluster=clu_all$cluster)
row_gap=cumsum(rle(clu_all$tissue)$lengths)
rownames(rowann_out)=rownames(MML_in_dt)

#Refine plotting parameters
colann <- data.frame(time=sub('.*-','',colnames(MML_in_dt)),tissue=sub('-.*','',colnames(MML_in_dt)),stringsAsFactors = F)
rownames(colann) <- colnames(MML_in_dt)
c1 <- mouse_color()
c2 <- brewer.pal(10,'Set3')
names(c2) <- 1:10
c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
names(c4) <- sort(unique(colann[,1]))
# tiff(paste0('../downstream/output/heatmap_acrosstissue/',n,'.tiff'),width=3000,height=3000,res=300)
# #png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
# pheatmap(mat,cluster_rows = F,annotation_row = rowann,cluster_cols = F,
#          annotation_col = colann,show_colnames = F,show_rownames = F,
#          gaps_row = row_gap,gaps_col = cumsum(rle(colann[,1])$lengths),
#          annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4,dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)))
# dev.off()

tiff('../downstream/output/mouse_analysis/clustering/all_sc_N17_ft_kmeans_10run_filtered_all_MML_non_scale.tiff',width=5000,height=5000,res=300)
#png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
pheatmap(MML_in_dt,cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
         annotation_col = colann,show_colnames = F,show_rownames = F,
         gaps_row = row_gap,gaps_col = cumsum(rle(colann[,2])$lengths),
         annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4
                                  #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)
         ))
dev.off()


# Non-liver ---------------------------------------------------------------

MML_in=readRDS('../downstream/input/mouse_analysis/MML_matrix_mouse_all_dedup_N2_all_regions.rds')
MML_in_dt=convert_GR(MML_in,direction='matrix')
colnames(MML_in_dt)=gsub('-all','',colnames(MML_in_dt))
cluster_region_out_fn='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_1/cluster_all_region_assignment_filtered_0_1.rds'
clu=readRDS(cluster_region_out_fn)
clu_all=do.call(rbind,lapply(clu,function(x) x[order(as.numeric(cluster))]))
MML_in_dt=MML_in_dt[clu_all$regions,!grepl('liver',colnames(MML_in_dt))]
rowann_out=data.frame(tissue_r=clu_all$tissue,cluster=clu_all$cluster)
row_gap=cumsum(rle(clu_all$tissue)$lengths)
rownames(rowann_out)=rownames(MML_in_dt)

#Refine plotting parameters
colann <- data.frame(time=sub('.*-','',colnames(MML_in_dt)),tissue=sub('-.*','',colnames(MML_in_dt)),stringsAsFactors = F)
rownames(colann) <- colnames(MML_in_dt)
c1 <- mouse_color()
c2 <- brewer.pal(10,'Set3')
names(c2) <- 1:10
c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
names(c4) <- sort(unique(colann[,1]))
# tiff(paste0('../downstream/output/heatmap_acrosstissue/',n,'.tiff'),width=3000,height=3000,res=300)
# #png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
# pheatmap(mat,cluster_rows = F,annotation_row = rowann,cluster_cols = F,
#          annotation_col = colann,show_colnames = F,show_rownames = F,
#          gaps_row = row_gap,gaps_col = cumsum(rle(colann[,1])$lengths),
#          annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4,dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)))
# dev.off()

tiff('../downstream/output/mouse_analysis/clustering/all_sc_N17_ft_kmeans_10run_filtered_all_MML_no_liver.tiff',width=5000,height=5000,res=300)
#png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
pheatmap(scalematrix(MML_in_dt),cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
         annotation_col = colann,show_colnames = F,show_rownames = F,
         gaps_row = row_gap,gaps_col = cumsum(rle(colann[,2])$lengths),
         annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4
                                  #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)
         ))
dev.off()

#MML only & Both
region_cat=readRDS('../downstream/output/mouse_analysis/correlation/UC_01_cutoff/tissue_out_N17_kmeans_10run_filtered_all_region.rds')
region_cat_MML=lapply(region_cat,function(x) x[region_type%in%c("Both","MML only")])
region_cat_MML=do.call(rbind,region_cat_MML)
MML_in=readRDS('../downstream/input/mouse_analysis/MML_matrix_mouse_all_dedup_N2_all_regions.rds')
MML_in_dt=convert_GR(MML_in,direction='matrix')
colnames(MML_in_dt)=gsub('-all','',colnames(MML_in_dt))
cluster_region_out_fn='../downstream/output/mouse_analysis/clustering/tissue_specific/UC_0_1/cluster_all_region_assignment_filtered_0_1.rds'
clu=readRDS(cluster_region_out_fn)
clu=lapply(clu,function(x) x[regions %in%region_cat_MML$region])
clu_all=do.call(rbind,lapply(clu,function(x) x[order(as.numeric(cluster))]))
MML_in_dt=MML_in_dt[clu_all$regions,]
rowann_out=data.frame(tissue_r=clu_all$tissue,cluster=clu_all$cluster)
row_gap=cumsum(rle(clu_all$tissue)$lengths)
rownames(rowann_out)=rownames(MML_in_dt)

#Refine plotting parameters
colann <- data.frame(time=sub('.*-','',colnames(MML_in_dt)),tissue=sub('-.*','',colnames(MML_in_dt)),stringsAsFactors = F)
rownames(colann) <- colnames(MML_in_dt)
c1 <- mouse_color()
c2 <- brewer.pal(10,'Set3')
names(c2) <- 1:10
c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
names(c4) <- sort(unique(colann[,1]))
# tiff(paste0('../downstream/output/heatmap_acrosstissue/',n,'.tiff'),width=3000,height=3000,res=300)
# #png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
# pheatmap(mat,cluster_rows = F,annotation_row = rowann,cluster_cols = F,
#          annotation_col = colann,show_colnames = F,show_rownames = F,
#          gaps_row = row_gap,gaps_col = cumsum(rle(colann[,1])$lengths),
#          annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4,dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)))
# dev.off()

tiff('../downstream/output/mouse_analysis/clustering/all_sc_N17_ft_kmeans_10run_filtered_all_MML_Both_MML_only.tiff',width=5000,height=5000,res=300)
#png(paste0('/dcl01/hongkai/data/zji4/ase/mouse/plot/heatmap/combine_nosubcluster/heatmap_acrosstissue/',n,'.png'),width = 800,height=800,res=300)
pheatmap(scalematrix(MML_in_dt),cluster_rows = F,annotation_row = rowann_out,cluster_cols = F,
         annotation_col = colann,show_colnames = F,show_rownames = F,
         gaps_row = row_gap,gaps_col = cumsum(rle(colann[,2])$lengths),
         annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4
                                  #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10)
         ))
dev.off()

