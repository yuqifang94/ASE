library(GenomicRanges)
library(dplyr)
library(glue)
library(ComplexHeatmap)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(matrixStats)


motif_analysis<-function(data_median_com_DNase,data_median_com_control,tissue,exclude=T){
	##exclude P0 time point and tissues not used
	sample_select_DNase=grepl(tissue,colnames(data_median_com_DNase))
	sample_select_control=grepl(tissue,colnames(data_median_com_control))
	if(exclude){
		sample_select_DNase=!sample_select_DNase
		sample_select_control=!sample_select_control
	}

	
	data_median_com_DNase <- data_median_com_DNase[,sample_select_DNase]
	data_median_com_control <- data_median_com_control[,sample_select_control]

	data_median_com_filter <- cbind(data_median_com_control,data_median_com_DNase)

	##get motif order based on the row mean difference between regulatory and non-regulatory DNA regions
	motif_NME_mean_DNase_combine_ave <- rowMeans(data_median_com_DNase, na.rm=TRUE)
	motif_NME_mean_control_combine_ave <- rowMeans(data_median_com_control, na.rm=TRUE)

	motif_order <- sort(motif_NME_mean_DNase_combine_ave - motif_NME_mean_control_combine_ave, decreasing = TRUE, index.return=TRUE)$x %>% names()

	##standardize the NME values across regulatory and non-regulatory DNA regions for each motif
	data_median_com_filter_sd <- apply(data_median_com_filter, 1, scale) %>% t()
	dimnames(data_median_com_filter_sd) <- dimnames(data_median_com_filter)
	colnames(data_median_com_filter_sd) <- sub("-all","",colnames(data_median_com_filter_sd))
	data_median_com_filter_sd <- data_median_com_filter_sd[motif_order,]

	##generate heatmap for all motifs
	tissue_group <- colnames(data_median_com_filter) %>% sub("-.*","",.)
	region_group <- rep(c("non-regulatory","regulatory"),each=46)

	##perform differential test between motif sites at non-regulatory DNA and regulatory DNA regions for each motif
	data_all <- data.frame(group = rep(c("non-regulatory","regulatory"), each = sum(sample_select_DNase)), t(data_median_com_filter))

	data_test_p <- lapply(colnames(data_all[,-1]), function(t){
	data_sub <- data_all[,c("group",t)]
		colnames(data_sub) <- c("group","NME")
		p_value <- wilcox.test(data_sub[data_sub$group=="non-regulatory","NME"],data_sub[data_sub$group=="regulatory","NME"],paired = TRUE)$p.value
	})

	data_test_p <- Reduce(c, data_test_p)
	data_test_FDR <- p.adjust(data_test_p, method="BH")

	data_all_stat <- colMeans(data_all[data_all$group == "regulatory", -1], na.rm=T) - colMeans(data_all[data_all$group == "non-regulatory", -1], na.rm=T)
	data_all_stat <- data.table(motif = names(data_all_stat),TF=sapply(names(data_all_stat), function(x) strsplit(x,"_")[[1]][2]),
																diff=data_all_stat, p_value=data_test_p, FDR=data_test_FDR)

	return(data_all_stat)
}

mouse_dir="../downstream/output/mouse_analysis/mouse_motif_Ken/"
mouse_non_tissue="brain|heart|liver|Lung|intestine|stomach|NT|kidney|P0"
mouse_median_com_DNase <- readRDS(paste0(mouse_dir,"mouse_motif_DNase_median.rds"))
mouse_median_com_DNase=mouse_median_com_DNase[,!grepl("P0",colnames(mouse_median_com_DNase))]
mouse_median_com_control <- readRDS(paste0(mouse_dir,"mouse_motif_control_median.rds"))
mouse_median_com_control=mouse_median_com_control[,!grepl("P0",colnames(mouse_median_com_control))]
mouse_not_shared_tissue=motif_analysis(mouse_median_com_DNase,mouse_median_com_control,tissue=mouse_non_tissue,exclude=T)
mouse_shared_tissue=motif_analysis(mouse_median_com_DNase,mouse_median_com_control,tissue="brain|heart|liver",exclude=F)

human_dir="../downstream/output/human_analysis/Ken_motif/homogeneous/"
human_median_com_DNase <- readRDS(paste0(human_dir,"human_motif_DNase_median_NME.rds"))
human_median_com_control <- readRDS(paste0(human_dir,"human_motif_control_median_NME.rds"))
human_tissue="brain|Brain|Ventricle|Aorta|Atrium|Liver"
human_not_shared_tissue=motif_analysis(human_median_com_DNase,human_median_com_control,tissue=human_tissue,exclude=T)
human_shared_tissue=motif_analysis(human_median_com_DNase,human_median_com_control,tissue=human_tissue,exclude=F)

sum(mouse_shared_tissue[FDR<=0.05]$motif%in%human_shared_tissue[FDR<=0.05]$motif)
sum(mouse_not_shared_tissue[FDR<=0.05]$motif%in%human_not_shared_tissue[FDR<=0.05]$motif)