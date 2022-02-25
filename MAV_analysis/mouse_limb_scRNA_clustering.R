#JHPCE only
#ml hdf5
#R CMD javareconf -e
source('mainFunctions_sub.R')
library(Seurat)
library(SeuratDisk)

#Convert("../downstream/data/mouseLimb/mouse_limb_scRNA_10x.h5ad", "../downstream/data/mouse_limb_scRNA_10x.h5seurat",assay="RNA")
#This is from: The changing mouse embryo transcriptome at whole tissue and single-cell resolution
#Trying the same method: 
#10x scRNA-seq clustering and t-SNE visualization
#Marker gene identification for C1 and 10x scRNA-seq data
#limb <- LoadH5Seurat("../downstream/data/mouseLimb/mouse_limb_scRNA_10x.h5seurat")#Need at least 200G
#Use SAVER to impute data
#limb_RNA_count=as.matrix(GetAssayData(object = limb, slot = "counts"))#dim = 43346 90637
#saveRDS(limb_RNA_count,"../downstream/data/mouseLimb/mouse_limb_scRNA_10x_rawCount.rds")
#See saver load
#It seems SAVER can only be on R/3.6.1 with glmnet v3.0
#They already annotated cell types
limb_cellType=fread("../downstream/data/mouse_scRNA/meta.tsv")
# Retrieve or set data in an expression matrix ('counts', 'data', and 'scale.data')
#limb_RNA_data=as.matrix(GetAssayData(object = limb, slot = "data"))#dim = 43346 90637
#limb_RNA_data=limb_RNA_data[(rowSums(limb_RNA_data!=0))>(ncol(limb_RNA_data)/1000),]#17530 90637
#row: gene, col: sample, the data is already filtered
#make a list for each cell type and get cell - expression matrix
#Reading in mouse data from Jason's saver

dirSAVER10x= '../downstream/data/mouseLimb/10x'

limb_cellType$merged_celltype=trimws(gsub(paste(1:10,collapse="|"),"",limb_cellType$cell_type))
limb_cellType$cellType_stage=paste0(limb_cellType$merged_celltype,"-E",limb_cellType$stage)
cellType_stage_uq=unique(limb_cellType$cellType_stage)
cellType_stageMAV<-lapply(cellType_stage_uq,function(x) {
    return(limb_RNA_data[,limb_cellType[cellType_stage==x]$index])

})
names(cellType_stageMAV)=cellType_stage_uq
#Get cell number of each cellType stage
limb_cellType_nCell=limb_cellType[,list(nCell=length(index)),by=list(merged_celltype,stage)][order(nCell,decreasing=T)]
#Find a cell type that have most nCell there
limb_cellType_nCellAvg=limb_cellType_nCell[,list(meanNCell=mean(nCell)),by=list(merged_celltype)]
#Go for Mesenchymal first
cellTypeStageMAV_Mesenchymal=cellType_stageMAV[grepl("Mesenchymal",names(cellType_stageMAV))]
#MAV calc
MAV_calc<-function(expr){
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  var=var[var!=0]
  gene_non0Var=names(var)
  m=m[gene_non0Var]
  lm=lm[gene_non0Var]

  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.table(gene=gene_non0Var,mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  return(res)
}
cellTypeStage_Mesenchymal_MAV=lapply(cellTypeStageMAV_Mesenchymal,MAV_calc)
cellTypeStage_Mesenchymal_MAV=lapply(names(cellTypeStage_Mesenchymal_MAV),function(x) {
                                    dtIn= cellTypeStage_Mesenchymal_MAV[[x]]
                                    dtIn$cellType_stage=x
                                    return(dtIn)
                                    
        })
cellTypeStage_Mesenchymal_MAV=do.call(rbind,cellTypeStage_Mesenchymal_MAV)
cellTypeStage_Mesenchymal_MAV$stage=gsub('.*-','',cellTypeStage_Mesenchymal_MAV$cellType)
cellTypeStage_Mesenchymal_MAV_sd=cellTypeStage_Mesenchymal_MAV[,list(sdMAV=sd(hypervar_logvar)),list(gene)]
cellTypeStage_Mesenchymal_MAV_sd[order(sdMAV,decreasing=T)][1:10]
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/highSdMAV_gene.pdf')
for (genes in cellTypeStage_Mesenchymal_MAV_sd[order(sdMAV,decreasing=T)][1:10]$gene){
    print(ggplot(cellTypeStage_Mesenchymal_MAV[gene==genes],aes(x=stage,y=hypervar_logvar))+geom_line()+geom_point()+ggtitle(genes))


}
dev.off()
cellTypeStage_Mesenchymal_MAV_dc=dcast.data.table(gene~stage,value.var ="hypervar_logvar",data=cellTypeStage_Mesenchymal_MAV)
NME_mm10=readRDS(NME_matrix_file)
NME_mm10=convert_GR(NME_mm10,direction="matrix")
cor_dt_filtered=readRDS(tissue_out_filtered_fn)
enhancer=readRDS(bin_enhancer_rds)#21441
NME_mm10_limb=NME_mm10[cor_dt_filtered$limb$region,grepl("limb",colnames(NME_mm10))]
NME_mm10_limb=as.data.table(NME_mm10_limb)
NME_mm10_limb$region=cor_dt_filtered$limb$region
NME_mm10_limb=melt.data.table(NME_mm10_limb,id.vars="region")
colnames(NME_mm10_limb)=c("region","stage","NME")
NME_mm10_limb$stage=gsub("limb-|-all","",NME_mm10_limb$stage)
#Get the annotated genes
olap=findOverlaps(convert_GR(NME_mm10_limb$region,direction="GR"),enhancer)
NME_mm10_limb_enhancer=NME_mm10_limb[queryHits(olap)]
NME_mm10_limb_enhancer$gene=enhancer$`Target Gene`[subjectHits(olap)]
NME_mm10_limb_enhancer$stage_gene=paste0(NME_mm10_limb_enhancer$stage,'-',NME_mm10_limb_enhancer$gene)
#relax about half day of the stage
cellTypeStage_Mesenchymal_MAV$stage_relax=paste0(gsub('.5',"",cellTypeStage_Mesenchymal_MAV$stage),".5")
cellTypeStage_Mesenchymal_MAV$stage_gene=paste0(cellTypeStage_Mesenchymal_MAV$stage_relax,'-',cellTypeStage_Mesenchymal_MAV$gene)
cellTypeStage_Mesenchymal_MAV[stage_gene%in%NME_mm10_limb_enhancer$stage_gene]#10554
NME_mm10_limb_enhancer$MAV=cellTypeStage_Mesenchymal_MAV[match(NME_mm10_limb_enhancer$stage_gene,stage_gene)]$hypervar_logvar
NME_mm10_limb_enhancer_MAV=NME_mm10_limb_enhancer[!is.na(MAV)]
NME_mm10_limb_enhancer_count=NME_mm10_limb_enhancer_MAV[,list(nStage=length(unique(stage))),by=list(region)]
NME_mm10_limb_enhancer_MAV_5stage=NME_mm10_limb_enhancer_MAV[region%in%NME_mm10_limb_enhancer_count[nStage>=5]$region]
NME_mm10_limb_enhancer_MAV_5stage_sdMAV=NME_mm10_limb_enhancer_MAV_5stage[,list(sdMAV=sd(MAV)),by=list(region)]
NME_mm10_limb_enhancer_cor=NME_mm10_limb_enhancer_MAV_5stage[,list(NME_MAV_cor=cor(NME,MAV,method="spearman"),gene=unique(gene)),
                by=list(region)]#MAV change median= 0.03

pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor.pdf')
print(ggplot(NME_mm10_limb_enhancer_cor,aes(x=NME_MAV_cor))+geom_density())
dev.off()
write.csv(NME_mm10_limb_enhancer_cor[,list(mean_cor=mean(NME_MAV_cor)),by=gene][order(mean_cor,decreasing=T)],
'../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_mean.csv',row.names=F)
#GO terms that are related to limb?
GO_out=readRDS(paste0(GO_01_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_01_enhancer.rds'))
GO_out_limb=GO_out$`NME only`$limb
GO_out_limb=do.call(rbind,lapply(GO_out_limb,function(x) x[[1]][FC>1.5&FDR>0.2]))
sig_term=c("negative regulation of actin filament bundle assembly",
            "negative regulation of stress fiber assembly")
           # "positive regulation of chondrocyte differentiation",
            #"establishment of epithelial cell polarity" )
sig_term_gene=unique(unlist(strsplit(GO_out_limb[Term %in% sig_term]$genes,';')))
NME_mm10_limb_enhancer_cor_sig_gene=NME_mm10_limb_enhancer_cor[gene%in%sig_term_gene][order(-NME_MAV_cor)]
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_vs_MAV.pdf')
for (i in 1:10){
    gene=NME_mm10_limb_enhancer_cor_sig_gene[i]$gene
    regionIn=NME_mm10_limb_enhancer_cor_sig_gene[i]$region
    correlation=NME_mm10_limb_enhancer_cor_sig_gene[i]$NME_MAV_cor
    GO_terms=paste(unique(GO_out_limb[Term %in% sig_term][grepl(gene,GO_out_limb[Term %in% sig_term]$genes)]$Term),collapse="\n")
    print(ggplot(NME_mm10_limb_enhancer[!is.na(MAV)&region==regionIn],aes(x=NME,y=MAV))+geom_point()+
            ggtitle(paste0(gene,": correlation=",correlation,"\n GO Terms:",GO_terms,"\n region:",regionIn)))
}
dev.off()
#Top GO terms for correlation > 0.5
NME_mm10_limb_enhancer_cor_mean=NME_mm10_limb_enhancer_cor[,list(mean_cor=mean(NME_MAV_cor)),by=gene][order(mean_cor,decreasing=T)]
high_cor_gene_GO=GO_run(NME_mm10_limb_enhancer_cor_mean[mean_cor>0.5]$gene,NME_mm10_limb_enhancer_cor_mean$gene,cluster=1)
high_cor_gene_GO[grepl("mesenchymal",Term)&FC>1.5]
mesen_term_gene=unique(unlist(strsplit(high_cor_gene_GO[grepl("mesenchymal",Term)&FC>1.5]$genes,';')))
NME_mm10_limb_enhancer_cor_mesen_gene=NME_mm10_limb_enhancer_cor[gene%in%mesen_term_gene][order(-NME_MAV_cor)]
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_vs_MAV_mesenchymal_related.pdf')
for (i in 1:23){
    gene=NME_mm10_limb_enhancer_cor_mesen_gene[i]$gene
    regionIn=NME_mm10_limb_enhancer_cor_mesen_gene[i]$region
    correlation=NME_mm10_limb_enhancer_cor_mesen_gene[i]$NME_MAV_cor
    GO_terms=paste(unique(high_cor_gene_GO[grepl("mesenchymal",Term)&FC>1.5][grepl(gene,high_cor_gene_GO[grepl("mesenchymal",Term)&FC>1.5]$genes)]$Term),collapse="\n")
    print(ggplot(NME_mm10_limb_enhancer[!is.na(MAV)&region==regionIn],aes(x=NME,y=MAV))+geom_point()+
            ggtitle(paste0(gene,": correlation=",correlation,"\n GO Terms:",GO_terms,"\n region:",regionIn)))
}
dev.off()
mm10CpG=getCpgSitesmm10()
 NME_mm10_limb_enhancer_cor_mesen_gene$NCG=countOverlaps(convert_GR(NME_mm10_limb_enhancer_cor_mesen_gene$region,direction="GR"),mm10CpG)
  NME_mm10_limb_enhancer_cor_sig_gene$NCG=countOverlaps(convert_GR(NME_mm10_limb_enhancer_cor_sig_gene$region,direction="GR"),mm10CpG)
