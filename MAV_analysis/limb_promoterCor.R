source('mainFunctions_sub.R')
library(Seurat)
library(SeuratDisk)

limb <- LoadH5Seurat("../downstream/data/mouseLimb/mouse_limb_scRNA_10x.h5seurat")#Need at least 200G
limb_RNA_data=as.matrix(GetAssayData(object = limb, slot = "data"))#dim = 43346 90637
limb_RNA_data=limb_RNA_data[(rowSums(limb_RNA_data!=0))>(ncol(limb_RNA_data)/1000),]#17530 90637
#Pre processing and MAV calc
limb_cellType=fread("../downstream/data/mouse_scRNA/meta.tsv")
limb_cellType$merged_celltype=trimws(gsub(paste(1:10,collapse="|"),"",limb_cellType$cell_type))
limb_cellType$cellType_stage=paste0(limb_cellType$merged_celltype,"-E",limb_cellType$stage)
limb_cellType=limb_cellType[stage!=13.5]
cellType_stage_uq=unique(limb_cellType$cellType_stage)
cellType_stageMAV<-lapply(cellType_stage_uq,function(x) {
    return(limb_RNA_data[,limb_cellType[cellType_stage==x]$index])

})
names(cellType_stageMAV)=cellType_stage_uq
#Get cell number of each cellType stage
limb_cellType_nCell=limb_cellType[,list(nCell=length(index)),by=list(merged_celltype,stage)][order(nCell,decreasing=T)]
#Find a cell type that have most nCell there
limb_cellType_nCellAvg=limb_cellType_nCell[,list(meanNCell=mean(nCell)),by=list(merged_celltype)]
#Mesenchymal
cellTypeStage_Mesenchymal_MAV=cellTypeMAVcalc(cellType_stageMAV,'Mesenchymal')

#Adding NME
NME_mm10=readRDS(NME_matrix_file)
NME_mm10=convert_GR(NME_mm10,direction="matrix")
cor_dt_filtered=readRDS(tissue_out_filtered_fn)
NME_mm10_limb=NME_mm10[cor_dt_filtered$limb$region,grepl("limb",colnames(NME_mm10))]
NME_mm10_limb=as.data.table(NME_mm10_limb)
NME_mm10_limb$region=cor_dt_filtered$limb$region
NME_mm10_limb=melt.data.table(NME_mm10_limb,id.vars="region")
colnames(NME_mm10_limb)=c("region","stage","NME")
NME_mm10_limb$stage=gsub("limb-|-all","",NME_mm10_limb$stage)
tss_mm10=get_mm10_tss()

plotHistCor<-function(NMEIn,MAV,tss){
    histAll=list()
    for (maxGap in c(1000,2000,3000,5000,10000,20000)){
        NMEInAnnotated=promAnno(NMEIn,tss,maxGap=maxGap)

        NMEInAnnotatedCor=corCalc(NMEInAnnotated,MAV,permute=F,seed=0,NMEVar="NME",MAVVar="hypervar_logvar")
        histAll[[as.character(maxGap)]]=ggplot(NMEInAnnotatedCor, aes(x=NME_MAV_cor))+geom_histogram(color="black", fill="white")+
                xlab("correlation")+
                ggtitle(paste0("correaltion at: ",maxGap," from TSS \n"))

    
    }
    return(histAll)

}
mm10Hist=plotHistCor(NME_mm10_limb,cellTypeStage_Mesenchymal_MAV,tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_mes_nearest.pdf')
    grid.arrange(grobs =mm10Hist, nrow = 3)
dev.off()
mm10HistAll=plotHistCor(NME_mm10_limb,MAV10x,tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_SAVER_nearest.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()
mm10HistAllnoSAVER=plotHistCor(NME_mm10_limb,cellType_stageOnlyMAV[stage!="E13.5"],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_noSAVER_nearest.pdf')
    grid.arrange(grobs =mm10HistAllnoSAVER, nrow = 3)
dev.off()







#Functions
encAnno<-function(NME,enhancer){
      #Get the annotated genes
    olap=findOverlaps(convert_GR(NME$region,direction="GR"),enhancer)
    NME_annotated=NME[queryHits(olap)]
    NME_annotated$gene=enhancer$`Target Gene`[subjectHits(olap)]
    return(NME_annotated)
}
promAnno<-function(NME,promoter,maxGap=2000){
      #Get the annotated genes
    #olap=findOverlaps(convert_GR(NME$region,direction="GR"),promoter,maxgap=maxGap)
    #NME_annotated=NME[queryHits(olap)]

    nearestDist=as.data.table(distanceToNearest(convert_GR(NME$region,direction="GR"),promoter))
    NME$dist=nearestDist$distance
    NME$gene=names(promoter)[nearestDist$subjectHits]
    return(NME[abs(dist)<=maxGap])
}
enhancerCorCalc<-function(NME,enhancer,MAV,permute=F,seed=0,NMEVar="NME",MAVVar="hypervar_logvar"){
    NME_annotated=encAnno(NME,enhancer)
    NME_annotated_cor=corCalc(NME_annotated,MAV,permute,seed,NMEVar,MAVVar)
    return(NME_annotated_cor)
}
corCalc<-function(NME_annotated,MAV,permute=F,seed=0,NMEVar="NME",MAVVar="hypervar_logvar"){
   MAV$hypervar_logvar=MAV[[MAVVar]]
    NME_annotated$NME=NME_annotated[[NMEVar]]   
    NME_annotated$stage_relax=gsub('\\.5',"",NME_annotated$stage)
    NME_annotated$stage_gene=paste0(NME_annotated$stage_relax,'-',NME_annotated$gene)
    #relax about half day of the stage
    MAV$stage_relax=gsub('\\.5',"",MAV$stage)
    MAV$stage_gene=paste0(MAV$stage_relax,'-',MAV$gene)
    MAV[stage_gene%in%NME_annotated$stage_gene]#10554
    NME_annotated$MAV=MAV[match(NME_annotated$stage_gene,stage_gene)]$hypervar_logvar
    NME_annotated_MAV=NME_annotated[!is.na(MAV)]
    NME_annotated_count=NME_annotated_MAV[,list(nStage=length(unique(stage))),by=list(region)]
    NME_annotated_MAV_5stage=NME_annotated_MAV[region%in%NME_annotated_count[nStage>=5]$region]
    #NME_annotated_MAV_5stage_sdMAV=NME_annotated_MAV_5stage[,list(sdMAV=sd(MAV)),by=list(region)]
    NME_annotated_MAV_5stage=NME_annotated_MAV_5stage[,list(NME=mean(NME),MAV=MAV),by=list(gene,stage_relax)]
    print(NME_annotated_MAV_5stage)
    if(permute){
        set.seed(seed)
            NME_annotated_cor=NME_annotated_MAV_5stage[,list(NME_MAV_cor=cor(NME,sample(MAV),method="pearson")),
                    by=list(gene)]#MAV change median= 0.03
    }
    else{
    NME_annotated_cor=NME_annotated_MAV_5stage[,list(NME_MAV_cor=cor(NME,MAV,method="pearson")),
                    by=list(gene)]#MAV change median= 0.03
    }
    return(NME_annotated_cor)
}