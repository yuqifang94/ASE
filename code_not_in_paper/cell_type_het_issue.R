source('mainFunctions_sub.R')

heatmapNME_markerGene<-function(NME_mm10_tissue,tissue,mm10_TSS,geneMarker,maxgapTSS,mouseAnalysis_dir="../downstream/output/mouse_analysis/"){
    
    mcols(NME_mm10_tissue)=mcols(NME_mm10_tissue)[,grepl(tissue,colnames(mcols(NME_mm10_tissue)))]
    #Promoter regions,2k

    
    if(maxgapTSS=="enhancer"){
        enhancer_bin=readRDS(bin_enhancer_rds)
        NME_mm10_tissue_gene=unique(subsetByOverlaps(NME_mm10_tissue,enhancer_bin[enhancer_bin$`Target Gene` ==geneMarker]))
    }else{
        NME_mm10_tissue_gene=unique(subsetByOverlaps(NME_mm10_tissue,mm10_TSS[geneMarker],maxgap=maxgapTSS))
    }
    if(length(NME_mm10_tissue_gene)>0){
    NME_mm10_heatmap_mt_raw=mcols(NME_mm10_tissue_gene)

    #scale rows
    scaleRow=TRUE
    if(scaleRow){
        NME_mm10_heatmap_mt=t(scale(t(NME_mm10_heatmap_mt_raw)))

    }

    fn=paste0(mouseAnalysis_dir,"cell_population/",tissue,"_",geneMarker,"_",maxgapTSS,".pdf")
    pheatmap(NME_mm10_heatmap_mt,cluster_rows =F,cluster_cols = F,
                show_colnames = T,show_rownames = F,
                filename=fn,display_numbers=round(as.matrix(NME_mm10_heatmap_mt_raw),digits=2),
                cellwidth=60,cellheight=25,annotation_legend = F,angle_col = "90",
                fontsize=10,legend = F)
            
    return(NME_mm10_heatmap_mt_raw)
    
    }else{
        cat("NME at ",maxgapTSS," of ", geneMarker," for ",tissue," does not exist\n")

    }
}
NME_mm10=readRDS(NME_matrix_file)
mm10_TSS=get_mm10_tss()
#Get tissue
tissue="heart"
maxgapTSS=2000
geneMarker="Ttn"
mouseAnalysis_dir="../downstream/output/mouse_analysis/"
#NME_mm10_tissue,tissue,mm10_TSS,geneMarker,maxgapTSS
colMeans(heatmapNME_markerGene(NME_mm10,tissue="heart",mm10_TSS=mm10_TSS,geneMarker="Ttn",maxgapTSS="enhancer"))

