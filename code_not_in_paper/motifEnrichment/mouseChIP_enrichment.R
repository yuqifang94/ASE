source('mainFunctions_sub.R')
ChIPOutDirMouse="../downstream/output/mouse_analysis/ChIP_analysis/"
preProcessChIP<-function(bedFile){
    factorIn=fread(bedFile,skip = 1,sep='\t')
    colnames(factorIn)=c('seqnames','start','end','metadata','log10qval','not_used','start2','end2','color')
    factorIn=factorIn[,list(seqnames,start,end,metadata,log10qval)]
    factorIn$metadata=gsub("%20","_",factorIn$metadata)
    factorIn$metadata=gsub("%3B","",factorIn$metadata)#87960
    factorIn$TF = unlist(lapply(strsplit(factorIn$metadata,";"), 
        function(x){
            TF = x[grepl("^Name=",x)]            
            TF = gsub("Name=","",TF)
            TF = gsub("\\_\\(\\@\\_.*","",TF)
            return(TF)
        }
    ))
    return(makeGRangesFromDataFrame(factorIn,keep.extra.columns=T))
}
#Generate Heatmap to check average NME and MML stratifying by tissue and developmental stage (sample)

heatmapGen<-function(datMatrix,ChIPIn,fn,minoverlap=200, control = T, selectedColumn =NULL,olapThreshold=100){
    if (length(selectedColumn)>0){
        mcols(datMatrix)=mcols(datMatrix)[,selectedColumn]

    }
    TFMatrix = lapply(unique(ChIPIn$TF),function(x){
                TFolap = subsetByOverlaps(datMatrix,ChIPIn[ChIPIn$TF == x],minoverlap=minoverlap)
                #cat("The overlap size for",x,"is",length(TFolap),"\n")
                return(convert_GR(TFolap,direction="matrix"))

    })
    names(TFMatrix)=unique(ChIPIn$TF)
    TFMatrix=TFMatrix[unlist(lapply(TFMatrix,function(x) nrow(x)>=olapThreshold))]
    TFMatrixMean=do.call(rbind,
                            lapply(names(TFMatrix),function(x){
                            matrixOut= t(as.matrix(colMeans(TFMatrix[[x]],na.rm=T)))
                            rownames(matrixOut) = x
                            return(matrixOut)

    }))
    TFMatrixMean=TFMatrixMean[order(rowMeans(TFMatrixMean),decreasing=T),]
    #Biotin can be used as reference
    controlList=c("Epitope_tags","Biotin","GFP")
    if (control){
    
    TFMatrixMean=rbind(TFMatrixMean[!rownames(TFMatrixMean)%in% controlList,],
                            TFMatrixMean[rownames(TFMatrixMean)%in%controlList,])
    }else{
         TFMatrixMean=rbind(TFMatrixMean[!rownames(TFMatrixMean)%in% controlList,])

    }
    colnames(TFMatrixMean)=gsub("-all","",colnames(TFMatrixMean))
    colann= data.frame(stage=gsub('.*-','',colnames(TFMatrixMean)),tissue=gsub('-.*','',colnames(TFMatrixMean)))
    rownames(colann)=colnames(TFMatrixMean)
    tissue_col=mouse_color()
    stage_col <- brewer.pal(length(unique(colann$stage)),'Set3')
    #Heatmap 
    pheatmap(TFMatrixMean,cluster_rows =F,cluster_cols = F,
                show_colnames = T,show_rownames = T,,border_color = NA,
                color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),
                breaks=seq(0,1,0.01),
                filename=fn,
                cellwidth=25,cellheight=10,annotation_legend = F,angle_col = "90",
                fontsize=10,legend = T,annotation_col = colann,
                annotation_colors = list(tissue=tissue_col,cluster=stage_col))
    return(list(TFMatrixMean=TFMatrixMean,TFMatrix=TFMatrix))
}
# Seleccted from Chip-atlas ---------------------------------------------
#embyro_heart: mm10->Embyro->Embyronic heart
#adult_heart: mm10->Cardiovascular->all
#embyro_all: mm10->Embyro->all
#Adult Neuro:  mm10->Neural->all
#embyro Heart: 20220120 download
mm10EmbyroChIP=preProcessChIP(paste0(chip_atlas_dir,'mm10_embyro_all_TF.bed'))
saveRDS(mm10EmbyroChIP,paste0(ChIPOutDirMouse,"mouseEmbyroChIP.rds"))
mm10EmbyroHeartChIP=preProcessChIP(paste0(chip_atlas_dir,'mm10_embyro_heart_allTF.bed'))
saveRDS(mm10EmbyroChIP,paste0(ChIPOutDirMouse,"mouseEmbyroHeartChIP.rds"))
mm10AdultHeartChIP=preProcessChIP(paste0(chip_atlas_dir,'mm10_heart_all_TF.bed'))
saveRDS(mm10EmbyroChIP,paste0(ChIPOutDirMouse,"mouseAdultHeartChIP.rds"))
mm10EmbyroBrainChIP=preProcessChIP(paste0(chip_atlas_dir,'mm10_embyro_brain_all_TF.bed'))
saveRDS(mm10EmbyroChIP,paste0(ChIPOutDirMouse,"mouseEmbyroBrainChIP.rds"))

mm10NME = readRDS(NME_matrix_file)#mm10
mm10MML = readRDS(MML_matrix_file)#mm10
#All embyro
NME_allEmbyro_TFMatrix=heatmapGen(mm10NME,mm10EmbyroChIP,paste0(ChIPOutDirMouse,"NME_mouse_heatmap_embyroTF.pdf"),control=F)
MML_allEmbyro_TFMatrix=heatmapGen(mm10MML,mm10EmbyroChIP,paste0(ChIPOutDirMouse,"MML_mouse_heatmap_embyroTF.pdf"),control=F)
#QC for histgram
pdf(paste0(ChIPOutDirMouse,'mm10_embyro_all_TF_QC_hist.pdf'))
    olapRegion=unlist(lapply(NME_allEmbyro_TFMatrix$TFMatrix,nrow))
    hist(ifelse(olapRegion>1000,1100,olapRegion),breaks=seq(0,1100,100),
            xlab="Number of overlapped regions", main="",xlim=c(0,1190))
            axis(1, at=seq(0,1100,100)-.5, labels=c(seq(0,1000,100),">1000"))
dev.off()
#Heart
NME_HeartEmbyro_HeartTFMatrix=heatmapGen(mm10NME,mm10EmbyroHeartChIP,paste0(ChIPOutDirMouse,"NME_mouseHeart_heatmap_embyroHeartTF.pdf"),control=F,selectedColumn=grepl("heart",colnames(mcols(mm10NME))))
MML_heartEmbyro_HeartTFMatrix=heatmapGen(mm10MML,mm10EmbyroHeartChIP,paste0(ChIPOutDirMouse,"MML_mouseHeart_heatmap_embyroHeartTF.pdf"),control=F,selectedColumn=grepl("heart",colnames(mcols(mm10NME))))
NME_HeartEmbyro_TFMatrix=heatmapGen(mm10NME,mm10EmbyroHeartChIP,paste0(ChIPOutDirMouse,"NME_mouse_heatmap_embyroHeartTF.pdf"),control=F)
MML_heartEmbyro_TFMatrix=heatmapGen(mm10MML,mm10EmbyroHeartChIP,paste0(ChIPOutDirMouse,"MML_mouse_heatmap_embyroHeartTF.pdf"),control=F)

#Brain
NME_BrainEmbyro_BrainTFMatrix=heatmapGen(mm10NME,mm10EmbyroBrainChIP,paste0(ChIPOutDirMouse,"NME_mouseBrain_heatmap_embyroBrainTF.pdf"),control=F,selectedColumn=grepl("brain",colnames(mcols(mm10NME))))
MML_BrainEmbyro_BrainTFMatrix=heatmapGen(mm10MML,mm10EmbyroBrainChIP,paste0(ChIPOutDirMouse,"MML_mouseBrain_heatmap_embyroBrainTF.pdf"),control=F,selectedColumn=grepl("brain",colnames(mcols(mm10NME))))
NME_BrainEmbyro_TFMatrix=heatmapGen(mm10NME,mm10EmbyroBrainChIP,paste0(ChIPOutDirMouse,"NME_mouse_heatmap_embyroBrainTF.pdf"),control=F)
MML_BrainEmbyro_TFMatrix=heatmapGen(mm10MML,mm10EmbyroBrainChIP,paste0(ChIPOutDirMouse,"MML_mouse_heatmap_embyroBrainTF.pdf"),control=F)
#Check NME on NME only
NMEOnly_human=fread("../downstream/output/human_analysis/motif_analysis/motif_prefer_high_NME_only_OMIM.csv")
MMLOnly_human=fread("../downstream/output/human_analysis/motif_analysis/motif_prefer_low_MML_only_OMIM.csv")
humanHighNMEOnly=unlist(lapply(rownames(NME_allEmbyro_TFMatrix$TFMatrixMean),function(x) any(grepl(x,NMEOnly_human$TF,ignore.case=T))))
humanLowMMLOnly=unlist(lapply(rownames(NME_allEmbyro_TFMatrix$TFMatrixMean),function(x) any(grepl(x,MMLOnly_human$TF,ignore.case=T))))

humanHighNMEOnly=unlist(lapply(rownames(NME_BrainEmbyro_BrainTFMatrix$TFMatrixMean),function(x) any(grepl(x,NMEOnly_human$TF,ignore.case=T))))
humanLowMMLOnly=unlist(lapply(rownames(NME_BrainEmbyro_BrainTFMatrix$TFMatrixMean),function(x) any(grepl(x,MMLOnly_human$TF,ignore.case=T))))

