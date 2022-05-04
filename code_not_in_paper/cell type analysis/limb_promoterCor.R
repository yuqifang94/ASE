source('mainFunctions_sub.R')
#Extract data, no SAVER
limb_data_raw=readRDS('../downstream/data/mouseLimb/limb_10x.rds')
#Mesenchymal
stage=unique(gsub('.*-|:.*','',colnames(limb_data_raw)))
cellTypeStage = lapply(stage,function(x) limb_data_raw[,grepl(x,colnames(limb_data_raw))])
names(cellTypeStage)=stage
cellTypeStage_Mesenchymal=lapply(cellTypeStage,function(x) x[,grepl("Mesenchymal",colnames(x))])
cellTypeStage_Mesenchymal_MAV=cellTypeMAVcalc(cellTypeStage_Mesenchymal,'*')
cellTypeStage_all_MAV=cellTypeMAVcalc(cellTypeStage,'*')
#Extract data, SAVER
mouseSAVER_dir='../downstream/data/mouseLimb/10xSAVER_run2/'
mouseSAVER_dir_Jason='../downstream/data/mouseLimb/10x/'
limb_RNA_data=list()
#,pattern="10xLimbSaver"
for (fn in dir(mouseSAVER_dir)){
    stage=gsub('10xLimbSaver_|.rds','',fn)
    SAVERIn=readRDS(paste0(mouseSAVER_dir,fn))
    #JasonSaver= readRDS(paste0(mouseSAVER_dir_Jason,gsub("E","",stage),'.rds'))
    #saverCellType=readRDS(paste0(mouseSAVER_dir,fn))
    #colnames(JasonSaver)=colnames(saverCellType)
    #limb_RNA_data[[stage]]=JasonSaver
    limb_RNA_data[[stage]]=SAVERIn
}

limb_RNA_data_MAV=cellTypeMAVcalc(limb_RNA_data,'*')
limb_RNA_data_mes=lapply(limb_RNA_data,function(x) x[,grepl("Mesenchymal",colnames(x))])
limb_RNA_data_MAV_mes=cellTypeMAVcalc(limb_RNA_data_mes,'*')


# #Comparing 10x from Jason and mine MAV
# mouse10xDir='../downstream/input/mouse_analysis/mouse_10x/'
# Jason_105=readRDS(paste0(mouse10xDir,'10.5.rds'))
# #In MARCC
# Jason_105=readRDS("../downstream/data/mouseLimb/10x/10.5.rds")
# mouse105=readRDS('../downstream/data/mouseLimb/10xSAVER_run/10xLimbSaver_E10.5.rds')
#diff=Jason_105-mouse105


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


mm10Hist=plotHistCor(NME_mm10_limb,cellTypeStage_Mesenchymal_MAV[stage!='E13.5'],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_NoSAVER_Mes.pdf')
    grid.arrange(grobs =mm10Hist, nrow = 3)
dev.off()
mm10Hist=plotHistCor(NME_mm10_limb,cellTypeStage_all_MAV[stage!='E13.5'],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_NoSAVER_all.pdf')
    grid.arrange(grobs =mm10Hist, nrow = 3)
dev.off()
mm10Hist=plotHistCor(NME_mm10_limb,limb_RNA_data_MAV_mes[stage!="E13.5"],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_mes_SAVER.pdf')
    grid.arrange(grobs =mm10Hist, nrow = 3)
dev.off()
#Need check remove E13.5, E13.0 or all
mm10HistAll=plotHistCor(NME_mm10_limb,limb_RNA_data_MAV,tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_SAVER_allStage.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()
mm10HistAll=plotHistCor(NME_mm10_limb,limb_RNA_data_MAV[stage!="E13.5"],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_SAVER_no135.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()
mm10HistAll=plotHistCor(NME_mm10_limb,limb_RNA_data_MAV[stage!="E13.0"],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_SAVER_no13.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()
#Mes
mm10HistAll=plotHistCor(NME_mm10_limb,limb_RNA_data_MAV_mes,tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_SAVER_allStage_mes.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()
mm10HistAll=plotHistCor(NME_mm10_limb,limb_RNA_data_MAV_mes[stage!="E13.5"],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_SAVER_no135_mes.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()
mm10HistAll=plotHistCor(NME_mm10_limb,limb_RNA_data_MAV_mes[stage!="E13.0"],tss_mm10)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_promoter_SAVER_no13_mes.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()
#Check SAVER enhancer again
enhancer=readRDS(bin_enhancer_rds)
NME_limb_enhancerCor_mesenchymal=enhancerCorCalc(NME_mm10_limb,enhancer,limb_RNA_data_MAV_mes)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_Mesenchymal_SAVER.pdf')
print(ggplot(NME_limb_enhancerCor_mesenchymal,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_mesenchymal$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()

NME_limb_enhancerCor_all=enhancerCorCalc(NME_mm10_limb,enhancer,limb_RNA_data_MAV)
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_All_SAVER.pdf')
print(ggplot(NME_limb_enhancerCor_all,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_all$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()

#Check SAVER enhancer again remove E13.5
enhancer=readRDS(bin_enhancer_rds)
NME_limb_enhancerCor_mesenchymal=enhancerCorCalc(NME_mm10_limb,enhancer,limb_RNA_data_MAV_mes[stage!='E13.5'])
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_Mesenchymal_SAVER_no135.pdf')
print(ggplot(NME_limb_enhancerCor_mesenchymal,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_mesenchymal$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()

NME_limb_enhancerCor_all=enhancerCorCalc(NME_mm10_limb,enhancer,limb_RNA_data_MAV[stage!='E13.5'])
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_All_SAVER_no135.pdf')
print(ggplot(NME_limb_enhancerCor_all,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_all$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()

#Check SAVER enhancer again remove E13.0
enhancer=readRDS(bin_enhancer_rds)
NME_limb_enhancerCor_mesenchymal=enhancerCorCalc(NME_mm10_limb,enhancer,limb_RNA_data_MAV_mes[stage!='E13'])
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_Mesenchymal_SAVER_no13.pdf')
print(ggplot(NME_limb_enhancerCor_mesenchymal,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_mesenchymal$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()

NME_limb_enhancerCor_all=enhancerCorCalc(NME_mm10_limb,enhancer,limb_RNA_data_MAV[stage!='E13'])
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_All_SAVER_no13.pdf')
print(ggplot(NME_limb_enhancerCor_all,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_all$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()

#After SAVER, do the bs, nonBS, SAVER, nonSAVER histogram analysis
#do dNME vs dMAV at enhancer and promoter
dNME_mm10_limb=diffCalc(NME_mm10_limb,"region",valueName="dNME",valueVar="NME")
limb_RNA_data_dMAV=diffCalc(limb_RNA_data_MAV[stage!="E13.5",list(gene,hypervar_logvar,stage)],"gene",valueName="dMAV",valueVar="hypervar_logvar")
limb_RNA_data_dMAV_mes=diffCalc(limb_RNA_data_MAV_mes[stage!="E13.5",list(gene,hypervar_logvar,stage)],"gene",valueName="dMAV",valueVar="hypervar_logvar")

mm10HistAll=plotHistCor(dNME_mm10_limb,limb_RNA_data_dMAV,tss_mm10,NMEVar='dNME',MAVVar='dMAV')
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/dNME_dMAV_cor_promoter_SAVER_allStage.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()

mm10HistAll=plotHistCor(dNME_mm10_limb,limb_RNA_data_dMAV_mes,tss_mm10,NMEVar='dNME',MAVVar='dMAV')
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/dNME_dMAV_cor_promoter_SAVER_mes_allStage.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()

#do dNME vs dMAV at enhancer and promoter no SAVER
dNME_mm10_limb=diffCalc(NME_mm10_limb,"region",valueName="dNME",valueVar="NME")
limb_RNA_data_dMAV=diffCalc(cellTypeStage_all_MAV[stage!="E13.5",list(gene,hypervar_logvar,stage)],"gene",valueName="dMAV",valueVar="hypervar_logvar")
limb_RNA_data_dMAV_mes=diffCalc(cellTypeStage_Mesenchymal_MAV[stage!="E13.5",list(gene,hypervar_logvar,stage)],"gene",valueName="dMAV",valueVar="hypervar_logvar")

mm10HistAll=plotHistCor(dNME_mm10_limb,limb_RNA_data_dMAV,tss_mm10,NMEVar='dNME',MAVVar='dMAV')
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/dNME_dMAV_cor_promoter_noSAVER_allStage.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()

mm10HistAll=plotHistCor(dNME_mm10_limb,limb_RNA_data_dMAV_mes,tss_mm10,NMEVar='dNME',MAVVar='dMAV')
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/dNME_dMAV_cor_promoter_noSAVER_mes_allStage.pdf')
    grid.arrange(grobs =mm10HistAll, nrow = 3)
dev.off()

enhancer=readRDS(bin_enhancer_rds)
NME_limb_enhancerCor_all=enhancerCorCalc(NME_mm10_limb,enhancer,cellTypeStage_all_MAV[stage!='E13.5'])
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/NME_MAV_cor_All_noSAVER_no135.pdf')
print(ggplot(NME_limb_enhancerCor_all,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_all$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()
dNME_limb_enhancerCor_all=enhancerCorCalc(dNME_mm10_limb,enhancer,limb_RNA_data_dMAV[stage!='E13.5'],NMEVar='dNME',MAVVar='dMAV')
pdf('../downstream/output/mouse_analysis/limb_cell_heterogenity/dNME_dMAV_cor_All_noSAVER_no135.pdf')
print(ggplot(dNME_limb_enhancerCor_all,aes(x=NME_MAV_cor))+geom_density())
hist(NME_limb_enhancerCor_all$NME_MAV_cor,xlab="correlation",main="NME MAV correlation at enhancer")
dev.off()
#Functions
plotHistCor<-function(NMEIn,MAV,tss,NMEVar="NME",MAVVar="hypervar_logvar"){
    histAll=list()
    for (maxGap in c(150,1000,2500,5000,10000,20000)){
        NMEInAnnotated=promAnno(NMEIn,tss,maxGap=maxGap,geneList=unique(MAV$gene),NMEVar=NMEVar)
        
        NMEInAnnotatedCor=corCalc(NMEInAnnotated,MAV,permute=F,seed=0,NMEVar="NME",MAVVar=MAVVar,meanCalc=F)
        histAll[[as.character(maxGap)]]=ggplot(NMEInAnnotatedCor, aes(x=NME_MAV_cor))+geom_histogram(color="black", fill="white")+
                xlab("correlation")+
                ggtitle(paste0("correlation at: ",maxGap," from TSS \n"))

    
    }
    return(histAll)

}
encAnno<-function(NME,enhancer){
      #Get the annotated genes
    olap=findOverlaps(convert_GR(NME$region,direction="GR"),enhancer)
    NME_annotated=NME[queryHits(olap)]
    NME_annotated$gene=enhancer$`Target Gene`[subjectHits(olap)]
    return(NME_annotated)
}
promAnno<-function(NME,promoter,maxGap=2000,geneList,NMEVar="NME"){
      #Get the annotated genes
     NME$NME = NME[[NMEVar]]
    promoter <- promoters(promoter,upstream=maxGap*2,downstream=maxGap)
    #promoter <- promoter[geneList]

    olap=findOverlaps(convert_GR(NME$region,direction="GR"),promoter)
    NME_annotated=NME[queryHits(olap)]
    NME_annotated$subjH= subjectHits(olap)
    NME_annotated = NME_annotated[,list(NME=mean(NME)),by=list(stage,subjH)]
    NME_annotated$gene=names(promoter)[NME_annotated$subjH]
    # nearestDist=as.data.table(distanceToNearest(convert_GR(NME$region,direction="GR"),promoter))
    # NME$dist=nearestDist$distance
    # NME$gene=names(promoter)[nearestDist$subjectHits]
    # return(NME[abs(dist)<=maxGap])
    print(NME_annotated)
    return(NME_annotated)
}
enhancerCorCalc<-function(NME,enhancer,MAV,permute=F,seed=0,NMEVar="NME",MAVVar="hypervar_logvar"){
    NME_annotated=encAnno(NME,enhancer)
    NME_annotated_cor=corCalc(NME_annotated,MAV,permute,seed,NMEVar,MAVVar)
    return(NME_annotated_cor)
}

corCalc<-function(NME_annotated,MAV,permute=F,seed=0,NMEVar="NME",MAVVar="hypervar_logvar",meanCalc=T){
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
   
    #NME_annotated_MAV_5stage_sdMAV=NME_annotated_MAV_5stage[,list(sdMAV=sd(MAV)),by=list(region)]
    if(meanCalc){
        NME_annotated_count=NME_annotated_MAV[,list(nStage=length(unique(stage))),by=list(region)]
        NME_annotated_MAV_5stage=NME_annotated_MAV[region%in%NME_annotated_count[nStage>=5]$region]
        NME_annotated_MAV_5stage=NME_annotated_MAV_5stage[,list(NME=mean(NME),MAV=mean(MAV)),by=list(gene,stage_relax)]
    }else{
        NME_annotated_count=NME_annotated_MAV[,list(nStage=length(unique(stage))),by=list(gene)]
        NME_annotated_MAV_5stage=NME_annotated_MAV[gene%in%NME_annotated_count[nStage>=5]$gene]
        NME_annotated_MAV_5stage=NME_annotated_MAV_5stage[,list(NME=mean(NME),MAV=mean(MAV)),by=list(gene,stage_relax)]
    
    }
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
#MAV calc
MAV_calc<-function(expr){
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
  expr <- expr[!grepl('^Rpl|^Rps',rownames(expr)),]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
   #hypervar_logvar<- resid(lm(logvar~bs(lm))) #This is bs
  res <- data.table(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar,gene=rownames(expr))
  return(res)
  #
 
  
}
#Go for Mesenchymal first
cellTypeMAVcalc<-function(limbDat,selectedCellType){
    cellTypeStageMAV=limbDat[grepl(selectedCellType,names(limbDat))]
    cellTypeStage_MAV=lapply(cellTypeStageMAV,MAV_calc)
    cellTypeStage_MAV=lapply(names(cellTypeStage_MAV),function(x) {
                                        dtIn= cellTypeStage_MAV[[x]]
                                        dtIn$cellType_stage=x
                                        return(dtIn)
                                        
            })
    cellTypeStage_MAV=do.call(rbind,cellTypeStage_MAV)
    cellTypeStage_MAV$stage=gsub('.*-','',cellTypeStage_MAV$cellType)
    return(cellTypeStage_MAV)
}
diffCalc<-function(dtIn,rnCol,valueName,valueVar){
    stage_diff=expand.grid(unique(dtIn$stage),unique(dtIn$stage))
    stage_diff=stage_diff[stage_diff$Var1!=stage_diff$Var2,]
    dtIn$id=dtIn[[rnCol]]
    dtIn=dcast.data.table(dtIn,id~stage,value.var=valueVar)
    rn=dtIn[,which(colnames(dtIn)=="id"),with=F][[1]]
    dtIn=as.matrix(dtIn)
    dtIn=dtIn[,which(colnames(dtIn)!="id")]
    rownames(dtIn)=rn
    
    colNamediff=c()
    dtOut=data.frame(matrix(,nrow=nrow(dtIn),ncol=0))
    for(i in 1:nrow(stage_diff)){
        stageDiff=paste0(stage_diff[i,1],'-',stage_diff[i,2])
        dtOut[[stageDiff]]=as.numeric(dtIn[,stage_diff[i,1]])-as.numeric(dtIn[,stage_diff[i,2]])


    }
    dtOut[[rnCol]]=rownames(dtIn)
    dtOut=as.data.table(dtOut)
    dtOut=melt.data.table(dtOut,id.vars=rnCol, variable.name = "stage",value.name = valueName)
    return(dtOut)

}