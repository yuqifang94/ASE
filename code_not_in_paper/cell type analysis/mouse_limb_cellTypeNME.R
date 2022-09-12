source('mainFunctions_sub.R')
limb_cellType=fread("../downstream/data/mouse_scRNA/meta.tsv")
limb_cellTypeProp=limb_cellType[,list(count=.N),by=list(stage,cell_type)]
limb_cellTypeProp=limb_cellTypeProp[,list(prop =count/sum(count),cell_type=cell_type),by=list(stage)]
limb_cellTypeEnt=limb_cellTypeProp[,list(NMEprop=-sum(prop*log2(prop))),by=list(stage)]
limb_cellTypeEnt$stage=gsub("\\.5|\\.0","",limb_cellTypeEnt$stage)
limb_cellTypeEnt=limb_cellTypeEnt[,list(NMEprop=mean(NMEprop)),by=list(stage)]
pdf("../downstream/output/mouse_analysis/cell_population/cellType_ent.pdf")
    ggplot(limb_cellTypeEnt,aes(x=as.numeric(stage),y=NMEprop))+geom_point()+theme_classic()
dev.off()
cellTypeCor<-function(datIn,statType,limb_cellTypeEnt,permute=F,seed=1){
    
    mcols(datIn)=mcols(datIn)[,grepl("limb",colnames(mcols(datIn)))]
    datIn=convert_GR(datIn,direction="DT")
    datIn=datIn[rowSums(is.na(datIn))==0]
    datIn=melt.data.table(datIn,id.vars= c("region"),variable.name = "stage", value.name = "datIn")
    datIn$stage=gsub("limb.E|.all","",datIn$stage)
    datIn$NMEprop=limb_cellTypeEnt[match(gsub("\\.5","",datIn$stage),stage)]$NMEprop
    regionFt=datIn[,list(NStage=.N),by=list(region)]
    totalStage=length(unique(datIn$stage))
    if(permute){
        set.seed(seed)
        datIn_NME_cor=datIn[region%in%regionFt[NStage==totalStage]$region,list(cellTypeNMECor=cor(scale(datIn)[,1],sample(scale(NMEprop)[,1]))),by=list(region)]
    }else{
        datIn_NME_cor=datIn[region%in%regionFt[NStage==totalStage]$region,list(cellTypeNMECor=cor(scale(datIn)[,1],scale(NMEprop)[,1])),by=list(region)]
    }
    datIn_NME_cor=datIn_NME_cor[!is.na(cellTypeNMECor)]
    saveRDS(datIn_NME_cor,paste0("../downstream/output/mouse_analysis/cell_population/",statType,"_cellTypeNME_cor_permute_",permute,".rds"))
    return(datIn)
}
NME_cellTypeCor=cellTypeCor(readRDS(NME_matrix_file),"NME",limb_cellTypeEnt)
MML_cellTypeCor=cellTypeCor(readRDS(MML_matrix_file),"MML",limb_cellTypeEnt)
NME_cellTypeCor=cellTypeCor(readRDS(NME_matrix_file),"NME",limb_cellTypeEnt,permute=T)
MML_cellTypeCor=cellTypeCor(readRDS(MML_matrix_file),"MML",limb_cellTypeEnt,permute=T)
#Analyze correlation in selected region
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
tissue_out_filtered=tissue_out_filtered$limb
plot_cor_dist<-function(statType,tissue_out_filtered){
    dat_cor_cellTypeNME=readRDS(paste0("../downstream/output/mouse_analysis/cell_population/",statType,"_cellTypeNME_cor_permute_FALSE.rds"))
    tissue_out_filtered$cellTypeNMECor=dat_cor_cellTypeNME[match(tissue_out_filtered$region,region)]$cellTypeNMECor
    tissue_out_filtered_plot=tissue_out_filtered[,list(region,region_type,cellTypeNMECor)]
    dat_cor_cellTypeNME_perm=readRDS(paste0("../downstream/output/mouse_analysis/cell_population/",statType,"_cellTypeNME_cor_permute_TRUE.rds"))[region%in%tissue_out_filtered$region]
    dat_cor_cellTypeNME_perm$region_type="random control"
    tissue_out_filtered_plot=rbind(tissue_out_filtered_plot,dat_cor_cellTypeNME_perm)
    pdf(paste0("../downstream/output/mouse_analysis/cell_population/cellType_ent_cor_",statType,".pdf"))
        print(ggplot(tissue_out_filtered_plot,aes(x=cellTypeNMECor,color=region_type))+geom_density()+xlab(statType))
    dev.off()
    #Two sided
    tissue_out_filtered_plot$empirP=1-ecdf(abs(tissue_out_filtered_plot[region_type=="random control"]$cellTypeNMECor))(abs(tissue_out_filtered_plot$cellTypeNMECor))
    tissue_out_filtered_plot_P=tissue_out_filtered_plot[region_type!="random control"]
    tissue_out_filtered_plot_P$FDR=p.adjust(tissue_out_filtered_plot_P$empirP,method='BH')#18214
    return(list(tissue_out_filtered_plot_P,tissue_out_filtered_plot))
}
NME_cor=plot_cor_dist("NME",tissue_out_filtered)#18214/72846 FDR<=0.1

#       region_type   medianCor
# 1:           Both  0.63114701
# 2:        Neither  0.63280089
# 3:       MML only  0.58768651
# 4:       NME only  0.55495401
# 5: random control -0.02416779
MML_cor=plot_cor_dist("MML",tissue_out_filtered)#28923/72846 
#       region_type    medianCor
# 1:           Both -0.734054611
# 2:        Neither -0.666352591
# 3:       MML only -0.903180416
# 4:       NME only -0.669475917
# 5: random control -0.006016005