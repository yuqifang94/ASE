setwd(".../../../../")
source('mainFunctions_sub.R')
library(tidyverse)
experimentList = read.table('../downstream/input/mouse_analysis/motif_analysis/chipatlas/experimentList.tab',sep="\t",fill=T)
experimentList = experimentList%>% unite(sampleMeta, 10:24,sep=';',remove=F)
experimentList[,paste0("V",10:24)]=NULL
experimentList$sampleMeta=gsub(";;+","",experimentList$sampleMeta)
colnames(experimentList)[1:9]=c("ExpID","GenomeAssembly","AntigenClass","Antigen","CelltypeClass","CellType","CellTypeDescription","ProcessLog","Title")
experimentList=as.data.table(experimentList)
experimentList = experimentList[AntigenClass =="TFs and others"&GenomeAssembly=="mm10"]
experimentList_embryos=experimentList[grepl("Embryonic|embryos",CellType,ignore.case = T)]
experimentList_embryos_table=experimentList_embryos[,list(numberTFs=length(unique(Antigen))),by=CellType]
write.csv(experimentList_embryos_table[order(numberTFs,decreasing=T)],"../downstream/output/mouse_analysis/ChIP_analysis/TFcount_cellType.csv",row.names = F)
experimentList_limb=experimentList_embryos[CellType=="Embryonic limb"]
experimentList_brain=experimentList_embryos[CellType=="Embryonic brains"]
experimentList_heart=experimentList_embryos[CellType=="Embryonic heart"]
intersect(unique(experimentList_limb$Antigen),unique(experimentList_brain$Antigen))
experimentList_embryos_highNMETF=experimentList_embryos[experimentList_embryos$Antigen%in%c("Foxg1","Cebpa","Jund","Meis1","Tal1","Tbx4","Rara","Sox10","Pbx1","Klf4", "Pknox1")]
lapply(unique(experimentList_embryos_highNMETF$Antigen),function(x) table(experimentList_embryos_highNMETF[Antigen==x]$CellType))

NMEOnly_human=fread("../downstream/output/human_analysis/motif_analysis/motif_prefer_high_NME_only_OMIM.csv")
MMLOnly_human=fread("../downstream/output/human_analysis/motif_analysis/motif_prefer_low_MML_only_OMIM.csv")
experimentList_embryos_NMEOnly=experimentList_embryos[unlist(lapply(Antigen,function(x) any(grepl(x,NMEOnly_human$TF,ignore.case = T))))]
experimentList_embryos_table_NMEOnly=experimentList_embryos_NMEOnly[,list(numberTFs=length(unique(Antigen))),by=CellType]
write.csv(experimentList_embryos_table_NMEOnly[order(numberTFs,decreasing=T)],"../downstream/output/mouse_analysis/ChIP_analysis/TFcount_cellType_NMEOnly.csv",row.names = F)

#TF and number of TF that are shared across different tissue
experimentList_embryos_NMEOnly_TF=experimentList_embryos_NMEOnly[,list(NCell=length(unique(CellType)),CellTypes=paste(unique(CellType),collapse=",")),by=Antigen]
write.csv(experimentList_embryos_NMEOnly_TF[order(NCell,decreasing=T)],"../downstream/output/mouse_analysis/ChIP_analysis/CellTypeCount_TF_NMEOnly.csv",row.names = F)
