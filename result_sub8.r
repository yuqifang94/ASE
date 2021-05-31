source('mainFunctions_sub.R')



#write all GO terms in the csv file
for (region_type in names(GO_out_all)){
  GO_in_rt=GO_out_all[[region_type]]
  for(ts in names(GO_in_rt)){
    GO_in_rt_ts=GO_in_rt[[ts]]
    GO_in_rt_ts_merge=do.call(rbind,lapply(GO_in_rt_ts,function(x) {
      x_terms=x$GO_out_cluster_all
      
      x_terms=x_terms[,.(GO.ID,Term,cluster,Annotated,Significant,Expected,FC,p_cond,genes)]
      x_terms$FDR=p.adjust(x_terms$p_cond,method="BH")
      
      return(x_terms)
      
    }))
    dir_out=paste0('../downstream/output/GO_sheets/Complete GO Terms/',region_type,'/')
    if(!dir.exists(dir_out)){dir.create(dir_out)}
    write.csv(GO_in_rt_ts_merge[p_cond<=0.1][order(FDR,decreasing=F)],paste0(dir_out,ts,'.csv'))
  }
}

GO_out_all=readRDS('../downstream/output/GO_out_all_dMML_dNME_0rm_FC_promoter.rds')
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
for(region_type in names(GO_out_all)){
  plot_GO_heatmap_all(tissue_all,GO_out_all[[region_type]],region_type=region_type,enc_type="promoter")
}



#Output the gene and terms


#dNME maxUC refers maxUC between continous time points

enhancer=import.bed('../downstream/input/mm9_ENCODE_putative_enc.bed')
mm9tomm10=import.chain('../downstream/input/mm9ToMm10.over.chain')
enhancer_mm10 = unlist(liftOver(enhancer, mm9tomm10))
enhancer_mm10=resize(enhancer_mm10,width=width(enhancer_mm10)+2000,fix="center")
hg19tomm10=import.chain('../downstream/input/hg19ToMm10.over.chain')
feDMR_in_gr_mm10=unlist(liftOver(feDMR_in_gr,hg19tomm10))
subsetByOverlaps(feDMR_in_gr_mm10,enhancer_mm10)
subsetByOverlaps(feDMR_in_gr,enhancer_mm10)



#Check proportion overlap: except NT more than 90% of the region type overlaps
tissue_out_filtered_olap=lapply(names(tissue_out_filtered),function(x){
  tissue_out_ts=tissue_out_filtered[[x]]
  tissue_out_ts$region_type_from_all=tissue_out_all[[x]]$region_type[match(tissue_out_ts$region,tissue_out_all[[x]]$region)]
  return(tissue_out_ts)
  
}

)
names(tissue_out_filtered_olap)=names(tissue_out_filtered)
saveRDS(tissue_out_filtered_olap,'../downstream/output/correlation/tissue_out_filtered_olap.rds')  


# Figures categorizing the regions ----------------------------------------
FDR_cutoff=0.2
tissue_out_filtered=readRDS('../downstream/output/correlation/tissue_out_filtered.rds')  
ts="heart"


dMML_dNME_cutoff_dt(tissue_out_filtered$heart)
ggplot()+xlim(c(-1,1))+ylim(c(-1,1))+
  geom_tile(data=tissue_out_filtered[[ts]],aes(x=dMML_cor,y=dNME_cor,stat = "density_2d_filled"),bins=200)+
  scale_fill_gradient(low = "light blue", high = "dark red")+geom_line(data=cutoff_dt,aes(x=dMML,y=dNME,group=region_type),color="black")


#change annotation cutoffs
enc_type="enhancer"
GO_out_all=readRDS(paste0('../downstream/output/GO_out_all_dMML_dNME_',enc_type,'.rds'))
pdf('../downstream/output/GO_FDR.pdf')
for(rtype in names(GO_out_all)){
  for(ts in names(GO_out_all[[rtype]])){
      for(clu in 1:10){
      clu_tissue_in=GO_out_all[[rtype]][[ts]][[clu]]$GO_out_cluster_all
      
      print(ggplot(clu_tissue_in,aes(x=Annotated,y=-log(as.numeric(classicFisher))))+geom_point(alpha=0.01)+xlim(c(0,500))+ylab("-log(p)")+
        ggtitle(paste0(rtype,':',ts,":",clu)))
      print(ggplot(clu_tissue_in,aes(x=Significant,y=-log(as.numeric(classicFisher))))+geom_point(alpha=0.01)+ylab("-log(p)")+
        geom_text(aes(x=0,y=0.5,label=paste0(round(sum(clu_tissue_in$Significant==0)/nrow(clu_tissue_in)*100,digits=1),"%")))+
        ggtitle(paste0(rtype,':',ts,":",clu)))
      }
  }
}
dev.off()


# repetitive element ------------------------------------------------------

