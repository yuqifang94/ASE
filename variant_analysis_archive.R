#Heatmap:

density_df=data.frame(density=factor(round(NME_ASM_het_sub$density,digits = 3)),diff=factor(round(NME_ASM_het_sub$diff,digits=2)))
density_df=count(density_df,c('density','diff'))
density_df_dc=dcast(density_df,diff~density,value.var = 'freq')
density_df_dc[is.na(density_df_dc)]=0

density_heatmap=density_df_dc[,2:ncol(density_df_dc)]
rownames(density_heatmap)=density_df_dc$diff
density_heatmap=scale(density_heatmap,center=FALSE)
heatmap.2(density_heatmap,scale='none',col = bluered(100),trace = "none", density.info = "none",dendrogram='none', Rowv=FALSE, Colv=FALSE)


#NME vs density
NME_ASM_allele=NME_allele_ASM_calc[[2]]
NME_ASM_allele_het=NME_ASM_allele[NME_ASM_allele$HetCpG]

#Higher density indicate higer NME
density_df=data.frame(density=round(log10(NME_ASM_allele_het$density),digits = 2),NME=NME_ASM_allele_het$Value)
density_agg=aggregate(density_df,by=list(density_df$density),FUN=median)
plot(density_agg$Group.1,density_agg$NME,xlab='log10(density)',ylab='allelic NME')
####Calculating enrichment in het CpG vs ASM####
NME_all_gr$density_bin=round(NME_all$density*2,digits = 2)/2
NME_all_gr_enrich = data.frame(density=unique(NME_all_gr$density_bin),enrichment=NA)
for (log_den in unique(NME_all_gr$density_bin)){
  NME_all_gr_enrich$enrichment[NME_all_gr_enrich$density==log_den]=ASM_het_enrichment(NME_all_gr[NME_all_gr$density_bin==log_den])$estimate
}
NME_all_gr_enrich=NME_all_gr_enrich[NME_all_gr_enrich$density<=0.1,]
plot(NME_all_gr_enrich$density,NME_all_gr_enrich$enrichment,xlab='CpG density',ylab='Odds Ratio')



