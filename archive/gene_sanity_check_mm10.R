source('mainFunctions_sub.R')
csv_test=fread('../downstream/input/mm10_cluster_all/forebrain.csv')
gtf <- fread('../downstream/input/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
genes <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
names(genes) <- gn
pro <- promoters(genes,upstream=2000,downstream=1000)
tss <- promoters(genes,upstream=0,downstream=1)

gene_check=nearest(convert_GR(csv_test$region),tss)
csv_test$gene_check=names(tss)[gene_check]
csv_test$dist_check=mcols(distanceToNearest(convert_GR(csv_test$region),tss))$distance
GO_out=list()
for(i in 1:10){
gl=csv_test[cluster==i&chromHMM_enhancer ==T]$gene_check
back=csv_test[chromHMM_enhancer ==T]$gene_check
GO_out[[i]]=GO_run(unique(gl),unique(back))
GO_out[[i]]$clu=i
}
GO_out=readRDS('../downstream/output/GO_out_limb.rds')
all_terms=unique(do.call(rbind,lapply(GO_out,function(x) x[FC>=1.5][1:5,1:8]))$Term)
GO_out_all=do.call(rbind,lapply(GO_out,function(x) x[Term %in%all_terms][,c(1:8,11)]))
GO_out_all$cluster=GO_out_all$clu
GO_in_main= dcast_matrix(GO_out_all,"FC")
GO_in_FDR= dcast_matrix(GO_out_all,"FDR")
GO_in_FDR[GO_in_FDR<=0.1]="*"
GO_in_FDR[GO_in_FDR>0.1]=""
GO_in_FDR=GO_in_FDR[rownames(GO_in_main),]
col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
c2 <- brewer.pal(10,'Set3')
names(c2) <- 1:10
breaksList = seq(-1, 1, by = 0.01)
colann= data.frame(cluster=as.character(1:10))
#pdf(paste0('../downstream/output/graphs/Figure6/GO_', tissue,'_',GO_anno,'_FC.pdf'),width=25,height=14)
pheatmap(scalematrix(GO_in_main),cluster_rows =F,cluster_cols = F,
         show_colnames = T,show_rownames = T,display_numbers=GO_in_FDR,border_color = NA,
         color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),
         filename=paste0('../downstream/output/graphs/Figure6/GO_chromHMM_another_anno_limb.pdf'),
         cellwidth=60,cellheight=25,annotation_colors = list(cluster=c2),annotation_col = colann, annotation_legend = F,
         fontsize=30,legend = F)
GO_run<-function(gl,back,ptcount=0){
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                                  annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")})
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  #sigres <- sigres[as.numeric(sigres$FDR) <= 0.1,]
  fc <- ((sigres[,"Significant"]+ptcount)/(sum(GOdata@allScores[GOdata@feasible]==1)+ptcount))/((sigres[,"Annotated"]+ptcount)/(sum(GOdata@feasible)+ptcount))
  sigres <- data.frame(sigres,FC=fc)
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
  sigres=as.data.table(sigres)
  siggene_forID=lapply(sigres$GO.ID,function(x,GOdata){
    gene=sigGenes(GOdata)[sigGenes(GOdata)%in%unlist(genesInTerm(GOdata, x))]
    gl_dt=data.table(rank=1:length(gl),gene=gl)
    mt=match(gl_dt$gene,gene)
    mt=mt[!is.na(mt)]
    #highest_rank=min(gl_dt$rank[gl_dt$gene %in% gene])
    highest_rank=NA
    #motif_in=paste(motif[motif%in%gene],collapse = ";")
    #if(length(motif_in)==0){motif_in=NA}
    return(list(paste(gene[mt],collapse =";"),highest_rank))
    
    
  },GOdata=GOdata)
  siggene=unlist(lapply(siggene_forID,function(x) x[[1]]))
  max_rank=unlist(lapply(siggene_forID,function(x) x[[2]]))
  #motif_in=unlist(lapply(siggene_forID,function(x) x[[3]]))
  if(nrow(sigres)>0){ 
    sigres$genes=siggene
    sigres$higest_ranks=max_rank
    # sigres$motif=motif_in
  }
  return(sigres)
}
