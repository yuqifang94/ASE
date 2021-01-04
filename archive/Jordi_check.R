#GO annotation
GO_run<-function(gl,back,motif){
  #GO object generation
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                                  annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")})
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  #Filter out the GO terms
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  #maybe remove Expected filter, ask Jason if it's a good idea add Expected filter, look for standard way, Expected is more strict, and redndent, keep one 
  sigres <- sigres[(as.numeric(sigres$FDR) <= 0.1)&(sigres$Expected>=5),] 
  #Calculate FC based
  fc <- (sigres[,"Significant"]/sum(GOdata@allScores[GOdata@feasible]==1))/(sigres[,"Annotated"]/sum(GOdata@feasible))
  sigres <- data.frame(sigres,FC=fc)
  #Order terms
  sigres <- sigres[order(sigres$FDR,-sigres$FC),]
  sigres=as.data.table(sigres)
  #Find significant genes with rank, GenTable function
  siggene_forID=lapply(sigres$GO.ID,function(x,GOdata){
    gene=sigGenes(GOdata)[sigGenes(GOdata)%in%unlist(genesInTerm(GOdata, x))]
    gl_dt=data.table(rank=1:length(gl),gene=gl)
    mt=match(gl_dt$gene,gene)
    mt=mt[!is.na(mt)]
    highest_rank=min(gl_dt$rank[gl_dt$gene %in% gene])
    
    motif_in=motif[motif%in%gene]
    if(length(motif_in)==0){motif_in=NA}
    return(list(paste(gene[mt],collapse =";"),highest_rank,paste(motif_in,collapse = ";")))
    
    
  },GOdata=GOdata)
  siggene=unlist(lapply(siggene_forID,function(x) x[[1]]))
  max_rank=unlist(lapply(siggene_forID,function(x) x[[2]]))
  motif_in=unlist(lapply(siggene_forID,function(x) x[[3]]))
  if(nrow(sigres)>0){ 
    sigres$genes=siggene
    sigres$higest_ranks=max_rank
    sigres$motif=motif_in
  }
  return(sigres[sigres$FC>=1.5])
}

#CSV files are from Jason

csv_files=dir('../downstream/input/mm10_cluster/',pattern="csv")
tissue=unique(sub(".csv*","",csv_files))
#jsd_enhancer=readRDS('../downstream/input/chromHMM_enhancer.rds')
for (ts in tissue){
  fn=paste0(ts,'.csv')
  #read in csv file for given tissue
  csv_in_ts=fread(paste0('../downstream/input/mm10_cluster/',fn))
  csv_in_ts=csv_in_ts[order(dNME_maxJSD_rank,decreasing = F)]
  # Getting enhancer based on choice
  if(enc_type=="chromHMM_enhancer"){csv_in_ts=csv_in_ts[csv_in_ts$chromHMM_enhancer]}else
    if(enc_type=="bin_enhancer"){
      
      csv_in_gr=convert_GR(csv_in_ts$region)
      mcols(csv_in_gr)=csv_in_ts
      olap=findOverlaps(csv_in_gr,enhancer)
      csv_in_gr=csv_in_gr[queryHits(olap)]
      csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
      csv_in_gr$distance=NA
      csv_in_ts=as.data.table(mcols(csv_in_gr))
    }else
      if(enc_type=="TSS"){csv_in_ts=csv_in_ts[abs(distance)<=500]}
  
  if(nrow(csv_in_ts)>1){
    #GO annotation for each cluster
    csv_out=fastDoCall("rbind",lapply(1:10,function(clu){
      sp=paste0(ts,'-',clu)
      csv_in_ts_clu=csv_in_ts[cluster==clu]
      #rank csv based on dNME_maxJSD
      csv_in_ts_clu=csv_in_ts_clu[order(dNME_maxJSD_rank,decreasing=F)]
      
      if(nrow(csv_in_ts_clu)>1){
        
        #GO annotation for each cluster
        GO_out_cluster=GO_run(csv_in_ts_clu$gene,unique(csv_in_ts$gene),motif=motif_all)
        if(any(GO_out_cluster$motif!="NA")){cat(ts,clu,"Have motif in GO\n")}
        write.csv(GO_out_cluster,
                  file=paste0('../downstream/output/mm10_result/',enc_type,'/cluster_GO/',sp,'_cluster_GO.csv'),row.names = F,quote = T)
        csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
        return(csv_in_ts_clu)
        
      }
      
    }))
    write.csv(csv_out,file=paste0('../downstream/output/mm10_result/',enc_type,'/all_gene_list/',ts,'_all.csv'))
  }
  
}


#calculating dNME,dMML

UC=readRDS('../downstream/output/uc_matrix_DNase.rds')
chromHMM=readRDS('../downstream/output/chromHMM_enhancer.rds')
mml <- readRDS('../downstream/output/mml_matrix_DNase.rds')
nme <- readRDS('../downstream/output/nme_matrix_DNase.rds')

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- GenomicFeatures::genes(txdb)
promoters <- promoters(genes,upstream=2000,downstream=1000)
promoters$gene_name=AnnotationDbi::select(Mus.musculus,key=as.character(promoters$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
#find time order
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
UC <- lapply(UC,function(i) {
  i <- i[rowSums(i > 0.1) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  i <- i[,order(match(colnames(i),timeorder))]
  i <- i[complete.cases(i),]
  return(i)
})

#dMML and dNME
dmml=list()
dnme=list()
tissues=unique(sub('-.*','',colnames(mml)))
for (n in tissues){
  dmml[[n]] <- sapply(colnames(UC[[n]]),function(i) {
    time <- strsplit(i,'-')
    sapply(time,function(time_in)
    {abs(mml[rownames(UC[[n]]),paste0(n,'-',time_in[1],'-all')]-mml[rownames(UC[[n]]),paste0(n,'-',time_in[2],'-all')])})
  })
  colnames(dmml[[n]])=colnames(UC[[n]])
  rownames(dmml[[n]])=rownames(UC[[n]])
  
  
  dnme[[n]] <- sapply(colnames(UC[[n]]),function(i) {
    time <- strsplit(i,'-')
    sapply(time,function(time_in)
    {abs(nme[rownames(UC[[n]]),paste0(n,'-',time_in[1],'-all')]-nme[rownames(UC[[n]]),paste0(n,'-',time_in[2],'-all')])})
    
  })
  colnames(dnme[[n]])=colnames(UC[[n]])
  rownames(dnme[[n]])=rownames(UC[[n]])
}

# quantile plot -----------------------------------------------------------
matrix_all=list()
for(sp in names(UC)){
  dmml_in=dmml[[sp]]
  dnme_in=dnme[[sp]]
  UC_in=UC[[sp]]
  colnames(dmml_in)=paste0(colnames(dmml_in),'_','dMML')
  colnames(dnme_in)=paste0(colnames(dnme_in),'_','dNME')
  colnames(UC_in)=paste0(colnames(UC_in),'_','UC')
  matrix_tissue=as.data.table(cbind(UC_in,dmml_in,dnme_in),keep.rownames = T)
  colnames(matrix_tissue)[1]="regions"
  olap_chromHMM=findOverlaps(convert_GR(matrix_tissue$regions),chromHMM[chromHMM$tissue==sp])
  olap_promoter=findOverlaps(convert_GR(matrix_tissue$regions),promoters)
  matrix_tissue$promoter=FALSE
  matrix_tissue$enhancer=FALSE
  matrix_tissue$tissue=sp
  matrix_tissue$promoter[queryHits(olap_promoter)]=TRUE
  matrix_tissue$enhancer[queryHits(olap_chromHMM)]=TRUE
  matrix_tissue=melt.data.table(matrix_tissue,id.vars = c("regions","promoter","enhancer","tissue"))
  matrix_tissue$stage=sub('_.*','',matrix_tissue$variable)
  matrix_tissue$stat=sub('.*_','',matrix_tissue$variable)
  matrix_tissue=dcast.data.table(matrix_tissue,tissue+enhancer+promoter+stage+regions~stat)
  
  
  matrix_all[[sp]]=matrix_tissue
}
matrix_all=fastDoCall('rbind',matrix_all)
#matrix_all=matrix_all[matrix_all$UC!=1]
#quantiles within enhancers vs in promoters? within each feature, removed 1
matrix_all$UC_quant=findInterval(matrix_all$UC,quantile(matrix_all$UC,prob=c(0,0.25,0.5,0.75),na.rm=T))

matrix_all$UC_quant[matrix_all$UC_quant==5]=4#5th quantile is the maximum number, move to 4th
quant_conv=c("Q1","Q2","Q3","Q4")
matrix_all$UC_quant=quant_conv[matrix_all$UC_quant]
matrix_all$region_type="Non-regulatory"
#duplicate those overlapped lines, ones is promoter and the other is enhancer
matrix_all$region_type[matrix_all$promoter]="Promoter"
matrix_all$region_type[matrix_all$enhancer]="Enhancer"
cat(paste0(round(table(matrix_all$region_type)/length(matrix_all$region_type)*100,digits=2),"%"),"\n")
# matrix_all_agg=matrix_all[,list(dMML=median(dMML_quant),dNME=median(dNME_quant),UC=median(UC),
#                                 dMML_top25=quantile(dMML,prob=0.75),dNME_top25=quantile(dNME,prob=0.75),
#                                 dMML_bottom25=quantile(dMML,prob=0.25),dNME_bottom25=quantile(dNME,prob=0.25)),by=list(UC_quant,region_type)]
pdf('../downstream/output/graphs/Figure5/Figure5B_dNME_dMML_enhancer_UC_quantile.pdf',width=3.5,height=3.5)
dNME_plot=ggplot(matrix_all[region_type!="Non-regulatory"],aes(x=UC_quant,y=dNME,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.75))
dMML_plot=ggplot(matrix_all[region_type!="Non-regulatory"],aes(x=UC_quant,y=dMML,fill=region_type))+geom_boxplot(outlier.shape = NA)+theme_glob+
  xlab('UC quantile')+theme(legend.position = "bottom",legend.title = element_blank())+ylim(c(0,0.3))
ggarrange(dNME_plot,dMML_plot,nrow=2,ncol=1,common.legend=T)
dev.off()

direction_calc_enriched_subj<-function(motif_loc,variant_all,gene_in,pval_cutoff=0.1){
  
  motif_direction_out=mclapply(gene_in,direction_enriched_sample,
                               variant_gene=variant_all,motif_gene_subj=motif_loc,pval_cutoff=pval_cutoff,mc.cores =1)
  
  
  
  return(do.call(rbind,motif_direction_out))
}

direction_enriched_sample<-function(tf,variant_gene,motif_gene_subj,pval_cutoff,nperm=0){
  motif_gene_subj=motif_gene_subj[motif_gene_subj$geneSymbol==tf]
  variant_gene=variant_gene[variant_gene$dNME_pval<=pval_cutoff]
  olap=findOverlaps(variant_gene,motif_gene_subj)
  variant_gene=variant_gene[queryHits(olap)]
  variant_gene$alleleDiff=motif_gene_subj$alleleDiff[subjectHits(olap)]
  variant_gene=variant_gene[!is.na(variant_gene$alleleDiff)]
  
  #alleleDiff is calculated use ref - alt, prefer low ent ones
  variant_gene$NME_diff=variant_gene$altNME-variant_gene$refNME
  variant_gene$MML_diff=variant_gene$altMML-variant_gene$refMML
  same_dir=sum(sign(variant_gene$alleleDiff)== sign(variant_gene$NME_diff),na.rm = TRUE)
  opposite_dir=sum(sign(variant_gene$alleleDiff)!= sign(variant_gene$NME_diff),na.rm = TRUE)

  total_data=same_dir+opposite_dir
  
  variant_gene_df=data.frame(alleleDiff=sign(variant_gene$alleleDiff),NME_diff=sign(variant_gene$NME_diff))
  len_x=nrow(variant_gene_df)
  if(nperm>0){
    same_dir_perm=replicate(nperm,
                            sum(sample(variant_gene_df$alleleDiff,len_x,replace = F)==sample(variant_gene_df$NME_diff,len_x,replace = F)))
    same_dir_perm_prob=same_dir_perm/total_data
  }else{same_dir_perm_prob=-1}
  if(same_dir >0 &opposite_dir>0){
    
    binom=binom.test(same_dir,(same_dir+opposite_dir),0.5)

    prob_binom=binom$estimate[[1]]
    if(nperm>0){binom.pval=sum(abs(same_dir_perm_prob-0.5)>=abs(prob_binom-0.5))/nperm}else{binom.pval=NA}
    return(data.frame(TF=unique(motif_gene_subj$geneSymbol),total_data=total_data,same_dir=same_dir,opposite_dir=opposite_dir,
                      binom.pval_perm=binom.pval,binom.pval=binom$p.value,prob=prob_binom,NSNP=length(variant_gene),stringsAsFactors = F))
  }
}

#Read in allele-agnositc model
read.agnostic<-function(file_in,GR_merge_in,allele_include=T){
  stat=toupper(sub('.*_','',sub('.bedGraph','',file_in)))
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    olap=findOverlaps(GR_merge_in,informME_in)
    cat('Percent overlap with dNME region:',length(unique(subjectHits(olap)))/length(informME_in)*100,'%\n')
    informME_in$score_original=informME_in$score
    informME_in$dMML=NA
    informME_in$dMML_pval=NA
    informME_in$dMML[subjectHits(olap)]=GR_merge_in$dMML[queryHits(olap)]
    informME_in$dMML_pval[subjectHits(olap)]=GR_merge_in$dMML_pval[queryHits(olap)]
    #add GR_merge data
    if(allele_include){
      #replace value instead of remove regions
      informME_in$score[subjectHits(olap)]=
        rowMeans(as.matrix(elementMetadata(GR_merge_in)[paste(stat,c('1','2'),sep='')]))[queryHits(olap)]
    }else( informME_in=informME_in[-subjectHits(olap)])
    informME_in=informME_in[!is.infinite(informME_in$score)]
    quant=c("0-25%","25%-50%","50%-75%","75%-100%")
    informME_in$quant_score=findInterval(informME_in$score,quantile(informME_in$score,prob=c(0,0.25,0.5,0.75),na.rm=T))
    informME_in$quant_score=quant[informME_in$quant_score]
    informME_in$Sample=unique(GR_merge_in$Sample)
    informME_in$hyper_var_fn=unique(GR_merge_in$hyper_var_fn)
    return(informME_in)
  }
}