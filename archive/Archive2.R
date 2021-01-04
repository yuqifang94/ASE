# # reading in allele-agnostic data -----------------------------------------
NME_in=GRanges()
MML_in=GRanges()
GR_merge=readRDS(GR_merge_file)
GR_merge$fn=paste(GR_merge$Subject,GR_merge$tissue,sep='_')
GR_merge$fn[GR_merge$Sample=="Adipose_Tissue_single - STL003"]="STL003_Adipose_Tissue_single"
in_dir='../downstream/data/allele_agnostic_hg19_cov5_3kb_FANTOM/'
#For the mean expression: log2 or not?
for (fn in unique(GR_merge$fn)){
  
  fn_mml=paste(in_dir,fn,'_phased_allele_agnostic_mml.bedGraph',sep='')
  fn_nme=paste(in_dir,fn,'_phased_allele_agnostic_nme.bedGraph',sep='')
  cat('Processing',fn,'\n')
  if(file.exists(fn_mml)&file.exists(fn_nme)){
    
    cat('reading in:',fn_mml,'\n')
    cat('reading in:',fn_nme,'\n')
    NME_in=c(NME_in,read.agnostic(fn_nme,GR_merge[GR_merge$fn==fn]))
    MML_in=c(MML_in,read.agnostic(fn_mml,GR_merge[GR_merge$fn==fn],stat='MML'))
  }
}
NME_in$NME=NME_in$score
MML_in$MML=MML_in$score
saveRDS(NME_in,NME_agnostic_file)
saveRDS(MML_in,MML_agnostic_file)
NME_in=readRDS(NME_agnostic_file)
MML_in=readRDS(MML_agnostic_file)
NME_in_matrix=agnostic_matrix_conversion(NME_in[NME_in$bioreplicate==2])#duplicated regions due to error
MML_in_matrix=agnostic_matrix_conversion(MML_in[MML_in$bioreplicate==2],'MML')#duplicated regions due to error
saveRDS(NME_in_matrix,'../downstream/output/NME_matrix_cov5.rds')
saveRDS(MML_in_matrix,'../downstream/output/MML_matrix_cov5.rds')


# Add GWAS to each SNP ----------------------------------------------------
traits_hg38 = makeCurrentGwascat()
ch = import.chain('../downstream/input/hg38ToHg19.over.chain')
seqlevelsStyle(traits_hg38) = 'UCSC'
traits_hg19 = liftOver(traits_hg38, ch)
traits_hg19=unlist(traits_hg19)
traits_hg19=lapply(unique(traits_hg19$`DISEASE/TRAIT`), function(x) traits_hg19[traits_hg19$`DISEASE/TRAIT`==x])
saveRDS(traits_hg19,'../downstream/output/variant_traits.rds')
# traits=as.data.frame(get_traits()@traits)
# tt1=proc.time()[[3]]
# variants=mclapply(traits$efo_id,trait_variant,traits=traits,mc.cores=18)
# proc.time()[[3]]-tt1

# Input age data ----------------------------------------------------------
# age_out=GRanges()
# for (i in 1:22){
#   age_in=as.data.frame(read.csv(paste('../downstream/data/Age_atlas/atlas.chr',i,'.csv.gz',sep=''),skip = 3))
#   colnames(age_in)[2]='chr'
#   colnames(age_in)[3]='start'
#   age_in$chr=paste('chr',age_in$chr,sep='')
#   age_in$end=age_in$start
#   age_out=c(age_out,makeGRangesFromDataFrame(age_in,keep.extra.columns = T))
# }
# saveRDS(age_out,'../downstream/input/age_atlas.rds')
# saveRDS(subsetByOverlaps(age_out,variant_HetCpG_meta),'../downstream/input/age_atlas_sub.rds')
# elementMetadata(age_out)=elementMetadata(age_out)[,c('VariantID','AlleleRef','AlleleAlt',
#                                                      'AgeMode_Jnt','AgeMean_Jnt','AgeMedian_Jnt','QualScore_Jnt')]
# saveRDS(age_out,'../downstream/input/age_atlas_sub_jnt.rds')
#Index of /pub/ftp.ensembl.org/release-81/variation/gvf/homo_sapiens
#genonmic_1k=import('../downstream/input/1000GENOMES-phase_3.gvf.gz',format='gvf')
# #genonmic_1k_sub=genomic_1k[geno$type=='SNV']
# genome_1k_variant=subsetByOverlaps(genonmic_1k_sub,do.call(c,variant_all),maxgap = 1)
# genome_1k_variant$MAF=unlist(mclapply(genome_1k_variant$global_minor_allele_frequency, 
#                                       function(x) as.numeric(strsplit(x,"\\|")[[1]][2]),mc.cores=24))
# saveRDS(genome_1k_variant,'genome_1k_variant.rds')


# # reading in motif result -------------------------------------------------
dir='../downstream/input/motif_all/JASPAR/'
for (i in 1:2){
  #check if overlapped result are same
  motif_in=readRDS(paste(dir,'motif_all_JASPAR_',i,'_default.rds',sep=''))
  if(!exists("motif_out")){motif_out=motif_in}else{
    olap=findOverlaps(motif_in,motif_out)
    print(length(olap))
    if(length(olap)>0){
      motif_in=motif_in[-queryHits(olap)]
    }
    motif_out=c(motif_out,motif_in)
  }
}
#
saveRDS(motif_out,'../downstream/output/motif_all_JASPAR_default.rds')

#For each region in the gff, Cont number of Het CPG in it
CpG_hg19=readRDS('../downstream/input/CpG_hg19.rds')
#gff regions with HetCpG inforamtion
non_gff=lapply(hetCpG_gff,function(x) which(x$N_nonhet != x$N))#on level of 20000-30000 about 5%
#This is modified for new run, no longer needs to find smaller regions
GR_merge=readRDS(GR_merge_file)
#It's might be better to do it in reversed way.
variant_HetCpG=readRDS(variant_HetCpG_file)
#GR_merge=readRDS(GR_merge_file)
# H1_hyper_var=as.data.frame(read.csv("../downstream/input/H1 data/H1_hypervar_result.csv"))
# #H1_hyper_var=H1_hyper_var[H1_hyper_var$gene_type=='miRNA',]
# H1_hyper_var$hypervar_logvar=H1_hyper_var$hypervar_log2
# saveRDS(H1_hyper_var,'../downstream/output/H1_hypervar.rds')




#mm10_GO Example


for(ts in names(csv_in_ts_out)){
  write_csv=csv_in_ts_out[[ts]]
  write_csv$correlation_nme_mml=abs(write_csv$correlation_nme_mml)
  write_csv=write_csv[write_csv$region%in%rownames(nme),]
  write.csv(write_csv,file=paste0('../downstream/output/dNME_GO/cluster/',ts,'_ranked.csv'))
}
for(sp in names(GO_out_dNME_only)){
  GO_write=GO_out_dNME_only[[sp]][GO_out_dNME_only[[sp]]$Expected>=5,]
  GO_write$sigGenes=unlist(lapply(GO_write$sigGenes,paste,collapse=','))
  GO_write$sigGenes=unlist(lapply(GO_write$sigGenes,function(x){
    paste(x[match(unique(csv_in_ts_clu_ft$gene),x)],collapse=',')
  }))
  write.csv(GO_write,file=paste0('../downstream/output/dNME_GO/',sp,'_dNME_only_GO_out_090.csv'),row.names = F,quote = F)
  
}


for(sp in names(GO_out_cluster)){
  GO_write=GO_out_cluster[[sp]][GO_out_cluster[[sp]]$Expected>=5,c("Term","FDR")]
  GO_write$sigGenes=unlist(lapply(GO_write$sigGenes,paste,collapse=','))
  write.table(GO_write,file=paste0('../downstream/output/dNME_GO/',sp,'_cluster_GO_out.csv'),row.names = F,col.names = FALSE,sep=',',quote = F)
  
}

#Plot horizontal barplot use top 3 GO terms in cluster with corresponding color
GO_out_cluster=readRDS('../downstream/output/dNME_GO/cluster_GO.rds')
cluster_col<- brewer.pal(10,'Set3')
for(sp in names(GO_out_cluster)){
  GO_plot=GO_out_cluster[[sp]][GO_out_cluster[[sp]]$Expected>=5,c("Term","classicFisher","FDR")][1:5,]
  GO_plot$logp=-log10(as.numeric(GO_plot$classicFisher))
  clu_n=sub('.*-','',sp)
  pdf(paste0('../downstream/output/dNME_GO/pdf/',sp,'_dNME_cluster_GO.pdf'),width=11,height=4)
  print(ggplot(GO_plot,aes(x=Term,y=logp))+geom_bar(stat="identity",fill=cluster_col[as.numeric(clu_n)])+coord_flip()+
          theme()+xlab('')+ylab('-log10(p)'))
  dev.off()
  Sys.sleep(1)
}

#get midbrain regions
midbrain_Jason=fread('../downstream/output/dNME_GO/cluster/midbrain_all_csv.csv')

midbrain_Jason$correlation_nme_mml=
  corfunc(nme[midbrain_Jason$region,sub('-.*','',colnames(nme))=='midbrain'],
          mml[midbrain_Jason$region,sub('-.*','',colnames(mml))=='midbrain'])
midbrain_Jason$correlation_nme_mml=abs(midbrain_Jason$correlation_nme_mml)
midbrain_Jason=midbrain_Jason[order(maxpair_rankratio)]
write.csv(midbrain_Jason[correlation_nme_mml<=0.2&(!grepl('Gm',gene))],'../downstream/output/dNME_GO/cluster/midbrain_all_csv_cor_filter.csv') 
#c('Lats2','Otx1','Zfp428','Polr3g','Aqp11','Fam43a','Dnmt3a','Tead3','Notch4','Trim32')
midbrain_Jason_select=midbrain_Jason[midbrain_Jason$gene=='Notch4']
midbrain_Jason_select=midbrain_Jason_select[order(midbrain_Jason_select$dNME_maxpair,decreasing=T)][1]
midbrain_select_df=lapply(midbrain_Jason_select$region,function(x){
  x_out=data.table(sample=sub('-all','',names(nme[x,])),nme=nme[x,],mml=mml[x,])
  print(cor(x_out$nme,x_out$mml))
  x_out$tissue=sub('-.*','',x_out$sample)
  x_out$stage=sub('.*-','',x_out$sample)
  x_out=x_out[tissue=='midbrain']
  x_out=melt(x_out,id.vars=c('tissue','sample','stage'),measure.vars=c('nme','mml'),value.name='value',variable.name='Stat')
  x_out$stage=sub('day0','P0',x_out$stage)
  x_out$stage=sub('day','E',x_out$stage)
  x_out$stage=sub('_','.',x_out$stage)
  x_out$stage=factor(x_out$stage,levels=c(paste0('E',10:16,'.5'),'P0'))
  x_out=x_out[x_out$stage!="P0"]
  return(x_out)
})
ggplot(midbrain_select_df[[1]],aes(x=stage,y=value,group=Stat,color=Stat))+geom_point()+geom_line(aes(linetype=Stat),size=1)+ylim(c(0,1))+
  theme(legend.position = "bottom")

nme_pt=ggplot(midbrain_select_df[[1]][midbrain_select_df[[1]]$Stat=='nme'],aes(x=stage,y=value,group=Stat))+
  geom_point(color='red')+geom_line(size=1,color='red')+ylim(c(0,1))+ggtitle('NME')+xlab('')+
  theme(plot.title = element_text(hjust = 0.5))
mml_pt=ggplot(midbrain_select_df[[1]][midbrain_select_df[[1]]$Stat=='mml'],aes(x=stage,y=value,group=Stat))+
  geom_point(color='blue')+geom_line(size=1,color='blue')+ylim(c(0,1))+ggtitle('MML')+xlab('')+
  theme(plot.title = element_text(hjust = 0.5))

library(ggpubr)
pdf("../downstream/output/dNME_GO/pdf/NOTCH4.pdf")
ggarrange(nme_pt, mml_pt, ncol = 1, nrow = 2)
dev.off()
#select regions and visualize using  Gviz
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
chr="chr17"
gen="mm10"
grtrack=GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene)
symbols <- unlist(mapIds(org.Mm.eg.db, genes(TxDb.Mmusculus.UCSC.mm10.knownGene)$gene_id, "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(grtrack) <- symbols[gene(grtrack)]
itrack <- IdeogramTrack(genome =gen, chromosome = chr)
gtrack <- GenomeAxisTrack()
enhancer_track=AnnotationTrack(chromHMM_enhancer[seqnames(chromHMM_enhancer)==chr&chromHMM_enhancer$tissue=='forebrain'],name="Enhancer",fill="red")
select_region=GRanges(seqnames=sub(':.*','',midbrain_Jason_select$region),
                      ranges=IRanges(start=as.numeric(sub('-.*','',sub('.*:','',midbrain_Jason_select$region))),
                                     end=as.numeric(sub('.*-','',midbrain_Jason_select$region))))
target_region=AnnotationTrack(select_region,name="example",fill="purple")
pdf("../downstream/output/dNME_GO/pdf/NOTCH4_gr.pdf")
plotTracks(list(itrack, gtrack,grtrack,enhancer_track,target_region),from = start(select_region)-5000, 
           to = end(select_region)+10000, transcriptAnnotation = "symbol")
dev.off()


csv_all_enhancer_gr=GRanges(seqnames=sub(':.*','',csv_all_enhancer$region),
                            IRanges(start=as.numeric(sub('-.*','',sub('.*:','',csv_all_enhancer$region))),
                                    end=as.numeric(sub('.*-','',csv_all_enhancer$region))))


#calculate allelic difference and plot
allele_calc_plot<-function(gr_allele_in,stat,outDir,picname){
  gr_diff_calc=allele_diff_merge(gr_allele_in)
  allele_plot(gr_diff_calc,stat,outDir,picname)
  return(gr_diff_calc)
}
#Calculate allelic difference *Check
allele_diff_merge<-function(allele_gr_in){
  
  allele_gr_in$CpGdiff=allele_gr_in$g1CG-allele_gr_in$g2CG
  sign=allele_gr_in$CpGdiff
  sign[sign!=0]=sign[sign!=0]/abs(sign[sign!=0])
  sign[sign==0]=1
  allele_gr_in$sign=sign
  allele_gr_in$diff_NME=(allele_gr_in$NME1-allele_gr_in$NME2)*sign
  allele_gr_in$diff_MML=(allele_gr_in$MML1-allele_gr_in$MML2)*sign
  allele_gr_in$HetCpG=allele_gr_in$CpGdiff!=0
  
  return(allele_gr_in)
}

GO_custom <-function(selected_gene,gmt,background_gene){
  GO_out=lapply(gmt,GO_calc,selected=selected_gene,background=background_gene)
  GO_out_name=names(GO_out)
  GO_out=do.call(rbind,GO_out)
  GO_out$GO=GO_out_name
  GO_out=GO_out[GO_out$OR>1&!is.infinite(GO_out$OR),]
  GO_out$FDR=p.adjust(GO_out$pvalue,method='BH')
  return(GO_out)
}

GO_calc<-function(GO,selected,background){
  non_select=background[-which(background%in% selected)]
  GO_selected=sum(selected %in% GO)
  non_GO_selected=length(selected)-GO_selected
  GO_non_selected=sum(non_select %in% GO)
  non_GO_non_selected=length(non_select)-GO_non_selected
  cont_table=matrix(c(GO_selected,non_GO_selected,GO_non_selected,non_GO_non_selected),nrow=2)
  OR=fisher.test(cont_table)
  return(data.frame(OR=OR$estimate[[1]],pvalue=OR$p.value,lowerCI=OR$conf.int[1],upperCI=OR$conf.int[2]))
}

#make plot
allele_plot<-function(gr_diff_calc,stat,outDir,picname){
  gr_diff=gr_diff_calc
  
  gr_diff$type='Non Het CpG'
  gr_diff$type[gr_diff$CpGdiff!=0]='Het CpG'
  
  #Use density plot
  plot_df=data.frame(elementMetadata(gr_diff)[,paste("diff",stat,sep='_')],gr_diff$type)
  colnames(plot_df)=c('Allele_difference','Allele_type')
  
  # 
  # diff_plot=ggplot(plot_df,aes(x=Allele_difference,color=Allele_type))+
  #   geom_density(alpha=0.6,size=0.5)+xlab(paste(stat, 'Difference'))+
  #   theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
  #   ggtitle(paste('Allilic difference at',picname))+ylim(0,3)
  
  # diff_plot_het_CpG=ggplot(plot_df[plot_df$Allele_type=='Het CpG',],aes(x=Allele_difference))+
  #   geom_density(alpha=0.6,fill='purple',color='black',size=1.5)+xlab(paste(stat, 'Difference'))+
  #   theme(legend.position="none",plot.title = element_text(hjust=0.5))+
  #   ggtitle(paste('Allilic difference (More CpG-Less CpG) at',picname, 'with Het CpG'))+ylim(0,2)
  # diff_plot=arrangeGrob(diff_plot_het_CpG,diff_plot_non_het,nrow=2,ncol=1)
  # pdf(paste(outDir,stat,' diff ',picname,'.pdf',sep=''))
  # print(diff_plot)
  # dev.off()
  #More CpG genome2-genome1
  allele_gr_ASM=gr_diff_calc
  allele_gr_ASM=allele_gr_ASM[allele_gr_ASM$CpGdiff!=0]
  print(allele_gr_ASM)
  #Figure 4A * check
  allele_het_df=rbind(data.frame(value=c(elementMetadata(allele_gr_ASM)[,paste(stat,'1',sep='')][allele_gr_ASM$CpGdiff>0],
                                         elementMetadata(allele_gr_ASM)[,paste(stat,'2',sep='')][allele_gr_ASM$CpGdiff<0]),type='More CpG'),
                      data.frame(value=c(elementMetadata(allele_gr_ASM)[,paste(stat,'2',sep='')][allele_gr_ASM$CpGdiff>0],
                                         elementMetadata(allele_gr_ASM)[,paste(stat,'1',sep='')][allele_gr_ASM$CpGdiff<0]),type='less CpG'))
  print(wilcox.test(allele_het_df$value[allele_het_df$type=='More CpG'],allele_het_df$value[allele_het_df$type=='less CpG']))
  dist_plot=ggplot(allele_het_df,aes(x=value,color=type))+
    geom_density(alpha=0.6,size=0.5)+xlab(stat)+ggtitle(paste(stat,'at',picname,'with different types of haplotypes'))+
    theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+ylim(0,4)+scale_color_manual(values=c("blue","red"))+
    geom_hline(yintercept=0, colour="white", size=1)
  pdf(paste(outDir,stat,' dist ',picname,'.pdf',sep=''))
  print(dist_plot)
  dev.off()
  # #reference distribution
  # allele_gr_ASM=gr_diff_calc
  # allele_het_df=rbind(data.frame(value=elementMetadata(allele_gr_ASM)[,paste(stat,'1',sep='')],type='Ref'),
  #                     data.frame(value=elementMetadata(allele_gr_ASM)[,paste(stat,'2',sep='')],type='Alt'))
  # 
  # dist_plot_ref=ggplot(allele_het_df,aes(x=value,fill=type))+
  #   geom_density(alpha=0.6)+xlab(stat)+ggtitle(paste(stat,'at',picname,'reference distribution'))+
  #   theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
  #   scale_fill_manual( values = c("blue","red"))+ylim(0,4)
  # pdf(paste(outDir,stat,' dist ref',picname,'.pdf',sep=''))
  # print(dist_plot_ref)
  # dev.off()
}

#GO analysis
GO_anno<-function(myInterestingGenes,geneNames,topNodes=200){
  #Define interesting gene list
  # myInterestingGenes=names(sort(table(dmr_anno[[1]]$symbol),decreasing=T))
  # myInterestingGenes=gene_loc_mono$gene
  #Make factor
  geneList=factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList)=geneNames
  #find annotation orgs annot=annFUN.org,mapping="org.Hs.eg.db"
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, mapping="org.Hs.eg.db",ID="symbol")
  resultF <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  #resultF.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  #resultF.weight <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  #GT=GenTable(GOdata,classicFisher=resultF,elimFisher=resultF.elim,weightFisher=resultF.weight,orderBy="classicFisher",topNodes=topNodes)
  GT=GenTable(GOdata,classicFisher=resultF,orderBy="classicFisher",topNodes=topNodes)
  #FC
  GT$FC=GT$Significant/GT$Expected
  library(GO.db)
  goterms <- Term(GOTERM)
  GT[,2]=unlist(goterms[GT[,1]])
  return(list(GT,GOdata))
}

#Enrichment analysis, use granges, find overlapped SNPs with HetCpG, redo it for significant ranges 
motif_enrich<-function(motif_gene_strong,variant_in,pval_cutoff=0.1,dist=500,cont_cutoff=2){
  variant_sig_out_df=data.frame(gene_name=c(),OR=c(),lower_CI=c(),upper_CI=c(),pval=c(),dMML=c(),dNME=c(),stringsAsFactors = F)
  #variant_sig_out=GRanges()
  for(i in unique(motif_gene_strong$geneSymbol)){
    SNP_in=motif_gene_strong[motif_gene_strong$geneSymbol==i]
    #Cont table: ASM? in motifBreak?
    if(dist==0){
      variant_SNP=variant_in
    }else{
      variant_SNP=subsetByOverlaps(variant_in,SNP_in,maxgap =dist)
    }
    #variant_SNP_sp=variant_in
    #Find significant variants
    variant_sig=variant_SNP[which(variant_SNP$dNME_pval<=pval_cutoff)]
    variant_non_sig=variant_SNP[which(variant_SNP$dNME_pval>pval_cutoff)]
    #Use CMH test here
    count_all=data.frame(subjet=NULL,ASM=NULL,motif=NULL,count=NULL)
    for (sp in unique(variant_sig$Sample)){
      variant_sig_sp=variant_sig[variant_sig$Sample==sp]
      variant_non_sig_sp=variant_non_sig[variant_non_sig$Sample==sp]
      ASM_break_sp=subsetByOverlaps(variant_sig_sp,SNP_in,type="equal")
      non_ASM_break_sp=subsetByOverlaps(variant_non_sig_sp,SNP_in,type="equal")
      
      ASM_break_n_sp=length(ASM_break_sp)
      non_ASM_break_n_sp=length(non_ASM_break_sp)
      ASM_non_break_n_sp=length(variant_sig_sp)-ASM_break_n_sp
      non_ASM_non_break_n_sp=length(variant_non_sig_sp)-non_ASM_break_n_sp
      if(ASM_break_n_sp>=cont_cutoff & ASM_non_break_n_sp>=cont_cutoff&non_ASM_non_break_n_sp>=cont_cutoff & non_ASM_break_n_sp>=cont_cutoff){
        cont_sp=data.frame(subject=factor(sp),
                           ASM=factor(c("ASM","ASM","NonASM","NonASM")),
                           feature=factor(c("Motif","NonMotif","Motif","NonMotif")),
                           count=c(ASM_break_n_sp,ASM_non_break_n_sp,non_ASM_break_n_sp,non_ASM_non_break_n_sp))
        count_all=rbind(count_all,cont_sp)
        
      }
    }
    #count contengency table
    # ASM_break=subsetByOverlaps(variant_sig,SNP_in,type="within")
    # non_ASM_break=subsetByOverlaps(variant_non_sig,SNP_in,type="within")
    # ASM_non_break_n=length(variant_sig)-length(ASM_break)
    # non_ASM_non_break_n=length(variant_non_sig)-length(non_ASM_break)
    # ASM_break_n=length(ASM_break)
    # non_ASM_break_n=length(non_ASM_break)
    #Construct contangency table
    # if(ASM_break_n>0 & ASM_non_break_n>0&non_ASM_non_break_n>0 & non_ASM_break_n>0){
    #   #Store the dataframe and Granges
    #   #print(matrix(c(ASM_break_n,ASM_non_break_n,non_ASM_break_n,non_ASM_non_break_n),nrow=2))
    #   # GR_length=length(GR_merge_sub)
    #   fs=fisher.test(matrix(c(ASM_break_n,ASM_non_break_n,non_ASM_break_n,non_ASM_non_break_n),nrow=2))
    #   #variant_sig_df=data.frame(gene_name=rep(i,GR_length),OR=rep(fs$estimate[[1]],GR_length),lower_CI=rep(fs$conf.int[1],GR_length),
    #   #                         upper_CI=rep(fs$conf.int[2]),pval=rep(fs$p.value),stringsAsFactors = F)
    #   #elementMetadata(GR_merge_sub)=cbind(elementMetadata(GR_merge_sub),variant_sig_df)
    #   #variant_sig_out=c(variant_sig_out,GR_merge_sub)
    #   #variant_sig_df_uq=unique(variant_sig_df)
    #Reformat for CMHtest in new function
    fs=CMH_test(count_all)
    
    # count_all$subject=factor(count_all$subject,levels=unique(count_all$subject))
    # count_all$motif=factor(count_all$motif,levels=unique(count_all$motif))
    # count_all$ASM=factor(count_all$ASM,levels=unique(count_all$ASM))
    # #return(count_all)
    # print(count_all)
    # if (nrow(count_all)>0){
    # CMH_table=xtabs(count~ASM+motif+subject,data=count_all)
    # if (length(dim(CMH_table)[3]!=0)){
    # 
    # if(dim(CMH_table)[3]>1){
    # fs=mantelhaen.test(CMH_table)
    # }
    # else if(dim(CMH_table)[3]==1){
    #   #print(as.matrix(CMH_table[,,1]))
    #   #print(count_all)
    #   fs=fisher.test(as.matrix(CMH_table[,,1]))
    #   
    # }
    #  
    if(length(fs)!=0){
      variant_sig_df=data.frame(gene_name=i,OR=fs$estimate[[1]],lower_CI=fs$conf.int[1],
                                upper_CI=fs$conf.int[2],pval=fs$p.value,stringsAsFactors = F)
      #variant_sig_df$dNME=sum(ASM_break$altNME)-sum(ASM_break$refNME)
      #variant_sig_df$dMML=sum(ASM_break$altMML)-sum(ASM_break$altMML)
      variant_sig_out_df=rbind(variant_sig_out_df,variant_sig_df)
    }
    #}   
    # }
  }
  variant_sig_out_df$qval=p.adjust(variant_sig_out_df$pval,method="BH")
  #variant_sig_out$qval=variant_sig_out_df$qval[match(variant_sig_out$gene_name,variant_sig_out_df$gene_name)]
  return(variant_sig_out_df)
}


# agnostic_gen ------------------------------------------------------------


# #Modify length
# start(TSS)=start(TSS)-gff_length
# end(TSS)=end(TSS)+gff_length+1
# TSS$gene=TSS$gene_symbol
# TSS$region_type='TSS'
# mcols(TSS)=mcols(TSS)[,c('region_type','gene')]
# FANTOM$region_type='enhancer'
# FANTOM$gene=FANTOM$name
# mcols(FANTOM)=mcols(FANTOM)[,c('region_type','gene')]
# TSS=c(TSS,FANTOM)
# TSS=sort(TSS)
# strand(TSS)="*"
# TSS[width(TSS)<1000]=resize(TSS[width(TSS)<1000],fix='center',1000)
# #Merge overlapped regions
# TSS=TSS[seqnames(TSS)%in%seqlevels(TSS)[1:22]]
# TSS_break=reduce(TSS)
# TSS_break=subdivideGRanges(TSS_break,250)
# olap=as.data.frame(findOverlaps(TSS_break,TSS))
# olap_agg=aggregate(olap$subjectHits,by=list(olap$queryHits),function(x) list(x))
# names(olap_agg)=c('qt','st')
# TSS_break$region_type[olap_agg$qt]=lapply(olap_agg$st,function(x) TSS$region_type[x])
# TSS_break$gene[olap_agg$qt]=lapply(olap_agg$st,function(x) TSS$gene[x])

#deal with those non-overlapped regions with size of 10001
# TSS_min_range=TSS[width(TSS)==(gff_length*2+1)]
# TSS_min_range_split=split_ranges_5k(TSS_min_range,n=gff_length*2,size=200)
# #FOr those TSS is larger than 10001
# TSS_larger_range=TSS[width(TSS)>(gff_length*2+1)]
# TSS_larger_range=TSS_larger_range[seqnames(TSS_larger_range)%in% seqlevels(TSS_larger_range)[1:24]]


# #Deal with each ranges
# TSS_larger_range_split=GRanges()
# for(i in 1:length(TSS_larger_range)){
#   TSS_larger_range_split=c(TSS_larger_range_split,split_ranges_larger_5k(TSS_larger_range[i]),size=200)
# }


#TSS_old=import.gff('../downstream/output/H1_allele_agnostic_analysis_old.gff')#All old regions are in current ones

#gff_example=import.gff('../downstream/data/gff_file/H1_het.cpelasm.gff')


# Finding examples -----------

#Not in Ecker's paper but functional
#select regions and visualize using  Gviz
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
chromHMM_enhancer=readRDS("../downstream/output/chromHMM_enhancer.rds")
bin_enhancer=readRDS("../downstream/output/bin_enhancer.rds")
#Prepare necessary tracks
grtrack=GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene)
symbols <- unlist(mapIds(org.Mm.eg.db, GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm10.knownGene)$gene_id, "SYMBOL", "ENTREZID", 
                         multiVals = "first"))
symbol(grtrack) <- symbols[gene(grtrack)]
grtrack=grtrack[!is.na(symbol(grtrack))]
gtrack <- GenomeAxisTrack()
plot_regions_mouse<-function(tissue,select_region,size_forward,size_backward,motif_annotation=AnnotationTrack()){
  #getting mml and nme
  mml_example=mml[select_region$region,sub('-.*','',colnames(mml))==tissue]
  nme_example=nme[select_region$region,names(mml_example)]
  dmml_example=dmml[[tissue]][select_region$region,]
  dnme_example=dnme[[tissue]][select_region$region,]
  UC_example=UC[[tissue]][select_region$region,]
  nme_mml_dt=data.table(stage=sub('-all','',sub(paste0(tissue,'-'),'',names(mml_example))),
                        MML=mml_example,NME=nme_example)
  nme_mml_dt=melt.data.table(nme_mml_dt[stage!="P0"],id.vars = "stage",variable.name = "Statistics")
  dnme_dmml_UC_dt=data.table(stage=names(UC_example),UC=UC_example,dMML=dmml_example,dNME=dnme_example)
  dnme_dmml_UC_dt=melt.data.table(dnme_dmml_UC_dt,id.vars = "stage",variable.name = "Statistics")
  #covert to GR
  select_region=convert_GR(select_region)
  target_region=AnnotationTrack(select_region,name="example",fill="purple")
  chr=as.character(seqnames(select_region))
  itrack <- IdeogramTrack(genome ="mm10", chromosome = chr)
  chromHMM_anno=AnnotationTrack(subsetByOverlaps(chromHMM_enhancer[chromHMM_enhancer$tissue==tissue],
                                                 select_region,maxgap = 1000),name="Enhancer chromHMM",fill="red")
  bin_enhancer_track=AnnotationTrack(subsetByOverlaps(bin_enhancer,select_region,maxgap = 1000),name="Enhancer ref",fill="pink")
  Gviz::plotTracks(list(itrack, gtrack,grtrack,bin_enhancer_track,chromHMM_anno,motif_annotation,target_region),from = start(select_region)-size_forward, 
                   to = end(select_region)+size_backward, transcriptAnnotation = "symbol")
  
  
  nme_mml_plot=ggplot(nme_mml_dt,aes(x=stage,y=value,group=Statistics,color=Statistics))+geom_line(size=1)+theme_glob+theme(legend.position = "bottom")+
    ylim(c(0,1))+scale_color_manual(values=c('red','blue'))+theme(axis.title.x=element_blank(),axis.title.y = element_blank())
  dNME_dMML_UC_plot=ggplot(dnme_dmml_UC_dt,aes(x=stage,y=value,group=Statistics,color=Statistics))+geom_line(size=1)+theme_glob+
    theme(legend.position = "bottom",axis.text.x = element_text(angle = 90),axis.title.x=element_blank(),axis.title.y = element_blank())+
    scale_color_manual(values=c('magenta','red','blue'))
  ggarrange(nme_mml_plot,dNME_dMML_UC_plot,nrow=1,ncol=2)
  
}

# Heart -------------------------------------------------------------------
tissue="heart"
#Select regions from Heart
heart_GO=fread('../downstream/output/mm10_result/bin_enhancer/all_gene_list/heart_all.csv')
heart_GO=heart_GO[(GO_result!="")&chromHMM_enhancer]
heart_GO=heart_GO[order(dNME_maxJSD,decreasing = T)]
DNase=readRDS('../downstream/input/mm10_DNase.rds')
#N>=3
heart_GON3=subsetByOverlaps(convert_GR(heart_GO$region),DNase[DNase$N>=5])
heart_GON3=paste0(seqnames(heart_GON3),':',start(heart_GON3),'-',end(heart_GON3))
heart_GO[region %in% heart_GON3,c("region","cluster","dNME_maxJSD")][1:10]
#min NME<=0.3?
#Tbx1: Epithelial Properties of the Second Heart Field
pdf('../downstream/output/graphs/example_Mouse/Popdc3.pdf',width=7,height=3.5)
plot_regions_mouse("heart",heart_GO[gene=="Popdc3","region"],25000,10000)
dev.off()
#Example motif BMP4, klf2: 
pdf('../downstream/output/graphs/example_Mouse/BMP4.pdf',width=7,height=3.5)
plot_regions_mouse("heart",heart_GO[gene=="Bmp4","region"],10000,80000)
dev.off()
#Wnt signaling through Dishevelled, Rac and JNK regulates dendritic development
pdf('../downstream/output/graphs/example_Mouse/Uncx.pdf',width=7,height=3.5)
plot_regions_mouse("forebrain",heart_GO[gene=="Wnt7b"&cluster==4,"region"],20000,10000)
dev.off()
#Example motif Gata4: 
pdf('../downstream/output/graphs/example_Mouse/Mef2a.pdf',width=7,height=3.5)
plot_regions_mouse("heart",heart_GO[gene=="Mef2a"&cluster==10,"region"],15000000,10000)
dev.off()

#Generate bed files from all csv files
gene_enhancer_out=GRanges()
for (sp in names(UC)){
  gene_enhancer_out=c(gene_enhancer_out,convert_GR(rownames(UC[[sp]])))
  
}
gene_enhancer_out=unique(gene_enhancer_out)
export.bed(resize(gene_enhancer_out,500,fix='center'),'../downstream/output/DNase_cluster.bed')

