source('mainFunctions_sub.R')
#Look at GO results
GO_out=readRDS(paste0(GO_01_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_01_enhancer.rds'))
UC_merge=readRDS(UC_merge_file)#Define all analyzed regions, were using UC_merge_max_loc_cluster01.rds,4626
cutoff_fn='01'
#Runnning
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
mapping="org.Mm.eg.db"
#prepare enhancer background gene list
uc_gr=lapply(UC_merge,function(x) rownames(x))
uc_gr=Reduce(intersect,uc_gr)
uc_gr=convert_GR(uc_gr)
enhancer=readRDS(bin_enhancer_rds)#21441
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)
rm(enhancer_bg)
rm(uc_gr)
geneList=rep(1,length(bg_enhancer))
names(geneList)=bg_enhancer
GOdata_all <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                                  annot = annFUN.org, mapping = mapping, ID = "Symbol")
tissue_specific_GO=data.table(tissue_all=tissue_all,
                              GO_grepl=c("facial|odontogenesis|chondrocyte|myoblast|oteoblast",
                                         "axon|neuron|oligodendrocyte|synapse|ganglion|brain|forebrain|telencephalon|gyrus",
                                         "heart|cardiac|valve|aorta|ventricle|ventricular",
                                         "postsynapse|nervous|cerebellar|neural|hindbrain",
                                         "cartilage|bone|chondrocyte|epithelial",
                                         "hemopoiesis|leukocyte|hemopothesis|lymphocyte|monosaccharide",
                                         "glial|spinal cord|diencephalon|collateral sprouting|neuron")
)
#For heart tissue
tissue_UC_diff=lapply(tissue_all,function(tissue){
    #Find UC in the tissue
    UC_tissue=UC_merge[[tissue]]
    #Annotate regions to genes 
    UC_tissue_olap=findOverlaps(convert_GR(rownames(UC_tissue)),enhancer)
    UC_tissue_enhancer=UC_tissue[queryHits(UC_tissue_olap),]
    UC_tissue_enhancer$gene="NA"
    UC_tissue_enhancer$gene=enhancer$`Target Gene`[subjectHits(UC_tissue_olap)]
    UC_tissue_enhancer=UC_tissue_enhancer[order(UC_tissue_enhancer$UC_max_pair,decreasing=T),]
    GO_tissue=GO_out$all[[tissue]]
    #Find GO output
    tissueRegion=fread(paste0(dir_out_cluster01,"/",tissue,".csv"))
    #Find significant terms
    GO_tissue_clu=do.call(rbind,lapply(GO_tissue,function(x) x$GO_out_cluster_all[FDR<=0.2&FC>=1.5]))
    #GO_tissue_clu=GO_tissue[[1]]$GO_out_cluster_all[FDR<=0.2&FC>=1.5]
    grep_term=tissue_specific_GO[tissue_all==tissue]$GO_grepl
    tissue_genes=unlist(genesInTerm(GOdata_all, GO_tissue_clu[grepl(grep_term,Term)]$GO.ID))
    write.csv(GO_tissue_clu[grepl(grep_term,Term)],paste0("../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/",tissue,"_specific_term.csv"))

    UC_tissue_enhancer=UC_tissue_enhancer[order(UC_tissue_enhancer$UC_max_pair,decreasing=T),]
    UC_tissue_enhancer$ts_GO=FALSE
    UC_tissue_enhancer[UC_tissue_enhancer$gene%in%tissue_genes,]$ts_GO=TRUE
    UC_tissue_enhancer_dt=as.data.table(UC_tissue_enhancer)
    UC_tissue_enhancer_dt$rankUC=1:nrow(UC_tissue_enhancer_dt)
    UC_tissue_enhancer_dt_gene=UC_tissue_enhancer_dt[,list(UC_max_gene=max(UC_max_pair),UC_mean_gene=mean(UC_max_pair)),by=list(gene,ts_GO)]
    UC_tissue_enhancer_dt_gene=UC_tissue_enhancer_dt_gene[order(UC_max_gene,decreasing=T)]
    UC_tissue_enhancer_dt_gene$tissue=tissue
    #UC_tissue_enhancer_dt_gene$rankUCGene=nrow(UC_tissue_enhancer_dt_gene)-rank(UC_tissue_enhancer_dt_gene$UC_max_gene)+1
    return(UC_tissue_enhancer_dt_gene)

})
saveRDS(tissue_UC_diff,'../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/tissue_UC_diff.rds')
#Calculate difference
tissue_UC_diff_summary=lapply(tissue_UC_diff,function(x){
    tstat_out=t.test(x[ts_GO==TRUE]$UC_max_gene,x[ts_GO==FALSE]$UC_max_gene,alternative="greater")
    return(data.table(diff_mean=-diff(tstat_out$estimate),p_value=tstat_out$p.value,tissue=unique(x$tissue)))

})
range(do.call(rbind,tissue_UC_diff_summary)$p_value)
pdf("../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/region_UC.pdf")
print(ggplot(UC_tissue_enhancer_dt,aes(x=rankUC,y=UC_max_pair,color=ts_GO))+geom_point(size=0.00001)+theme(axis.text.x=element_blank()))
dev.off()

pdf("../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/gene_UC.pdf")
print(ggplot(UC_tissue_enhancer_dt_gene,aes(x=rankUCGene,y=UC_max_gene,color=ts_GO))+geom_point(size=0.00001)+theme(axis.text.x=element_blank()))
dev.off()
pdf("../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/gene_UC_ts.pdf")
print(ggplot(UC_tissue_enhancer_dt_gene,aes(x=ts_GO,y=UC_max_gene))+geom_boxplot()+xlab("gene related to tissue")+ylab("UC"))
dev.off()
pdf("../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/region_UC_ts.pdf")
print(ggplot(UC_tissue_enhancer_dt,aes(x=ts_GO,y=UC_max_pair))+geom_boxplot(outlier.shape=NA)+xlab("region associated with gene related to tissue")+ylab("UC")+ylim(c(0,0.15)))
dev.off()
pdf("../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/region_UC_ts_clu1.pdf")
print(ggplot(UC_tissue_enhancer_dt,aes(x=ts_GO,y=`UC-heart-E10.5-E11.5-all`))+geom_boxplot(outlier.shape=NA)+xlab("region associated with gene related to tissue")+ylab("UC")+ylim(c(0,0.15)))
dev.off()

UC_tissue_enhancer_dt_gene=UC_tissue_enhancer_dt[,list(UC_max_gene=max(UC-heart-E10.5-E11.5-all),UC_mean_gene=mean(UC-heart-E10.5-E11.5-all)),by=list(gene,ts_GO)]
UC_tissue_enhancer_dt_gene=UC_tissue_enhancer_dt_gene[order(UC_max_gene,decreasing=T)]
UC_tissue_enhancer_dt_gene$rankUCGene=nrow(UC_tissue_enhancer_dt_gene)-rank(UC_tissue_enhancer_dt_gene$UC_max_gene)+1
pdf("../downstream/output/mouse_analysis/GO_analysis/GO_enrichment/gene_UC_ts_clu1.pdf")
print(ggplot(UC_tissue_enhancer_dt_gene,aes(x=ts_GO,y=UC_max_gene))+geom_boxplot(outlier.shape=NA)+xlab("region associated with gene related to tissue")+ylab("UC")+ylim(c(0,0.15)))
dev.off()