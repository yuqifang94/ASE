# Currently Not in use ----------------------------------------------------
# 
# 
# #collapase variants
# variants_collapase<-function(varsDiff){
#   variants <- paste(as.character(varsDiff$REF),as.character(unlist(varsDiff$ALT)),sep="-")
#   # Combine same variants (ALT/REF -> REF/ALT)
#   variants[variants %in% c("A-C","C-A")] <- "A-C"
#   variants[variants %in% c("A-G","G-A")] <- "A-G"
#   variants[variants %in% c("A-T","T-A")] <- "A-T"
#   variants[variants %in% c("C-G","G-C")] <- "C-G"
#   variants[variants %in% c("C-T","T-C")] <- "C-T"
#   variants[variants %in% c("G-T","T-G")] <- "G-T"
#   return(variants)
# }
# 
# #Merge allele agnositc data 
# allele_agnostic_merge<-function(GR_in,nme_in,mml_in,pval_cutoff=0.1){
#   NME_in=import.bedGraph(nme_in)
#   if(all(seqlevels(NME_in)==gsub('chr','',seqlevels(NME_in)))){seqlevels(NME_in)=paste('chr',seqlevels(NME_in),sep='')}
#   MML_in=import.bedGraph(mml_in)
#   if(all(seqlevels(MML_in)==gsub('chr','',seqlevels(MML_in)))){seqlevels(MML_in)=paste('chr',seqlevels(MML_in),sep='')}
#   olap=findOverlaps(NME_in,GR_in)
#   agnostic_diff_df=data.frame(NME=NME_in$score[queryHits(olap)],dNME=GR_in$dNME[subjectHits(olap)],
#                               dNME_pval=GR_in$dNME_pval[subjectHits(olap)],dMML=GR_in$dMML[subjectHits(olap)],
#                               dMML_pval=GR_in$dMML_pval[subjectHits(olap)],MML=MML_in$score[queryHits(olap)])
#   agnostic_diff_df$dNME_ASM=agnostic_diff_df$dNME_pval<=pval_cutoff
#   agnostic_diff_df$dMML_ASM=agnostic_diff_df$dMML_pval<=pval_cutoff
#   agnostic_diff_df$dNME_ASM_non_dMML=agnostic_diff_df$dMML_pval>pval_cutoff & agnostic_diff_df$dNME_pval<=pval_cutoff
#   agnostic_diff_df$sample=unique(GR_in$Sample)
#   return(agnostic_diff_df)
#   
# }
# 
# #Processing RNA-seq data
# RNA_seq_process<-function(dir_in,fds,condition_rep=3){
#   files=paste(dir_in,fds,"/t_data.ctab",sep="")
#   tmp=read_tsv(files[1])
#   tx2gene <- tmp[, c("t_name", "gene_name")]
#   txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)
#   txi[1:3]=lapply(txi[1:3],function(x) {colnames(x)=fds 
#   return(x)})
#   sampleTable <- data.frame(condition = factor(rep(c("genome1", "genome2"), each =condition_rep)))
#   rownames(sampleTable) <- colnames(txi$counts)
#   ####Perfrom differential RNA analysis
#   dds_RNA<- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
#   if(ncol(as.data.frame(assay(dds_RNA)))==2){
#     res_RNA<-as.data.frame(cpm(assay(dds_RNA)))
#     res_RNA$gene_name=rownames(res_RNA)
#     res_RNA=res_RNA[(res_RNA[,1]!=0&res_RNA[,2]!=0),]
#     res_RNA=res_RNA[(res_RNA[,1]>10&res_RNA[,2]>10),]
#     res_RNA$log2FoldChange=log2((res_RNA[,2])/(res_RNA[,1]))
#   }else{
#     dds_RNA<-DESeq(dds_RNA)
#     res_RNA<-results(dds_RNA,name= "condition_genome2_vs_genome1")
#   }
#   return(res_RNA)
# }
# 
# #Tissue to germlayer
tissue_to_germlayer<-function(GR_input){
  GR_input$germlayer=NA
  tissue_ectoderm=c("foreskin_keratinocyte_paired",
                    "foreskin_melanocyte_paired",
                    "ectoderm_paired",
                    "brain_cerebellum_tissue_paired",
                    "brain_germinal_matrix_tissue_paired",
                    "Brain_substantia_nigra_paired",
                    "Brain_Hippocampus_middle_paired" )
  tissue_mesoderm=c("mesoderm_23_paired","Adipose_single",
                    "Left_Ventricle_single","Psoas_Muscle_single" ,
                    "Right_Ventricle_single","Right_Atrium_single","Spleen_single",
                    "Adrenal_Gland_single","Aorta_single","Ovary_single")
  tissue_endoderm=c("Small_Intestine_single","Lung_single","endoerm_27_paired",
                    "Bladder_single" ,"Gastric_single", "Sigmoid_Colon_single",
                    "Thymus_single","Esophagus_single", "Pancreas_single" ,"Liver_single")
  tissue_ESC=c("rep1","rep2","merged","42_embryonic_stem_cell_single" , "stem_27_undifferentiated_paired",'ESC')
  GR_input$germlayer[GR_input$tissue %in% tissue_ectoderm]='ectoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_mesoderm]='mesoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_endoderm]='endoderm'
  GR_input$germlayer[GR_input$tissue %in% tissue_ESC]='ESC'
  return(GR_input)
}
# 
# 



# 
# 
# 
# OR_VMR<-function(NME_dat,vmr,percent,NME_quant='quant_score'){
#   NME_dat$VMR=FALSE
#   olap=findOverlaps(NME_dat,vmr)
#   NME_dat$VMR[unique(queryHits(olap))]=TRUE
#   NME_VMR=sum((elementMetadata(NME_dat)[[NME_quant]]%in%percent)&NME_dat$VMR)
#   nonNME_VMR=sum((!elementMetadata(NME_dat)[[NME_quant]]%in%percent)&NME_dat$VMR)
#   nonNME_nonVMR=sum((!elementMetadata(NME_dat)[[NME_quant]]%in%percent)&!NME_dat$VMR)
#   NME_nonVMR=sum((elementMetadata(NME_dat)[[NME_quant]]%in%percent)&!NME_dat$VMR)
#   #print(matrix(c(NME_VMR,nonNME_VMR,NME_nonVMR,nonNME_nonVMR),nrow = 2))
#   fisher.test(matrix(c(NME_VMR,nonNME_VMR,NME_nonVMR,nonNME_nonVMR),nrow = 2))
# }
# 
# get_traits_GWAS<-function(variant_in,trait_gr_in,pval_cutoff=0.1,count_cutoff=3,stat='dNME_pval',CMH=FALSE,maxgap=500,ncores=15){
#   traits_ls=mclapply(trait_gr_in,get_traits_GWAS_all_trait,
#                    variant_in=variant_in,pval_cutoff=pval_cutoff,count_cutoff=count_cutoff,stat=stat,CMH=CMH,maxgap=maxgap,
#                    mc.cores=ncores)
#   return(traits_ls)
#   #list(do.call(rbind,lapply(traits_ls,function(x)x[[1]])),do.call(rbind,lapply(traits_ls,function(x)x[[2]])))
# }
# get_traits_GWAS_all_trait<-function(trait_gr,variant_in,pval_cutoff,count_cutoff,stat,CMH,maxgap=500){
#   #trait_gr=trait_gr[trait_gr$trait%in%trait]
#   trait=unique(trait_gr$`DISEASE/TRAIT`)
#   OR_output=data.frame()
#   CMH_df=data.frame()
#   dNME_sig=variant_in[elementMetadata(variant_in)[,stat]<=pval_cutoff]
#   dNME_non_sig=variant_in[elementMetadata(variant_in)[,stat]>pval_cutoff]
#   dNME_traits_gr=findOverlaps(dNME_sig,trait_gr,maxgap = maxgap)
#   # dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
#   # non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
#   # dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
#   # non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
#   dNME_trait=length(unique(queryHits(dNME_traits_gr)))
#   dNME_non_trait=length(dNME_sig)-dNME_trait
#   non_dNME_trait=length(subsetByOverlaps(dNME_non_sig,trait_gr,maxgap = maxgap))
#   non_dNME_non_trait=length(dNME_non_sig)-non_dNME_trait
#   trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
#   journal=paste(unique(trait_gr$JOURNAL),collapse = ',')
#   trait_gr_out=trait_gr[subjectHits(dNME_traits_gr)]
#   mcols(trait_gr_out)=mcols(trait_gr_out)[c('MAPPED_GENE','STRONGEST SNP-RISK ALLELE','DISEASE/TRAIT')]
#   names(mcols(trait_gr_out))=c('genes','risk allele','traits')
#   rm(dNME_sig)
#   rm(dNME_non_sig)
#   rm(trait_gr)
#   gc()
#   if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
#     
#     CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
#                       count=trait_count,subject='All')
#     cat('processing trait',trait,'\n')
#     
#     if(length(CMH_df)>0){
#       if (CMH){return(CMH_df)}
#       else{
#         
#         OR=CMH_test(CMH_df)
#         OR_output=cbind(t(CMH_df$count),data.frame(trait=trait,OR=OR$estimate,
#                                                    p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
#         colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
#         OR_output$journal=journal
#         
#         return(list(OR_output,trait_gr_out))
#       }
#     }
#   }
# }
# #Modify this to adapt trait analysis
# get_traits_GWAS_all_trait_single<-function(variant_in,trait_gr,pval_cutoff,count_cutoff,stat,CMH,maxgap){
#  
#   OR_output=data.frame()
#   CMH_df=data.frame()
#   dNME_sig=variant_in[elementMetadata(variant_in)[,stat]<=pval_cutoff]
#   dNME_non_sig=variant_in[elementMetadata(variant_in)[,stat]>pval_cutoff]
#   
#   # dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
#   # non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
#   # dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
#   # non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
#   dNME_trait=length(subsetByOverlaps(dNME_sig,trait_gr,maxgap = maxgap))
#   dNME_non_trait=length(dNME_sig)-dNME_trait
#   non_dNME_trait=length(subsetByOverlaps(dNME_non_sig,trait_gr,maxgap = maxgap))
#   non_dNME_non_trait=length(dNME_non_sig)-non_dNME_trait
#   trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
#   rm(dNME_sig)
#   rm(dNME_non_sig)
#   rm(trait_gr)
#   gc()
#   if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
#     
#     CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
#                       count=trait_count,subject='All')
# 
#     
#     if(length(CMH_df)>0){
#       if (CMH){return(CMH_df)}
#       else{
#         
#         OR=CMH_test(CMH_df)
#         OR_output=cbind(t(CMH_df$count),data.frame(OR=OR$estimate,
#                                                    p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
#         colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
#         OR_output$journal=paste(unique(trait_gr$JOURNAL),collapse = ',')
#         return(OR_output)
#       }
#     }
#   }
# }
# 
# 
# 
# get_traits_GWAS_trait<-function(trait,variant_in,pval_cutoff,count_cutoff,stat,CMH){
#   traits_sp_ls=lapply(unique(variant_in$germlayer),get_traits_GWAS_sp_trait,trait=trait,
#                       variant_in=variant_in,pval_cutoff=pval_cutoff,count_cutoff=count_cutoff,stat=stat,CMH=CMH)
#   traits_sp=do.call(rbind,traits_sp_ls)
#   if(CMH&length(traits_sp)>0){
#     OR=CMH_test(traits_sp)
#     OR_output=data.frame(trait=trait,OR=OR$estimate,p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2])
#     return(OR_output)
#   }else if(!CMH){
#     return(traits_sp)
#   }
# 
# 
# }
# 
# 
# 
# 
# get_traits_GWAS_sp_trait<-function(sp,trait,variant_in,pval_cutoff,count_cutoff,stat,CMH){
#  
#       OR_output=data.frame()
#       CMH_df=data.frame()
#       variant_in_sp=variant_in[variant_in$germlayer==sp&!is.na(variant_in$trait),]
#  
#       dNME_sig=variant_in_sp[,stat]<=pval_cutoff
#       dNME_non_sig=variant_in_sp[,stat]>pval_cutoff
#     
#       dNME_trait=sum(unlist(variant_in_sp$trait[dNME_sig])==trait)
#       non_dNME_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])==trait)
#       dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_sig])!=trait)
#       non_dNME_non_trait=sum(unlist(variant_in_sp$trait[dNME_non_sig])!=trait)
#       trait_count=c(dNME_trait,non_dNME_trait,dNME_non_trait,non_dNME_non_trait)
#       if(all(!is.na(trait_count))&all(trait_count>=count_cutoff)){
#        
#         CMH_df=data.frame(ASM=c('Yes','Yes','No','No'),feature=c('Yes','No','Yes','No'),
#                                        count=trait_count,subject=sp)
#         cat('processing',sp,'with trait',trait,'\n')
# 
#       if(length(CMH_df)>0){
#         if (CMH){return(CMH_df)}
#         else{
#           OR=CMH_test(CMH_df)
#           OR_output=cbind(t(CMH_df$count),data.frame(subject=sp,trait=trait,OR=OR$estimate,
#                                                      p_value=OR$p.value,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
#           colnames(OR_output)[1:4]=c('dNME_trait','non_dNME_trait','dNME_non_trait','non_dNME_non_trait')
#           return(OR_output)
#         }
#       }
#      }
# }

# 
# 
# ASM_het_enrichment<-function(gr_in){
#   ASM_het=sum(gr_in$ASM=='Yes' & gr_in$HetCpG)
#   ASM_non_het=sum(gr_in$ASM=='Yes' & !gr_in$HetCpG)
#   non_ASM_non_het=sum(gr_in$ASM=='No' & !gr_in$HetCpG)
#   non_ASM_het=sum(gr_in$ASM=='No' & gr_in$HetCpG)
#   cont_table=matrix(c(ASM_het,ASM_non_het,non_ASM_het,non_ASM_non_het),nrow=2)
#   print(cont_table)
#   fisher.test(cont_table)
# }
# testEnrichmentFeature <- function(dataGR,featureGR){
#   
#   # Find ranges overlapping with feature
#   olaps <- findOverlaps(dataGR,featureGR,type="any",select="all")
#   indFeature <- queryHits(olaps)
#   featureData <- dataGR[indFeature]
#   complementaryData <- dataGR[-indFeature]
#   
#   # Enrichment of dMML-HASM in feature
#   featureDmml <- featureData[featureData$Statistic=="dMML"]
#   complementaryDmml <- complementaryData[complementaryData$Statistic=="dMML"]
#   contTableDmml <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
#   rownames(contTableDmml) <- c("Feature","Complementary")
#   contTableDmml[1,]$ASM <- sum(featureDmml$ASM=="Yes")
#   contTableDmml[1,]$nonASM <- sum(featureDmml$ASM=="No")
#   contTableDmml[2,]$ASM <- sum(complementaryDmml$ASM=="Yes")
#   contTableDmml[2,]$nonASM <- sum(complementaryDmml$ASM=="No")
#   dmmlFisher <- fisher.test(contTableDmml)
#   
#   # Enrichment of dNME-HASM in feature
#   featureDnme <- featureData[featureData$Statistic=="dNME"]
#   complementaryDnme <- complementaryData[complementaryData$Statistic=="dNME"]
#   contTableDnme <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
#   rownames(contTableDnme) <- c("Feature","Complementary")
#   contTableDnme[1,]$ASM <- sum(featureDnme$ASM=="Yes")
#   contTableDnme[1,]$nonASM <- sum(featureDnme$ASM=="No")
#   contTableDnme[2,]$ASM <- sum(complementaryDnme$ASM=="Yes")
#   contTableDnme[2,]$nonASM <- sum(complementaryDnme$ASM=="No")
#   dnmeFisher <- fisher.test(contTableDnme)
#   
#   # Enrichment of UC-HASM in feature
#   featureUc <- featureData[featureData$Statistic=="UC"]
#   complementaryUc <- complementaryData[complementaryData$Statistic=="UC"]
#   contTableUc <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
#   rownames(contTableUc) <- c("Feature","Complementary")
#   contTableUc[1,]$ASM <- sum(featureUc$ASM=="Yes")
#   contTableUc[1,]$nonASM <- sum(featureUc$ASM=="No")
#   contTableUc[2,]$ASM <- sum(complementaryUc$ASM=="Yes")
#   contTableUc[2,]$nonASM <- sum(complementaryUc$ASM=="No")
#   ucFisher <- fisher.test(contTableUc)
#   
#   # Return list of Fisher's test
#   return(list(dmmlFisher,dnmeFisher,ucFisher))
#   
# }
# 
# getGeneralFeats_mm9 <- function(CpGdir,enhancerDir='',chrsOfInterest=paste("chr",1:21,sep="")){
#   
#   # Features included
#   featureNickNames <- c("genome-wide","CpG island","CpG shore","CpG shelf","CpG open sea",
#                         "gene body","exon","intron","intergenic")
#   
#   # Define list of feature GRs
#   outGR <- GRangesList()
#   GRtemp <- unlist(tileGenome(seqinfo(Mus),ntile=1))
#   
#   outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
#   #Redefining CpG islands using hidden Markov models 
#   CpG_all <- readRDS(paste(CpGdir,"CpG_hg19.rds",sep=""))
#   CpG_all<-setGenomeLengths(CpG_all)
#   #Could also use UCSC genome browser CpG file
#   cpg_islands <- readRDS(paste(CpGdir,"cpg_islands_hg19.rds",sep=""))
#   cpg_islands<-subsetByOverlaps(CpG_all,cpg_islands)
#   outGR[["CpG island"]] <- setGenomeLengths(cpg_islands)
#   
#   # extract the shore defined by 2000 bp upstream and downstream of cpg islands
#   shore1 <- flank(cpg_islands, 2000)
#   shore2 <- flank(cpg_islands,2000,FALSE)
#   shore1_2 <- reduce(c(shore1,shore2))
#   
#   # extract the features (ranges) that are present in shores only and not in
#   # cpg_islands (ie., shores not overlapping islands)
#   cpgi_shores <- setdiff(shore1_2, cpg_islands)
#   olap=findOverlaps(CpG_all,cpgi_shores)
#   cpgi_shores<-subsetByOverlaps(CpG_all,cpgi_shores)
#   outGR[["CpG shore"]] <- setGenomeLengths(cpgi_shores)
#   
#   # extract the shore defined by 4000 bp upstream and downstream of cpg islands
#   shelves1 <- flank(cpg_islands, 4000)
#   shelves2 <- flank(cpg_islands,4000,FALSE)
#   shelves1_2 <- reduce(c(shelves1,shelves2))
#   
#   # create a set of ranges consisting CpG Islands, Shores
#   island_shores <- c(cpg_islands,cpgi_shores)
#   
#   # extract the features (ranges) that are present in shelves only
#   # and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
#   cpgi_shelves <- setdiff(shelves1_2, island_shores)
#   cpgi_shelves<-subsetByOverlaps(CpG_all,cpgi_shelves)
#   outGR[["CpG shelf"]] <- setGenomeLengths(cpgi_shelves)
#   
#   # Open sea
#   open_sea <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]],outGR[["CpG shelf"]]))
#   open_sea<-subsetByOverlaps(CpG_all,open_sea)
#   outGR[["CpG open sea"]] <- setGenomeLengths(open_sea)
#   
#   # Enhancers 
#   #enhancers <- import.bed(paste(enhancerDir,"enhancers.bed",sep=""))[,c()]
#   
#   #outGR[["enhancer"]] <- setGenomeLengths(enhancers)
#   
#   # Other generic features
#   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#   genes <- genes(txdb)
#   outGR[["gene body"]] <- setGenomeLengths(genes)
#   exons <- exons(txdb)
#   outGR[["exon"]] <- setGenomeLengths(exons[,c()])
#   introns <- intronicParts(txdb)
#   outGR[["intron"]] <- setGenomeLengths(introns[,c()])
#   intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
#   outGR[["intergenic"]] <- setGenomeLengths(intergenic)
#   #Use annotation hub for TSS
#   proms <- promoters(genes,upstream=2000,downstream=2000)
#   outGR[["promoter"]] <- setGenomeLengths(proms)
#   TSS<-promoters(genes,upstream=0,downstream=0)
#   outGR[["TSS"]] <- setGenomeLengths(TSS)
#   # Gene name mapping
#   geneBodyNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["gene body"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
#   outGR[["gene body"]]$gene_name <- geneBodyNameMap$SYMBOL
#   promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["promoter"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
#   outGR[["promoter"]]$gene_name <- promNameMap$SYMBOL
#   promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["TSS"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
#   outGR[["TSS"]]$gene_name <- promNameMap$SYMBOL
#   # Return
#   return(outGR)
#   
# }

# 
# #Get CpG density for each chromosome
# getCpgdensH19 <- function(chrsOfInterest=paste("chr",1:22,sep="")){
#   # Obtain all CpG sites
#   cgs <- lapply(chrsOfInterest, function(x) start(matchPattern("CG", Hsapiens[[x]])))
#   cgs_dist=lapply(cgs,function(x) x[2:length(x)]-x[1:length(x)-1])
#   cgs_df=data.frame(CG_number=unlist(lapply(cgs,length)),
#                     CG_dist=unlist(lapply(cgs_dist,mean)),
#                     total_length=unlist(lapply(chrsOfInterest,function(x) length(Hsapiens[[x]]))))
#   cgs_df$CG_density=cgs_df$CG_number/cgs_df$total_length
#   return(list(cgs_df,cgs_dist))
# }
# 
# 
# 
# # Function to plot a GR using Gviz
# plotGR <- function(CpGdir,enhancerDir,GR,startHight,highSize=500,reverseStrand=FALSE,chr="chr11",lim=c(2010000,2022500)){
#   
#   # Get genome
#   gen <- "hg19"
# 
#   # GRanges to intersect with and keep the relevant data
#   windowGR <- GRanges(seqnames=chr,ranges=IRanges(start=lim[1],end=lim[2]),strand="*")
#   
#   # Subset GR
#   GR <- subsetByOverlaps(GR,windowGR)
# 
#   # Create data tracks
#   dmmlTrack <- DataTrack(GR[GR$Statistic=="dMML",c("Value")],name="dMML")
#   dnmeTrack <- DataTrack(GR[GR$Statistic=="dNME",c("Value")],name="dNME")
#   ucTrack <- DataTrack(GR[GR$Statistic=="UC",c("Value")],name="UC")
#   
#   # Gene track 1
#   bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="grch37.ensembl.org",
#                 path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
#   biomTrack <- BiomartGeneRegionTrack(genome=gen,chromosome=chr,start=lim[1],end=lim[2],name="ENSEMBL",biomart=bm)
# 
#   # Add CpG island annotation track
#   genomicFeatures <- getGeneralFeats(CpGdir,enhancerDir)
#   cpgIslands <- genomicFeatures[["CpG island"]]
#   cpgIslands <- subsetByOverlaps(cpgIslands,windowGR)
#   islandTrack <- AnnotationTrack(cpgIslands,name="CpG islands")
#   
#   # Chromosome information tracks
#   gtrack <- GenomeAxisTrack()
#   itrack <- IdeogramTrack(genome=gen,chromosome=chr)
#   
#   # Highlight
#   ht <- HighlightTrack(trackList=list(dmmlTrack,dnmeTrack,ucTrack),start=startHight,width=highSize,chromosome=chr)
#   
#   # Return plot
#   plotTracks(list(itrack,gtrack,biomTrack,islandTrack,ht),from=lim[1],to=lim[2],
#              transcriptAnnotation="symbol",type=c("gradient"),stacking="squish",reverseStrand=reverseStrand,
#              collapseTranscripts = "meta")
#   
# }
# # Function to assign genome
# assignGenomeFile <- function(row) {
#   
#   genome = NA
#   if((row[1]=="REF")&(row[2] %in% c("0|1","0/1"))){
#     genome = 1
#   } else if((row[1]=="ALT")&(row[2] %in% c("0|1","0/1"))){
#     genome = 2
#   } else if((row[1]=="REF")&(row[2] %in% c("1|0","1/0"))) {
#     genome = 2
#   } else if((row[1]=="ALT")&(row[2] %in% c("1|0","1/0"))){
#     genome = 1
#   } else {
#     # nothing
#   }
#   
#   return(genome)
# }
# 
# # Function that returns bool if variant results in a CpG site
# hetCpgSite <- function(row) {
#   
#   # Check if we have C or G, otherwise return false
#   if(row[3] %in% c("C","G")){
#     # continue
#   } else {
#     return(FALSE)
#   }
#   
#   # Initialize output
#   hetCpg <- FALSE
#   context <- toString(getSeq(Hsapiens,row[1],start=as.numeric(row[2])-1,end=as.numeric(row[2])+1,strand="+"))
#   if((row[3]=="C")&(substr(context,3,3)=="G")){
#     hetCpg <- TRUE
#   } else if((row[3]=="G")&(substr(context,1,1)=="C")){
#     hetCpg <- TRUE
#   } else {
#     # nothing
#   }
#   
#   # Return binary vector
#   return(hetCpg)
#   
# }
# 
# #Find min width, no longer needed for new output
# gff_min<-function(op_more,olap_sub,qt,gffwid){
#   gff_idx=olap_sub[which(qt==op_more)]
#   return(gff_idx[which.min(gffwid[gff_idx])])
# }
# 
# ###End of calculation
# 
# 
# 
# #generate and export bed file
# 
# ASM_bed_gen_sp<-function(gr_ASM,gr_allele,sp,outdir){
#   gr_ASM=gr_ASM[gr_ASM$CpGdiff!=0]
#   gr_ASM=gr_ASM[order(gr_ASM$diff,decreasing = TRUE)]
#   #do it for each sample
#   #sp='Gastric - STL001'
#   sp_ASM=gr_ASM[gr_ASM$Sample==sp]#stats for this subject
#   sp_allele=gr_allele[gr_allele$Sample==sp & gr_allele$HetCpG]
#   sp_bed=granges(sp_ASM)
#   sp_bed$diff=sp_ASM$diff
#   olap1=findOverlaps(sp_allele[sp_allele$Genome==1],sp_bed,type='equal')
#   olap2=findOverlaps(sp_allele[sp_allele$Genome==2],sp_bed,type='equal')
#   sp_bed$A1[subjectHits(olap1)]=sp_allele[sp_allele$Genome==1]$Value[queryHits(olap1)]
#   sp_bed$A2[subjectHits(olap2)]=sp_allele[sp_allele$Genome==2]$Value[queryHits(olap2)]
#   #Export bed file
#   export_bed(sp_bed,'diff',paste(outdir,sp,'_hetASM_diff.bedGraph',sep=''))
#   export_bed(sp_bed,'A1',paste(outdir,sp,'_hetASM_A1.bedGraph',sep=''))
#   export_bed(sp_bed,'A2',paste(outdir,sp,'_hetASM_A2.bedGraph',sep=''))
# }
# export_bed<-function(gr_in,dat,out_name){
#   bed_out=granges(gr_in)
#   bed_out$score=elementMetadata(gr_in)[,dat]
#   export(bed_out,out_name,format='bedGraph')
# }
# 
# 
# 
# 
# #Generate motif from variant file
# 
# #From vcf file, extract het CpG information
# extractmotif<-function(vcfDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
#   cat('Processing subject:', sub,'\n')
#   tt1=proc.time()[[3]]
#   #genomeGr <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
#   #genomeGr <- setGenomeLengths(genomeGr)
#   vcf <- readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19")
#   gt <- as.vector(geno(vcf)$GT)
#   vcf <- rowRanges(vcf)
#   vcf$GT <- gt
#   vcf$snpId <- paste(sub,seq(1:length(vcf)),sep="-")
#   # Keep only relevant variables
#   vcf <- vcf[,c("REF","ALT","GT","snpId")]
#   vcf$REF <- as.character(vcf$REF)
#   vcf$ALT <- as.character(unlist(vcf$ALT))
#   names(vcf)=NULL
#   # Delete labels
#   vcf=motif_df(vcf)
#   vcf$sub=sub
#   cat('Done processing',sub,'using:', proc.time()[[3]]-tt1,'\n')
#   return(vcf)
#   
# }
# #ID motif from vcf file
# motif_df<-function(var_gr){
#   #Get 11 nucleotide location, 3 nucleotide 4-6, 5 nulceotide 3-7, 7 nucleotide 2-8
#   nucleo_11_alt<-nucleo_11_X<-nucleo_11_ref<-as.matrix(getSeq(Hsapiens,resize(var_gr,11,fix='center')))
#   nucleo_11_X[,6]='N'
#   nucleo_11_alt[,6]=var_gr$ALT
#   #Running too long. using matrix operation
#   #cat matrix into 1 string
#   nucleo_11_ref=apply(nucleo_11_ref,1,function(x) paste(x,collapse=''))
#   nucleo_11_X=apply(nucleo_11_X,1,function(x) paste(x,collapse=''))
#   nucleo_11_alt=apply(nucleo_11_alt,1,function(x) paste(x,collapse=''))
#   #Get 3,5,7,9 nucleotide location
#   nucleo_9_ref=substr(nucleo_11_ref,start=2,stop=10)
#   nucleo_7_ref=substr(nucleo_11_ref,start=3,stop=9)
#   nucleo_5_ref=substr(nucleo_11_ref,start=4,stop=8)
#   nucleo_3_ref=substr(nucleo_11_ref,start=5,stop=7)
#   
#   nucleo_9_alt=substr(nucleo_11_alt,start=2,stop=10)
#   nucleo_7_alt=substr(nucleo_11_alt,start=3,stop=9)
#   nucleo_5_alt=substr(nucleo_11_alt,start=4,stop=8)
#   nucleo_3_alt=substr(nucleo_11_alt,start=5,stop=7)
#   
#   nucleo_9_X=substr(nucleo_11_X,start=2,stop=10)
#   nucleo_7_X=substr(nucleo_11_X,start=3,stop=9)
#   nucleo_5_X=substr(nucleo_11_X,start=4,stop=8)
#   nucleo_3_X=substr(nucleo_11_X,start=5,stop=7)
#   
#   #Assign value to var_gr
#   var_df=data.frame(nucleo_11_ref=nucleo_11_ref,
#   nucleo_9_ref=nucleo_9_ref,
#   nucleo_7_ref=nucleo_7_ref,
#   nucleo_5_ref=nucleo_5_ref,
#   nucleo_3_ref= nucleo_3_ref,
#   
#   nucleo_11_alt=nucleo_11_alt,
#   nucleo_9_alt=nucleo_9_alt,
#   nucleo_7_alt=nucleo_7_alt,
#   nucleo_5_alt=nucleo_5_alt,
#   nucleo_3_alt=nucleo_3_alt,
#   
#   nucleo_11_X=nucleo_11_X,
#   nucleo_9_X=nucleo_9_X,
#   nucleo_7_X=nucleo_7_X,
#   nucleo_5_X=nucleo_5_X,
#   nucleo_3_X=nucleo_3_X,stringsAsFactors = FALSE)
#   var_gr=makeGRangesFromDataFrame(cbind(var_gr,var_df),keep.extra.columns = TRUE)
#   return(var_gr)
# }
# 
# #For each sample, add dMML and dNME information
# extract_diff_values<-function(sp,diff,variant){
#   #Get subject information for sp
#   subj= strsplit(sp,' - ')[[1]][2]
#   variant=variant[[subj]]
#   #For this sample, extract the loci
#   #dMML
#   dMML=diff[diff$Statistic=='dMML']
#   olap_dMML=findOverlaps(variant,dMML,type='within')
#   outGR=variant[queryHits(olap_dMML)]
#   outGR$dMML=dMML$Value[subjectHits(olap_dMML)]
#   outGR$dMML_pval=dMML$pvalue[subjectHits(olap_dMML)]
#   #dNME
#   dNME=diff[diff$Statistic=='dNME']
#   olap_dNME=findOverlaps(outGR,dNME,type='within')
#   outGR$dNME=NA
#   outGR$dNME_pval=NA
#   outGR[queryHits(olap_dNME)]$dNME=dNME$Value[subjectHits(olap_dNME)]
#   outGR[queryHits(olap_dNME)]$dNME_pval=dNME$pvalue[subjectHits(olap_dNME)]
#   outGR$Sample=sp
#   return(outGR)
# }
# 
# 
# #Reshape each motif and keep extra column
# reshape_sample_variant<-function(variant){
#   variant=as.data.frame(variant)
#   melt_id_var=colnames(variant)
#   melt_id_var=melt_id_var[-which(melt_id_var=="REF" | melt_id_var=="ALT")]
#  
#   # Melt: origin = ref/alt, variant= ATCG
#   variant <- melt(data=variant,id.vars=melt_id_var,
#                 measure.vars=c("REF","ALT"),
#                 variable.name="Origin",value.name="Variant")
#   #Assign genome
#   variant$Genome <- apply(variant[,c("Origin","GT")],1,assignGenomeFile)
#   #Reshape ref columns
#   ref_id=c('nucleo_11_ref','nucleo_9_ref','nucleo_7_ref','nucleo_5_ref','nucleo_3_ref')
#   alt_id=c('nucleo_11_alt','nucleo_9_alt','nucleo_7_alt','nucleo_5_alt','nucleo_3_alt')
#   new_id=c('nucleo_11','nucleo_9','nucleo_7','nucleo_5','nucleo_3')
#   #reassign ID
#   variant[,new_id]=NA
#   variant[variant$Genome==1,new_id]=variant[variant$Genome==1,ref_id]
#   variant[variant$Genome==2,new_id]=variant[variant$Genome==2,alt_id]
#   variant=variant[,-which(colnames(variant)%in%c(ref_id,alt_id))]
#   #Make Granges
#   variant=makeGRangesFromDataFrame(variant,keep.extra.columns = TRUE)
#   return(variant)
# }
# assignGenomeFile <- function(row) {
#   
#   genome = NA
#   if((row[1]=="REF")&(row[2] %in% c("0|1","0/1"))){
#     genome = 1
#   } else if((row[1]=="ALT")&(row[2] %in% c("0|1","0/1"))){
#     genome = 2
#   } else if((row[1]=="REF")&(row[2] %in% c("1|0","1/0"))) {
#     genome = 2
#   } else if((row[1]=="ALT")&(row[2] %in% c("1|0","1/0"))){
#     genome = 1
#   } else {
#     # nothing
#   }
#   
#   return(genome)
# }
# #Extract allele values for each allele
# extract_allele_value<-function(outGR,cpelAllele){
#   # Cross resulting GR with MML of genome 1
#   cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="1")]
#   olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
#   outGR$MML1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
#   outGR$MML1[outGR$Genome=="2"] <- NA
#   
#   # Cross resulting GR with MML of genome 2
#   cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="2")]
#   olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
#   outGR$MML2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
#   outGR$MML2[outGR$Genome=="1"] <- NA
#   
#   # Consolidate MML1 and MML2 columns into single column
#   outGR$MML <- rowSums(cbind(outGR$MML1,outGR$MML2), na.rm=T)
#   outGR$MML1 <- outGR$MML2 <- NULL
#   
#   # Cross resulting GR with NME of genome 1
#   cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="1")]
#   olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
#   outGR$NME1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
#   outGR$NME1[outGR$Genome=="2"] <- NA
#   
#   # Cross resulting GR with NME of genome 2
#   cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="2")]
#   olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
#   outGR$NME2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
#   outGR$NME2[outGR$Genome=="1"] <- NA
#   
#   # Consolidate NME1 and NME2 columns into single column
#   outGR$NME <- rowSums(cbind(outGR$NME1,outGR$NME2), na.rm=T)
#   outGR$NME1 <- outGR$NME2 <- NULL
# 
#   # Return
#   return(outGR)
# }
# #Enrichment test
# motif_enrichment<-function(GR_allele,pval_cutoff=0.1,p_stat,motif,motif_type){#pstat either dNME_pval or dMML_pval
#   GR_allele=as.data.frame(GR_allele)
#   GR_allele$ASM=NA
#   GR_allele$ASM[GR_allele[,p_stat]<=pval_cutoff] =  TRUE
#   GR_allele$ASM[GR_allele[,p_stat]>pval_cutoff] = FALSE
#   ASM_motif=sum(GR_allele$ASM & GR_allele[,motif_type]==motif)
#   ASM_not_motif=sum(GR_allele$ASM & !GR_allele[,motif_type]==motif)
#   notASM_motif=sum(!GR_allele$ASM & GR_allele[,motif_type]==motif)
#   notASM_notmotif=sum(!GR_allele$ASM & !GR_allele[,motif_type]==motif)
#   cont_table=matrix(c(ASM_motif,ASM_not_motif,notASM_motif,notASM_notmotif),nrow=2,byrow=TRUE)
#   return(fisher.test(cont_table))
# }
# 
# # #for given sample
# # variant_meta_sp<-function(variant_subj,GR_st){
# #   st=unique(GR_st$Statistic)
# #   gr_out_st = GRanges()
# #   for (statistics in st){
# #     gr_out_st=c(gr_out_st,variant_meta_sp_st(variant_subj,GR_st[GR_st$Statistic==statistics]))
# #     
# #   }
# #   return(gr_out_st)
# # }
# 
# #Enrichment of variants
# variant_enrich<-function(variant_in,cutoff=0.1){
#   variant_in$ASM=FALSE
#   variant_in$ASM[variant_in$dNME_pval<=cutoff]=TRUE
#   cont_table=matrix(c(sum(variant_in$variant & variant_in$ASM),
#                       sum(variant_in$variant & !variant_in$ASM),
#                       sum(!variant_in$variant & variant_in$ASM),
#                       sum(!variant_in$variant & !variant_in$ASM)),
#                     nrow=2,byrow=TRUE)
#   colnames(cont_table)=c('ASM','Not ASM')
#   rownames(cont_table)=c('Vairant','Not variant')
#   print(cont_table)
#   return(fisher.test(cont_table))
# }
# #Subset for SNP-containing ranges
# SNP_conmtaining_hap<-function(gr_in,variant_in){
#   SNP_sub=GRanges()
#   SNP_not=GRanges()
#   for (subj in unique(gr_in$Subject)){
#     olap=findOverlaps(gr_in[gr_in$Subject==subj],variant_in[[subj]])
#     SNP_sub=c(SNP_sub,subsetByOverlaps(gr_in[gr_in$Subject==subj], variant_in[[subj]]))
#     SNP_not=c(SNP_not,gr_in[gr_in$Subject==subj][-queryHits(olap)])
#   }
# return(list(SNP_containing=SNP_sub,Non_SNP_containing=SNP_not))
# }
# #ASM_het_CpG_enrichment
# ASM_het_enrich<-function(gr_in,title){
# OR_df=data.frame(sp=NULL,OR=NULL,lower_CI=NULL,upper_CI=NULL)
# for (sp in unique(gr_in$Sample)){
#   OR=ASM_het_enrichment(gr_in[gr_in$Sample==sp])
#   OR_df=rbind(OR_df,data.frame(sp=sp,subjects=strsplit(sp,' - ')[[1]][2],
#             OR=OR$estimate,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
# }
# theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
# ggplot(OR_df,aes(x=sp,y=OR,fill=subjects)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 13)+
#   ggtitle(title)+xlab('Sample name')+ylab('Odds Ratio')+
#   geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
#                 position=position_dodge(.9))+  theme_bar
# 
# }
# 
# #Plot global distribution of density vs CG type
# plot_density<-function(CG_df,ylab,title,ylim,xlab='CG type'){
# ggplot(CG_df,aes(x=CG_type,y=density,fill=CG_type))+#scale_fill_manual(values = c("blue","red"))+
#   geom_boxplot(outlier.shape = NA)+xlab(xlab)+ylab(ylab)+ylim(ylim)+
#   theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
#   ggtitle(title)+theme(legend.title = element_blank())
# }
# #CpG feature enrichment plot
# genome_feature_plot<-function(gr,feature,Stats,title,ylim=c(0,2),pval_cutoff=0.1){
#   groups=unique(gr$Group)
#   gr$ASM=NA
#   gr$ASM[which(elementMetadata(gr)[,paste(Stats,"pval",sep="_")]<=pval_cutoff)]="Yes"
#   gr$ASM[which(elementMetadata(gr)[,paste(Stats,"pval",sep="_")]>pval_cutoff)]="No"
#   OR_df=data.frame(group=groups,OR=0,lower_CI=0,upper_CI=0)
#   for(gp in groups){
#     OR_out=testEnrichmentFeature_stat(gr[gr$Group==gp],feature)[[2]]
#     OR_df$OR[OR_df$group==gp]=OR_out$estimate
#     OR_df$lower_CI[OR_df$group==gp]=OR_out$conf.int[1]
#     OR_df$upper_CI[OR_df$group==gp]=OR_out$conf.int[2]
#   }  
#   print(OR_df)
#   theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
#   ggplot(OR_df,aes(x=group,y=OR,fill=group)) + ylim(ylim)+
#     geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
#     geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
#                   position=position_dodge(.9))+theme_bar
#     
# }
# 
# #calculating distance
# 
# multiple_reduce<-function(qt_idx,hits){
# 
#   hits_red=hits[queryHits(hits)==qt_idx]
# 
#   idx_min=which.min(elementMetadata(hits_red)$distance)
#   return(c(subjectHits(hits_red)[idx_min],elementMetadata(hits_red)$distance[idx_min]))
#   
# }
# compare_dist<-function(x){
#   if(any(is.na(x))){
#     return(x[!is.na(x)])
#   }else{
#     return(x[which.min(abs(x))])
#   }
#   
# }
# compare_dist_gene<-function(x){
#  
#   if(any(is.na(x[1:2]))){
#     return(x[which(!is.na(x[1:2]))+2])
#   }else{
#     return(x[which.min(abs(as.numeric(x[1:2])))+2])
#   }
#   
#   
# }
# #Distance to given granges
# gr_distance<-function(gr_in,gr_feature,xlab,main,ylim){
#   gr_dist=GRanges()
#   for (subj in unique(gr_in$Sample)){
#     gr_dist=c(gr_dist,dist_calc(gr_in[gr_in$Sample == subj],gr_feature))
#               }
#   gr_dist$dist_round=gr_dist$dist_round/1000
#   #gr_all_close=gr_dist[abs(gr_dist$dist_round)<=50]
#   gr_all_close=gr_dist[abs(gr_dist$dist_round)<=15]
#   gr_count=table(gr_all_close$dist_round)
#   gr_plot_df=data.frame(dist=as.numeric(names(gr_count)),percent_ASM=gr_count/length(gr_all_close))
#   plot(gr_plot_df$dist,gr_plot_df$percent_ASM.Freq,pch=1,cex=0.8,ylab='Proportion of ASM',xlab=xlab,main=main,ylim=ylim,xlim=c(-15,15))
#   lines(gr_plot_df$dist,gr_plot_df$percent_ASM.Freq,lwd=1.5)
#   abline(h=mean(gr_plot_df$percent_ASM.Freq),lty=2,lwd=4)
#   return(gr_all_close)
# }
# 
# readEnhancer <- function(enhancerDir){
#   
#   #Function to convert colnames to granges
#   load(paste(enhancerDir,"enhancers_intersect.RData",sep=""))
#   enhancer_gr_all <- do.call('c',lapply(rownames(max_states),rownames2Granges))
#   
#   # Return
#   return(enhancer_gr_all)
# }
# #Find the overlap event in ASM
# olap_ASM<-function(GR_in){
#   dMML=GR_in[GR_in$Statistic=='dMML']
#   dNME=GR_in[GR_in$Statistic=='dNME']
#   UC=GR_in[GR_in$Statistic=='UC']
#   dMML_dNME=length(subsetByOverlaps(dMML,dNME))
#   dMML_UC=length(subsetByOverlaps(dMML,UC))
#   dNME_UC=length(subsetByOverlaps(dMML,dNME))
#   olap3=length(subsetByOverlaps(subsetByOverlaps(UC,dNME),dMML))
#   return(data.frame(dMML_dNME=dMML_dNME-olap3,dMML_UC=dMML_UC-olap3,
#                     dNME_UC=dNME_UC-olap3,olap3=olap3,dMML=length(dMML),
#                     dNME=length(dNME),UC=length(UC),
#                     dMML_nonolap=length(dMML)-dMML_dNME-dMML_UC-olap3,
#                     dNME_nonolap=length(dNME)-dMML_dNME-dNME_UC-olap3,
#                     UC_nonolap=length(UC)-dNME_UC-dMML_UC-olap3,
#                     sample=unique(GR_in$Sample),subject=unique(GR_in$Subject)))
# }
# 
# #Calculate enrichment of each variants in ASM
# variants_OR<-function(variant_gr,variant,statistics,cutoff=0.1){
#   variant_gr$pvalue=elementMetadata(variant_gr)[,paste(statistics,'pval',sep="_")]
#   invariant=variant_gr[variant_gr$variants==variant]
#   nonvariant=variant_gr[variant_gr$variants!=variant]
#   variant_ASM=sum(invariant$pvalue<=cutoff)
#   variant_nonASM=sum(invariant$pvalue>cutoff)
#   nonvariant_ASM=sum(nonvariant$pvalue<=cutoff)
#   nonvariant_nonASM=sum(nonvariant$pvalue>cutoff)
#   cont_table=matrix(c(variant_ASM,variant_nonASM,nonvariant_ASM,nonvariant_nonASM),nrow=2)
#   return(fisher.test(cont_table))
# }
# #calculate odds ratio of trinucleotide
# tri_nucleo_OR<-function(gr,tri,stat_type,cutoff=0.05){
#   gr$pvalue=elementMetadata(gr)[,paste(stat_type,'pval',sep='_')]
#   tri_gr=gr[gr$mask_tri==tri]
#   nontri=gr[gr$mask_tri!=tri]
#   tri_ASM=sum(tri_gr$pvalue<=cutoff)
#   tri_nonASM=sum(tri_gr$pvalue>cutoff)
#   nontri_ASM=sum(nontri$pvalue<=cutoff)
#   nontri_nonASM=sum(nontri$pvalue>cutoff)
#   cont_table=matrix(c(tri_ASM,tri_nonASM,nontri_ASM,nontri_nonASM),byrow = T,nrow=2)
#   return(fisher.test(cont_table))
# }
# #Merge 3 stats for each sample
# 
# 
# density_plot_hyper<-function(GR_in,genes_hypervaribility_in,genomic_features,quant_dat="hyper_var_promoter"){
#   #$hyper_var is used to calculating quantile
#   GR_in=dist_calc(GR_in,genomic_features$TSS,k_round=100)
#   plot_out=data.frame(NME=GR_in$mean_NME,dNME=GR_in$dNME,dist=GR_in$dist,dNME_pval=GR_in$dNME_pval,
#                       MML=GR_in$mean_MML,dMML=GR_in$dMML,dMML_pval=GR_in$dMML_pval,
#                       hyper_var=elementMetadata(GR_in)[,quant_dat])
#   plot_out$quant=findInterval(plot_out$hyper_var,quantile(plot_out$hyper_var,prob=c(0,0.25,0.5,0.75),na.rm=T))
#   quant=c("0-25%","25%-50%","50%-75%","75%-100%")
#   plot_out$quant=quant[plot_out$quant]
#   plot_out=plot_out[!is.na(plot_out$quant),]
#   print(table(plot_out$quant))
#   #within 1k of promoter
#   NME=ggplot(plot_out,aes(x=NME,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
#     ggtitle("NME distribution")
#   dNME=ggplot(plot_out,aes(x=dNME,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
#     ggtitle("dNME distribution")
#   print(grid.arrange(NME,dNME,nrow=1))
#   MML=ggplot(plot_out,aes(x=MML,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
#     ggtitle("MML distribution")
#   dMML=ggplot(plot_out,aes(x=dMML,color=quant))+geom_density(size=1)+theme(legend.position="bottom")+
#     ggtitle("dMML distribution")
#   print(grid.arrange(MML,dMML,nrow=1))
#   dNME_dot=ggplot(plot_out,aes(x=hyper_var,y=dNME))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
#     ggtitle("dNME vs hypervaribility")
#   NME_dot=ggplot(plot_out,aes(x=hyper_var,y=NME))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
#     ggtitle("NME vs hypervaribility")
#   print(grid.arrange(NME_dot,dNME_dot,nrow=1))
#   dMML_dot=ggplot(plot_out,aes(x=hyper_var,y=dMML))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
#     ggtitle("dMML vs hypervaribility")
#   MML_dot=ggplot(plot_out,aes(x=hyper_var,y=MML))+geom_smooth(size=1)+geom_point(alpha=0.05)+theme(legend.position="bottom")+
#     ggtitle("MML vs hypervaribility")
#   print(grid.arrange(MML_dot,dMML_dot,nrow=1))
#   #print(table(NME_plot_out[abs(NME_plot_out$dist)<=1000,"quant"]))
#   dist_NME=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=NME,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("NME distribution vs distance to TSS")
#   dist_dNME=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=dNME,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("dNME distribution vs distance to TSS")
#   print(grid.arrange(dist_NME,dist_dNME,nrow=1))
#   dist_MML=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=MML,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("MML distribution vs distance to TSS")
#   dist_dMML=ggplot(plot_out[abs(plot_out$dist)<=5000,],aes(x=dist,y=dMML,color=quant))+geom_smooth(size=1,na.rm=TRUE,method="loess",se=TRUE)+theme(legend.position="bottom")+ggtitle("dMML distribution vs distance to TSS")
#   print(grid.arrange(dist_MML,dist_dMML,nrow=1))
# }  
# 
# 
# #Look for gene being bind by those TF
# motif_family_enrich<-function(hits,background,family_anno,pse_count=25){
#   family_anno$Family=family_anno$Family
#   bg_non_hit=background[!background %in% hits]
#   df_out=data.frame(OR=NULL,pvalue=NULL,lowerCI=NULL,upperCI=NULL,family=NULL,hits_in=NULL)
#   for(fam in unique(family_anno$Family)){
#     
#     family=family_anno$Name[family_anno$Family==fam]
#     
#     hit_in_family<-sum(hits %in% family)+pse_count
#     hit_not_family<-sum(!hits %in% family)+pse_count
#     bg_in_family<-sum(bg_non_hit %in% family)+pse_count
#     bg_not_family<-sum(!bg_non_hit %in% family)+pse_count
# 
#     cont_table=matrix(c(hit_in_family,hit_not_family,bg_in_family,bg_not_family),nrow=2)
#     fs_pse=fisher.test(cont_table)
#     fs=fisher.test(cont_table-(pse_count-1))
#     df_out=rbind(df_out,data.frame(OR=fs_pse$estimate[[1]],pvalue=fs$p.value,
#                                    lowerCI=fs$conf.int[[1]],upperCI=fs$conf.int[[2]],family=fam,
#                                    hits_in=hit_in_family-pse_count,hit_not_family=hit_not_family-pse_count,
#                                    bg_in_family=bg_in_family-pse_count,bg_not_family=bg_not_family-pse_count,
#                                    variance=fs_pse$conf.int[[2]]-fs_pse$conf.int[[1]]))
#   }
#   df_out=df_out[df_out$hits_in>1,]
#   df_out$qvalue=p.adjust(df_out$pvalue,method='BH')
#   return(df_out[order(df_out$pvalue,decreasing=F),])
# }
# 
# split_ranges<-function(TSS_in,size=350){
#   n=width(TSS_in)
#   n=n-size
#   start_TSS_in=start(TSS_in)
#   strand_in=strand(TSS_in)
#   chr_in=seqnames(TSS_in)
#   #gene_name_TSS=TSS_in$gene_symbol
#   #i=1
#   for (num in seq(size,n,size)){
#     start_TSS_in=c(start_TSS_in,start(TSS_in)+num)
#     strand_in=c(strand_in,strand(TSS_in))
#     chr_in=c(chr_in,seqnames(TSS_in))
#     #gene_name_TSS=c(gene_name_TSS,paste(TSS_in$gene_symbol,i,sep='-'))
#      #i=i+1
#   }
#   
#   TSS_in_break=data.frame(seqnames=chr_in,start=start_TSS_in,end=start_TSS_in+size-1,strand=strand_in)
#   TSS_in_break=makeGRangesFromDataFrame(TSS_in_break[which(TSS_in_break$start>0),],keep.extra.columns = T)
#   TSS_in_break=TSS_in_break[seqnames(TSS_in_break)%in% seqlevels(TSS_in_break)[1:24]]#select only auto+xy
#   return(TSS_in_break)
# }
# #Need to fix it, previous setting =200
# # split_ranges_larger_5k<-function(TSS_in,size=350){
# #   
# #   TSS_out=split_ranges_5k(TSS_in,))
# #   return(TSS_out)
# # }
# 
# make_gr_dnase<-function(x){
#   x=strsplit(x,':')[[1]]
#   range_sp=strsplit(x[2],"-")[[1]]
#   return(makeGRangesFromDataFrame(data.frame(seqnames=x[1],start=range_sp[1],end=range_sp[2])))
# }
# 
# stat_hyper_enrich<-function(GR_merge_sp,stat='NME'){
#   stat_in_quant=quantile(c(elementMetadata(GR_merge)[,paste(stat,'1',sep='')],elementMetadata(GR_merge)[,paste(stat,'2',sep='')]),prob=0.75)
#   GR_merge_sp$hyper_var_genes=GR_merge_sp$hyper_var_promoter>=GR_merge_sp$hyper_var_upper
#   GR_merge_sp=GR_merge_sp[!is.na(GR_merge_sp$hyper_var_genes)]
#   # (elementMetadata(GR_merge_sp)[,paste(stat,'1',sep='')]>=stat_in_quant|
#   #     elementMetadata(GR_merge_sp)[,paste(stat,'2',sep='')]>=stat_in_quant)
#   ASM_hypervar=sum(GR_merge_sp$dNME_pval<=0.1 &GR_merge_sp$hyper_var_genes)
#   ASM_not_hypervar=sum(GR_merge_sp$dNME_pval<=0.1&!GR_merge_sp$hyper_var_genes)
#   
#   notASM_hypervar=sum(GR_merge_sp$dNME_pval>0.1  &GR_merge_sp$hyper_var_genes)
#   notASM_not_hypervar=sum(GR_merge_sp$dNME_pval>0.1 &
#                             !GR_merge_sp$hyper_var_genes)
#   cont_table=matrix(c(ASM_hypervar,ASM_not_hypervar,notASM_hypervar,notASM_not_hypervar),nrow=2)
#  
#   table_CMH=data.frame(subject=factor(tissue),
#                        ASM=factor(c("ASM","ASM","NonASM","NonASM")),
#                        feature=factor(c("Hypervar","NonHypervar","Hypervar","NonHypervar")),
#                        count=c(ASM_hypervar,ASM_not_hypervar,notASM_hypervar,notASM_not_hypervar))
#   print(table_CMH)
# 
#     if(all(table_CMH$count>1)){
#       print(fisher.test(cont_table))
#       return(table_CMH)
#     }
#   
# }
# 
# 
# RNA_df<-function(GR_in,RNA_in){
#   #Checking overlapping with H1
#   cat('Number of genes covered:',sum(GR_in$genes_promoter %in% rownames(RNA_in)),'\n')
#   rna_asm_GR=rownames(RNA_in)[rownames(RNA_in) %in% GR_in$genes_promoter]
#   #Make a summary df for those genes, GR_merge have most dNME
#   rna_asm_hyper=data.table(dNME_promo=GR_in$dNME[GR_in$genes_promoter %in% rna_asm_GR])
#   rna_asm_hyper$dNME_pval=GR_in$dNME_pval[GR_in$genes_promoter %in% rna_asm_GR]
#   rna_asm_hyper$genes=GR_in$genes_promoter[GR_in$genes_promoter %in% rna_asm_GR]
#   rna_asm_hyper$ASE_log2FC=RNA_in$log2FoldChange[match(rna_asm_hyper$genes,rownames(RNA_in))]
#   rna_asm_hyper$dMML_pval=GR_in$dMML_pval[GR_in$genes_promoter %in% rna_asm_GR]
#   rna_asm_hyper$dMML=GR_in$dMML[GR_in$genes_promoter %in% rna_asm_GR]
#   rna_asm_hyper$dMML_relative=(GR_in$MML2-GR_in$MML1)[GR_in$genes_promoter %in% rna_asm_GR]
#   rna_asm_hyper$dNME=GR_in$dNME[GR_in$genes_promoter %in% rna_asm_GR]
#   rna_asm_hyper$dNME_relative=(GR_in$NME2-GR_in$NME1)[GR_in$genes_promoter %in% rna_asm_GR]
#   return(rna_asm_hyper)
# }
# 
# 
# DEVG_analysis<-function(GR_in,tissue,stat_var='hyper_var_TSS'){
#   if(!is.na(tissue)){GR_in=GR_in[GR_in$tissue==tissue]}
#   GR_in$hyper_var_genes=elementMetadata(GR_in)[[stat_var]]>=GR_in$hyper_var_upper
#   NME_quant=quantile(c(GR_in$NME1,GR_in$NME2),prob=0.75)
#   high_NME_two_hyprvarible=GR_in[which(GR_in$NME1>=NME_quant&GR_in$NME2>=NME_quant& GR_in$hyper_var_genes)]
#   high_NME_one_hyprvarible=GR_in[which((GR_in$NME1>=NME_quant|GR_in$NME2>=NME_quant)
#                                        &(GR_in$dNME_pval<=0.01)&
#                                          GR_in$hyper_var_genes)]
#   low_NME_dNME=GR_in[which(GR_in$NME1<NME_quant&GR_in$NME2<NME_quant&
#                              GR_in$dNME_pval<=pval_cutoff& GR_in$hyper_var_genes)]
#   low_dNME_non_dNME=GR_in[which(GR_in$NME1<NME_quant&GR_in$NME2<NME_quant&
#                                   GR_in$dNME_pval>pval_cutoff& GR_in$hyper_var_genes)]
#   if(!is.na(tissue)){
#     return(data.frame(high_NME_two_hyprvarible=length(high_NME_two_hyprvarible),
#                       high_NME_one_hyprvarible=length(high_NME_one_hyprvarible),
#                       low_NME_dNME=length(low_NME_dNME),
#                       low_dNME_non_dNME=length(low_dNME_non_dNME),
#                       tissue=tissue))
#   }else{
#     
#     return(list(high_NME_two_hyprvarible=high_NME_two_hyprvarible,
#                 high_NME_one_hyprvarible=high_NME_one_hyprvarible,
#                 low_NME_dNME=low_NME_dNME,low_dNME_non_dNME=low_dNME_non_dNME))
#   }
# }
# 
# 
# motif_reformat<-function(motif_in){
#   motif_in$TF=gsub('\\(var.2\\)','',motif_in$TF)
#   motif_in$TF=gsub('\\(var.3\\)','',motif_in$TF)
#   motif_in=unlist(strsplit(motif_in$TF,"::"))
#   return(motif_in)
# }
# motif_pref<-function(variant_in,motif_gene,motif_dir_in){
#   olap=findOverlaps(variant_in,motif_gene)
#   variant_in$alleleDiff=NA
#   variant_in$alleleDiff[queryHits(olap)]=motif_gene$alleleDiff[subjectHits(olap)]
#   motif_sig=motif_gene[motif_gene$geneSymbol %in% motif_dir_in$TF[motif_dir_in$qvalue<=0.1]]
#   
#   variant_in=subsetByOverlaps(variant_in,motif_sig)
#   #variant_in=variant_in[variant_in$dNME_pval<=pval_cutoff]
#   variant_in=variant_in[sign(variant_in$NME2-variant_in$NME1) == sign(variant_in$alleleDiff)]
#   return(variant_in)
# }

# trait_variant=function(x,traits) {
#   variant=NA
#   
#   variant=tryCatch(get_variants(efo_id = x),error=function(e) NA)
#   if(!is.na(variant)){
#     if(length(variant@variants$variant_id)>0){
#       cat('Processing',traits$trait[traits$efo_id==x],'\n')
#       
#       var_out=as.data.frame(variant@variants)
#       var_out=var_out[!is.na(var_out$chromosome_position)&!is.na(var_out$chromosome_name),]
#       if(length(var_out$chromosome_position)>0&length(var_out$chromosome_name)>0&length(var_out$chromosome_name)==length(var_out$chromosome_position)){
#         var_gr=makeGRangesFromDataFrame(var_out,
#                                         seqnames.field = 'chromosome_name',start.field = 'chromosome_position',end.field = 'chromosome_position')
#         seqlevels(var_gr)=paste('chr',seqlevels(var_gr),sep='')
#         var_gr$efo_id=traits$trait[traits$efo_id== x]
#         var_gr$rsid=variant@variants$variant_id
#         #check genomic version
#         snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
#         var_out=snpsById(snps,variant_id[1])
#         if(length(subsetByOverlaps(var_gr[1],var_out))==0){
#           
#         }
#     return(var_out)
#       }
#     }
#   }
# }
# 
# agnostic_matrix_conversion<-function(gr_in,stat='NME'){
#   gr_out=granges(unique(gr_in))
#   olap=findOverlaps(gr_in,gr_out,type='equal')
#   stat_in_df=elementMetadata(gr_in[queryHits(olap)])[c(stat,'Sample')]
#   stat_in_df$idx=NA
#   stat_in_df$idx[queryHits(olap)]=subjectHits(olap)
#   stat_in_df=as.data.table(stat_in_df)
#  
#   stat_in_df_stat=dcast.data.table(data=stat_in_df,formula=idx~Sample,value.var = stat,fun.aggregate=mean)#remove agg.fun for new run
#   gr_out=gr_out[stat_in_df_stat$idx]
#   mcols(gr_out)=stat_in_df_stat[,-1]
#   return(gr_out)
#  
# }
# PCA_df_prep<-function(UC_in_matrix_sub){
#   UC_in_matrix_sub_df=as.data.frame(t(as.matrix(mcols(UC_in_matrix_sub))))
#   UC_in_matrix_sub_PCA=prcomp(UC_in_matrix_sub_df, scale. = TRUE)
#   UC_in_matrix_sub_df$sample=unlist(lapply(strsplit(rownames(UC_in_matrix_sub_df),'-'),function(x) paste(x[-length(x)],collapse = '-')))
#  
#   return(list(UC_in_matrix_sub_PCA,UC_in_matrix_sub_df))
# }
# 
# col_fun <- colorRampPalette(
#   c(
#     ##rev(brewer.pal(8,"RdYlBu"))[1:4],
#     rev(brewer.pal(8,"RdYlBu"))[1:2],
#     #"white",
#     rev(brewer.pal(8,"RdYlBu"))[5:8]
#   )
# )
# 
# #UC
# read.agnostic.mouse.uc<-function(file_in,matrix=FALSE,fileter_N=1,gff_in=NA){
#   
#   cat('processing:',file_in,'\n')
#   informME_in=import.bedGraph(file_in)
#   if(length(informME_in)>0){
#     colnames(elementMetadata(informME_in))=c('score','N','K')
#     if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
#     #fit  bedGraph reads, import.bedGraph will remove 1 from start
#     start(informME_in)=start(informME_in)-1
#     #process file name
#     file_in=strsplit(file_in,'\\/')[[1]]
#     file_in=file_in[length(file_in)]
#     comp= strsplit(strsplit(file_in,'_uc.bedGraph')[[1]],'-vs-')[[1]]
#     strain=unlist(lapply(strsplit(comp,'_'),function(x) x[1]))
#     #if  contain BL6DBA, use ref is BL6DBA
#     strain=ifelse('BL6DBA'%in%strain,'BL6DBA','mm10')
#     comp=unlist(lapply(strsplit(comp,'_'),function(x) paste(x[-1],collapse = '_')))
#     comp=comp[comp!='']
#     comp_stage=unlist(lapply(comp,function(x) {x_split=strsplit(x,'_')[[1]]
#     x_split=x_split[-length(x_split)][-1]
#     x_split=paste(x_split,collapse = '_')
#     return(x_split)}))
#     tissue1=strsplit(comp[1],'_')[[1]][1]
#     tissue2=strsplit(comp[2],'_')[[1]][1]
#     #if BL6DBA, the 1st comp_stage is empty
#     
#     comp_stage=gsub('_5','.5',comp_stage)
#     comp_stage=gsub('day','E',comp_stage)
#     comp_stage=gsub('E0','P0',comp_stage)
#     replicate=strsplit(comp[1],'_')[[1]][length(strsplit(comp[1],'_')[[1]])]
#     replicate=gsub('merged','',replicate)
#     informME_in$Sample=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
#     informME_in=informME_in[informME_in$N>=fileter_N]
#     cat('Minimum N:',min(informME_in$N),'\n')
#     #informME_in$Ref=strain
#     if(matrix){
#       informME_in_dt=as.data.table(mcols(informME_in))[,c("score","Sample")]
#       colnames(informME_in_dt)=c("UC","Sample")
#       informME_in_dt$UC=as.numeric(informME_in_dt$UC)
#       informME_in_dt$region=paste0(seqnames(informME_in),":",start(informME_in),"-",end(informME_in))
#       informME_in_dt=informME_in_dt[match(gff_in,region),"UC"]
#       colnames(informME_in_dt)=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
#       return(informME_in_dt)
#     }
#     else{return(informME_in)}
#   }
# }
# 
# #Get CpG sites from hg19

trait_overlap<-function(trait,NME_in_gr_mean_ts,variant_trait,high_NME){
  
  NME_trait=subsetByOverlaps(NME_in_gr_mean_ts,variant_trait[variant_trait$`DISEASE/TRAIT`==trait])
  DNase_high=sum(NME_trait$NME_mean>=high_NME&NME_trait$DNase=="DNase")
  DNase_total=sum(NME_trait$DNase=="DNase")
  # control_high=sum(NME_trait$NME_mean>=high_NME&NME_trait$DNase=="control")
  # control_total=sum(NME_trait$DNase=="control")
  DNase_mean=mean(NME_trait[NME_trait$DNase=="DNase"]$NME_mean)
  #control_mean_quantile=mean(NME_ecdf(NME_trait[NME_trait$DNase=="control"]$NME_mean))
  return(data.table(trait=trait,high_NME=DNase_high,total=DNase_total,
                    #control_high=control_high,control_total=control_total,
                    #DNase_mean_quantile=DNase_mean_quantile,
                    DNase_mean=DNase_mean
  ))
  #control_mean_quantile=control_mean_quantile))
  
  
}