# Function to set genome and chromosome lengths to a GR object
setGenomeLengths <- function(GR,chrsOfInterest=paste("chr",1:22,sep="")){
  # Get genome info
  
  genome.seqinfo <- Seqinfo(genome="hg19")
  genome.seqinfo <- genome.seqinfo[chrsOfInterest]
  GR <- GR[seqnames(GR) %in% chrsOfInterest]
  genome(GR) <- genome(genome.seqinfo)
  seqlevels(GR) <- seqlevels(genome.seqinfo)
  seqlengths(GR) <- seqlengths(genome.seqinfo)
  
  return(GR)
}

#Function to read in single GR object:
read.diffGR<-function(subjects,tissues,gff_in,inDir,cutoff=0.05,chrsOfInterest=paste("chr",1:22,sep="")){
  #Initialization
  GRs=GRanges()
  # dmml
  filename=paste(inDir,subjects,"_",tissues,"_phased_dmml_pvals.bedGraph",sep="")
  
  GR.in=import.ASMbed(subjects,tissues,gff_in,filename,pvalue=TRUE,'dMML')
  GRs <- append(GRs,GR.in)
  # dnme
  filename=paste(inDir,subjects,"_",tissues,"_phased_dnme_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,gff_in,filename,pvalue=TRUE,'dNME')
  GRs <- append(GRs,GR.in)
  # uc
  filename=paste(inDir,subjects,"_",tissues,"_phased_uc_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,gff_in,filename,pvalue=TRUE,'UC')
  GRs <- append(GRs,GR.in)
  #Check ASM
  GRs <- GRs[!is.na(GRs$pvalue)]
  GRs$ASM <- NA
  GRs[GRs$pvalue<=cutoff]$ASM <- "Yes"
  GRs[GRs$pvalue>cutoff]$ASM <- "No"
  # Add sample field
  GRs$Sample <- paste(tissues,"-",GRs$Subject)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  return(GRs)
}
#file_ends can be 
#c('dmml_pvals,dnme_pvals,uc_pvals')
#c('mml1,mml2,nme1,nme2')
read.alleleGR<-function(subjects,tissues,gff_in,inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  #MML1
  GRs=GRanges()
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,gff_in,filename,pvalue=FALSE,'MML')
  GR.in$Genome="1"
  GRs <- append(GRs, GR.in)
  #MML2
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,gff_in,filename,pvalue=FALSE,'MML')
  GR.in$Genome="2"
  GRs <- append(GRs, GR.in)
  #NME1
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,gff_in,filename,pvalue=FALSE,'NME')
  GR.in$Genome="1"
  GRs <- append(GRs, GR.in)
  #NME2
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,gff_in,filename,pvalue=FALSE,'NME')
  GR.in$Genome="2"
  GRs <- append(GRs, GR.in)
  # Add sample field
  GRs$Sample <- paste(tissues,"-",subjects)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  GRs$K=GRs$'NA.1'
  GRs$'NA.1'=NULL
  return(GRs)
}
#Read in each bed file
import.ASMbed<-function(subjects,tissues,gff_in,filename,pvalue=TRUE,Statistic,chrsOfInterest=paste("chr",1:22,sep="")){
  GR <- import.bedGraph(filename)
  if (!is.na(gff_in)){
    gff_in=gff_in[gff_in$Subject==subjects]
    #resize to gff ranges
    olap_gff=findOverlaps(GR,gff_in,type='within')
    GR_out=gff_in[subjectHits(olap_gff)]
    olap_gff=findOverlaps(GR,GR_out,type='within')
  }else{
    GR_out=GRanges(GR)
    olap_gff=findOverlaps(GR,GR_out)
    }
  GR_out$ASM=NULL
  GR_out$Data=NULL
  GR_out$Subject <- subjects
  GR_out$Statistic <- Statistic
  GR_out$Value[subjectHits(olap_gff)]<- GR$score[queryHits(olap_gff)]
  if(pvalue){GR_out$pvalue[subjectHits(olap_gff)] <- as.numeric(GR[queryHits(olap_gff)]$NA.)}
  else{GR_out$N[subjectHits(olap_gff)] <- GR$NA.1[queryHits(olap_gff)]} #not diff
  GR_out <- GR_out[!duplicated(GR[,c()])]
  GR_out <- resize(GR_out, width(GR_out) + 1, fix="end")
  GR_out <- resize(GR_out, width(GR_out) + 1, fix="start")
  GR_out <- setGenomeLengths(GR_out,chrsOfInterest=chrsOfInterest)
  return(GR_out)
}
#read in all sample tissue diff
import.subject<-function(inDir,gff_in,calc='diff'){
  #for calc: diff -> dMML etc, allele -> NME etc
  # H1
  H1_subject <- rep("H1",3)
  H1_subject_labels <- rep("CL3",3)
  H1_tissues <- c(
    "rep1","rep2","merged"
  )             
  H1_tissue_labels <- c(
    "Embryonic Stem Cell", "Embryonic Stem Cell", "Embryonic Stem Cell"
  )
  H1_gtex_labels <- c(
    "", "Embryonic Stem Cell", "Embryonic Stem Cell"
  )
  H1_diff_labels <- c(
    "Undifferentiated","Undifferentiated","Undifferentiated"
  )
  # GM12878
  GM12878_subject <- rep("GM12878",3)
  GM12878_subject_labels <- rep("CL4",3)
  GM12878_tissues <- c(
    "1","2","merged"
  )             
  GM12878_tissue_labels <- c(
    "EBV B-cell","EBV B-cell","EBV B-cell"
  )
  GM12878_gtex_labels <- c(
    "","",""
  )
  GM12878_diff_labels <- c(
    "","",""
  )
  # H9
  H9_subject <- rep("H9",1)
  H9_subject_labels <- rep("CL1",1)
  H9_tissues <- c(
    "42_embryonic_stem_cell_single"
  )             
  H9_tissue_labels <- c(
    "Embryonic Stem Cell"
  )
  H9_gtex_labels <- c(
    ""
  )
  H9_diff_labels <- c(
    "Undifferentiated"
  )
  
  # HUES64
  HUES64_subject <- rep("HUES64",4)
  HUES64_subject_labels <- rep("CL2",4)
  HUES64_tissues <- c(
    "stem_27_undifferentiated_paired",
    "ectoderm_paired",
    "endoerm_27_paired",
    "mesoderm_23_paired"
  )
  HUES64_tissue_labels <- c(
    "Embyonic Stem Cell",
    "Ectoderm",
    "Endoderm",
    "Mesoderm"
  )
  HUES64_gtex_labels <- c("","","","")
  HUES64_diff_labels <- c(
    "Undifferentiated",
    "Semidifferentiated",
    "Semidifferentiated",
    "Semidifferentiated"
  )
  
  # skin03
  skin03_subject <- rep("skin03",2)
  skin03_subject_labels <- rep("D11",2)
  skin03_tissues <- c(
    "foreskin_keratinocyte_paired",
    "foreskin_melanocyte_paired"
  )
  skin03_tissue_labels <- c(
    "Foreskin Keratinocyte",
    "Foreskin Melanocyte"
  )
  skin03_gtex_labels <- c("","")
  skin03_diff_labels <- rep("Differentiated",length(skin03_subject))
  
  # STL001
  stl001_subject <- rep("STL001",11)
  stl001_subject_labels <- rep("D5",11)
  stl001_tissues <- c(
    "Adipose_single",
    "Small_Intestine_single",
    "Bladder_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Lung_single",
    "Psoas_Muscle_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Spleen_single",
    "Thymus_single"
  )
  stl001_tissue_labels <- c(
    "Adipose",
    "Small Intestine",
    "Bladder",
    "Gastric",
    "Left Ventricle",
    "Lung",
    "Psoas Muscle",
    "Right Ventricle",
    "Sigmoid Colon",
    "Spleen",
    "Thymus"
  )
  stl001_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "",
    "",
    "",
    "Heart_Left_Ventricle",
    "Lung",
    "",
    "",
    "Colon_Transverse",
    "",
    ""
  )
  stl001_diff_labels <- rep("Differentiated",length(stl001_subject))
  
  # STL002
  stl002_subject <- rep("STL002",11)
  stl002_subject_labels <- rep("D6",11)
  stl002_tissues <- c(
    "Adipose_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Lung_single",
    "Ovary_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Small_Intestine_single",
    "Spleen_single"
  )
  stl002_tissue_labels <- c(
    "Adipose",
    "Adrenal Gland",
    "Aorta",
    "Esophagus",
    "Gastric",
    "Lung",
    "Ovary",
    "Pancreas",
    "Psoas Muscle",
    "Small Intestine",
    "Spleen"
  )
  stl002_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Esophagus_Mucosa",
    "",
    "Lung",
    "Uterus",
    "Pancreas",
    "",
    "",
    ""
  )
  stl002_diff_labels <- rep("Differentiated",length(stl002_subject))
  
  # STL003
  stl003_subject <- rep("STL003",13)
  stl003_subject_labels <- rep("D7",13)
  stl003_tissues <- c(
    "Adipose_Tissue_single",
    "Adrenal_Gland_single",
    "Aorta_single",
    "Esophagus_single",
    "Gastric_single",
    "Left_Ventricle_single",
    "Pancreas_single",
    "Psoas_Muscle_single",
    "Right_Atrium_single",
    "Right_Ventricle_single",
    "Sigmoid_Colon_single",
    "Small_Intestine_single",
    "Spleen_single"
  )
  stl003_tissue_labels <- c(
    "Adipose",
    "Adrenal Gland",
    "Aorta",
    "Esophagus",
    "Gastric",
    "Left Ventricle",
    "Pancreas",
    "Psoas Muscle",
    "Right Atrium",
    "Right Ventricle",
    "Sigmoid Colon",
    "Small Intestine",
    "Spleen"
  )
  stl003_gtex_labels <- c(
    "Adipose_Subcutaneous",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Esophagus",
    "",
    "Heart_Left_Ventricle",
    "Pancreas",
    "",
    "",
    "",
    "Colon_Transverse",
    "",
    ""
  )
  stl003_diff_labels <- rep("Differentiated",length(stl003_subject))
  
  # STL011
  stl011_subject <- rep("STL011",1)
  stl011_subject_labels <- rep("D8",1)
  stl011_tissues <- c(
    "Liver_single"
  )
  stl011_tissue_labels <- c(
    "Liver"
  )
  stl011_gtex_labels <- c(
    "Liver"
  )
  stl011_diff_labels <- rep("Differentiated",length(stl011_subject))
  
  # Create single vectors
  subjects <- c(H9_subject,HUES64_subject,skin03_subject,stl001_subject,
                stl002_subject,stl003_subject,stl011_subject,H1_subject,GM12878_subject)
  tissues <- c(H9_tissues,HUES64_tissues,skin03_tissues,stl001_tissues,
               stl002_tissues,stl003_tissues,stl011_tissues,H1_tissues,GM12878_tissues)
  subject_labels <- c(H9_subject_labels,HUES64_subject_labels,skin03_subject_labels,
                      stl001_subject_labels,stl002_subject_labels,stl003_subject_labels,
                      stl011_subject_labels,H1_subject_labels,GM12878_subject_labels)
  tissue_labels <- c(H9_tissue_labels,HUES64_tissue_labels,skin03_tissue_labels,
                     stl001_tissue_labels,stl002_tissue_labels,stl003_tissue_labels,
                     stl011_tissue_labels,H1_tissue_labels,GM12878_tissue_labels)
  gtex_labels <- c(H9_gtex_labels,HUES64_gtex_labels,skin03_gtex_labels,
                   stl001_gtex_labels,stl002_gtex_labels,stl003_gtex_labels,
                   stl011_gtex_labels,H1_gtex_labels,GM12878_gtex_labels)
  diff_labels <- c(H9_diff_labels,HUES64_diff_labels,skin03_diff_labels,stl001_diff_labels,
                   stl002_diff_labels,stl003_diff_labels,stl011_diff_labels,H1_diff_labels,GM12878_diff_labels)
  GRs=GRanges()
  for (i in 1:length(subjects)) {
    # Print sample being loaded
    print(paste("Loading sample:",subjects[i],tissues[i]))
    if (calc=='diff'){
      GR.in=read.diffGR(subjects[i],tissues[i],gff_in,inDir,cutoff=0.05)
    }else if(calc=='allele'){
      GR.in=read.alleleGR(subjects[i],tissues[i],gff_in,inDir)
    }else {cat('Wrong calc \n')}
    GR.in$SubjectLabel <- subject_labels[i]
    GR.in$Tissue <- tissue_labels[i]
    GR.in$GTEx <- gtex_labels[i]
    GR.in$State <- diff_labels[i]
    GRs=append(GRs,GR.in)
  }
  return(GRs)
  
}
testEnrichmentFeature_stat<-function(dataGR,featureGR,statistic){
  # Find ranges overlapping with feature
  olaps <- findOverlaps(dataGR,featureGR,type="any",select="all",maxgap = 200)

  indFeature <- queryHits(olaps)
  featureData <- dataGR[indFeature]
  complementaryData <- dataGR[-indFeature]
  
  # Enrichment of in feature
  featurestatistic <- featureData[featureData$Statistic==statistic]
  complementarystatistic <- complementaryData[complementaryData$Statistic==statistic]
  contTablestatistic <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTablestatistic) <- c("Feature","Complementary")
  contTablestatistic[1,]$ASM <- sum(featurestatistic$ASM=="Yes")
  contTablestatistic[1,]$nonASM <- sum(featurestatistic$ASM=="No")
  contTablestatistic[2,]$ASM <- sum(complementarystatistic$ASM=="Yes")
  contTablestatistic[2,]$nonASM <- sum(complementarystatistic$ASM=="No")
  print(contTablestatistic)
   fisher.test(contTablestatistic)
  
}

ASM_het_enrichment<-function(gr_in){
  ASM_het=sum(gr_in$ASM=='Yes' & gr_in$HetCpG)
  ASM_non_het=sum(gr_in$ASM=='Yes' & !gr_in$HetCpG)
  non_ASM_non_het=sum(gr_in$ASM=='No' & !gr_in$HetCpG)
  non_ASM_het=sum(gr_in$ASM=='No' & gr_in$HetCpG)
  cont_table=matrix(c(ASM_het,ASM_non_het,non_ASM_het,non_ASM_non_het),nrow=2)
  print(cont_table)
  fisher.test(cont_table)
}
testEnrichmentFeature <- function(dataGR,featureGR){
  
  # Find ranges overlapping with feature
  olaps <- findOverlaps(dataGR,featureGR,type="any",select="all")
  indFeature <- queryHits(olaps)
  featureData <- dataGR[indFeature]
  complementaryData <- dataGR[-indFeature]
  
  # Enrichment of dMML-HASM in feature
  featureDmml <- featureData[featureData$Statistic=="dMML"]
  complementaryDmml <- complementaryData[complementaryData$Statistic=="dMML"]
  contTableDmml <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableDmml) <- c("Feature","Complementary")
  contTableDmml[1,]$ASM <- sum(featureDmml$ASM=="Yes")
  contTableDmml[1,]$nonASM <- sum(featureDmml$ASM=="No")
  contTableDmml[2,]$ASM <- sum(complementaryDmml$ASM=="Yes")
  contTableDmml[2,]$nonASM <- sum(complementaryDmml$ASM=="No")
  dmmlFisher <- fisher.test(contTableDmml)
  
  # Enrichment of dNME-HASM in feature
  featureDnme <- featureData[featureData$Statistic=="dNME"]
  complementaryDnme <- complementaryData[complementaryData$Statistic=="dNME"]
  contTableDnme <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableDnme) <- c("Feature","Complementary")
  contTableDnme[1,]$ASM <- sum(featureDnme$ASM=="Yes")
  contTableDnme[1,]$nonASM <- sum(featureDnme$ASM=="No")
  contTableDnme[2,]$ASM <- sum(complementaryDnme$ASM=="Yes")
  contTableDnme[2,]$nonASM <- sum(complementaryDnme$ASM=="No")
  dnmeFisher <- fisher.test(contTableDnme)
  
  # Enrichment of UC-HASM in feature
  featureUc <- featureData[featureData$Statistic=="UC"]
  complementaryUc <- complementaryData[complementaryData$Statistic=="UC"]
  contTableUc <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTableUc) <- c("Feature","Complementary")
  contTableUc[1,]$ASM <- sum(featureUc$ASM=="Yes")
  contTableUc[1,]$nonASM <- sum(featureUc$ASM=="No")
  contTableUc[2,]$ASM <- sum(complementaryUc$ASM=="Yes")
  contTableUc[2,]$nonASM <- sum(complementaryUc$ASM=="No")
  ucFisher <- fisher.test(contTableUc)
  
  # Return list of Fisher's test
  return(list(dmmlFisher,dnmeFisher,ucFisher))
  
}
#Get features
getGeneralFeats_CpG <- function(CpGdir,enhancerDir='',chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Features included
  featureNickNames <- c("genome-wide","CpG island","CpG shore","CpG shelf","CpG open sea",
                        "gene body","exon","intron","intergenic")
  
  # Define list of feature GRs
  outGR <- GRangesList()
  GRtemp <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  
  outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
  #Redefining CpG islands using hidden Markov models 
  CpG_all <- readRDS(paste(CpGdir,"CpG_hg19.rds",sep=""))
  CpG_all<-setGenomeLengths(CpG_all)
  #Could also use UCSC genome browser CpG file
  cpg_islands <- readRDS(paste(CpGdir,"cpg_islands_hg19.rds",sep=""))
  cpg_islands<-subsetByOverlaps(CpG_all,cpg_islands)
  outGR[["CpG island"]] <- setGenomeLengths(cpg_islands)
  
  # extract the shore defined by 2000 bp upstream and downstream of cpg islands
  shore1 <- flank(cpg_islands, 2000)
  shore2 <- flank(cpg_islands,2000,FALSE)
  shore1_2 <- reduce(c(shore1,shore2))
  
  # extract the features (ranges) that are present in shores only and not in
  # cpg_islands (ie., shores not overlapping islands)
  cpgi_shores <- setdiff(shore1_2, cpg_islands)
  olap=findOverlaps(CpG_all,cpgi_shores)
  cpgi_shores<-subsetByOverlaps(CpG_all,cpgi_shores)
  outGR[["CpG shore"]] <- setGenomeLengths(cpgi_shores)
  
  # extract the shore defined by 4000 bp upstream and downstream of cpg islands
  shelves1 <- flank(cpg_islands, 4000)
  shelves2 <- flank(cpg_islands,4000,FALSE)
  shelves1_2 <- reduce(c(shelves1,shelves2))
  
  # create a set of ranges consisting CpG Islands, Shores
  island_shores <- c(cpg_islands,cpgi_shores)
  
  # extract the features (ranges) that are present in shelves only
  # and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
  cpgi_shelves <- setdiff(shelves1_2, island_shores)
  cpgi_shelves<-subsetByOverlaps(CpG_all,cpgi_shelves)
  outGR[["CpG shelf"]] <- setGenomeLengths(cpgi_shelves)
  
  # Open sea
  open_sea <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]],outGR[["CpG shelf"]]))
  open_sea<-subsetByOverlaps(CpG_all,open_sea)
  outGR[["CpG open sea"]] <- setGenomeLengths(open_sea)
  
  # Enhancers 
  #enhancers <- import.bed(paste(enhancerDir,"enhancers.bed",sep=""))[,c()]
  
  #outGR[["enhancer"]] <- setGenomeLengths(enhancers)
  
  # Other generic features
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- genes(txdb)
  outGR[["gene body"]] <- setGenomeLengths(genes)
  exons <- exons(txdb)
  outGR[["exon"]] <- setGenomeLengths(exons[,c()])
  introns <- intronicParts(txdb)
  outGR[["intron"]] <- setGenomeLengths(introns[,c()])
  intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
  outGR[["intergenic"]] <- setGenomeLengths(intergenic)
  #Use annotation hub for TSS
  proms <- promoters(genes,upstream=2000,downstream=2000)
  outGR[["promoter"]] <- setGenomeLengths(proms)
  TSS<-promoters(genes,upstream=0,downstream=0)
  outGR[["TSS"]] <- setGenomeLengths(TSS)
  # Gene name mapping
  geneBodyNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["gene body"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["gene body"]]$gene_name <- geneBodyNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["promoter"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["promoter"]]$gene_name <- promNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["TSS"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["TSS"]]$gene_name <- promNameMap$SYMBOL
  # Return
  return(outGR)
  
}
getGeneralFeats_mm9 <- function(CpGdir,enhancerDir='',chrsOfInterest=paste("chr",1:21,sep="")){
  
  # Features included
  featureNickNames <- c("genome-wide","CpG island","CpG shore","CpG shelf","CpG open sea",
                        "gene body","exon","intron","intergenic")
  
  # Define list of feature GRs
  outGR <- GRangesList()
  GRtemp <- unlist(tileGenome(seqinfo(Mus),ntile=1))
  
  outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
  #Redefining CpG islands using hidden Markov models 
  CpG_all <- readRDS(paste(CpGdir,"CpG_hg19.rds",sep=""))
  CpG_all<-setGenomeLengths(CpG_all)
  #Could also use UCSC genome browser CpG file
  cpg_islands <- readRDS(paste(CpGdir,"cpg_islands_hg19.rds",sep=""))
  cpg_islands<-subsetByOverlaps(CpG_all,cpg_islands)
  outGR[["CpG island"]] <- setGenomeLengths(cpg_islands)
  
  # extract the shore defined by 2000 bp upstream and downstream of cpg islands
  shore1 <- flank(cpg_islands, 2000)
  shore2 <- flank(cpg_islands,2000,FALSE)
  shore1_2 <- reduce(c(shore1,shore2))
  
  # extract the features (ranges) that are present in shores only and not in
  # cpg_islands (ie., shores not overlapping islands)
  cpgi_shores <- setdiff(shore1_2, cpg_islands)
  olap=findOverlaps(CpG_all,cpgi_shores)
  cpgi_shores<-subsetByOverlaps(CpG_all,cpgi_shores)
  outGR[["CpG shore"]] <- setGenomeLengths(cpgi_shores)
  
  # extract the shore defined by 4000 bp upstream and downstream of cpg islands
  shelves1 <- flank(cpg_islands, 4000)
  shelves2 <- flank(cpg_islands,4000,FALSE)
  shelves1_2 <- reduce(c(shelves1,shelves2))
  
  # create a set of ranges consisting CpG Islands, Shores
  island_shores <- c(cpg_islands,cpgi_shores)
  
  # extract the features (ranges) that are present in shelves only
  # and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
  cpgi_shelves <- setdiff(shelves1_2, island_shores)
  cpgi_shelves<-subsetByOverlaps(CpG_all,cpgi_shelves)
  outGR[["CpG shelf"]] <- setGenomeLengths(cpgi_shelves)
  
  # Open sea
  open_sea <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]],outGR[["CpG shelf"]]))
  open_sea<-subsetByOverlaps(CpG_all,open_sea)
  outGR[["CpG open sea"]] <- setGenomeLengths(open_sea)
  
  # Enhancers 
  #enhancers <- import.bed(paste(enhancerDir,"enhancers.bed",sep=""))[,c()]
  
  #outGR[["enhancer"]] <- setGenomeLengths(enhancers)
  
  # Other generic features
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- genes(txdb)
  outGR[["gene body"]] <- setGenomeLengths(genes)
  exons <- exons(txdb)
  outGR[["exon"]] <- setGenomeLengths(exons[,c()])
  introns <- intronicParts(txdb)
  outGR[["intron"]] <- setGenomeLengths(introns[,c()])
  intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
  outGR[["intergenic"]] <- setGenomeLengths(intergenic)
  #Use annotation hub for TSS
  proms <- promoters(genes,upstream=2000,downstream=2000)
  outGR[["promoter"]] <- setGenomeLengths(proms)
  TSS<-promoters(genes,upstream=0,downstream=0)
  outGR[["TSS"]] <- setGenomeLengths(TSS)
  # Gene name mapping
  geneBodyNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["gene body"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["gene body"]]$gene_name <- geneBodyNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["promoter"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["promoter"]]$gene_name <- promNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["TSS"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["TSS"]]$gene_name <- promNameMap$SYMBOL
  # Return
  return(outGR)
  
}
#Get CpG sites from hg19
getCpgSitesH19 <- function(chrsOfInterest=paste("chr",1:22,sep="")){
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cpgr <- do.call(c,lapply(1:length(chrsOfInterest), function(x) GRanges(names(Hsapiens)[x],IRanges(cgs[[x]],width=2))))
  # Set genome and seqlengths
  cpgr <- setGenomeLengths(cpgr)
  # Return CpG site GR
  return(cpgr)
}

#Get CpG density for each chromosome
getCpgdensH19 <- function(chrsOfInterest=paste("chr",1:22,sep="")){
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cgs_dist=lapply(cgs,function(x) x[2:length(x)]-x[1:length(x)-1])
  cgs_df=data.frame(CG_number=unlist(lapply(cgs,length)),
                    CG_dist=unlist(lapply(cgs_dist,mean)),
                    total_length=unlist(lapply(chrsOfInterest,function(x) length(Hsapiens[[x]]))))
  cgs_df$CG_density=cgs_df$CG_number/cgs_df$total_length
  return(list(cgs_df,cgs_dist))
}
#Read in gff file
readAllGff <- function(inDir,subjects,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Loop over all GFF files
  outGR <- GRanges()
  #subjects <- c("H9","HUES64","skin03","HuFGM02","STL001","STL002","STL003","STL011")
  for (sub in subjects) {
    # Import GFF file
    cat('importing',sub,'\n')
    tmpGR <- readSubGff(inDir,sub,chrsOfInterest)
    # Append to output GR
    outGR <- append(outGR,tmpGR)
  }
  
  # Add data and ASM fields
  outGR$ASM <- 0
  outGR$Data <- 0
  
  # Return GR with all haplotypes
  return(outGR)
  
}
#Read in gff file for each subject
readSubGff <- function(inDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Import GFF
  outGR <- import.gff(paste(inDir,sub,"_het.cpelasm.gff",sep=""))
  
  # Retain a the required subset
  outGR <- outGR[,c("N","CpGs","hetCpGg1","hetCpGg2")]
  
  # Add subject metadata column
  outGR$Subject <- sub
  
  # Make N numeric and filter
  outGR$N <- as.numeric(outGR$N)
  
  # Add genome info 
  outGR <- setGenomeLengths(outGR,chrsOfInterest=chrsOfInterest)
  
  # Return GR with all haplotypes
  return(outGR)
  
}

# Function to plot a GR using Gviz
plotGR <- function(CpGdir,enhancerDir,GR,startHight,highSize=500,reverseStrand=FALSE,chr="chr11",lim=c(2010000,2022500)){
  
  # Get genome
  gen <- "hg19"

  # GRanges to intersect with and keep the relevant data
  windowGR <- GRanges(seqnames=chr,ranges=IRanges(start=lim[1],end=lim[2]),strand="*")
  
  # Subset GR
  GR <- subsetByOverlaps(GR,windowGR)

  # Create data tracks
  dmmlTrack <- DataTrack(GR[GR$Statistic=="dMML",c("Value")],name="dMML")
  dnmeTrack <- DataTrack(GR[GR$Statistic=="dNME",c("Value")],name="dNME")
  ucTrack <- DataTrack(GR[GR$Statistic=="UC",c("Value")],name="UC")
  
  # Gene track 1
  bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="grch37.ensembl.org",
                path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
  biomTrack <- BiomartGeneRegionTrack(genome=gen,chromosome=chr,start=lim[1],end=lim[2],name="ENSEMBL",biomart=bm)

  # Add CpG island annotation track
  genomicFeatures <- getGeneralFeats(CpGdir,enhancerDir)
  cpgIslands <- genomicFeatures[["CpG island"]]
  cpgIslands <- subsetByOverlaps(cpgIslands,windowGR)
  islandTrack <- AnnotationTrack(cpgIslands,name="CpG islands")
  
  # Chromosome information tracks
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome=gen,chromosome=chr)
  
  # Highlight
  ht <- HighlightTrack(trackList=list(dmmlTrack,dnmeTrack,ucTrack),start=startHight,width=highSize,chromosome=chr)
  
  # Return plot
  plotTracks(list(itrack,gtrack,biomTrack,islandTrack,ht),from=lim[1],to=lim[2],
             transcriptAnnotation="symbol",type=c("gradient"),stacking="squish",reverseStrand=reverseStrand,
             collapseTranscripts = "meta")
  
}
# Function to assign genome
assignGenomeFile <- function(row) {
  
  genome = NA
  if((row[1]=="REF")&(row[2] %in% c("0|1","0/1"))){
    genome = 1
  } else if((row[1]=="ALT")&(row[2] %in% c("0|1","0/1"))){
    genome = 2
  } else if((row[1]=="REF")&(row[2] %in% c("1|0","1/0"))) {
    genome = 2
  } else if((row[1]=="ALT")&(row[2] %in% c("1|0","1/0"))){
    genome = 1
  } else {
    # nothing
  }
  
  return(genome)
}

# Function that returns bool if variant results in a CpG site
hetCpgSite <- function(row) {
  
  # Check if we have C or G, otherwise return false
  if(row[3] %in% c("C","G")){
    # continue
  } else {
    return(FALSE)
  }
  
  # Initialize output
  hetCpg <- FALSE
  context <- toString(getSeq(Hsapiens,row[1],start=as.numeric(row[2])-1,end=as.numeric(row[2])+1,strand="+"))
  if((row[3]=="C")&(substr(context,3,3)=="G")){
    hetCpg <- TRUE
  } else if((row[3]=="G")&(substr(context,1,1)=="C")){
    hetCpg <- TRUE
  } else {
    # nothing
  }
  
  # Return binary vector
  return(hetCpg)
  
}

###Here most current version of Het CpG analysis

#From vcf file, extract het CpG information
extractHetCpG<-function(vcfDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  cat('Processing subject:', sub,'\n')
  tt1=proc.time()[[3]]
  #genomeGr <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  #genomeGr <- setGenomeLengths(genomeGr)
  vcf <- readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19")
  gt <- as.vector(geno(vcf)$GT)
  vcf <- rowRanges(vcf)
  vcf$GT <- gt
  vcf$snpId <- paste(sub,seq(1:length(vcf)),sep="-")
  # Keep only relevant variables
  vcf <- vcf[,c("REF","ALT","GT","snpId")]
  vcf$REF <- as.character(vcf$REF)
  vcf$ALT <- as.character(unlist(vcf$ALT))
  names(vcf)=NULL
  # Delete labels
  vcf=het_CpG_df(vcf)
  vcf$sub=sub
  cat('Done processing',sub,'using:', proc.time()[[3]]-tt1,'\n')
  return(vcf)
  
}
#ID heterozygous CpG from vcf file
het_CpG_df<-function(var_gr){
  #Get plus one location
  plus_loc=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,1)))
  minus_loc=unlist(getSeq(Hsapiens,GenomicRanges::shift(var_gr,-1)))
  #get dinucleotide for ref, alt, plus and minus
  var_gr$REF_plus=paste_nucleotide(var_gr$REF,plus_loc,'plus')
  var_gr$REF_minus=paste_nucleotide(var_gr$REF,minus_loc,'minus')
  var_gr$ALT_plus=paste_nucleotide(var_gr$ALT,plus_loc,'plus')
  var_gr$ALT_minus=paste_nucleotide(var_gr$ALT,minus_loc,'minus')
  #get trinucleotide
  var_gr$REF_tri=paste_trinucleotide(minus_loc,var_gr$REF,plus_loc)
  var_gr$ALT_tri=paste_trinucleotide(minus_loc,var_gr$ALT,plus_loc)
  #check if heterogygouze: note rowSum =2 have trinucleotide form CGG with ref =G alt =C
  var_gr$npmCG=rowSums(as.data.frame(var_gr)[,c('REF_plus','REF_minus','ALT_plus','ALT_minus')]=='CG')
  var_gr$HetCpg=var_gr$npmCG>0
  return(var_gr)
}
#Get dinucleotide 
paste_nucleotide<-function(ref,seq,direction){
  df=data.frame(ref=ref,seq=unlist(seq))
  if(direction=='plus'){return(paste(df$ref,df$seq,sep=''))}
  else if(direction=='minus'){return(paste(df$seq,df$ref,sep=''))}
  else{print('wrong direction')}
}
#Get trinucleotide
paste_trinucleotide<-function(minus,ref,plus){
  df=data.frame(minus=unlist(minus),ref=ref,plus=unlist(plus))
  return(paste(df$minus,df$ref,df$plus,sep=''))
}

#Count number of Heterozygous CpG at each gff region
gff_hetCpG_count<-function(sub,gff_in,vcf_in,CpG){
  #Read in gff file
  gff_sub=gff_in[gff_in$Subject==sub]
  #Read in vcf het CpG information
  vcf_sub=vcf_in[[sub]]
  #Filter het CpG regions
  vcf_sub_het=vcf_sub[vcf_sub$HetCpg]
  #For each gff region count overlaps
  #gff_het_count=countOverlaps(gff_sub,vcf_sub_het)
  #cat(sub,':',sum(gff_het_count>1)/length(gff_het_count)*100, '% all gff region have more than 1 het CpG\n', sep='')
  #cat(sub,':',sum(gff_het_count>1)/sum(gff_het_count>0)*100, '% all het CpG region have more than 1 het CpG\n', sep='')
  #cat(sub,':',sum(gff_het_count>0)/length(gff_het_count)*100, '% gff region have het CpG\n', sep='')
  olap=findOverlaps(vcf_sub_het,gff_sub,type='within',select='all')
  #Count number of CG
  vcf_sub_het$refCG=(vcf_sub_het$REF_plus=='CG')+(vcf_sub_het$REF_minus=='CG')
  vcf_sub_het$altCG=(vcf_sub_het$ALT_plus=='CG')+(vcf_sub_het$ALT_minus=='CG')
  #Add Het CpG information
  gff_sub$HetCpG=FALSE
  gff_sub[subjectHits(olap)]$HetCpG=TRUE
  #refCG, SNP have CG in ref 
  #altCG SNP have CG in alt
  gff_sub$refCG=0
  gff_sub$altCG=0
  #Count number of Het CpG here
  df_sub=data.frame(qt=subjectHits(olap),refCG=vcf_sub_het$refCG[queryHits(olap)],altCG=vcf_sub_het$altCG[queryHits(olap)])
  agg_sub=aggregate(df_sub,by=list(df_sub$qt),FUN=sum)
  gff_sub$refCG[agg_sub$Group.1]=agg_sub$refCG
  gff_sub$altCG[agg_sub$Group.1]=agg_sub$altCG
  gff_sub$N_nonhet=countOverlaps(gff_sub,CpG)-countOverlaps(gff_sub,vcf_sub_het[vcf_sub_het$refCG>0])
  gff_sub$N_hg19=countOverlaps(gff_sub,CpG)
  return(gff_sub)
}
#Count number of hetCG at each allele GR from result, add gff size, and N

#Put hetCpG count into each sample in gr_allele
hetCGallele_sub<-function(sub,gr_allele,gff,CpG,vcf_in){
  cat('Analyzing',sub,'\n')
  t1=proc.time()[[3]]
  sub_vcf=vcf_in[[sub]]
  het_vcf=sub_vcf[sub_vcf$HetCpg]
  sub_allele=gr_allele[gr_allele$Subject==sub]
  sub_het=gff[[sub]]
  #gr_allele got resized with start and end +1, use +2 to include equal start & end
  sub_het=resize(sub_het, width(sub_het) + 4, fix="center")
  sub_ref=sub_allele[sub_allele$Genome==1]
  sub_alt=sub_allele[ sub_allele$Genome==2]
  sub_ref=allele_hetCG(sub_ref,sub_het,'refCG')
  sub_alt=allele_hetCG(sub_alt,sub_het,'altCG')
  gr_out=c(sub_ref,sub_alt)
  gr_out=GR_resize(gr_out,CpG,het_vcf,gene_size=500,sub)
  cat('Finish analyzing',sub,proc.time()[[3]]-t1,'\n')
  return(gr_out)
}
#Resize for each subject
GR_resize_sub<-function(sub,gr_allele,CpG,vcf_in){
  cat(paste('Analyzing',sub,'\n',sep=' '))
  tt1=proc.time()[3]
  sub_vcf=vcf_in[[sub]]
  het_vcf=sub_vcf[sub_vcf$HetCpg]
  gr_out=gr_allele[gr_allele$Subject==sub]
  gr_out=GR_resize(gr_out,CpG,het_vcf,sub)
  #saveRDS(gr_out,paste('../downstream/temp/gr_resize_sub2',sub,'.rds',sep=''))
  cat(paste('Finishing analyzing',sub,'in',proc.time()[3]-tt1,'\n',sep=' '))
  return(gr_out)
}

#for each type of CG find gff overlap, add width and N
allele_hetCG<-function(allele,gff,cgtype){
  #Find overlap region
  olap=findOverlaps(allele,gff,type='within',select="all")
  #NAs to be solved: check size of na
  olap_na=which(is.na(olap))
  cat("size of NA in olap:",length(olap_na),'\n')
  if (length(olap_na)!=0){
    allele=allele[-olap_na]
    olap=olap[-olap_na]
  }
  #separate 1 overlap and more than 1 overlap
  qt=queryHits(olap)
  qt_table=table(qt)
  olap_sub=subjectHits(olap)
  #Find overlap regions =1
  olap_1=as.numeric(names(qt_table[qt_table==1]))
  #for unique data, gff may have more than 1 times subsetted
  olap_gff_1=olap_sub[which(qt%in%olap_1)]
  allele=allele_addMeta(allele,olap_1,olap_gff_1,gff,cgtype)
  #find overlap regions >1
  olap_more=as.numeric(names(qt_table[qt_table>1]))
  #for data more than 1 hits, use smaller region in gff
  gff_out=unlist(lapply(olap_more,gff_min,olap_sub=olap_sub,qt=qt,gffwid=width(gff)))
  #merge unique ranges
  allele=allele_addMeta(allele,olap_more,gff_out,gff,cgtype)
  
  return(allele)
}
#Find min width
gff_min<-function(op_more,olap_sub,qt,gffwid){
  gff_idx=olap_sub[which(qt==op_more)]
  return(gff_idx[which.min(gffwid[gff_idx])])
}
#add metadata to allele GR
allele_addMeta<-function(allele,qh,sh,gff,cgtype){
  allele$CpGallele[qh]=elementMetadata(gff)[sh,cgtype]
  allele$N[qh]=gff$N[sh]
  allele$HetCpG[qh]=gff$HetCpG[sh]
  allele$N_hg19[qh]=gff$N_hg19[sh]
  allele$gff_size[qh]=width(gff[sh])
  return(allele)
}
###End of calculation

#Add ASM information to allele values:

add_ASM<-function(allele_GR,ASM_GR){
  allele_GR_out=GRanges()
  for (sp in unique(allele_GR$Sample)){
    ASM_sub=ASM_GR[ASM_GR$Sample==sp]
    allele_sub=allele_GR[allele_GR$Sample==sp]
    allele_sub$ASM=NA
    allele_sub$pval=NA
    olap=findOverlaps(allele_sub,ASM_sub)
    allele_sub$ASM[queryHits(olap)]=ASM_sub$ASM[subjectHits(olap)]
    allele_sub$pval[queryHits(olap)]=ASM_sub$pvalue[subjectHits(olap)]
    allele_GR_out=c(allele_GR_out,allele_sub)
  }
  return(allele_GR_out) 
}

#calculate allelic difference and plot
allele_calc_plot<-function(gr_allele_in,stat,outDir,picname){
  gr_diff_calc=allele_diff(gr_allele_in)
  allele_plot(gr_diff_calc,stat,outDir,picname)
  return(gr_diff_calc)
}
#make plot
allele_plot<-function(gr_diff_calc,stat,outDir,picname){
  
  gr_diff=gr_diff_calc[[1]]
  gr_diff$type='Non Het CpG'
  gr_diff$type[gr_diff$CGcount_diff!=0]='Het CpG'
  
  #Use density plot
  plot_df=data.frame(gr_diff$diff,gr_diff$type)
  colnames(plot_df)=c('Allele_difference','Allele_type')
  
  diff_plot_non_het=ggplot(plot_df[plot_df$Allele_type=='Non Het CpG',],aes(x=Allele_difference))+
    geom_density(alpha=0.6,fill='purple',color='black',size=1.5)+xlab(paste(stat, 'Difference'))+
    theme(legend.position="none",plot.title = element_text(hjust=0.5))+
    ggtitle(paste('Allilic difference at',picname,'with no Het-CpG'))+ylim(0,3)
  
  diff_plot_het_CpG=ggplot(plot_df[plot_df$Allele_type=='Het CpG',],aes(x=Allele_difference))+
    geom_density(alpha=0.6,fill='purple',color='black',size=1.5)+xlab(paste(stat, 'Difference'))+
    theme(legend.position="none",plot.title = element_text(hjust=0.5))+
    ggtitle(paste('Allilic difference (More CpG-Less CpG) at',picname, 'with Het CpG'))+ylim(0,3)
  diff_plot=arrangeGrob(diff_plot_het_CpG,diff_plot_non_het,nrow=2,ncol=1)
  ggsave(paste(outDir,stat,' diff ',picname,'.png',sep=''),diff_plot,width=8)
  #More CpG genome2-genome1
  allele_gr_ASM=gr_diff_calc[[2]]
  allele_gr_ASM=allele_gr_ASM[allele_gr_ASM$alleleCpG!=0]
  allele_het_df=rbind(data.frame(value=allele_gr_ASM$Value[allele_gr_ASM$CpGstat=='More'],type='More CpG'),
                      data.frame(value=allele_gr_ASM$Value[allele_gr_ASM$CpGstat=='Less'],type='Less CpG'))
  
  dist_plot=ggplot(allele_het_df,aes(x=value,color=type))+
    geom_density(alpha=0.6,size=1)+xlab(stat)+ggtitle(paste(stat,'at',picname,'with different types of haplotypes'))+
    theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+ylim(0,4)+scale_color_manual(values=c("blue","red"))+
    geom_hline(yintercept=0, colour="white", size=1)
  ggsave(paste(outDir,stat,' dist ',picname,'.png'),dist_plot)
  #reference distribution
  allele_gr_ASM=gr_diff_calc[[2]]
  allele_het_df=rbind(data.frame(value=allele_gr_ASM$Value[allele_gr_ASM$Genome==1],type='Ref'),
                      data.frame(value=allele_gr_ASM$Value[allele_gr_ASM$Genome==2],type='Alt'))
  
  dist_plot_ref=ggplot(allele_het_df,aes(x=value,fill=type))+
    geom_density(alpha=0.6)+xlab(stat)+ggtitle(paste(stat,'at',picname,'reference distribution'))+
    theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
    scale_fill_manual( values = c("blue","red"))+ylim(0,4)
  ggsave(paste(outDir,stat,' dist ref ',picname,'.png'),dist_plot_ref)
}

#Calculate allelic difference

allele_diff<-function(allele_gr_in){
  out_gr=GRanges()
  out_allele=GRanges()
  for(sample in unique(allele_gr_in$Sample)){
    cat('Processing:',sample,'\n')
    #Sort genome
    genome1=sort(allele_gr_in[allele_gr_in$Genome==1 & allele_gr_in$Sample==sample])
    genome2=sort(allele_gr_in[allele_gr_in$Genome==2 & allele_gr_in$Sample==sample])
    #Mroe CpG or less CpG
    genome2$CpGstat=NA
    genome1$CpGstat=NA
    cat('Identical check:',all(genome1$HetCpG==genome2$HetCpG),'\n')
    gr=granges(genome1)
    gr$N=genome1$N
    #always ref - alt, but need to be more CpG-less, if ref > alt, need ref-ale, then g1-g2 >0
    # else if ref<alt, need alt-ref, then g1-g2<0
    sign=(genome1$CpGallele-genome2$CpGallele)
    #Assign allele stat based on sign
    #For genome 2=alt, sign >0 mean ref > alt, then the alt is less
    genome2$CpGstat[sign!=0]='More'
    genome2$CpGstat[sign>0]='Less'
    genome2$alleleCpG=sign
    #Assign allele stat based on sign
    #For genome 1=ref, sign >0 mean ref > alt, then the alt is more
    genome1$alleleCpG=sign
    genome1$CpGstat[sign!=0]='Less'
    genome1$CpGstat[sign>0]='More'
    #Put sign here
    gr$CpGdiff=sign
    #Normalize sign
    sign[sign!=0]=sign[sign!=0]/abs(sign[sign!=0])
    sign[sign==0]=1
    #Use normalized sign to apply sign to the genome value
    #all is ref - alt, if sign ==1, then nCpG ref > alt, keep same
    #if sign ==-1, then nCpg alt > ref, then apply negative sign here because it should be alt -ref
    gr$diff=(genome1$Value-genome2$Value)*sign
    gr$mean=(genome1$Value+genome2$Value)/2
    #add meta information
    gr$Sample=sample
    gr$Statistic=paste('d',unique(genome1$Statistic),sep='')
    gr$HetCpG=genome1$HetCpG
    gr$pval=genome1$pval
    gr$ASM=genome1$ASM
    gr$Subject =genome1$Subject 
    #Calculating density
    gr$density =genome1$CGcount_hg19_extend
    #gr$density_diff=(genome1$density-genome2$density)*sign
    gr$TpA_CpG_diff=(genome1$TpA_CpG-genome2$TpA_CpG)*sign
    gr$CGcount_diff=(genome1$CGcount_allele_extend-genome2$CGcount_allele_extend)*sign
    gr$CGcount_hg19_extend=genome1$CGcount_hg19_extend
    gr$CG_count_extend=genome1$CG_count_extend
    gr$AT_count_extend=genome1$AT_count_extend
    gr$TpA_count_extend=genome1$TpA_count_extend
    gr$CpG_count_extend=genome1$CG_hg19_extend
    gr$gff_size_extend=genome1$gff_size_extend
    
    #Put allele information and gr information
    out_allele=c(out_allele,genome1,genome2)
    out_gr=c(out_gr,gr)
  }
  return(list(out_gr,out_allele))
}

#generate and export bed file

ASM_bed_gen_sp<-function(gr_ASM,gr_allele,sp,outdir){
  gr_ASM=gr_ASM[gr_ASM$CpGdiff!=0]
  gr_ASM=gr_ASM[order(gr_ASM$diff,decreasing = TRUE)]
  #do it for each sample
  #sp='Gastric - STL001'
  sp_ASM=gr_ASM[gr_ASM$Sample==sp]#stats for this subject
  sp_allele=gr_allele[gr_allele$Sample==sp & gr_allele$HetCpG]
  sp_bed=granges(sp_ASM)
  sp_bed$diff=sp_ASM$diff
  olap1=findOverlaps(sp_allele[sp_allele$Genome==1],sp_bed,type='equal')
  olap2=findOverlaps(sp_allele[sp_allele$Genome==2],sp_bed,type='equal')
  sp_bed$A1[subjectHits(olap1)]=sp_allele[sp_allele$Genome==1]$Value[queryHits(olap1)]
  sp_bed$A2[subjectHits(olap2)]=sp_allele[sp_allele$Genome==2]$Value[queryHits(olap2)]
  #Export bed file
  export_bed(sp_bed,'diff',paste(outdir,sp,'_hetASM_diff.bedGraph',sep=''))
  export_bed(sp_bed,'A1',paste(outdir,sp,'_hetASM_A1.bedGraph',sep=''))
  export_bed(sp_bed,'A2',paste(outdir,sp,'_hetASM_A2.bedGraph',sep=''))
}
export_bed<-function(gr_in,dat,out_name){
  bed_out=granges(gr_in)
  bed_out$score=elementMetadata(gr_in)[,dat]
  export(bed_out,out_name,format='bedGraph')
}

#Generate density from GR allele, at least need 200 bp for density
GR_resize<-function(GR.in,CpG_sites,hetCpG,sub,gene_size=500){
  ##From definitino of CpG island, use 200 bp regions
  GR.extend=resize(GR.in,width=width(GR.in)+gene_size,fix='center')
  GR.in$CG_hg19_extend=countOverlaps(GR.extend,CpG_sites)
  GR.in$CG_nonhet_extend=GR.in$CG_hg19_extend-countOverlaps(GR.extend,hetCpG[hetCpG$REF_plus=='CG'|hetCpG$REF_minus=='CG'])+
    countOverlaps(GR.extend,hetCpG[hetCpG$ALT_plus=='CG'|hetCpG$ALT_minus=='CG'])
  #Count CpG in genome 1
  GR.in$CG_allele_extend=NA
  GR.in$CG_allele_extend[GR.in$Genome==1]=GR.in$CG_nonhet_extend[GR.in$Genome==1]+countOverlaps(GR.extend[GR.extend$Genome==1],hetCpG[hetCpG$REF_plus=='CG'|hetCpG$REF_minus=='CG'])
  GR.in$CG_allele_extend[GR.in$Genome==2]=GR.in$CG_nonhet_extend[GR.in$Genome==2]+countOverlaps(GR.extend[GR.extend$Genome==2],hetCpG[hetCpG$ALT_plus=='CG'|hetCpG$ALT_minus=='CG'])
  GR.in$CG_nonhet_extend=NULL
  GR.in$gff_size_extend=width(GR.extend)
  #calculate Odds ratio for expected CG vs actual CG
  #Expected CG number C * number G/total length
  # Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes". Journal of Molecular Biology.
  #Wiki, actual: ((number of C + Number of G)/2)^2/length of genomics Normalized CpG content, whole genome ~25%
  #Check this command
  gr_seq=getSeq(Hsapiens,GR.extend,as.character=T)
  GR.in$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
  GR.in$TpA_count_extend=do.call('c',lapply(gr_seq,function(x){countPattern('TA',x)}))
  GR.in$AT_count_extend=do.call('c',lapply(gr_seq,function(x){countPattern('T',x)}))+do.call('c',lapply(gr_seq,function(x){countPattern('A',x)}))
  GR.in$CG_count_extend=do.call('c',lapply(gr_seq,function(x){countPattern('C',x)}))+do.call('c',lapply(gr_seq,function(x){countPattern('G',x)}))
  GR.in$CGcount_hg19_extend=GR.in$CG_hg19_extend/GR.in$CGcont_exp
  GR.in$CGcount_allele_extend=GR.in$CG_allele_extend/GR.in$CGcont_exp
  saveRDS(GR.in,paste('../downstream/temp/gr_resize_sub2',sub,'.rds',sep=''))
  return(GR.in)
}

#CountPattern and return data.frame
countCGOR<-function(x){ #x=input seq
  #calculate Odds ratio for expected CG vs actual CG
  #Expected CG number C * number G/total length
  # Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes". Journal of Molecular Biology.
  #Wiki, actual: ((number of C + Number of G)/2)^2/length of genomics Normalized CpG content, whole genome ~25%
  NC=countPattern('C',x)
  NG=countPattern('G',x)
  #NCG=countPattern('CG',x)
  CG_exp=NC*NG/nchar(x) #PC*PG*length
  CG_exp_norm=((NC+NG)/2)^2/nchar(x) #assuming PC=PG = (NC+NG)/2/length
  return(CG_exp_norm)
}

#Generate motif from variant file

#From vcf file, extract het CpG information
extractmotif<-function(vcfDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  cat('Processing subject:', sub,'\n')
  tt1=proc.time()[[3]]
  #genomeGr <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  #genomeGr <- setGenomeLengths(genomeGr)
  vcf <- readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19")
  gt <- as.vector(geno(vcf)$GT)
  vcf <- rowRanges(vcf)
  vcf$GT <- gt
  vcf$snpId <- paste(sub,seq(1:length(vcf)),sep="-")
  # Keep only relevant variables
  vcf <- vcf[,c("REF","ALT","GT","snpId")]
  vcf$REF <- as.character(vcf$REF)
  vcf$ALT <- as.character(unlist(vcf$ALT))
  names(vcf)=NULL
  # Delete labels
  vcf=motif_df(vcf)
  vcf$sub=sub
  cat('Done processing',sub,'using:', proc.time()[[3]]-tt1,'\n')
  return(vcf)
  
}
#ID motif from vcf file
motif_df<-function(var_gr){
  #Get 11 nucleotide location, 3 nucleotide 4-6, 5 nulceotide 3-7, 7 nucleotide 2-8
  nucleo_11_alt<-nucleo_11_X<-nucleo_11_ref<-as.matrix(getSeq(Hsapiens,resize(var_gr,11,fix='center')))
  nucleo_11_X[,6]='N'
  nucleo_11_alt[,6]=var_gr$ALT
  #Running too long. using matrix operation
  #cat matrix into 1 string
  nucleo_11_ref=apply(nucleo_11_ref,1,function(x) paste(x,collapse=''))
  nucleo_11_X=apply(nucleo_11_X,1,function(x) paste(x,collapse=''))
  nucleo_11_alt=apply(nucleo_11_alt,1,function(x) paste(x,collapse=''))
  #Get 3,5,7,9 nucleotide location
  nucleo_9_ref=substr(nucleo_11_ref,start=2,stop=10)
  nucleo_7_ref=substr(nucleo_11_ref,start=3,stop=9)
  nucleo_5_ref=substr(nucleo_11_ref,start=4,stop=8)
  nucleo_3_ref=substr(nucleo_11_ref,start=5,stop=7)
  
  nucleo_9_alt=substr(nucleo_11_alt,start=2,stop=10)
  nucleo_7_alt=substr(nucleo_11_alt,start=3,stop=9)
  nucleo_5_alt=substr(nucleo_11_alt,start=4,stop=8)
  nucleo_3_alt=substr(nucleo_11_alt,start=5,stop=7)
  
  nucleo_9_X=substr(nucleo_11_X,start=2,stop=10)
  nucleo_7_X=substr(nucleo_11_X,start=3,stop=9)
  nucleo_5_X=substr(nucleo_11_X,start=4,stop=8)
  nucleo_3_X=substr(nucleo_11_X,start=5,stop=7)
  
  #Assign value to var_gr
  var_df=data.frame(nucleo_11_ref=nucleo_11_ref,
  nucleo_9_ref=nucleo_9_ref,
  nucleo_7_ref=nucleo_7_ref,
  nucleo_5_ref=nucleo_5_ref,
  nucleo_3_ref= nucleo_3_ref,
  
  nucleo_11_alt=nucleo_11_alt,
  nucleo_9_alt=nucleo_9_alt,
  nucleo_7_alt=nucleo_7_alt,
  nucleo_5_alt=nucleo_5_alt,
  nucleo_3_alt=nucleo_3_alt,
  
  nucleo_11_X=nucleo_11_X,
  nucleo_9_X=nucleo_9_X,
  nucleo_7_X=nucleo_7_X,
  nucleo_5_X=nucleo_5_X,
  nucleo_3_X=nucleo_3_X,stringsAsFactors = FALSE)
  var_gr=makeGRangesFromDataFrame(cbind(var_gr,var_df),keep.extra.columns = TRUE)
  return(var_gr)
}

#For each sample, add dMML and dNME information
extract_diff_values<-function(sp,diff,variant){
  #Get subject information for sp
  subj= strsplit(sp,' - ')[[1]][2]
  variant=variant[[subj]]
  #For this sample, extract the loci
  #dMML
  dMML=diff[diff$Statistic=='dMML']
  olap_dMML=findOverlaps(variant,dMML,type='within')
  outGR=variant[queryHits(olap_dMML)]
  outGR$dMML=dMML$Value[subjectHits(olap_dMML)]
  outGR$dMML_pval=dMML$pvalue[subjectHits(olap_dMML)]
  #dNME
  dNME=diff[diff$Statistic=='dNME']
  olap_dNME=findOverlaps(outGR,dNME,type='within')
  outGR$dNME=NA
  outGR$dNME_pval=NA
  outGR[queryHits(olap_dNME)]$dNME=dNME$Value[subjectHits(olap_dNME)]
  outGR[queryHits(olap_dNME)]$dNME_pval=dNME$pvalue[subjectHits(olap_dNME)]
  outGR$Sample=sp
  return(outGR)
}


#Reshape each motif and keep extra column
reshape_sample_variant<-function(variant){
  variant=as.data.frame(variant)
  melt_id_var=colnames(variant)
  melt_id_var=melt_id_var[-which(melt_id_var=="REF" | melt_id_var=="ALT")]
 
  # Melt: origin = ref/alt, variant= ATCG
  variant <- melt(data=variant,id.vars=melt_id_var,
                measure.vars=c("REF","ALT"),
                variable.name="Origin",value.name="Variant")
  #Assign genome
  variant$Genome <- apply(variant[,c("Origin","GT")],1,assignGenomeFile)
  #Reshape ref columns
  ref_id=c('nucleo_11_ref','nucleo_9_ref','nucleo_7_ref','nucleo_5_ref','nucleo_3_ref')
  alt_id=c('nucleo_11_alt','nucleo_9_alt','nucleo_7_alt','nucleo_5_alt','nucleo_3_alt')
  new_id=c('nucleo_11','nucleo_9','nucleo_7','nucleo_5','nucleo_3')
  #reassign ID
  variant[,new_id]=NA
  variant[variant$Genome==1,new_id]=variant[variant$Genome==1,ref_id]
  variant[variant$Genome==2,new_id]=variant[variant$Genome==2,alt_id]
  variant=variant[,-which(colnames(variant)%in%c(ref_id,alt_id))]
  #Make Granges
  variant=makeGRangesFromDataFrame(variant,keep.extra.columns = TRUE)
  return(variant)
}
assignGenomeFile <- function(row) {
  
  genome = NA
  if((row[1]=="REF")&(row[2] %in% c("0|1","0/1"))){
    genome = 1
  } else if((row[1]=="ALT")&(row[2] %in% c("0|1","0/1"))){
    genome = 2
  } else if((row[1]=="REF")&(row[2] %in% c("1|0","1/0"))) {
    genome = 2
  } else if((row[1]=="ALT")&(row[2] %in% c("1|0","1/0"))){
    genome = 1
  } else {
    # nothing
  }
  
  return(genome)
}
#Extract allele values for each allele
extract_allele_value<-function(outGR,cpelAllele){
  # Cross resulting GR with MML of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with MML of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="MML")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$MML2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$MML2[outGR$Genome=="1"] <- NA
  
  # Consolidate MML1 and MML2 columns into single column
  outGR$MML <- rowSums(cbind(outGR$MML1,outGR$MML2), na.rm=T)
  outGR$MML1 <- outGR$MML2 <- NULL
  
  # Cross resulting GR with NME of genome 1
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="1")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME1[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME1[outGR$Genome=="2"] <- NA
  
  # Cross resulting GR with NME of genome 2
  cpelAlleleTmp <- cpelAllele[(cpelAllele$Statistic=="NME")&(cpelAllele$Genome=="2")]
  olaps <- findOverlaps(outGR,cpelAlleleTmp,type="within",select="all")
  outGR$NME2[queryHits(olaps)] <- cpelAlleleTmp$Value[subjectHits(olaps)]
  outGR$NME2[outGR$Genome=="1"] <- NA
  
  # Consolidate NME1 and NME2 columns into single column
  outGR$NME <- rowSums(cbind(outGR$NME1,outGR$NME2), na.rm=T)
  outGR$NME1 <- outGR$NME2 <- NULL

  # Return
  return(outGR)
}
#Enrichment test
motif_enrichment<-function(GR_allele,pval_cutoff=0.1,p_stat,motif,motif_type){#pstat either dNME_pval or dMML_pval
  GR_allele=as.data.frame(GR_allele)
  GR_allele$ASM=NA
  GR_allele$ASM[GR_allele[,p_stat]<=pval_cutoff] =  TRUE
  GR_allele$ASM[GR_allele[,p_stat]>pval_cutoff] = FALSE
  ASM_motif=sum(GR_allele$ASM & GR_allele[,motif_type]==motif)
  ASM_not_motif=sum(GR_allele$ASM & !GR_allele[,motif_type]==motif)
  notASM_motif=sum(!GR_allele$ASM & GR_allele[,motif_type]==motif)
  notASM_notmotif=sum(!GR_allele$ASM & !GR_allele[,motif_type]==motif)
  cont_table=matrix(c(ASM_motif,ASM_not_motif,notASM_motif,notASM_notmotif),nrow=2,byrow=TRUE)
  return(fisher.test(cont_table))
}
#Give each SNP an ASM information for each subject
variant_meta_subj<- function(subj,variant_in,GR_in){ #variant_in for each subject, GR_in for each subject
  cat('Processing',subj,'\n')
  GR_in_subj=GR_in[GR_in$Subject==subj]
  variant_in_subj=variant_in[[subj]]
  sp=unique(GR_in_subj$Sample)
  gr_out_sp = GRanges()
  for (sps in sp){
    gr_out_sp=c(gr_out_sp,variant_meta_sp(variant_in_subj,GR_in_subj[GR_in_subj$Sample==sps]))
    
  }
  return(gr_out_sp)
  
  
}
#for given sample
variant_meta_sp<-function(variant_subj,GR_st){
  st=unique(GR_st$Statistic)
  gr_out_st = GRanges()
  for (statistics in st){
    gr_out_st=c(gr_out_st,variant_meta_sp_st(variant_subj,GR_st[GR_st$Statistic==statistics]))
    
  }
  return(gr_out_st)
}
#for given stat: need to debug
variant_meta_sp_st <-function(variant_subj,GR_sp){
  olap=findOverlaps(variant_subj,GR_sp,type = 'within')
  #Find the overlap region
  gr_out=variant_subj[queryHits(olap)]
  olap=findOverlaps(gr_out,GR_sp,type = 'within')
  gr_out$pvalue[queryHits(olap)]=GR_sp$pvalue[subjectHits(olap)]
  gr_out$Value[queryHits(olap)]=GR_sp$Value[subjectHits(olap)]
  gr_out$Statistic[queryHits(olap)]=GR_sp$Statistic[subjectHits(olap)]
  gr_out$Sample[queryHits(olap)]=GR_sp$Sample[subjectHits(olap)]
  gr_out$Subject[queryHits(olap)]=GR_sp$Subject[subjectHits(olap)]
  gr_out$HetCpG=FALSE
  gr_out$HetCpG=((gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & !(gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG')) |
    (!(gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & (gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG'))
  return(gr_out)
}
#Enrichment of variants
variant_enrich<-function(variant_in,cutoff=0.1){
  variant_in$ASM=FALSE
  variant_in$ASM[variant_in$dNME_pval<=cutoff]=TRUE
  cont_table=matrix(c(sum(variant_in$variant & variant_in$ASM),
                      sum(variant_in$variant & !variant_in$ASM),
                      sum(!variant_in$variant & variant_in$ASM),
                      sum(!variant_in$variant & !variant_in$ASM)),
                    nrow=2,byrow=TRUE)
  colnames(cont_table)=c('ASM','Not ASM')
  rownames(cont_table)=c('Vairant','Not variant')
  print(cont_table)
  return(fisher.test(cont_table))
}
#Subset for SNP-containing ranges
SNP_conmtaining_hap<-function(gr_in,variant_in){
  SNP_sub=GRanges()
  SNP_not=GRanges()
  for (subj in unique(gr_in$Subject)){
    olap=findOverlaps(gr_in[gr_in$Subject==subj],variant_in[[subj]])
    SNP_sub=c(SNP_sub,subsetByOverlaps(gr_in[gr_in$Subject==subj], variant_in[[subj]]))
    SNP_not=c(SNP_not,gr_in[gr_in$Subject==subj][-queryHits(olap)])
  }
return(list(SNP_containing=SNP_sub,Non_SNP_containing=SNP_not))
}
#ASM_het_CpG_enrichment
ASM_het_enrich<-function(gr_in,title){
OR_df=data.frame(sp=NULL,OR=NULL,lower_CI=NULL,upper_CI=NULL)
for (sp in unique(gr_in$Sample)){
  OR=ASM_het_enrichment(gr_in[gr_in$Sample==sp])
  OR_df=rbind(OR_df,data.frame(sp=sp,subjects=strsplit(sp,' - ')[[1]][2],
            OR=OR$estimate,lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
}
theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
ggplot(OR_df,aes(x=sp,y=OR,fill=subjects)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0, 13)+
  ggtitle(title)+xlab('Sample name')+ylab('Odds Ratio')+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                position=position_dodge(.9))+  theme_bar

}

#Plot global distribution of density vs CG type
plot_density<-function(CG_df,ylab,title,ylim,xlab='CG type'){
ggplot(CG_df,aes(x=CG_type,y=density,fill=CG_type))+#scale_fill_manual(values = c("blue","red"))+
  geom_boxplot(outlier.shape = NA)+xlab(xlab)+ylab(ylab)+ylim(ylim)+
  theme(legend.position="bottom",plot.title = element_text(hjust=0.5))+
  ggtitle(title)+theme(legend.title = element_blank())
}
#CpG feature enrichment plot
genome_feature_plot<-function(gr,feature,Stats,title,ylim=c(0,2)){
  subjects=unique(gr$Subject)
  OR_df=data.frame(subject=subjects,OR=0,lower_CI=0,upper_CI=0)
  for(subj in subjects){
    OR_out=testEnrichmentFeature_stat(gr[gr$Subject==subj],feature,Stats)
    OR_df$OR[OR_df$subject==subj]=OR_out$estimate
    OR_df$lower_CI[OR_df$subject==subj]=OR_out$conf.int[1]
    OR_df$upper_CI[OR_df$subject==subj]=OR_out$conf.int[2]
  }  
  print(OR_df)
  theme_bar=theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="bottom")
  ggplot(OR_df,aes(x=subject,y=OR,fill=subject)) + ylim(ylim)+
    geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle(title)+
    geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,
                  position=position_dodge(.9))+theme_bar
    
}

#Distance to given granges
gr_distance<-function(gr_in,gr_feature,xlab,main,ylim){
  #gene following gr
  prec.tss.gr=precede(gr_in,gr_feature,select='all',ignore.strand=TRUE)
  gr_in.prec=gr_in[queryHits(prec.tss.gr)]
  #Positive distance for gene follows gr
  gr_in.prec$distance=distance(gr_in.prec,gr_feature[subjectHits(prec.tss.gr)])
  #gene ahead of gr
  follow.tss.gr=follow(gr_in,gr_feature,select='all',ignore.strand=TRUE)
  gr_in.follow=gr_in[queryHits(follow.tss.gr)]
  #Positive distance for gene follows gr
  gr_in.follow$distance=-distance(gr_in.follow,gr_feature[subjectHits(follow.tss.gr)])
  #combine 2 distance
  gr_all=c(gr_in.prec,gr_in.follow)
  #Round and calculate proportion
  gr_all$dist_round=round(gr_all$distance/200)*200
  gr_all_close=gr_all[abs(gr_all$dist_round)<=5000]
  #gr_all_count=table(gr_all_close$dist_round)
  gr_count=table(gr_all_close$dist_round)
  gr_plot_df=data.frame(dist=as.numeric(names(gr_count)),percent_ASM=gr_count/length(gr_all_close))
  plot(gr_plot_df$dist,gr_plot_df$percent_ASM.Freq,pch=1,cex=0.8,ylab='Proportion of ASM',xlab=xlab,main=main,ylim=ylim)
  lines(gr_plot_df$dist,gr_plot_df$percent_ASM.Freq,lwd=1.5)
  abline(h=mean(gr_plot_df$percent_ASM.Freq),lty=2,lwd=4)
}
readEnhancer <- function(enhancerDir){
  
  #Function to convert colnames to granges
  load(paste(enhancerDir,"enhancers_intersect.RData",sep=""))
  enhancer_gr_all <- do.call('c',lapply(rownames(max_states),rownames2Granges))
  
  # Return
  return(enhancer_gr_all)
}
#Find the overlap event in ASM
olap_ASM<-function(GR_in){
  dMML=GR_in[GR_in$Statistic=='dMML']
  dNME=GR_in[GR_in$Statistic=='dNME']
  UC=GR_in[GR_in$Statistic=='UC']
  dMML_dNME=length(subsetByOverlaps(dMML,dNME))
  dMML_UC=length(subsetByOverlaps(dMML,UC))
  dNME_UC=length(subsetByOverlaps(dMML,dNME))
  olap3=length(subsetByOverlaps(subsetByOverlaps(UC,dNME),dMML))
  return(data.frame(dMML_dNME=dMML_dNME-olap3,dMML_UC=dMML_UC-olap3,
                    dNME_UC=dNME_UC-olap3,olap3=olap3,dMML=length(dMML),
                    dNME=length(dNME),UC=length(UC),
                    dMML_nonolap=length(dMML)-dMML_dNME-dMML_UC-olap3,
                    dNME_nonolap=length(dNME)-dMML_dNME-dNME_UC-olap3,
                    UC_nonolap=length(UC)-dNME_UC-dMML_UC-olap3,
                    sample=unique(GR_in$Sample),subject=unique(GR_in$Subject)))
}
#collapase variants
variants_collapase<-function(varsDiff){
  variants <- paste(as.character(varsDiff$REF),as.character(unlist(varsDiff$ALT)),sep="-")
  # Combine same variants (ALT/REF -> REF/ALT)
  variants[variants %in% c("A-C","C-A")] <- "A-C"
  variants[variants %in% c("A-G","G-A")] <- "A-G"
  variants[variants %in% c("A-T","T-A")] <- "A-T"
  variants[variants %in% c("C-G","G-C")] <- "C-G"
  variants[variants %in% c("C-T","T-C")] <- "C-T"
  variants[variants %in% c("G-T","T-G")] <- "G-T"
  return(variants)
}
#Calculate enrichment of each variants in ASM
variants_OR<-function(variant_gr,variant,cutoff=0.05){
  invariant=variant_gr[variant_gr$variants==variant]
  nonvariant=variant_gr[variant_gr$variants!=variant]
  variant_ASM=sum(invariant$pvalue<=cutoff)
  variant_nonASM=sum(invariant$pvalue>cutoff)
  nonvariant_ASM=sum(nonvariant$pvalue<=cutoff)
  nonvariant_nonASM=sum(nonvariant$pvalue>cutoff)
  cont_table=matrix(c(variant_ASM,variant_nonASM,nonvariant_ASM,nonvariant_nonASM),nrow=2)
  return(fisher.test(cont_table))
}
#calculate odds ratio of trinucleotide
tri_nucleo_OR<-function(gr,tri,cutoff=0.05){
  tri_gr=gr[gr$mask_tri==tri]
  nontri=gr[gr$mask_tri!=tri]
  tri_ASM=sum(tri_gr$pvalue<=cutoff)
  tri_nonASM=sum(tri_gr$pvalue>cutoff)
  nontri_ASM=sum(nontri$pvalue<=cutoff)
  nontri_nonASM=sum(nontri$pvalue>cutoff)
  cont_table=matrix(c(tri_ASM,tri_nonASM,nontri_ASM,nontri_nonASM),byrow = T,nrow=2)
  return(fisher.test(cont_table))
}
#Merge 3 stats for each sample
stat_merge<-function(gr_in){
  dMML=gr_in[gr_in$Statistic=="dMML"]
  dNME=gr_in[gr_in$Statistic=="dNME"]
  UC=gr_in[gr_in$Statistic=="UC"]
  gr=granges(gr_in)
  olap=findOverlaps(dMML,dNME)
  gr$dMML=dMML$Value
  gr$dMML_pval=dMML$pvalue
  gr$Sample=gr_in$Sample
  gr$dNME=dNME$Value
  gr$dNME_pval=dNME$pvalue
  if (length(UC$Value)>0){gr$UC=UC$Value}
  if (length(UC$pvalue)>0){gr$UC=UC$pvalue}
  gr$Sample=gr_in$Sample
  return(gr)
}

