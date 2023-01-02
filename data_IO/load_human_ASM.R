rm(list=ls())
source("mainFunctions_sub.R")
#read in all sample tissue diff
import.subject<-function(inDir,calc='diff'){
  #for calc: diff -> dMML etc, allele -> NME etc
  bed_in=dir(inDir,pattern="bedGraph")
  sample_in=unique(sub('_phased.*','',bed_in))
  GRs=GRanges()
  for (sp in sample_in) {
    subjects=sub('_.*','',sp)
    tissues=sub(paste0(subjects,"_"),'',sp)
    # Print sample being loaded
    cat("Loading:",subjects,'-',tissues,'\n')
    if (calc=='diff'){
      GR.in=read.diffGR(subjects,tissues,inDir)
    }else if(calc=='allele'){
      GR.in=read.alleleGR(subjects,tissues,inDir)
    }else {cat('Wrong calc \n')}
    GRs=c(GRs,GR.in)
  }
  return(GRs)
  
}
#Function to read in single GR object:
read.diffGR<-function(subjects,tissues,inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  #Initialization
  GRs=GRanges()
  #Make sure the inputs are unique
  # dmml
  filename=paste(inDir,subjects,"_",tissues,"_phased_tmml_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dMML')
  
  GRs <- c(GRs,GR.in)
  # dnme
  filename=paste(inDir,subjects,"_",tissues,"_phased_tnme_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dNME')
  GRs <- c(GRs,GR.in)
  # uc
  filename=paste(inDir,subjects,"_",tissues,"_phased_tpdm_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'UC')
  GRs <- c(GRs,GR.in)
  #Check if pvalue available
  GRs <- GRs[!is.na(GRs$pvalue)]
  # Add sample field
  GRs$Sample <- paste(tissues,"-",GRs$Subject)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  return(GRs)
}
#file_ends can be 
#c('dmml_pvals,dnme_pvals,uc_pvals')
#c('mml1,mml2,nme1,nme2')
read.alleleGR<-function(subjects,tissues,inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  #MML1
  GRs=GRanges()
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'MML')
  GR.in$Genome="1"
  GRs <- c(GRs, GR.in)
  #MML2
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'MML')
  GR.in$Genome="2"
  GRs <- c(GRs, GR.in)
  #NME1
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="1"
  GRs <- c(GRs, GR.in)
  #NME2
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="2"
  GRs <- c(GRs, GR.in)
  # Add sample field
  GRs$Sample <- paste(tissues,"-",subjects)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  GRs$K=GRs$'NA.1'
  GRs$'NA.1'=NULL
  return(GRs)
}
#Read in each bed file, for new method, no need to resize
import.ASMbed<-function(subjects,tissues,filename,pvalue=TRUE,Statistic,chrsOfInterest=paste("chr",1:22,sep="")){
  GR <- import.bedGraph(filename)
  #fit  bedGraph reads, import.bedGraph will remove 1 from start
  #Check if all files are 0 based or 1 based? Check on genome browser, UCSC: check SNP location (0 based)
  start(GR)=start(GR)-1
  GR_out=chr_check(GR)
  GR_out$ASM=NULL
  GR_out$Data=NULL
  GR_out$Subject <- subjects
  GR_out$tissue<-tissues
  GR_out$Statistic <- Statistic
  GR_out$Value<- GR$score
  #for differential analysis, make sure the GR_out are unique
  if(pvalue){
    GR_out$pvalue <- as.numeric(GR$NA.)
    GR_out <- GR_out[!duplicated(GR[,c()])]
    }
  else
    {GR_out$N <- GR$NA.} #not diff for new samples
  return(GR_out)
}
# reading in stat for each sample -------------------------------
#Note here we're using coverage cutoff=5 and boundary check == true
GR=import.subject('../downstream/data/bedGraph_diff/')
saveRDS(GR,GR_file)
GR_allele=import.subject('../downstream/data/bedGraph_allele/',calc='allele')
saveRDS(GR_allele,GR_allele_file)

gff_in=readRDS(gff_in_file)

# Sanity check output regions not in original gff file --------------------
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003",
           "STL011","H1","HuFGM02","112","149","150")
for(subj in subjects){
  cat(paste(subj,':\n'))
  cat(length(subsetByOverlaps(GR_allele[GR_allele$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-
        length(GR_allele[GR_allele$Subject==subj]),'\n')
  cat(length(subsetByOverlaps(GR[GR$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-length(GR[GR$Subject==subj]),'\n')
}
#Result all 0
