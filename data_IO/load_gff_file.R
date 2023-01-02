rm(list=ls())
source("mainFunctions_sub.R")
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
#Read in gff file
readAllGff <- function(inDir,subjects,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Loop over all GFF files
  outGR <- GRanges()
  #subjects <- c("H9","HUES64","skin03","HuFGM02","STL001","STL002","STL003","STL011")
  for (sub in subjects) {
    # Import GFF file
    cat('importing',sub,'\n')
    tmpGR <- readSubGff(inDir,sub,chrsOfInterest)
    start(tmpGR)=start(tmpGR)
    # Append to output GR
    outGR <- c(outGR,tmpGR)
  }
  
  # Return GR with all haplotypes
  return(outGR)
  
}
# reading in analyzed regions for each sample -----------------------------
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
saveRDS(gff_in,gff_in_file)
