rm(list=ls())
source("mainFunctions_sub.R")
# reading in vcf files ----------------------------------------------------
#Linux for converting vcf to vcf.gz in ../downstream/data/vcfFiles/:
#for fn in *.vcf; do bgzip -c $fn > $fn.gz; done
#From vcf file, extract het CpG information
extractHetCpG<-function(vcfDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  cat('Processing subject:', sub,'\n')
  tt1=proc.time()[[3]]
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
  var_gr=chr_check(var_gr)
  plus_loc=as.character(Views(Hsapiens,GenomicRanges::shift(var_gr,1)))
  minus_loc=as.character(Views(Hsapiens,GenomicRanges::shift(var_gr,-1)))
  #get dinucleotide for ref, alt, plus and minus, find some examples region randomly: check if those match
  var_gr$REF_plus=paste0(var_gr$REF,plus_loc)
  var_gr$REF_minus=paste0(minus_loc,var_gr$REF)
  var_gr$ALT_plus=paste0(var_gr$ALT,plus_loc)
  var_gr$ALT_minus=paste0(minus_loc,var_gr$ALT)
  #get trinucleotide
  var_gr$REF_tri=paste0(minus_loc,var_gr$REF,plus_loc)
  var_gr$ALT_tri=paste0(minus_loc,var_gr$ALT,plus_loc)
  #Add CG content for genome1 and genome2 based on GT
  var_gr$genome1_plus=NA
  var_gr$genome1_minus=NA
  var_gr$genome1_tri=NA
  var_gr$genome2_plus=NA
  var_gr$genome2_minus=NA
  var_gr$genome2_tri=NA
  ####Fix the issue with genome file order does not equal to the ref/alt order, calculate the CG or het CG in genome 1 based on actual ref/alt order
  #Genome1
  var_gr$genome1_plus[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_plus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome1_minus[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_minus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome1_tri[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_tri[var_gr$GT %in% c("0/1","0|1")]
  
  var_gr$genome1_plus[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_plus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome1_minus[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_minus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome1_tri[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_tri[var_gr$GT %in% c("1/0","1|0")]
  #Genome2
  var_gr$genome2_plus[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_plus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome2_minus[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_minus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome2_tri[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_tri[var_gr$GT %in% c("0/1","0|1")]
  
  var_gr$genome2_plus[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_plus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome2_minus[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_minus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome2_tri[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_tri[var_gr$GT %in% c("1/0","1|0")]
  return(var_gr)
}


variant_HetCpG=lapply(subjects,function(x) extractHetCpG('../downstream/data/vcfFiles/',x)) 
names(variant_HetCpG)=subjects
saveRDS(variant_HetCpG,variant_HetCpG_file)