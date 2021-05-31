#' Draw methylation circle plot without SNP
#'
#' Draws CpG site methylation status as points, in reads containing a specific
#' CpG site. Generates one plot per bam file.
#'
#' @param cpgsite GRanges object containing a single CpG site location of
#'   interest
#' @param bamFile bismark bam file path
#' @param refFile fasta reference file path
#' @param pointSize Size of methylation circles. Default = 3.
#' @param dame (optional) GRanges object containing a region to plot
#' @param order Whether reads should be sorted by methylation status. Default=
#' False.
#' @param sampleName Plot title.
#' @param sampleReads Whether a subset of reads should be plotted.
#'   Default = FALSE.
#' @param numReads Number of reads to plot, if sampleReads is TRUE. Default = 20
#' @return Plot
#' @examples
#' DATA_PATH_DIR <- system.file('extdata', '.', package = 'DAMEfinder')
#' get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
#' bam_files <- get_data_path('NORM1_chr19_trim.bam')
#' sample_names <- 'NORM1'
#' #reference_file
#' suppressPackageStartupMessages({library(BSgenome.Hsapiens.UCSC.hg19)})
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' seqnames(genome) <- gsub("chr","",seqnames(genome))
#' dna <- DNAStringSet(genome[[19]], use.names = TRUE)
#' names(dna) <- 19
#'
#' cpg <- GenomicRanges::GRanges(19, IRanges::IRanges(292082, width = 1))
#' methyl_circle_plotCpG(cpgsite = len,
#'  bamFile = bam_files,
#'  refFile = dna)
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @import ggplot2
#'
#' @export
source('DAME/split_reads.R')
library(stringr)
methyl_circle_plotCpG <- function(cpgsite = cpgsite, bamFile = bamFile, 
    pointSize = 3, refFile = refFile, dame = NULL, order = FALSE, 
    sampleName = NULL, sampleReads = FALSE, numReads = 20,softclip=5,lwd=0.5) {
    #Here adjust for single end
    alns.pairs <- GenomicAlignments::readGAlignmentsList (bamFile, 
        param = Rsamtools::ScanBamParam(tag = c("MD", "XM", "XR", 
            "XG"), which = cpgsite), use.names = TRUE)
    
    
    if (sampleReads) {
        ran.names <- sample(names(alns.pairs), numReads)
        alns.pairs <- alns.pairs[names(alns.pairs) %in% ran.names]
    }
    alns.pairs_rf=list()
    for(pn in names(alns.pairs)){
      for(i in 1:length(alns.pairs[[pn]])){
        alns.pairs_rf[[paste0(pn,'.r',i)]]=alns.pairs[[pn]][i]
        
      }
      
      
    }
    alns.pairs=alns.pairs_rf
    alns.pairs_rf=NULL
    
    alns <- unlist(alns.pairs)
    
    #### get reference and CpG
    #### positions####-------------------------------------
    
    # Get limits for plotting
    
    if (is.null(dame)) {
        left <- min(start(alns))
        right <- max(end(alns))
        window <- GRanges(seqnames(cpgsite), IRanges(left, right))
    } else {
        window <- dame
        left <- start(dame)
        right <- end(dame)
    }
    
    message("Reading reference")
    # open reference seq to get correct CpG locations within that
    # reference chunk
    if (typeof(refFile) == "character") {
        fa <- open(Rsamtools::FaFile(refFile, index = paste0(refFile, 
            ".fai")))
        dna <- Rsamtools::scanFa(fa, param = window)
    } else {
      print(window)
        dna <- refFile[window]
    }
    
    cgsite <- stringr::str_locate_all(dna, "CG")[[1]][, 1]  #also look at GpCs?
    
    if (length(cgsite) < 1) {
        stop("No CpG sites in these reads")
    }
    
    mepos <- cgsite/Biostrings::nchar(dna)  #location of CpG sites
    
    ##### Use MD tag from bam to extract methylation status for each
    ##### site detected above ####----------------------------
    message("Getting meth state per read-pair")

    conversion <- vapply(alns.pairs, function(x) {
        
        # change C locations for G locations if reference context is
        # different
        if (S4Vectors::mcols(x)$XG[1] == "GA") {
            cgsite <- cgsite + 1
        }
        
        # fuse pair info in one row
        #conversion <- matrix(NA, 1, length(cgsite))
        
        # a <- mcols(x)$MD
        read.start <- start(x) - left + 1-(softclip)  #start of read
      
            
            # Get the locations from the specific read: Count along the
            # numerical MDtag until you reach the actual CpG start
            #top and bottom 5bp is softclipped
            this.mepos <- cgsite - read.start + 1
            
            # Fill in conversion table with meth state for each read 21
            # is unmethylated 19 is methylated 0 is not in read
            #z is unmet, Z is met
            return(c("U","M",NA)[match(substring(mcols(x)$XM,this.mepos,this.mepos), c("z","Z","") )])
      
        
    }, character(length(cgsite)))
    
    # Remove reads without CpG sites
    #print(conversion)
    rem <- colSums(is.na(conversion)) >= dim(conversion)[1]
    conversion <- conversion[, !rem]
    
    alns.pairs <- alns.pairs[names(alns.pairs) %in% colnames(conversion)]
    
    #### Get CpG of interest in window
    #### ####------------------------------------
    
    cpg.start <- start(cpgsite) - left + 1
    
    #### plot
    #### ####-------------------------------------------------------------
    
    message("Plotting")
    #Extend read 1 bp to start to include the reverse compliment read, they're 1 bp offf
    xstart <- vapply(names(alns.pairs), function(x) {
        start(alns.pairs[[x]]) - left + 1-softclip-1
    }, FUN.VALUE = double(1))
    xstart[xstart<0]=0
    xend <- vapply(names(alns.pairs), function(x) {
        end(alns.pairs[[x]]) - left + 1+softclip
    }, FUN.VALUE = double(1))
    xend[xend>width(dame)]=width(dame)
    start_order=names(xstart)[order(xstart,decreasing=F)]
    #Looking for adjacent reads
    
  
    # data for points
    d <- data.table(CpG = rep(cgsite, length(alns.pairs)), 
        reads = rep(seq(from = 1, 
        to = length(alns.pairs), by = 1), each = length(cgsite)), 
        value = as.vector(conversion[,start_order]),
        read_name=rep(start_order, each = length(cgsite)))
    
    d2 <- data.table(xstart=xstart[d$read_name], xend=xend[d$read_name], reads=d$reads,CpG=d$CpG )
    d2$line_s=as.numeric(NA)
    d$line_s=as.numeric(NA)
    #Due to the nature of previous sorting, look from end to begining
    line_num=1
    for(read_in in 1:max(d2$reads)){
      if(is.na(unique(d2[reads==read_in]$line_s))){
        #Find the next start greater than it's end
        x_end=unique(d2[reads==read_in]$xend)
        #Try look for unassigned ones
        same_line=sort(unique(d2[xstart>x_end&is.na(line_s)]$reads),decreasing=F)
        
        if(length(same_line)>0){
          sl_all=c(read_in)
          for(sl in same_line){
            
            if (unique(d2[reads==sl]$xstart)>x_end){
              
              sl_all=c(sl_all,sl)
              x_end=unique(d2[reads==sl]$xend)
            }
            
          }
          # #3 senarios: 1. all NA, assign new ones, 2. some NA, assign existed ones (should not appear)3. all not NA, skip
          # line_s_unq=is.na(d2[reads%in%sl_all]$line_s)
       
          #Check if there's already reads in the same line
          #If there's not add new line
          
          # if(any(line_s_unq)){
            # reads_NA=unique(d2[reads%in%sl_all][line_s_unq]$reads)
            d2[reads%in%sl_all]$line_s=line_num
            d[reads%in%sl_all]$line_s=line_num
            line_num=line_num+1
          # }
        }else{
          #If there's no line num here, add one
          
          
            d2[reads==read_in]$line_s=line_num
  
            d[reads==read_in]$line_s=line_num

            line_num=line_num+1
          }
          
      }
      
      
    }
    
   print(width(dame))
    print(pointSize)
    ggplot() + 
        scale_shape_identity() + 
        theme_void() + 
        geom_segment(data = d2, aes(x = xstart, y = -line_s, xend = xend,
            yend = -line_s), colour = "grey", size = lwd) +
        geom_point(data = d, aes(x = CpG, y = -line_s, group=value,fill =value,color=value),size = pointSize,stroke =0) +
      scale_fill_manual(values=c("U"="blue","M"="red","NA"=NA))+
      scale_color_manual(values=c("U"="blue","M"="red","NA"=NA))+
        # geom_point(aes(x = cpg.start, y = 0), shape = 24, size = 3, 
        #     fill = "green") + 
        guides(color = FALSE) + 
       xlim(c(0,width(dame)))+
       # ggtitle(sampleName)+
      theme(legend.position = "none")
    
}

