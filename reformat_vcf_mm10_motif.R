library(data.table)
 library(rtracklayer)
 library(BSgenome.Mmusculus.UCSC.mm10)
vcf_in=fread('../downstream/output/mm10_SNP.txt',header=F,stringsAsFactors=F,sep='\t')
vcf_in=vcf_in[,c(2,3,5)]
vcf_in=makeGRangesFromDataFrame(vcf_in,seqnames.field="V2",start.field="V3",end.field="V3",keep.extra.columns=T)
vcf_in$REF=sub('/.*','',vcf_in$V5)
vcf_in$ALT=sub('.*/','',vcf_in$V5)
vcf_in$V5=NULL
vcf_in$REF_BS=as.character(Views(Mmusculus,vcf_in))
attributes(vcf_in)$genome.package="BSgenome.Mmusculus.UCSC.mm10"
#Similar SNP size
bedGraph_in=import.bedGraph('/home-4/yfang27@jhu.edu/scratch_feinberg/yfang/allele_agnostic_mouse_all/CPEL/CPEL_agnostic/cpelasm/cpelasm/BL6DBA_Epiblast_merged_paired_phased_tnme_pvals.bedGraph')
vcf_in=subsetByOverlaps(vcf_in,bedGraph_in,maxgap=100)
chrsOfInterest=paste("chr",c(1:19,"X","Y"),sep="")
mm10<-getBSgenome("mm10")
genome.seqinfo <- seqinfo(mm10)
genome.seqinfo <- genome.seqinfo[chrsOfInterest]
vcf_in <-vcf_in[seqnames(vcf_in) %in% chrsOfInterest]
genome(vcf_in) <- genome(genome.seqinfo)
seqlevels(vcf_in) <- seqlevels(genome.seqinfo)
seqlengths(vcf_in) <- seqlengths(genome.seqinfo)
saveRDS(vcf_in,'../downstream/output/SNP_DBA.rds')