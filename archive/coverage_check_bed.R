library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(exomeCopy)
library(Mus.musculus)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#human
bed_out=GRanges()
for (fn in dir(pattern='*.cov')){
  bed_in=import.bedGraph(fn)
  sample=gsub('_phased.all.sorted.bam.agnostic.cov','',fn)
  bed_in$coverage=bed_in$NA.2
  bed_in$Sample=sample
  bed_in=bed_in[seqlevels(bed_in)%in%1:22]
  mcols(bed_in)=mcols(bed_in)[,c('coverage','Sample')]
  bed_out=c(bed_out,bed_in)
  
  
}
#mouse use all the CpG sites
chrs_mm10 <- names(Mmusculus)[1:19]#2276796
cgs_mm10 <- lapply(chrs_mm10, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr_mm10 <- do.call(c, lapply(1:19, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs_mm10[[x]], width = 2)))) #use first location
export.bed(cpgr_mm10,'../downstream/output/mm10_all_CpG.bed')
#reading in resulted bed file

region_in=readRDS('mm10_PRC.rds')
region_in=readRDS('mm10_DNase.rds')
blacklist=import.bed('mm10.blacklist.bed.gz')

import_bed=mclapply(dir(pattern=c('sort.dup.bam.agnostic.cov')),function(fn){
  tt1=proc.time()[[3]]
  cat('start processing:',fn,'\n')
  bed_in=import.bedGraph(fn)
  # sample=gsub('_merge.sort.bam.agnostic.cov','',fn)
  # sample=gsub('_merged.sort.bam.agnostic.cov','',fn)
  sample=gsub('.sort.dup.bam.agnostic.cov','',fn)
  sample=gsub('merged','',sample)
  bed_in$coverage=bed_in$NA.2
  bed_in$Sample=sample
  #bed_in=bed_in[seqlevels(bed_in)%in%paste('chr',1:22,sep='')]
  mcols(bed_in)=mcols(bed_in)[,c('coverage','Sample')]
 #  # jpeg(paste("cov_hist_all/",sample,".jpg",sep=''))
 #  # hist(log10(bed_in$coverage),xlab='coverage',main=sample)
 #  # dev.off()
 #  bed_in=subsetByOverlaps(bed_in,region_in)
 #  #bed_in_olap=findOverlaps(bed_in,blacklist)
 # # bed_in=bed_in[-queryHits(bed_in_olap)]
 #  bed_out=rbind(bed_out,data.frame(sample=sample,cov_mean=mean(bed_in$coverage),cov_sd=sd(bed_in$coverage)))
  cat('finish processing',sample,'in',proc.time()[[3]]-tt1, '\n')
  return(bed_in)
  
},mc.cores=24)

saveRDS(import_bed,'coverage_ls_dedup.rds')
coverage_summary=mclapply(import_bed,function(x) data.frame(sample=unique(x$Sample),mean=mean(x$coverage),sd=sd(x$coverage)),mc.cores=24)
region_in=readRDS('mm10_DNase.rds')
coverage_summary_DNase=mclapply(import_bed,function(x) {x=subsetByOverlaps(x,region_in)
                              data.frame(sample=unique(x$Sample),mean=mean(x$coverage),sd=sd(x$coverage))},mc.cores=24)
region_in=readRDS('mm10_PRC.rds')
coverage_summary_PRC=mclapply(import_bed,function(x) {x=subsetByOverlaps(x,region_in)
data.frame(sample=unique(x$Sample),mean=mean(x$coverage),sd=sd(x$coverage))},mc.cores=24)


hg19_gr_out=readRDS('../downstream/output/hg19_gr_out.rds')
seqlevels(hg19_gr_out)=paste('chr',seqlevels(hg19_gr_out),sep='')
genomic_features=readRDS(genomic_features_file)
hg19_gr_out_islands=subsetByOverlaps(hg19_gr_out,genomic_features$`CpG island`)#29
hg19_gr_out_shores=subsetByOverlaps(hg19_gr_out,genomic_features$`CpG shore`)#34.29
hg19_gr_out_open_sea=subsetByOverlaps(hg19_gr_out,genomic_features$`CpG open sea`)#38
hg19_gr_out_shelf=subsetByOverlaps(hg19_gr_out,genomic_features$`CpG shelf`)#33.5
cov_mean_out=data.frame()
for (sp in unique(hg19_gr_out$Sample)){
 cov_mean_out=rbind(cov_mean_out,
                    data.frame(islands=mean(hg19_gr_out_islands$coverage[hg19_gr_out_islands$coverage>0&hg19_gr_out_islands$Sample==sp]),#30
                    shores=mean(hg19_gr_out_shores$coverage[hg19_gr_out_shores$coverage>0&hg19_gr_out_shores$Sample==sp]),#35
                    shelves=mean(hg19_gr_out_shelf$coverage[hg19_gr_out_shelf$coverage>0&hg19_gr_out_shelf$Sample==sp]),#35
                    open_seas=mean(hg19_gr_out_open_sea$coverage[hg19_gr_out_open_sea$coverage>0&hg19_gr_out_open_sea$Sample==sp])#40
 ))
}

coverage_feature_df=rbind(data.frame(coverage=hg19_gr_out_islands$coverage,feature='CpG islands'),
                          data.frame(coverage=hg19_gr_out_shores$coverage,feature='CpG shores'),
                          data.frame(coverage=hg19_gr_out_shelf$coverage,feature='CpG shelf'),
                          data.frame(coverage=hg19_gr_out_open_sea$coverage,feature='CpG open sea'))
                          
jpeg('../downstream/output/genomic_feature.jpg')                          
boxplot(coverage~feature,data=coverage_feature_df, main="",
        xlab="Genomic features", ylab="CpG coverage",outline=F)                        
                          
dev.off()
