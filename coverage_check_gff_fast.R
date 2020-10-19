library(bsseq)
library(rtracklayer)
library(Rsamtools)
library(rlist)

bam_gr<-function(aln){
  aln_sub=aln[c('qname','rname','strand','pos')]
  aln_sub=as.data.frame(aln_sub)
  aln_sub$start=aln_sub$pos
  aln_sub$end=aln[['qwidth']]+aln_sub$start
  aln_sub$strand='*'
  colnames(aln_sub)=c('qname','seqnames','strand', 'pos','start','end')
  makeGRangesFromDataFrame(aln_sub)
}
cov_bam<-function(CpG_gff,bam){countOverlaps(CpG_gff,bam)}
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003","STL011","H1","HuFGM02","112","149","150")
gff_in=readRDS('../downstream/input/hetCpG_gff_new.rds')
#mouse
p_tag=ScanBamParam(what=c('qname',"rname", "strand", "pos", "qwidth"))

coverage_sample_out=list()
#coverage_sample_out=readRDS('.../downstream/coverag_sample_out_temp.rds')
bam_dir='../../../../../../allele_specific/bamfiles'
for (subj in subjects){
  files_in=unique(gsub('.bai','',dir(paste(bam_dir,subj,sep='/'))))
  for(fn in files_in){
    bam_in_fn=paste(bam_dir,subj,fn,sep='/')
    tt1=proc.time()[[3]]
    cat('Reading in: ',bam_in_fn,'\n')
    bam_in=scanBam(BamFile(bam_in_fn),param=p_tag)
    bam_in=bam_gr(bam_in[[1]])
    seqlevels(bam_in)=paste('chr',seqlevels(bam_in),sep='')
    coverage_sample_gr=gff_in[[subj]]
    coverage_sample_gr$coverage=NA
    coverage_sample_gr$bam=fn
    coverage_sample_gr$coverage=countOverlaps(coverage_sample_gr,bam_in)
    coverage_sample_out[[fn]]=coverage_sample_gr
    cat('Finish ',fn,'in',proc.time()[[3]]-tt1,'\n')
  }
}
###########Coverage#################
coverage_all=readRDS('../downstream/output/coverage_sample_fast.rds')
strsplit(names(coverage_all),'\\.')[[1]]
sp_name=unlist(lapply(strsplit(names(coverage_all),'\\.'),function(x) x[[1]]))
cov_sample_out=GRanges()
for (sp in unique(sp_name)){
  cov_in=coverage_all[which(sp_name==sp)]
  gr_out=cov_in[[1]]
  gr_out$bam=NA
  gr_out$Sample=sp
  gr_out$coverage=NA
  cov_in_name=unlist(lapply(strsplit(names(cov_in),'\\.'),function(x) x[[3]]))
  gr_out$genome1_cov=cov_in[[which(cov_in_name=='genome1')]]$coverage
  gr_out$genome2_cov=cov_in[[which(cov_in_name=='genome2')]]$coverage
  cov_sample_out=c(cov_sample_out,gr_out)
}
saveRDS(cov_sample_out,'../downstream/output/coverage_sample_fast.rds')
cov_sample_out=readRDS('../downstream/output/coverage_sample_fast.rds')
GR_merge=readRDS("../downstream/output/gr_merge_run3_promter_enhancer_density.rds")
coverage_cutoff=data.frame(cutoff=1:20,N=NA)

coverage_cutoff_agg=data.frame()
coverage_cutoff_agg_all=data.frame()
for (cutoff in coverage_cutoff$cutoff){
  coverage_cutoff$N[coverage_cutoff$cutoff==cutoff]=sum(cov_sample_out$genome1_cov>=cutoff & cov_sample_out$genome2_cov>=cutoff)
  coverage_cutoff_sp=data.frame(regions_cov=cov_sample_out$genome1_cov>=cutoff & cov_sample_out$genome2_cov>=cutoff,
                                sp=cov_sample_out$Sample)
  coverage_cutoff_sp_agg=aggregate(coverage_cutoff_sp$regions_cov,by=list(coverage_cutoff_sp$sp), sum)
  colnames(coverage_cutoff_sp_agg)=c('Sample','N')
  coverage_cutoff_sp_agg$cutoff=cutoff
  coverage_cutoff_agg=rbind(coverage_cutoff_agg,coverage_cutoff_sp_agg)
  coverage_cutoff_sp$cutoff=cutoff
  #coverage_cutoff_agg_all=rbind(coverage_cutoff_agg_all,coverage_cutoff_sp[coverage_cutoff_sp$regions_cov,])
}
ratio_raw=length(GR_merge)/coverage_cutoff$N[5]
coverage_cutoff_agg$N_estimate=coverage_cutoff_sp_agg$N*ratio_raw

coverage_cutoff$type='raw'
coverage_cutoff=rbind(coverage_cutoff,coverage_cutoff)
coverage_cutoff$N[21:nrow(coverage_cutoff)]=coverage_cutoff$N[1:20]*(length(GR_merge)/coverage_cutoff$N[4])
coverage_cutoff$type[21:nrow(coverage_cutoff)]='estimate'
ggplot(coverage_cutoff,aes(x=cutoff,y=N,color=type))+geom_line()+
  xlab('coverage cutoff')+ylab('number of regions')+geom_text(aes(label=round(N)),check_overlap = TRUE)
jpeg('coverage.jpeg')
ggplot(coverage_cutoff_agg,aes(x=cutoff,y=N,color=Sample))+geom_point()+
  xlab('coverage cutoff')+theme(legend.position = 'bottom')
dev.off()

coverage_cutoff_agg$N[coverage_cutoff_agg$Sample=='STL003_Adrenal_Gland_single_phased'&coverage_cutoff_agg$cutoff==10]/
  coverage_cutoff_agg$N[coverage_cutoff_agg$Sample=='STL003_Adrenal_Gland_single_phased'&coverage_cutoff_agg$cutoff==8]
