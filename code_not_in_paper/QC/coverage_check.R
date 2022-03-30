library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(exomeCopy)
library(Mus.musculus)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#mouse use all the CpG sites
chrs_mm10 <- names(Mmusculus)[1:19]#2276796
cgs_mm10 <- lapply(chrs_mm10, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr_mm10 <- do.call(c, lapply(1:19, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs_mm10[[x]], width = 2)))) #use both location
export.bed(cpgr_mm10,'../downstream/output/mouse_analysis/QC/mm10_all_CpG.bed')
#In bam file folder: 

#for fn in mm10*all.sort.dup.bam; do sbatch coverage_calc.sh mm10_all_CpG.sort.bed $fn coverage_file_all_CG/$fn.agnostic.cov; done


# Human coverage ----------------------------------------------------------

chrs <- names(Hsapiens)[1:25]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:25, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2)))) #use first location
seqlevels(cpgr)=gsub('chr','',seqlevels(cpgr))
export.bed(sort(cpgr),'../downstream/output/human_analysis/QC/hg19_all_CpG.bed')
#In coverage folder: /ibox/afeinbe2/yfang/human_analysis_bam_all/coverage
#for fn in ../bam_all/*phased.all.sorted.bam; do sbatch coverage_calc.sh hg19_all_CpG.bed $fn "${fn/\.\.\/bam_all\//}".agnostic.cov; done
#Need to rerun human coverage analysis
bed_out=data.table()
cov_dir='../downstream/data/coverage_human/'
for (fn in dir(cov_dir,pattern='_phased.sort.all.bam.agnostic.cov')){
  bed_in=fread(paste0(cov_dir,fn))[,list(V1,V2,V3,V7)]
  colnames(bed_in)=c("chr","start","end","coverage")
  if(nrow(bed_in)>0){
  sample=gsub('_phased.sort.all.bam.agnostic.cov','',fn)

  bed_in$Sample=sample
  bed_in=bed_in[chr%in%1:22]
  bed_out=rbind(bed_out,bed_in)
  }
  else{
    cat("Need redo coverage analysis on:", fn,"\n")

  }
  
}
bed_out$region=paste0(bed_out$chr,bed_out$start,bed_out$end)
saveRDS(bed_out,'../downstream/output/human_analysis/QC/human_coverage.rds')
#STL001_Psoas_Muscle_single_phased.all.sorted.bam.agnostic.cov
#STL001_Right_Ventricle_single_phased.all.sorted.bam.agnostic.cov
#Mouse coverage
cov_dir='../downstream/data/coverage_mouse/'
bed_out=GRanges()
for (fn in dir(cov_dir,pattern='*.cov')){
  bed_in=import.bedGraph(paste0(cov_dir,fn))
  sample=gsub('_all.sort.dup.bam.agnostic.cov','',fn)
  bed_in$coverage=bed_in$NA.2
  bed_in$Sample=sample
  bed_in=bed_in[seqlevels(bed_in)%in%paste0('chr',1:19)]
  mcols(bed_in)=mcols(bed_in)[,c('coverage','Sample')]
  bed_out=c(bed_out,bed_in)
  
}
saveRDS(bed_out,'../downstream/output/mouse_analysis/mouse_coverage.rds')
#Mbias
#for fn in ../bam_all/{149,150,HuFGM02,H1}*.bam; do echo sbatch mbias.sh $fn $PWD ~/data/yfang/referenceGenome/hg19_Arioc/; done
#ASM coverage
#for fn in ../bam_asm/*.bam; do sbatch coverage_calc.sh hg19_all_CpG.bed $fn "${fn/\.\.\/bam_asm\//}".agnostic.cov; donebed_out=GRanges()
cov_dir='../downstream/data/coverage_human/'
#Obtain all CpGs
CpGs = import.bed('../downstream/output/human_analysis/QC/hg19_all_CpG.bed')
bed_out=data.table(CpGs=convert_GR(CpGs,direction="DT")$region)
for (fn in dir(cov_dir,pattern='sort.genome')){
  bed_in=import.bedGraph(paste0(cov_dir,fn))
  sample=gsub('merged_sort.*|_phased.sort.*','',fn)
  genome = unlist(strsplit(fn,"\\."))
  genome=genome[grepl("genome",genome)]
  sample = paste0(sample,gsub("genome","_",genome))
  cat("Processing: ",fn,"\n")
  bed_in$coverage=bed_in$NA.2
  bed_in=bed_in[seqlevels(bed_in)%in%1:22]
  bed_in=convert_GR(bed_in,direction="DT")
  bed_out[[sample]]=bed_in[match(bed_out$CpGs,region)]$coverage
  
}
saveRDS(bed_out,'../downstream/output/human_analysis/QC/human_coverage_asm.rds')
#Estimate average coverage per region
gff_in=readRDS(gff_in_file)
bed_out_gr=CpGs
seqlevels(bed_out_gr)=paste0("chr",seqlevels(bed_out_gr))
olap_gff=findOverlaps(bed_out_gr,gff_in,minoverlap=2)
gff_in_dt=convert_GR(gff_in,direction="DT")
bed_out=bed_out[queryHits(olap_gff)]
gff_in_dt=gff_in_dt[subjectHits(olap_gff)]
bed_out_na=bed_out
bed_out_na[bed_out_na==0] = NA
bed_out_na=bed_out_na[,!(names(bed_out_na)%in% c("Subject_gff","regions","CpGs")|grepl("GM12878",names(bed_out_na))),with=F]
gff_in_dt$mean_cov=rowMeans(bed_out_na,na.rm=T)
#Sum coverage with no zero terms
gff_in_dt_mean=gff_in_dt[,list(mean_cov_mean=mean(mean_cov,na.rm=T)),by=list(region,N,Subject)]

saveRDS(gff_in_dt_mean,'../downstream/output/human_analysis/QC/human_coverage_asm_gff.rds')
gff_in_dt_mean=gff_in_dt_mean[,list(mean_cov_mean=mean(mean_cov_mean,na.rm=T)),by=list(region,N)]
gff_in_dt_mean_plot = rbind(data.table(N_CpG = gff_in_dt_mean[mean_cov_mean>=5]$N,cov_cutoff = "5"),
                            data.table(N_CpG = gff_in_dt_mean[mean_cov_mean>=10]$N,cov_cutoff = "10"))
gff_in_dt_mean_plot=gff_in_dt_mean_plot[,.N,by=list(N_CpG,cov_cutoff)]        
gff_in_dt_mean_plot[,freq:=N/sum(N),by=list(cov_cutoff)]        
gff_in_dt_mean_plot$N_CpG = factor(gff_in_dt_mean_plot$N_CpG,levels = unique(sort(gff_in_dt_mean_plot$N_CpG)))
pdf("../downstream/output/human_analysis/QC/N_dist_coverage.pdf",width=10)
ggplot(gff_in_dt_mean_plot,aes(x=N_CpG,y=freq,fill=cov_cutoff))+geom_bar(stat="identity",position=position_dodge())+
  ylab("frequency")+xlab("N CpG")+ theme(legend.position="bottom") + ylim(c(0,0.3))+
  guides(fill=guide_legend(title="Coverage cutoff"))
dev.off()
