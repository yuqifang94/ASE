source('mainFunctions_sub.R')
dirIn='../downstream/data/coverage_human/'
CovOut=data.table()
for (fn in dir(dirIn,pattern=".cov")){
    bedIn = fread(paste0(dirIn,fn))
    if(nrow(bedIn)>0){
    colnames(bedIn)=c("chr","start","end","V4","V5","V6","count")
    sampleIn = gsub('\\..*|_paired|_single|_phased',"",fn)
    cov = bedIn[count!=0,list(coverage=mean(count))]

    genome = unlist(strsplit(fn,'\\.'))
    genome = genome[which(genome %in% c("all","genome1","genome2"))]
    cov$sample = sampleIn
    cov$genome = genome


    CovOut=rbind(CovOut,cov)
    }
}
saveRDS(CovOut,'../downstream/output/human_analysis/coverage.rds')

#bedtools coverage -f 1 -sorted -counts -a CpG_hg19.bed -b ../H1/bam/H1_ESC_paired_phased.sort.genome1.bam > H1_ESC_paired_phased.sort.genome1.bam.cov
#bedtools coverage -f 1 -sorted -counts -a CpG_hg19.bed -b ../H1/bam/H1_ESC_paired_phased.sort.genome2.bam > H1_ESC_paired_phased.sort.genome2.bam.cov