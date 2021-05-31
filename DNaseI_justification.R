source('mainFunctions_sub.R')
DNAase=readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp.rds')
enhancer=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')
gtf=fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table=F)
promoter_in=gtf <- gtf[gtf[,3]=='gene',]
type <- sub('\".*','',sub('.*gene_type \"','',gtf[,9]))
gtf <- gtf[type=='protein_coding',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
names(gr) <- gn
promoters <- promoters(gr,upstream=0,downstream=1)

length(subsetByOverlaps(enhancer,DNAase))/length(enhancer)
length(subsetByOverlaps(promoters,DNAase,maxgap = 2000))/length(promoters)
