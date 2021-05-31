source('mainFunctions_sub.R')
#Get enhancer and promoter regions
enhancer=readRDS('../downstream/output/bin_enhancer.rds')
gtf=fread('../downstream/input/grcm38.gtf',data.table=F)
promoter_in=gtf <- gtf[gtf[,3]=='gene',]
type <- sub('\".*','',sub('.*gene_type \"','',gtf[,9]))
gtf <- gtf[type=='protein_coding',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
names(gr) <- gn
tss <- promoters(gr,upstream=2000,downstream=2000)
enhancer_ZEB2=enhancer[enhancer$`Target Gene`=="Zeb2"]
promoter_Zeb2=tss["Zeb2"]
UC_merge=readRDS('../downstream/output/UC_merge_max_loc_ft_N17.rds')
UC_merge=UC_merge$forebrain
olap=findOverlaps(convert_GR(rownames(UC_merge)),c(enhancer_ZEB2,promoter_Zeb2))
UC_merge_Zeb2=UC_merge[queryHits(olap),]
UC_merge_Zeb2_gr=convert_GR(rownames(UC_merge_Zeb2))
mcols(UC_merge_Zeb2_gr)=UC_merge_Zeb2

#Get continous ones
col_n=colnames(mcols(UC_merge_Zeb2_gr))
mcols(UC_merge_Zeb2_gr)=mcols(UC_merge_Zeb2_gr)[,grepl(paste(c("E10.5-E11.5",
                                                               "E11.5-E12.5",
                                                               "E12.5-E13.5",
                                                               "E13.5-E14.5",
                                                               "E14.5-E15.5",
                                                               "E15.5-E16.5"),collapse="|"),
                                                             col_n)]
olap=findOverlaps(UC_merge_Zeb2_gr,c(enhancer_ZEB2))
UC_merge_Zeb2_gr$region_type="promoter"
UC_merge_Zeb2_gr$region_type[queryHits(olap)]="enhancer"
UC_merge_Zeb2_gr=UC_merge_Zeb2_gr[order(rowMaxs(as.matrix(mcols(UC_merge_Zeb2_gr)[,grepl("UC",colnames(mcols(UC_merge_Zeb2_gr)))])),decreasing=T)]
csv_in=fread('../downstream/input/ts_cluster_0_1/forebrain.csv')
#Look for visualization method
UC_merge_Zeb2_gr$tissue_specific=FALSE
UC_merge_Zeb2_gr$tissue_specific[convert_GR(UC_merge_Zeb2_gr,direction="DT")$region %in% csv_in$region]=TRUE
UC_merge_Zeb2_gr_sub=UC_merge_Zeb2_gr[rowMaxs(as.matrix(mcols(UC_merge_Zeb2_gr)[,grepl("UC",colnames(mcols(UC_merge_Zeb2_gr)))]))>=0.1]
