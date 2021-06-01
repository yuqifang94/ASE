# counting ASM ------------------------------------------------------------

#Examples:
#dNME example1:
subsetByOverlaps(GR_merge[GR_merge$Sample=="ectoderm_paired - HUES64"],GRanges(seqnames="chr14",IRanges(start=104552150,end=104552495)))
#dNM example2:
subsetByOverlaps(GR_merge[GR_merge$Sample=="Psoas_Muscle_single - STL003"],GRanges(seqnames="chr17",IRanges(start=33750066,end=33770266)))
#dMML_example1
subsetByOverlaps(GR_merge[GR_merge$Sample=="stem_27_undifferentiated_paired - HUES64"],GRanges(seqnames="chr11",IRanges(start=2720817,end=2721033)))
#dMML_example2
subsetByOverlaps(GR_merge[GR_merge$Sample=="endoerm_27_paired - HUES64"],GRanges(seqnames="chr20",IRanges(start=32308087,end=32308287)))
#Find examples for dNME
GR_merge_dNME=GR_merge[GR_merge$N>=4&GR_merge$dNME_pval<=pval_cutoff]
GR_merge_dNME=GR_merge_dNME[order(GR_merge_dNME$dNME,decreasing=T)]
GR_merge_dNME_dt=convert_GR(GR_merge_dNME,direction="DT")
#Example 1: dNME
#chr8:3,770,843-3,771,336
#HUES64 stem
#chr14:104,566,670-104,566,870
#chr4:7,308,694-7,309,119
#chr12:125,583,449-125,583,649
#chr11:103,242,820-103,243,199

#Examples mpre dNME
#chr6:20,481,578-20,481,925, edoderm
#     chr2 100438919-100439371,ectoderm  
#chr3   76733036-76733603 ectoderm
#   chr1   18993859-18994059 , endoderm
# chr9     4023395-4023595   esc

Imprinted_Genes <- as.data.frame(read_excel("../downstream/input/human_analysis/imprinting_ASE/Imprinted Genes.xlsx"))
GR_merge_dMML=GR_merge[GR_merge$genes_promoter%in%Imprinted_Genes$Gene&GR_merge$dMML_pval<=pval_cutoff]
#chr6 144329267-144329363 149
#chr6:144,329,572-144,329,772 150
#chr6 144329572-144329772

#Mouse analysis
NME_in=readRDS('../downstream/input/mouse_analysis/NME_agnostic_mouse_all_merged.rds')
MML_in=readRDS('../downstream/input/mouse_analysis/MML_agnostic_mouse_all_merged.rds')
gtf=fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table=F)
promoter_in=gtf <- gtf[gtf[,3]=='gene',]
type <- sub('\".*','',sub('.*gene_type \"','',gtf[,9]))
gtf <- gtf[type=='protein_coding',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gene_body <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
gene_body$gene=gn
olap=findOverlaps(MML_in,gene_body)
MML_in$gene=NA
MML_in$gene[queryHits(olap)]=gene_body$gene
olap=findOverlaps(NME_in,gene_body)
NME_in$gene=NA
NME_in$gene[queryHits(olap)]=gene_body$gene

MML_in_gata4=subsetByOverlaps(MML_in,gene_body[gene_body$gene=="Gata4"])
MML_in_gata4=sort(MML_in_gata4)

NME_in_gata4=subsetByOverlaps(NME_in,gene_body[gene_body$gene=="Gata4"])
NME_in_gata4=sort(NME_in_gata4)
MML_in_gata4[MML_in_gata4$N>=4&MML_in_gata4$MML>=0.96&MML_in_gata4$tissue%in%c("heart","forebrain","limb")]#Forebrain E16.5: chr14 63212488-63212737

NME_in_gata4[NME_in_gata4$N>=4&NME_in_gata4$NME>=0.97&NME_in_gata4$tissue%in%c("heart","forebrain","limb")]#forebrain E11.5
