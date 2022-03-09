source('mainFunctions_sub.R')
hg19_CpG=getCpgSitesH19()
export.bed(hg19_CpG,"../downstream/output/human_analysis/CpG_hg19.bed")
mm10_CpG=getCpgSitesmm10()
export.bed(mm10_CpG,"../downstream/output/human_analysis/mm10_all_CpG.sort.bed")

