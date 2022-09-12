source('mainFunctions_sub.R')
NME_motif=fread("../downstream/output/human_analysis/Ken_motif/diff_test_TFBS_human_NME.csv")[order(abs(diff),decreasing=T)]
MML_motif=fread("../downstream/output/human_analysis/Ken_motif/diff_test_TFBS_human_MML.csv")[order(abs(diff),decreasing=T)]
intersect(NME_motif$motif,MML_motif$motif)
intersect(NME_motif$motif[1:100],MML_motif$motif[1:100])#No intersect