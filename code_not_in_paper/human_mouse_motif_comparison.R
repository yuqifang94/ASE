source('mainFunctions_sub.R')
mouse_motif_dir="../downstream/output/mouse_analysis/mouse_motif_Ken/"
mouse_motif_DNase=readRDS(paste0(mouse_motif_dir,"mouse_motif_DNase_median.rds"))
mouse_motif_control=readRDS(paste0(mouse_motif_dir,"mouse_motif_control_median.rds"))
mouse_motif_diff=mouse_motif_DNase-mouse_motif_control
mouse_motif_diff=mouse_motif_diff[,grepl("brain|heart|liver|limb|EFP",colnames(mouse_motif_diff))&(!grepl("P0",colnames(mouse_motif_diff)))]
human_motif_dir="../downstream/output/human_analysis/Ken_motif/homogeneous/"
human_motif_DNase=readRDS(paste0(human_motif_dir,"human_motif_DNase_median_NME.rds"))
human_motif_control=readRDS(paste0(human_motif_dir,"human_motif_control_median_NME.rds"))
human_motif_diff=human_motif_DNase-human_motif_control

mouse_brain=mouse_motif_diff[,grepl("brain",colnames(mouse_motif_diff))]
human_brain=human_motif_diff[,grepl("brain|Brain",colnames(human_motif_diff))]

length(intersect(names(which(rowMeans(mouse_brain)>0)),names(which(rowMeans(human_brain)>0))))#560

mouse_heart=mouse_motif_diff[,grepl("heart",colnames(mouse_motif_diff))]
human_heart=human_motif_diff[,grepl("Ventricle|Aorta|Atrium",colnames(human_motif_diff))]
length(intersect(names(which(rowMeans(mouse_heart)>0)),names(which(rowMeans(human_heart)>0))))#550


mouse_liver=mouse_motif_diff[,grepl("liver",colnames(mouse_motif_diff))]
human_liver=human_motif_diff[,grepl("Liver",colnames(human_motif_diff))]
length(intersect(names(which(rowMeans(mouse_liver)>0)),names(which(human_liver>0))))#496

mouse_intestine=mouse_motif_diff[,grepl("intestine",colnames(mouse_motif_diff))]
human_intestine=human_motif_diff[,grepl("Intestine|Colon",colnames(human_motif_diff))]
length(intersect(names(which(rowMeans(mouse_intestine)>0)),names(which(rowMeans(human_intestine)>0))))#580

mouse_stomach=mouse_motif_diff[,grepl("stomach",colnames(mouse_motif_diff))]
human_stomach=human_motif_diff[,grepl("Gastric|Esophagus",colnames(human_motif_diff))]
length(intersect(names(which(rowMeans(mouse_stomach)>0)),names(which(rowMeans(human_stomach)>0))))#581

mouse_lung=mouse_motif_diff[,grepl("Lung",colnames(mouse_motif_diff))]
human_lung=human_motif_diff[,grepl("Lung",colnames(human_motif_diff))]
length(intersect(names(which(rowMeans(mouse_lung)>0)),names(which(rowMeans(human_lung)>0))))#566

c(
    length(intersect(names(which(rowMeans(mouse_brain)>0)),names(which(rowMeans(human_brain)>0)))),
    length(intersect(names(which(rowMeans(mouse_heart)>0)),names(which(rowMeans(human_heart)>0)))),
    length(intersect(names(which(rowMeans(mouse_liver)>0)),names(which(human_liver>0))))
)/736

Reduce(intersect,list(
    intersect(names(which(rowMeans(mouse_brain)>0)),names(which(rowMeans(human_brain)>0))),
    intersect(names(which(rowMeans(mouse_heart)>0)),names(which(rowMeans(human_heart)>0))),
    intersect(names(which(rowMeans(mouse_liver)>0)),names(which(human_liver>0)))
))
mouse_non_matched=colnames(mouse_motif_diff)[!grepl("brain|heart|liver|Lung|intestine|stomach|NT|kidney|P0",colnames(mouse_motif_diff))]#12/46
human_non_matched=colnames(human_motif_diff)[!grepl("brain|Brain|Ventricle|Aorta|Atrium|Liver",colnames(human_motif_diff))]#36

mouse_motif_diff_non_match=mouse_motif_diff[,mouse_non_matched]
human_motif_diff_non_match=human_motif_diff[,human_non_matched]
intersect(names(which(rowMeans(mouse_motif_diff_non_match)>0)),names(which(rowMeans(human_motif_diff_non_match)>0)))#573

mouse_matched=colnames(mouse_motif_diff)[grepl("brain|heart|liver|Lung|intestine|stomach|NT|kidney|P0",colnames(mouse_motif_diff))]#12/46
human_matched=colnames(human_motif_diff)[grepl("brain|Brain|Ventricle|Aorta|Atrium|Liver",colnames(human_motif_diff))]#36


mouse_Ken=fread("../downstream/output/mouse_analysis/mouse_motif_Ken/mouse_motif_NME.csv",skip=1,header=T)
human_Ken=fread("../downstream/output/mouse_analysis/mouse_motif_Ken/human_NME_motif.csv",skip=1,header=T)
colnames(human_Ken)[1]="TF"
mouse_motif_diff_match=mouse_motif_diff[,mouse_matched]
human_motif_diff_match=human_motif_diff[,human_matched]
intersect(names(which(rowMeans(mouse_motif_diff_match)>0)),names(which(rowMeans(human_motif_diff_match)>0)))#559

