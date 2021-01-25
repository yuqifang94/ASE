source('mainFunctions_sub.R')

#Get the tissue x region
#get nme and mml
mml <- readRDS('../downstream/input/MML_matrix_mouse_all_dedup_N2.rds')
nme <- readRDS('../downstream/input/NME_matrix_mouse_all_dedup_N2.rds')
region_annotation<- readRDS("../downstream/output/mm10_allele_agnostic_analysis.rds")#Regions without CpG are filtered out
uc <- readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
bar_region_number<-function(DNase_num,control_num){
  sum_table=data.table(region_type=c('DNase','control'),
                       count_region=c(DNase_num,control_num))
 print(ggplot(sum_table, aes(x=region_type,y=count_region,fill=region_type))+
   geom_bar(stat='identity')+geom_text(label=sum_table$count_region,position =  position_dodge(0.9)))
}
bar_region_number(sum(region_annotation$region_type=="DNase"),sum(region_annotation$region_type=="control"))
region_annotation<-region_annotation[region_annotation$N>=2]
bar_region_number(sum(region_annotation$region_type=="DNase"),sum(region_annotation$region_type=="control"))
region_annotation_gr <- convert_GR(region_annotation,direction = 'DT')
dnase_gr <- region_annotation_gr[region_type=='DNase']$region
control_gr <- region_annotation_gr[region_type=='control']$region
#getting dnase regions
uc<- sapply(uc,function(gr) {
  grv <- paste0(as.character(seqnames(gr)),':',start(gr),'-',end(gr))
  am <- as.matrix(mcols(gr))
  rownames(am) <- grv
  am <- am[,!grepl('P0',colnames(am))]
  am <- am[complete.cases(am),]
})
uc <- uc[sapply(uc,ncol) > 6]
#jsd only
#names(uc) <- sapply(uc,function(i) sub('_.*','',sub('_all-all','',colnames(i)[1])))
names(uc) <- sapply(uc,function(i) sub('-.*','',colnames(i)[1]))
for (i in 1:length(uc)) colnames(uc[[i]]) <- gsub('-all','',gsub(paste0(names(uc)[[i]],'-'),'',colnames(uc[[i]])))
#Do it for NME and MML
grv <- paste0(as.character(seqnames(mml)),':',start(mml),'-',end(mml))
grv2 <- paste0(as.character(seqnames(nme)),':',start(nme),'-',end(nme))
cat('NME,MML identical check:',identical(grv,grv2),'\n')
mml <- as.matrix(mcols(mml))
nme <- as.matrix(mcols(nme))
rownames(mml) <- rownames(nme) <- grv
mml_DNase=mml[rownames(mml) %in% dnase_gr,]
nme_DNase=nme[rownames(nme) %in% dnase_gr,]
mml_control=mml[rownames(mml) %in% control_gr,]
nme_control=nme[rownames(nme) %in% control_gr,]
uc_DNase=sapply(uc,function(gr){
  gr[rownames(gr) %in% dnase_gr,]
  
})
uc_control=sapply(uc,function(gr){
  gr[rownames(gr) %in% control_gr,]
  
})
for(ts in names(uc_DNase)){
pdf(paste0('../downstream/output/uc_issue/',ts,'_DNase_control.pdf'))
  bar_region_number(nrow(uc_DNase[[ts]]),nrow(uc_control[[ts]]))
  dev.off()
}
for(ts in names(uc_DNase)){
png(paste0('../downstream/output/uc_issue/',ts,'_DNase_control_qq.png'))
qqplot(as.numeric(uc_control[[ts]]),as.numeric(uc_DNase[[ts]]))
abline(a=0,b=1)
dev.off()
}
saveRDS(mml_DNase,file='../downstream/output/mml_matrix_DNase.rds')
saveRDS(nme_DNase,file='../downstream/output/nme_matrix_DNase.rds')
saveRDS(uc_DNase,file='../downstream/output/uc_matrix_DNase.rds')

saveRDS(mml_control,file='../downstream/output/mml_matrix_control.rds')
saveRDS(nme_control,file='../downstream/output/nme_matrix_control.rds')
saveRDS(uc_control,file='../downstream/output/uc_matrix_control.rds')

t_null=ecdf(as.numeric(uc_control$heart))
t_obs=as.data.table(uc_DNase$heart)
t_obs$region=rownames(uc_DNase$heart)
t_obs=melt.data.table(t_obs,id.vars='region')
t_obs$pvalue=1-t_null(t_obs$value)
t_obs$FDR=p.adjust(t_obs$pvalue,method='BH')

#https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-8-8