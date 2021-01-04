source('mainFunctions_sub.R')
corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}
#Get the tissue x region
#get nme and mml
mml <- readRDS('../downstream/input/MML_matrix_mouse_all_dedup_N2.rds')
nme <- readRDS('../downstream/input/NME_matrix_mouse_all_dedup_N2.rds')
dnase<- readRDS('../downstream/input/mm10_DNase.rds')
uc <- readRDS('../downstream/input/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
dnase<-dnase[dnase$N>=2]
dnase_gr <- paste0(as.character(seqnames(dnase)),':',start(dnase),'-',end(dnase))
dnase_gr <- dnase_gr[mcols(dnase)$region_type=='DNase']
#getting dnase regions
uc <- sapply(uc,function(gr) {
  grv <- paste0(as.character(seqnames(gr)),':',start(gr),'-',end(gr))
  am <- as.matrix(mcols(gr))
  rownames(am) <- grv
  am <- am[rownames(am) %in% dnase_gr,]
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
mml=mml[rownames(mml) %in% dnase_gr,]
nme=nme[rownames(nme) %in% dnase_gr,]
saveRDS(mml,file='../downstream/output/mml_matrix_DNase.rds')
saveRDS(nme,file='../downstream/output/nme_matrix_DNase.rds')
saveRDS(uc,file='../downstream/output/uc_matrix_DNase.rds')

