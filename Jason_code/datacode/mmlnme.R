suppressMessages(library(GenomicRanges))
mml <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/raw/MML_matrix_mouse_all_dedup_N2.rds')
nme <- readRDS('/dcl01/hongkai/data/zji4/ase/mouse/data/raw/NME_matrix_mouse_all_dedup_N2.rds')
grv <- paste0(as.character(seqnames(mml)),':',start(mml),'-',end(mml))
grv2 <- paste0(as.character(seqnames(nme)),':',start(nme),'-',end(nme))
mml <- as.matrix(mcols(mml))
nme <- as.matrix(mcols(nme))
rownames(mml) <- rownames(nme) <- grv
saveRDS(mml,file='/dcl01/hongkai/data/zji4/ase/mouse/data/proc/mml.rds')
saveRDS(nme,file='/dcl01/hongkai/data/zji4/ase/mouse/data/proc/nme.rds')

