UC_NT_selected=UC_raw$NT[dMML_perm_in[region%in% csv_in$region&value>=0.95]$region,]

UC_NT_selected=UC_raw$NT[dmml_perm[value>=0.9]$region,]
mml=readRDS('../downstream/output/mml_matrix_DNase.rds')
mml_NT=mml[rownames(UC_NT_selected),grep("NT",colnames(mml))]
colnames(mml_NT)=gsub('-all','',gsub('NT-','',colnames(mml_NT)))

dmml_NT=do.call(cbind,lapply(colnames(UC_NT_selected),function(x){
  diff=matrix(abs(mml_NT[,gsub('.*-','',x)]-mml_NT[,gsub('-.*','',x)]),ncol=1)
  colnames(diff)=x
  return(diff)
  }))

nme=readRDS('../downstream/output/nme_matrix_DNase.rds')
nme_NT=nme[rownames(UC_NT_selected),grep("NT",colnames(nme))]
colnames(nme_NT)=gsub('-all','',gsub('NT-','',colnames(nme_NT)))

dnme_NT=do.call(cbind,lapply(colnames(UC_NT_selected),function(x){
  diff=matrix(abs(nme_NT[,gsub('.*-','',x)]-nme_NT[,gsub('-.*','',x)]),ncol=1)
  colnames(diff)=x
  return(diff)
}))
adj_timepoints=c("E11.5-E12.5","E12.5-E13.5","E13.5-E14.5","E14.5-E15.5")
#Due to low number of time points, we may drop it because it's heard to generate permute
pheatmap(scalematrix(UC_NT_selected),show_rownames = F)
