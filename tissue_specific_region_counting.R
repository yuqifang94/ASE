 source('mainFunctions_sub.R')
 UC_in=readRDS(UC_in_matrix_cluster_file)
 cut=0.1
 #Non tissue-specific high UC
  aid <- sapply(names(UC_in),function(i) {
    names(which(rowSums(UC_in[[i]] > cut) > 0))
  })  
 ts_aid <- sapply(names(UC_in),function(i) {
    #This is for one and only one
    return(setdiff(aid[[i]],unlist(aid[names(aid)!=i])))
    
})
length(unique(unlist(aid)))
length(unique(unlist(ts_aid)))
table(table(unlist(aid)))
#    1      2      3      4      5      6      7
# 588189 281624 192486 202206 278883 324387 218109
#How do they overlap enhancers
