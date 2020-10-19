rm(list=ls())
source("mainFunctions_sub.R")
library(Gmisc)
#Define ggplot theme

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
#Density analysis
GR_merge=readRDS(GR_merge_file)
#####Subsetting by DNase region will not give us enough power to do it

###Reading in data
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
motif_gene=fastDoCall('c',lapply(dir('../downstream/input/JASPAR_out/'),function(x){
  cat(x,'\n')
  readRDS(paste0('../downstream/input/JASPAR_out/',x))
  
}))
saveRDS(motif_gene,motif_gene_file)
motif_gene <- readRDS(motif_gene_file)
motif_dir=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1)
motif_dir_N1=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta[variant_HetCpG_meta$N>=1],
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1)
motif_dir$qval_binom=p.adjust(motif_dir$binom.pval,method="BH")
human_mono_motif_TSV=as.data.frame(read.table('../downstream/input/JASPAR_human_redundant_2020.csv',sep=',',header=T,stringsAsFactors = F))
more_ent_enrich=motif_family_enrich(unique(motif_dir$TF[motif_dir$qval_binom<=0.1 & motif_dir$prob>0.5]),
                                    unique(motif_gene$geneSymbol),human_mono_motif_TSV)
more_less_ent_enrich=motif_family_enrich(unique(motif_dir$TF[motif_dir$qval_binom<=0.1 & motif_dir$prob<0.5]),
                                    unique(motif_gene$geneSymbol),human_mono_motif_TSV)
write.csv(unique(motif_dir$TF[motif_dir$qval_binom<=0.1 & motif_dir$prob>0.5],
                 '../downstream/output/graphs/tableS1_motif_prefer_ent.csv'))
write.csv(more_ent_enrich,'../downstream/output/graphs/table1_motif_family_prefer_ent.csv')
write.csv(more_less_ent_enrich,'../downstream/output/graphs/table2_motif_family_not_prefer_ent.csv')
OMIM=read.csv('../downstream/input/genemap2.txt',fill = T,sep='\t',header=T,stringsAsFactors =F,na.strings = "NA",skip=3)
OMIM$Gene.Symbols=strsplit(as.character(OMIM$Gene.Symbols),', ')

OMIM_annotation<-function(motif_in,OMIM){
motif_in=motif_in[order(motif_in$prob,decreasing=T),]
motif_in$OMIM=NA
for(tf in unique(motif_in$TF)){
    tf_in=gsub('\\(var.2\\)','',tf)
    tf_in=gsub('\\(var.3\\)','',tf_in)
    
    tf_in=unlist(strsplit(tf_in,'::'))
    #print(tf_in)
    OMIM_disease=OMIM$Phenotypes[which(unlist(lapply(OMIM$Gene.Symbols,function(x) any(x%in% tf_in))))]
    if(length(OMIM_disease)>0){motif_in$OMIM[motif_in$TF==tf]=as.list(OMIM_disease)}
  }
  motif_in$OMIM=unlist(lapply(motif_in$OMIM,function(x) paste(x,collapase=',')))
  return(motif_in)
}
motif_ent_OMIM=OMIM_annotation(motif_dir[motif_dir$qval_binom<=0.1 & motif_dir$prob>0.5,],OMIM)
write.csv(motif_ent_OMIM[c('TF','total_data','same_dir','opposite_dir','prob','qval_binom','OMIM')],
          '../downstream/output/graphs/tableS1_motif_prefer_ent_OMIM.csv')
motif_non_ent_OMIM=OMIM_annotation(motif_dir[motif_dir$qval_binom<=0.1 & motif_dir$prob<0.5,],OMIM)
write.csv(motif_non_ent_OMIM[c('TF','total_data','same_dir','opposite_dir','prob','qval_binom','OMIM')],
          '../downstream/output/graphs/table3_motif_not_prefer_ent_OMIM.csv')


