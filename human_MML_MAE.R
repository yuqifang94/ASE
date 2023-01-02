source('mainFunctions_sub.R')
#MAE enrich
MAE_enrich<-function(GR_merge,pval_cutoff,genes='genes_promoter',stat='dMML_pval',MAE=MAE){
  #GR_merge=GR_merge[!is.na(GR_merge$genes_promoter)]
  GR_merge=elementMetadata(GR_merge)
  stat_gene=GR_merge[[genes]][GR_merge[[stat]]<=pval_cutoff]
  non_stat_gene=GR_merge[[genes]][GR_merge[[stat]]>pval_cutoff]
  stat_MAE=sum(unlist(lapply(stat_gene,function(x) any(x %in% MAE))))
  nonstat_MAE=sum(unlist(lapply(non_stat_gene,function(x) any(x %in% MAE)))) 
  nonstat_nonMAE=sum(!unlist(lapply(non_stat_gene,function(x) any(x %in% MAE))))
  stat_nonMAE=sum(!unlist(lapply(stat_gene,function(x) any(x %in% MAE))))
  print(matrix(c(stat_MAE,stat_nonMAE,nonstat_MAE,nonstat_nonMAE),nrow=2))
  ft=fisher.test(matrix(c(stat_MAE,stat_nonMAE,nonstat_MAE,nonstat_nonMAE),nrow=2))
  return(data.frame(OR=ft$estimate,pvalue=ft$p.value,lowerCI= ft$conf.int[1],upperCI=ft$conf.int[2]))
}
# Finding the overlap between monoallelic expressed gene ------------------
GR_merge=readRDS(GR_merge_file)
MAE_BAE_data_Gimelbrant <- as.data.frame(read_excel("../downstream/input/human_analysis/imprinting_ASE/MAE_BAE_data_Gimelbrant.xlsx"),stringsAsFactors=F)
#Enrichment in MAE dMML , not quiet enriched
MAE=MAE_BAE_data_Gimelbrant$Gene[ MAE_BAE_data_Gimelbrant$`MAE=1_BAE=0`==1]
#Find GR overlap the MAE:NME gene promoter
MAE_GR_dNME=GR_merge[GR_merge$dNME_pval<=pval_cutoff][GR_merge$genes_promoter[GR_merge$dNME_pval<=pval_cutoff] %in% MAE]
#write(unique(MAE_GR_dNME$genes_promoter),"../downstream/output/dNME_sig_gene_body_MAE.txt",pval_cutoff=0.1)
#dMML enrichment promoter
MAE_enrich(GR_merge[GR_merge$genes_promoter%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_promoter',stat='dMML_pval',MAE=MAE)

#dNME enrichment promoter
MAE_enrich(GR_merge[GR_merge$genes_promoter%in% MAE_BAE_data_Gimelbrant$Gene],
           pval_cutoff=0.1,genes='genes_promoter',stat='dNME_pval',MAE=MAE)


# #dNME enrichment body
# MAE_enrich(GR_merge[GR_merge$genes_body%in% MAE_BAE_data_Gimelbrant$Gene],
#            pval_cutoff=0.1,genes='genes_body',stat='dNME_pval',MAE=MAE)
# #dMML enrichment body
# MAE_enrich(GR_merge[GR_merge$genes_body%in% MAE_BAE_data_Gimelbrant$Gene],
#            pval_cutoff=0.1,genes='genes_body',stat='dMML_pval',MAE=MAE)

#Imprinted genes
Imprinted_Genes <- as.data.frame(read_excel("../downstream/input/human_analysis/imprinting_ASE/Imprinted Genes.xlsx"))
imrinting_out=rbind(cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',Imprinted_Genes$Gene),
      data.table(imprinting="All imprinted\ngenes")),
      
      cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',
           Imprinted_Genes$Gene[Imprinted_Genes[,5]=='Paternal']),
           data.table(imprinting="Paternally\nexpressed")),

      cbind(MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dMML_pval',
           Imprinted_Genes$Gene[Imprinted_Genes[,5]=='Maternal']), data.table(imprinting="Maternally\nexpressed"))
)
imrinting_out$imprinting=factor(imrinting_out$imprinting,levels=imrinting_out$imprinting)
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.text.x=element_text(size=9),
                                 axis.text.y=element_text(size=9))
pdf('../downstream/output/human_analysis/imprinting/dMML_imprinting_OR.pdf',width=4.75,height=4.75)
ggplot(imrinting_out,aes(x=imprinting,y=OR))+  
  geom_bar(stat="identity", color="light blue",fill="light blue", position=position_dodge()) +
  geom_errorbar(aes(ymin=lowerCI  , ymax=upperCI), width=.2,position=position_dodge(.9))+
  geom_text(aes(label=round(OR,digit=1),y=upperCI),vjust=-0.5)+xlab("")+ylab("Enrichment")+
  theme_glob

dev.off()

#This is for dNME
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',Imprinted_Genes$Gene)
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',
           Imprinted_Genes$Gene[Imprinted_Genes[,5]=='Paternal'])
MAE_enrich(GR_merge[!is.na(GR_merge$genes_promoter)],0.1,'genes_promoter','dNME_pval',
           Imprinted_Genes$Gene[Imprinted_Genes[,5]=='Maternal'])
