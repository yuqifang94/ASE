source('mainFunctions_sub.R')
variant_HetCpG_meta_dt=readRDS(variant_HetCpG_meta_dt_uq_file)
variant_HetCpG_meta_dt$CpG_change=gsub('CG','CpG',variant_HetCpG_meta_dt$CpG_change)
#Generating Figure 4B and calculating OR for dNME
#Param initialization & color theme
SNP_all=list()
SNP_het=list()
SNP_box=list()
sig_v=0.75
sig_h_pos=-0.05
sig_h_neg=1.2
color_theme=c(rainbow(length(unique(variant_HetCpG_meta_dt$SNP))))
variant_SNP_tri=data.table()
variant_SNP_tri_out=list()
names(color_theme)=unique(variant_HetCpG_meta_dt$SNP)
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.text.x=element_text(size=7),
                                 axis.text.y=element_text(size=7))
#MML
for (sn in unique(variant_HetCpG_meta_dt$SNP)){
  #OR calculation
  for(tri in unique(variant_HetCpG_meta_dt[SNP==sn]$tri_SNP_unique)){
    variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dMML_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
    variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
    
      variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
      
      
   
    variant_SNP_tri=rbind(variant_SNP_tri,variant_SNP_tri_OR)
  }
  #get HetCpG
  variant_SNP_tri$CpG_change=factor(variant_SNP_tri$CpG_change,levels=c('Lose CpG','No CpG change'))
  variant_SNP_tri=variant_SNP_tri[order(OR,decreasing=F)]
 
  variant_SNP_tri$SNP=gsub('->','\u2794',variant_SNP_tri$SNP)
  variant_SNP_tri$SNP=factor(variant_SNP_tri$SNP,levels = variant_SNP_tri$SNP)
  
  variant_SNP_tri$FDR=p.adjust(variant_SNP_tri$pvalue,method='BH')
  variant_SNP_tri$significant=add.significance.stars(variant_SNP_tri$FDR, cutoffs = c(0.05, 0.01, 0.001))
  #Plotting
   SNP_het[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=log(OR),fill=CpG_change))+geom_bar(stat="identity")+ylab('')+xlab("")+
    geom_errorbar(aes(ymin=log(lowerCI), ymax=log(upperCI)), width=.4,position=position_dodge(.9),size=0.25)+ggtitle(gsub('->',' \u2794 ',sn))+#ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+
    theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
     ylim(c(-1.5,1.5))+
     scale_fill_manual(values=c("No CpG change"="grey","Lose CpG"="light blue"))+
     #geom_text(data=variant_SNP_tri[OR>1],aes(label=significant,y=log(upperCI)*1),vjust =sig_v,hjust=sig_h_pos)+
     #geom_text(data=variant_SNP_tri[OR<1],aes(label=significant,y=log(lowerCI)*1),vjust =sig_v,hjust=sig_h_neg)+
     coord_flip()
   

   variant_SNP_tri_out[[sn]]=variant_SNP_tri
   variant_SNP_tri=data.table()

}
saveRDS(SNP_het,'../downstream/output/human_analysis/CpG_density/SNP_het_CpG_dMML.rds')
saveRDS(variant_SNP_tri_out,'../downstream/output/human_analysis/CpG_density/variant_SNP_tri_out_CpG_dMML.rds')

#Getting png files with mono-spaced font in windows setting, use the variant_SNP_tri_out
variant_SNP_tri_out=readRDS('../downstream/output/human_analysis/CpG_density/variant_SNP_tri_out_CpG_dMML.rds')
library(extrafont)
loadfonts()
SNP_het=readRDS('../downstream/output/human_analysis/CpG_density/SNP_het_CpG_dMML.rds')
png('../downstream/output/human_analysis/CpG_density/variant_OR_tri3_two_cat_greater_CG_bg_rev_hg19_dMML.png',
    width=7,height=7,units='in',res=1080, family = 'Consolas')
#SNP_het=SNP_het[c("C>G", names(SNP_het)[names(SNP_het)!="C>G"])]
ggarrange(plotlist=lapply(SNP_het,function(x) x+ ylab("log(Odds Ratio)")+theme( axis.title.x=element_text(hjust=0.5,size=16,face="bold"))), 
          nrow=2,ncol=2,common.legend = T,legend="top")
dev.off()

#Calculate OR for SNPs gaining CG, numbers in text
OR_calc(variant_HetCpG_meta_dt[dMML_pval<=pval_cutoff],'Lose CpG',"CpG_change")
