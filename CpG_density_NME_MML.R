rm(list=ls())
source("mainFunctions_sub.R")
#Define ggplot theme

theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=20,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=20,face="bold"),
                 axis.text.x=element_text(size=12),
                 axis.text.y=element_text(size=12))
                # text=element_text(family="Space Mono"))

# Preprocess SNP files to get unique SNP and trinucleotide -----------------------------------------------------
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
olap=findOverlaps(variant_HetCpG_meta,genomic_features$`CpG island`)
variant_HetCpG_meta_dt=convert_GR(variant_HetCpG_meta,direction='DT')
variant_HetCpG_meta_dt$SNP=apply(variant_HetCpG_meta_dt[,list(REF,ALT)],1,function(x) paste(x,collapse = '->'))
variant_HetCpG_meta_dt$tri_SNP=paste0(variant_HetCpG_meta_dt$REF_tri,'->',variant_HetCpG_meta_dt$ALT_tri)
saveRDS(variant_HetCpG_meta_dt,variant_HetCpG_meta_dt_file)
single_SNP_unique=unique_mutation(unique(variant_HetCpG_meta_dt$SNP))
#manually order unique mutation
single_SNP_unique[single_SNP_unique=="T->C"]="C->T"
single_SNP_unique[single_SNP_unique=="T->G"]="G->T"
single_SNP_unique[single_SNP_unique=="A->G"]="G->A"
single_SNP_unique[single_SNP_unique=="G->C"]="C->G"
single_SNP_unique[single_SNP_unique=="A->T"]="T->A"
single_SNP_unique[single_SNP_unique=="A->C"]="C->A"
mutation_tri_unique=unique_mutation(unique(variant_HetCpG_meta_dt$tri_SNP))
# #correct order based on single SNP mutation
for(i in 1:length(mutation_tri_unique)){
  tri=mutation_tri_unique[i]

  if(!tri_to_SNP(tri)%in%single_SNP_unique){
    tri_1=gsub('->.*','',tri)
    tri_2=gsub('.*->','',tri)
    mutation_tri_unique[i]=paste0(tri_2,'->',tri_1)


  }

}
#Make a hashtable-like to easy access data
mutation_tri_unique_dt=data.table(raw=names(mutation_tri_unique),unique_fw=mutation_tri_unique)
#Merge reverse complement
single_SNP_unique[single_SNP_unique=="G->A"]="C->T"
single_SNP_unique[single_SNP_unique=="G->T"]="C->A"
variant_HetCpG_meta_dt$SNP=single_SNP_unique[variant_HetCpG_meta_dt$SNP]
mutation_tri_unique_dt$rev_comp=unlist(lapply(mutation_tri_unique_dt$unique_fw,reverse_comp_SNP))
mutation_tri_unique_dt$rev_comp_inv=paste0(gsub('.*->','',mutation_tri_unique_dt$rev_comp),'->',
                                           gsub('->.*','',mutation_tri_unique_dt$rev_comp))
mutation_tri_unique_dt$rev_comp_unique="NA"
for(i in 1:nrow(mutation_tri_unique_dt)){
  #Only look for trinucleotide in the single SNP catogry
  tri_in=mutation_tri_unique_dt$unique_fw[i]
  if(tri_to_SNP(tri_in)%in%single_SNP_unique){
    mutation_tri_unique_dt[(rev_comp==tri_in|unique_fw==tri_in|rev_comp_inv==tri_in)&rev_comp_unique=="NA"]$rev_comp_unique=tri_in
    
    
  }
  
  
}
#Should be 52 unique ones
mutation_tri_unique=mutation_tri_unique_dt$rev_comp_unique
names(mutation_tri_unique)=mutation_tri_unique_dt$raw
#Check how many of them have CG in right, should be 0 since the trinucleotide is ordered
gainCG_idx=which(grepl("CG",sub('.*->','',mutation_tri_unique))&!grepl("CG",sub('->.*','',mutation_tri_unique)))#Should only be 0

# calculate the relative NME using the order of SNP -----------------------

variant_HetCpG_meta_dt$tri_SNP_unique=mutation_tri_unique[variant_HetCpG_meta_dt$tri_SNP]
variant_HetCpG_meta_dt$dNME_relative=as.numeric(NA)
#Use minus sign to get NME right-NME left, .e.g (NME in GTG- NME in GCG) for (GCG->GTG)
variant_HetCpG_meta_dt$dNME_relative=-variant_HetCpG_meta_dt[,list(dNME_relative=dNME_relative_calc(genome1_tri,genome2_tri,NME1,NME2,tri_SNP,tri_SNP_unique,SNP)), 
                                                            by = seq_len(nrow(variant_HetCpG_meta_dt))]$dNME_relative
  

#Convert everything to gain CG, id SNPs with CG changes
variant_HetCpG_meta_dt$CpG_change='Gain CG'
variant_HetCpG_meta_dt[((grepl('CG',REF_tri )) & (grepl('CG',ALT_tri)))|(!grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri))]$CpG_change='No CG change'
saveRDS(variant_HetCpG_meta_dt,variant_HetCpG_meta_dt_uq_file)
variant_HetCpG_meta_dt=readRDS(variant_HetCpG_meta_dt_uq_file)

#Generating Figure 4B and calcualting OR for dNME
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
for (sn in unique(variant_HetCpG_meta_dt$SNP)){
  #OR calculation
  for(tri in unique(variant_HetCpG_meta_dt[SNP==sn]$tri_SNP_unique)){
    variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
    variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
    
      variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
      
      
   
    variant_SNP_tri=rbind(variant_SNP_tri,variant_SNP_tri_OR)
  }
  #get HetCpG
  variant_SNP_tri$CpG_change=factor(variant_SNP_tri$CpG_change,levels=c('Gain CG','No CG change'))
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

     scale_fill_manual(values=c("No CG change"="grey","Gain CG"="light blue"))+
     geom_text(data=variant_SNP_tri[OR>1],aes(label=significant,y=log(upperCI)*1),vjust =sig_v,hjust=sig_h_pos)+
     geom_text(data=variant_SNP_tri[OR<1],aes(label=significant,y=log(lowerCI)*1),vjust =sig_v,hjust=sig_h_neg)+
     coord_flip()
   

   variant_SNP_tri_out[[sn]]=variant_SNP_tri
   variant_SNP_tri=data.table()

}
#Getting png files with mono-spaced font
library(extrafont)
loadfonts()
png('../downstream/output/graphs/Figure2/Figure-4B-variant_OR_tri3_two_cat_greater_CG_bg_rev.png',
    width=7,height=7,units='in',res=1080, family = 'Consolas')
#SNP_het=SNP_het[c("C>G", names(SNP_het)[names(SNP_het)!="C>G"])]
ggarrange(plotlist=SNP_het, nrow=2,ncol=2,common.legend = T,legend="bottom",label.x = 'Odds ratio')
dev.off()

#Calculate OR for SNPs gaining CG, numbers in text
OR_calc(variant_HetCpG_meta_dt[dNME_pval<=pval_cutoff],'Gain CG',"CpG_change")

# Density analysis using allele-specific way using regions------------------------------
GR_merge=readRDS(GR_merge_file)

GR_merge$CpGdiff=GR_merge$g1CG-GR_merge$g2CG
#CpG density vs dNME,here uses an extended density 
GR_merge$dNME_relative=GR_merge$NME1-GR_merge$NME2
GR_merge$dMML_relative=GR_merge$MML1-GR_merge$MML2
GR_merge_dt=convert_GR(GR_merge,'DT')
#density difference difference
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
#Correlation 
cor.test(GR_merge_dt[dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$dNME_relative, 
         GR_merge_dt[dNME_pval<=pval_cutoff&dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$density_diff)

#Figure 3A
#Catogrizing regions
GR_merge_dt$CpG_stat="No difference"
GR_merge_dt[CpGdiff!=0]$CpG_stat="With CpG difference"
GR_merge_dt$CpG_stat=factor(GR_merge_dt$CpG_stat,levels = c("With CpG difference","No difference"))
#Alwasy using allele with more CG minus alleles with less CG
GR_merge_dt$dNME_relative_more_less=GR_merge_dt$dNME_relative
GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative_more_less=GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative*sign(GR_merge_dt[GR_merge_dt$CpGdiff!=0]$CpGdiff)
t.test(GR_merge_dt[CpGdiff!=0&dNME_pval<=pval_cutoff]$dNME_relative_more_less,alternative="less")
#Figure 3C
pdf('../downstream/output/graphs/Figure3/Figure3A_CpG_number_NME.pdf',width=7,height=7)
ggplot(GR_merge_dt[dNME_pval<=pval_cutoff],aes(y=dNME_relative_more_less,x=CpG_stat,fill=CpG_stat))+
  geom_violin()+xlab("")+
  theme_glob+ylab('relative dNME')+theme(legend.position = "none")
dev.off()
#Line plot for dNME
GR_merge_dt_sig_density_diff=GR_merge_dt[dNME_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_sig_density_diff$density_difference_quantile=ecdf(GR_merge_dt_sig_density_diff$density_diff)(GR_merge_dt_sig_density_diff$density_diff)
pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dNME_ratio.pdf',width=7,height=7)
ggplot(GR_merge_dt_sig_density_diff,aes(x=density_difference_quantile,y=dNME_relative))+geom_smooth(fill='light blue')+
  xlab("CpG density ratio quantile")+ylab("relative dNME")+
  theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
#Line plot for dMML
GR_merge_dt_sig_density_diff=GR_merge_dt[dMML_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_sig_density_diff$density_difference_quantile=ecdf(GR_merge_dt_sig_density_diff$density_diff)(GR_merge_dt_sig_density_diff$density_diff)
pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dMML_ratio.pdf',width=7,height=7)
ggplot(GR_merge_dt_sig_density_diff,aes(x=density_difference_quantile,y=dMML_relative))+geom_smooth(fill='light blue')+
  xlab("CpG density ratio quantile")+ylab("relative dMML")+
  theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#allele-agnotic density ---------------------------------------
NME_in=readRDS(NME_agnostic_file)
CpG_hg19=readRDS(CpG_hg19_file)
NME_in$CG_hg19=countOverlaps(NME_in,CpG_hg19)
NME_in_gr=unique(granges(NME_in))
gr_seq=getSeq(Hsapiens,NME_in_gr,as.character=T)
NME_in_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
NME_in_olap=findOverlaps(NME_in,NME_in_gr,type='equal')
NME_in$CGcont_exp[queryHits(NME_in_olap)]=NME_in_gr$CGcont_exp[subjectHits(NME_in_olap)]
NME_in$density=NME_in$CG_hg19/NME_in$CGcont_exp
NME_in=readRDS(NME_agnostic_file)
NME_in=NME_in[seqnames(NME_in)%in%paste0("chr",1:22)]
cor.test(NME_in$density,NME_in$NME,method='pearson')
#Making boxplot of this with different interval
NME_in$density_quant=findInterval(NME_in$density,seq(0,1,0.1))
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
NME_in$density_quant=factor(quant_conv[NME_in$density_quant],levels=quant_conv)
#PLotting boxplot
pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot_CG_exp.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in)),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()   

#Feature enrichment of low NME regions
genomic_features=readRDS(genomic_features_file)
#Figure S5
olap_islands=findOverlaps(NME_in,genomic_features$`CpG island`)
olap_shores=findOverlaps(NME_in,genomic_features$`CpG shore`)
olap_shelf=findOverlaps(NME_in,genomic_features$`CpG shelf`)
olap_open_sea=findOverlaps(NME_in,genomic_features$`CpG open sea`)

olap_islands=findOverlaps(MML_in,genomic_features$`CpG island`)

CpG_density_NME=rbind(data.table(NME=NME_in$NME[queryHits(olap_islands)],feature='islands'),
                      data.table(NME=NME_in$NME[queryHits(olap_shores)],feature='shores'),
                      data.table(NME=NME_in$NME[queryHits(olap_shelf)],feature='shelf'),
                      data.table(NME=NME_in$NME[queryHits(olap_open_sea)],feature='open sea'))

# In mouse context --------------------------------------------------------
mml=readRDS(MML_matrix_file)
nme=readRDS(NME_matrix_file)
mm10_CpG=getCpgSitesmm10()
nme$CG_mm10=countOverlaps(nme,mm10_CpG)
mml$CG_mm10=countOverlaps(mml,mm10_CpG)
density_gr=unique(c(granges(nme),granges(mml)))
gr_seq=getSeq(Mmusculus,density_gr,as.character=T)
density_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
#the regions are the same for mml and nme so use same gr files
saveRDS(density_gr,CG_density_mouse)
nme_olap=findOverlaps(nme,density_gr,type='equal')
nme$CGcont_exp[queryHits(nme_olap)]=density_gr$CGcont_exp[subjectHits(nme_olap)]
nme$density=nme$CG_mm10/nme$CGcont_exp
mml_olap=findOverlaps(mml,density_gr,type='equal')
mml$CGcont_exp[queryHits(mml_olap)]=density_gr$CGcont_exp[subjectHits(mml_olap)]
mml$density=mml$CG_mm10/mml$CGcont_exp
nme=nme[seqnames(nme) %in% paste0("chr",1:19)]
mml=mml[seqnames(mml) %in% paste0("chr",1:19)]
nme$density_quant=findInterval(nme$density,seq(0,1,0.1))
mml$density_quant=findInterval(mml$density,seq(0,1,0.1))
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
nme$density_quant=factor(quant_conv[nme$density_quant],levels=quant_conv)

mml$density_quant=factor(quant_conv[mml$density_quant],levels=quant_conv)
density_mouse_calc<-function(gr_in,stat_name="NME"){
  gr_in=mcols(gr_in)
  gr_in=melt.data.table(as.data.table(gr_in),id.vars = c("CG_mm10","CGcont_exp","density"),variable.name = "Sample",value.name="stat_in")
  gr_in$density_quant=findInterval(gr_in$density,seq(0,1,0.1))
  #NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
  quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
  gr_in$density_quant=factor(quant_conv[gr_in$density_quant],levels=quant_conv)
  print(ggplot(as.data.frame(gr_in),aes(x=density_quant, y=stat_in))+
          ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
          ylab(stat_name)+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
}
pdf('../downstream/output/graphs/Figure3/mouse_NME_density_boxplot.pdf',width=3.5,height=3.5)
density_mouse_calc(nme,stat_name="NME")
dev.off()