rm(list=ls())
source("mainFunctions_sub.R")
#Define ggplot theme

theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=20,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=20,face="bold"),
                 axis.text.x=element_text(size=12),
                 axis.text.y=element_text(size=12))
                # text=element_text(family="Space Mono"))

# Check SNP frequency -----------------------------------------------------
genomic_features=readRDS(genomic_features_file)
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
olap=findOverlaps(variant_HetCpG_meta,genomic_features$`CpG island`)
variant_HetCpG_meta$CpG_island=FALSE
variant_HetCpG_meta$CpG_island[queryHits(olap)]=TRUE
variant_HetCpG_meta_dt=convert_GR(variant_HetCpG_meta,direction='DT')

#variant_HetCpG_meta$SNP=apply(variant_HetCpG_meta[,list(REF,ALT)],1,function(x) paste(x,collapse = '>'))
#variant_HetCpG_meta$tri_SNP=paste0(variant_HetCpG_meta$REF_tri,'>',variant_HetCpG_meta$ALT_tri)

variant_HetCpG_meta_dt$SNP=apply(variant_HetCpG_meta_dt[,list(REF,ALT)],1,function(x) paste(x,collapse = '->'))
variant_HetCpG_meta_dt$tri_SNP=paste0(variant_HetCpG_meta_dt$REF_tri,'->',variant_HetCpG_meta_dt$ALT_tri)
saveRDS(variant_HetCpG_meta_dt,'../downstream/output/human_analysis/CPEL_outputs/variant_HetCpG_meta_dt.rds')
#Finding unique trinucleotide and mutation
variant_HetCpG_meta_dt=readRDS('../downstream/output/human_analysis/CPEL_outputs/variant_HetCpG_meta_dt.rds')
#find unique ones
unique_mutation<-function(mutation_in){
  mutation_in=data.table(all_SNP=mutation_in)
  mutation_in$unique_class="NA"
  for(i in 1:nrow(mutation_in)){
    SNP_in=mutation_in$all_SNP[i]
    SNP_in_rev=paste0(gsub('.*->','',SNP_in),'->',gsub('->.*','',SNP_in))
    SNP_in2=c(SNP_in,SNP_in_rev)
    if(all(mutation_in[all_SNP%in%SNP_in2]$unique_class=="NA")){
         mutation_in[all_SNP%in%SNP_in2]$unique_class=SNP_in
    }
  
  }
  mutation_out=mutation_in$unique_class
  names(mutation_out)=mutation_in$all_SNP
  return(mutation_out)
}
tri_to_SNP<-function(tri){
  tri_l=gsub('->.*','',tri)
  tri_r=gsub('.*->','',tri)
  return(paste0(substring(tri_l, 2, 2),'->',substring(tri_r, 2, 2)))
}
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
mutation_tri_unique_dt=data.table(raw=names(mutation_tri_unique),unique_fw=mutation_tri_unique)
reverse_comp_SNP<-function(x){
  return(paste0(reverseComplement(DNAString(gsub('->.*','',x))),'->',
                reverseComplement(DNAString(gsub('.*->','',x)))))
  
  
}

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

gainCG_idx=which(grepl("CG",sub('.*->','',mutation_tri_unique))&!grepl("CG",sub('->.*','',mutation_tri_unique)))#Should only be 0

#mutation_tri_unique[gainCG_idx]=paste0(sub('.*->','', mutation_tri_unique[gainCG_idx]),'->',sub('->.*','', mutation_tri_unique[gainCG_idx]))
#Here uses genome1_tri and NME1 

variant_HetCpG_meta_dt$tri_SNP_unique=mutation_tri_unique[variant_HetCpG_meta_dt$tri_SNP]
variant_HetCpG_meta_dt$dNME_relative=as.numeric(NA)
dNME_relative_calc<-function(genome1_tri,genome2_tri,NME1,NME2,tri_SNP,tri_SNP_unique,SNP) {
  #Check if trinucleotide is in tri_SNP_unique

  genome1_tri_rev=as.character(reverseComplement(DNAString(genome1_tri)))
  genome2_tri_rev=as.character(reverseComplement(DNAString(genome2_tri)))
  
  if(all(c(genome1_tri,genome2_tri) %in% unlist(strsplit(tri_SNP_unique,'->')))){
    return(c(NME1,NME2)[which(c(genome1_tri,genome2_tri)==gsub('->.*','',tri_SNP_unique))]-
             c(NME1,NME2)[which(c(genome1_tri,genome2_tri)==gsub('.*->','',tri_SNP_unique))])
  }  else 
    if(all(c(genome1_tri_rev,genome2_tri_rev)%in%unlist(strsplit(tri_SNP_unique,'->')))){
      return(c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri_rev)==gsub('->.*','',tri_SNP_unique))]-
               c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri_rev)==gsub('.*->','',tri_SNP_unique))])
    }else 
      if(all(c(genome1_tri_rev,genome2_tri)%in%unlist(strsplit(tri_SNP_unique,'->')))){
        return(c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri)==gsub('->.*','',tri_SNP_unique))]-
                 c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri)==gsub('.*->','',tri_SNP_unique))])
      }else 
        if(all(c(genome1_tri,genome2_tri_rev)%in%unlist(strsplit(tri_SNP_unique,'->')))){
          return(c(NME1,NME2)[which(c(genome1_tri,genome2_tri_rev)==gsub('->.*','',tri_SNP_unique))]-
                   c(NME1,NME2)[which(c(genome1_tri,genome2_tri_rev)==gsub('.*->','',tri_SNP_unique))])
        }
  
}
#Use minus sign to get NME2-NME1
variant_HetCpG_meta_dt$dNME_relative=-variant_HetCpG_meta_dt[,list(dNME_relative=dNME_relative_calc(genome1_tri,genome2_tri,NME1,NME2,tri_SNP,tri_SNP_unique,SNP)), 
                                                            by = seq_len(nrow(variant_HetCpG_meta_dt))]$dNME_relative
  
#Use minus sign to get NME2-NME1 
variant_HetCpG_meta_dt$dNME_relative=-variant_HetCpG_meta_dt$dNME_relative
# # #Lumping all
# variant_HetCpG_meta_dt$CpG_change="No change"
# variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri ))]$CpG_change='Lose CG'
# variant_HetCpG_meta_dt[(!grepl('CG',REF_tri )) & (grepl('CG',ALT_tri ))]$CpG_change='Gain CG'
# variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (grepl('CG',ALT_tri))]$CpG_change='No change'

# #Calculate relative dNME
# variant_HetCpG_meta_dt$dNME_relative=variant_HetCpG_meta_dt$refNME-variant_HetCpG_meta_dt$altNME
# #convert all to less CG - more CG, in this case, a minus sign is added in the ones "Lose CG" since in that case altNME has fewer CG than refNME
# variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dNME_relative=-variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dNME_relative
#Convert everything to gain CG
variant_HetCpG_meta_dt$CpG_change='Gain CG'
variant_HetCpG_meta_dt[((grepl('CG',REF_tri )) & (grepl('CG',ALT_tri)))|(!grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri))]$CpG_change='No CG change'
saveRDS(variant_HetCpG_meta_dt,'../downstream/output/human_analysis/CPEL_outputs/variant_HetCpG_meta_dt_relative_dNME2uq.rds')
variant_HetCpG_meta_dt=readRDS('../downstream/output/human_analysis/CPEL_outputs/variant_HetCpG_meta_dt_relative_dNME2uq.rds')
#substr(variant_HetCpG_meta_dt$tri_SNP,2,2)="X"
#For dNME
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
  for(tri in unique(variant_HetCpG_meta_dt[SNP==sn]$tri_SNP_unique)){
    variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
    variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
    #For C->G SNP, change background
    # if(sn != 'C->G'){
    #   #Write in method
    #   variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
    # variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
    # }else{
    #   #Separate G->C in SNP and C->G SNP
    #   #G->C SNP
    #   tri_l=gsub('->.*','',tri)
    #   tri_r=gsub('.*->','',tri)
    #   if(paste0(substring(tri_l, 2, 2),'->',substring(tri_r, 2, 2))=='G->C'){
    #     
    #     variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[paste0(substring(gsub('->.*','',tri_SNP_unique), 2, 2),'->',substring(gsub('.*->','',tri_SNP_unique), 2, 2))=='G->C' &
    #                                                         dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
    #     
    #   }else{
    #     variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[paste0(substring(gsub('->.*','',tri_SNP_unique), 2, 2),'->',substring(gsub('.*->','',tri_SNP_unique), 2, 2))=='C->G' &
    #                                                         dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
    #     
    #     
    #   }
      variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
      
      
    #}
    variant_SNP_tri=rbind(variant_SNP_tri,variant_SNP_tri_OR)
  }
  #get HetCpG
  variant_SNP_tri$CpG_change=factor(variant_SNP_tri$CpG_change,levels=c('Gain CG','No CG change'))
  variant_SNP_tri=variant_SNP_tri[order(OR,decreasing=F)]
  # variant_SNP_tri$SNP=paste0(gsub('.*->','',variant_SNP_tri$SNP),'->',gsub('->.*','',variant_SNP_tri$SNP))
  # sn=  paste0(gsub('.*->','',sn),'->',gsub('->.*','',sn))
  variant_SNP_tri$SNP=gsub('->','\u2794',variant_SNP_tri$SNP)
  variant_SNP_tri$SNP=factor(variant_SNP_tri$SNP,levels = variant_SNP_tri$SNP)
  
  variant_SNP_tri$FDR=p.adjust(variant_SNP_tri$pvalue,method='BH')
  variant_SNP_tri$significant=add.significance.stars(variant_SNP_tri$FDR, cutoffs = c(0.05, 0.01, 0.001))
  
  # SNP_all[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR))+geom_bar(fill=color_theme[sn],stat="identity")+ylab('Enrichment')+
  #   geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.1,position=position_dodge(.9))+ggtitle(sn)+xlab("")+
  # theme_glob+theme(legend.position = "none")+
  #   #theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #   ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+geom_text(aes(label=significant,y=upperCI*1+sig_h),vjust = sig_v,hjust=sig_h)+coord_flip()
   
   SNP_het[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=log(OR),fill=CpG_change))+geom_bar(stat="identity")+ylab('')+xlab("")+
    geom_errorbar(aes(ymin=log(lowerCI), ymax=log(upperCI)), width=.4,position=position_dodge(.9),size=0.25)+ggtitle(gsub('->',' \u2794 ',sn))+#ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+
    theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
     ylim(c(-1.5,1.5))+
     #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
     scale_fill_manual(values=c("No CG change"="grey","Gain CG"="light blue"))+
     geom_text(data=variant_SNP_tri[OR>1],aes(label=significant,y=log(upperCI)*1),vjust =sig_v,hjust=sig_h_pos)+
     geom_text(data=variant_SNP_tri[OR<1],aes(label=significant,y=log(lowerCI)*1),vjust =sig_v,hjust=sig_h_neg)+
     coord_flip()
   
   # SNP_box[[sn]]=ggplot(variant_HetCpG_meta_dt[SNP==sn & dNME_pval<=pval_cutoff],aes(x=reorder(tri_SNP_unique, -NME_relative, FUN = median) ,y=NME_relative,fill=CpG_change))+
   #   geom_boxplot()+
   #   ylab('refNME-altNME')+theme_glob+theme(legend.position = "bottom")+ ylim(c(-1,1))+ 
   #   scale_fill_manual(values=c("No CG change"="grey","Gain CG"="light blue"))+
   #   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
   #   xlab('')+coord_flip()
   variant_SNP_tri_out[[sn]]=variant_SNP_tri
   variant_SNP_tri=data.table()

}
# pdf('../downstream/output/graphs/Figure3/variant_OR_tri2.pdf',width=14,height=14)
# SNP_all=SNP_all[c("C>G", names(SNP_all)[names(SNP_all)!="C>G"])]
# ggarrange(plotlist=SNP_all, nrow=4,ncol=3)
# dev.off()
#After import_font, use loadfont(device='pdf')
# font_import('C:/Users/vince/Downloads/dejavu-fonts-ttf-2.37')
# loadfonts(device='pdf')
#cairo_pdf('../downstream/output/graphs/Figure2/Figure-4B-variant_OR_tri3_two_cat_greater_CG_bg.pdf',width=10,height=7,family = "DejaVu Sans")
png('../downstream/output/graphs/Figure2/Figure-4B-variant_OR_tri3_two_cat_greater_CG_bg_rev.png',width=7,height=7,units='in',res=1080, family = 'Consolas')
#SNP_het=SNP_het[c("C>G", names(SNP_het)[names(SNP_het)!="C>G"])]
ggarrange(plotlist=SNP_het, nrow=2,ncol=2,common.legend = T,legend="bottom",label.x = 'Odds ratio')
dev.off()
# pdf('../downstream/output/graphs/Figure3/variant_box_tri.pdf',width=14,height=14)
# SNP_box=SNP_box[c("C>G", names(SNP_box)[names(SNP_box)!="C>G"])]
# ggarrange(plotlist=SNP_box, nrow=4,ncol=3,common.legend = T,legend="bottom")
# dev.off()
#Two catogries do not need this since they're complmentary
OR_all_SNP_change=data.table()

for(sn in unique(variant_HetCpG_meta_dt$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc(variant_HetCpG_meta_dt[dNME_pval<=pval_cutoff],sn,"CpG_change"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_2cat.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

#For dMML

#Calculate relative dMML
variant_HetCpG_meta_dt=readRDS('../downstream/output/human_analysis/CPEL_outputs/variant_HetCpG_meta_dt.rds')
single_SNP_unique=unique_mutation(unique(variant_HetCpG_meta_dt$SNP))
variant_HetCpG_meta_dt$SNP=single_SNP_unique[variant_HetCpG_meta_dt$SNP]
mutation_tri_unique=unique_mutation(unique(variant_HetCpG_meta_dt$tri_SNP))
#Finding the ones losing CG
gainCG_idx=which(!grepl("CG",sub('.*-','',mutation_tri_unique))&grepl("CG",sub('-.*','',mutation_tri_unique)))
mutation_tri_unique[gainCG_idx]=paste0(sub('.*-','', mutation_tri_unique[gainCG_idx]),'-',sub('-.*','', mutation_tri_unique[gainCG_idx]))
variant_HetCpG_meta_dt$tri_SNP_unique=mutation_tri_unique[variant_HetCpG_meta_dt$tri_SNP]
#Lumping all
variant_HetCpG_meta_dt$CpG_change="No change"
variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri ))]$CpG_change='Lose CG'
variant_HetCpG_meta_dt[(!grepl('CG',REF_tri )) & (grepl('CG',ALT_tri ))]$CpG_change='Gain CG'
variant_HetCpG_meta_dt[(grepl('CG',REF_tri )) & (grepl('CG',ALT_tri))]$CpG_change='No change'
variant_HetCpG_meta_dt$dMML_relative=variant_HetCpG_meta_dt$refMML-variant_HetCpG_meta_dt$altMML
#convert all to less CG - more CG, in this case, a minus sign is added in the ones "Lose CG" since in that case altNME has fewer CG than refNME
variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dMML_relative=-variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$dMML_relative
#Convert everything to gain CG
variant_HetCpG_meta_dt[CpG_change=="Lose CG"]$CpG_change='Gain CG'
SNP_all=list()
SNP_het=list()
SNP_box=list()
sig_v=0.75
sig_h=-0.05
color_theme=c(rainbow(length(unique(variant_HetCpG_meta_dt$SNP))))
variant_SNP_tri=data.table()
names(color_theme)=unique(variant_HetCpG_meta_dt$SNP)
for (sn in unique(variant_HetCpG_meta_dt$SNP)){
  for(tri in unique(variant_HetCpG_meta_dt[SNP==sn]$tri_SNP_unique)){
    variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dMML_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff,stat_in="MML")
    variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
    variant_SNP_tri=rbind(variant_SNP_tri,variant_SNP_tri_OR)
  }
  #get HetCpG
  variant_SNP_tri$CpG_change=factor(variant_SNP_tri$CpG_change,levels=c('Gain CG','No change','Lose CG'))
  variant_SNP_tri=variant_SNP_tri[order(OR,decreasing=F)]
  variant_SNP_tri$SNP=factor(variant_SNP_tri$SNP,levels = variant_SNP_tri$SNP)
  
  variant_SNP_tri$FDR=p.adjust(variant_SNP_tri$pvalue,method='BH')
  variant_SNP_tri$significant=add.significance.stars(variant_SNP_tri$FDR, cutoffs = c(0.05, 0.01, 0.001))
  variant_SNP_tri=variant_SNP_tri[!is.infinite(OR)]
  SNP_all[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR))+geom_bar(fill=color_theme[sn],stat="identity")+ylab('Enrichment')+
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.1,position=position_dodge(.9))+ggtitle(sn)+xlab("")+
    theme_glob+theme(legend.position = "none")+
    #theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+geom_text(aes(label=significant,y=upperCI*1+sig_h),vjust = sig_v,hjust=sig_h)+coord_flip()
  
  SNP_het[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=OR,fill=CpG_change))+geom_bar(stat="identity")+ylab('Enrichment')+xlab("")+
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.4,position=position_dodge(.9),size=0.25)+ggtitle(sn)+ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+
    theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+geom_text(aes(label=significant,y=upperCI*1),vjust =sig_v,hjust=sig_h)+coord_flip()
  

  
  
  variant_SNP_tri=data.table()
}
# pdf('../downstream/output/graphs/Figure3/variant_OR_tri2.pdf',width=14,height=14)
# SNP_all=SNP_all[c("C>G", names(SNP_all)[names(SNP_all)!="C>G"])]
# ggarrange(plotlist=SNP_all, nrow=4,ncol=3)
# dev.off()

pdf('../downstream/output/graphs/Figure2/Figure-S3-variant_OR_tri3_dMML_cat.pdf',width=10,height=7)
#SNP_het=SNP_het[c("C>G", names(SNP_het)[names(SNP_het)!="C>G"])]
ggarrange(plotlist=SNP_het, nrow=2,ncol=3,common.legend = T,legend="bottom")
dev.off()
# pdf('../downstream/output/graphs/Figure3/variant_box_tri.pdf',width=14,height=14)
# SNP_box=SNP_box[c("C>G", names(SNP_box)[names(SNP_box)!="C>G"])]
# ggarrange(plotlist=SNP_box, nrow=4,ncol=3,common.legend = T,legend="bottom")
# dev.off()
#All
OR_all_SNP_change=data.table()

for(sn in unique(variant_HetCpG_meta_dt$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc(variant_HetCpG_meta_dt[dMML_pval<=pval_cutoff],sn,"CpG_change",stat_in="MML"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_dMML.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

#CpG island: -0.189
OR_all_SNP_change=data.table()

for(sn in unique( variant_HetCpG_meta_dt[CpG_island==TRUE]$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc( variant_HetCpG_meta_dt[CpG_island==TRUE][dMML_pval<=pval_cutoff],sn,"CpG_change",stat_in="MML"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_dMML_island.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

#CpG island: 0.333
OR_all_SNP_change=data.table()

for(sn in unique(variant_HetCpG_meta_dt[CpG_island==FALSE]$CpG_change)){
  OR_all_SNP_change=rbind(OR_all_SNP_change,OR_calc(variant_HetCpG_meta_dt[CpG_island==FALSE][dMML_pval<=pval_cutoff],sn,"CpG_change",stat_in="MML"))
  
}
OR_all_SNP_change=OR_all_SNP_change[order(OR,decreasing=F)]
OR_all_SNP_change$SNP=factor(OR_all_SNP_change$SNP,levels=OR_all_SNP_change$SNP)
OR_all_SNP_change$FDR=p.adjust(OR_all_SNP_change$pvalue,method='BH')
OR_all_SNP_change$sig=add.significance.stars(OR_all_SNP_change$FDR, cutoffs = c(0.05, 0.01, 0.001))
OR_all_SNP_change[SNP=="Gain CG"]$OR/OR_all_SNP_change[SNP=="Lose CG"]$OR
pdf('../downstream/output/graphs/Figure2/variant_OR_tri_all_dMML_non_island.pdf',width=7,height=7)
ggplot(OR_all_SNP_change,aes(x=SNP,y=OR,fill=SNP))+geom_bar(stat="identity")+ylab('OR (REF > ALT)')+xlab("")+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=.2,position=position_dodge(.9))+ylim(c(0,2.5))+
  theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=sig,y=upperCI*1),vjust = -0.5)+ scale_fill_manual(values=c("Lose CG"="red","No change"="grey","Gain CG"="blue"))+ggtitle("")+
  coord_flip()
dev.off()

# Density analysis using allele-specific way ------------------------------


GR_merge=readRDS(GR_merge_file)

GR_merge$CpGdiff=GR_merge$g1CG-GR_merge$g2CG

#plot_CpG_number(GR_merge[GR_merge$dNME_pval<=pval_cutoff])

#CpG density vs dNME,here uses an extended density 
GR_merge$dNME_relative=GR_merge$NME1-GR_merge$NME2
GR_merge$dMML_relative=GR_merge$MML1-GR_merge$MML2
olap=findOverlaps(GR_merge,genomic_features$`CpG island`)
GR_merge_dt=convert_GR(GR_merge,'DT')
#abs difference
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
#GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1/CG_allele_extend_g2)]
#dNME
#ratio

cor.test(GR_merge_dt[dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$dNME_relative, 
         GR_merge_dt[dNME_pval<=pval_cutoff&dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$density_diff)


GR_merge_dt$CpG_stat="No difference"
GR_merge_dt[CpGdiff!=0]$CpG_stat="With CpG difference"
GR_merge_dt$CpG_stat=factor(GR_merge_dt$CpG_stat,levels = c("With CpG difference","No difference"))
#Make this one better
#Figure 3A
GR_merge_dt$dNME_relative_more_less=GR_merge_dt$dNME_relative
GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative_more_less=GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative*sign(GR_merge_dt[GR_merge_dt$CpGdiff!=0]$CpGdiff)
t.test(GR_merge_dt[CpGdiff!=0&dNME_pval<=pval_cutoff]$dNME_relative_more_less,alternative="less")
pdf('../downstream/output/graphs/Figure3/Figure3A_CpG_number_NME.pdf',width=7,height=7)
ggplot(GR_merge_dt[dNME_pval<=pval_cutoff],aes(y=dNME_relative_more_less,x=CpG_stat,fill=CpG_stat))+
  geom_violin()+xlab("")+
  theme_glob+ylab('relative dNME')+theme(legend.position = "none")

dev.off()
GR_merge_dt_sig_density_diff=GR_merge_dt[dNME_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_sig_density_diff$density_difference_quantile=ecdf(GR_merge_dt_sig_density_diff$density_diff)(GR_merge_dt_sig_density_diff$density_diff)
pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dNME_ratio.pdf',width=7,height=7)

ggplot(GR_merge_dt_sig_density_diff,aes(x=density_difference_quantile,y=dNME_relative))+geom_smooth(fill='light blue')+
  xlab("CpG density ratio quantile")+ylab("relative dNME")+
  #stat_summary(fun=median, geom="point")+
  theme_glob+
  #ylim(c(-0.5,0.3))+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# ggplot(GR_merge_dt_sig_density_diff,aes(x=as.factor(round(density_diff,digits = 2)),y=dNME_relative))+geom_violin(fill='light blue')+
#   xlab("CpG density ratio")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
#   theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
GR_merge_dt_sig_density_diff=GR_merge_dt[dMML_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_sig_density_diff$density_difference_quantile=ecdf(GR_merge_dt_sig_density_diff$density_diff)(GR_merge_dt_sig_density_diff$density_diff)
pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dMML_ratio.pdf',width=7,height=7)

ggplot(GR_merge_dt_sig_density_diff,aes(x=density_difference_quantile,y=dMML_relative))+geom_smooth(fill='light blue')+
  xlab("CpG density ratio quantile")+ylab("relative dMML")+
  #stat_summary(fun=median, geom="point")+
  theme_glob+
 # ylim(c(-0.5,0.3))+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# ggplot(GR_merge_dt_sig_density_diff,aes(x=as.factor(round(density_diff,digits = 2)),y=dNME_relative))+geom_violin(fill='light blue')+
#   xlab("CpG density ratio")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
#   theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dNME_quantile.pdf',width=7,height=7)
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
GR_merge_dt_S4=GR_merge_dt[dNME_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_S4$density_diff=findInterval(GR_merge_dt_S4$density_diff,quantile(GR_merge_dt_S4$density_diff,prob=seq(0,0.9,0.1)))
GR_merge_dt_S4$density_diff=factor(paste0(GR_merge_dt_S4$density_diff*10,'%'),levels=paste0( seq(0,100,10),"%"))
ggplot(GR_merge_dt_S4,aes(x=density_diff,y=dNME_relative))+geom_violin(fill='light blue')+
  xlab("CpG density difference quantile")+ylab("relative dNME")+stat_summary(fun=median, geom="point")+theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#dMML

#0.214, significant
cor.test(GR_merge_dt[dMML_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$dMML_relative, 
         GR_merge_dt[dMML_pval<=pval_cutoff&dMML_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$density_diff)

GR_merge_dt$CpG_stat="No difference"
GR_merge_dt[CpGdiff!=0]$CpG_stat="With CpG difference"
GR_merge_dt$CpG_stat=factor(GR_merge_dt$CpG_stat,levels = c("With CpG difference","No difference"))
#Make this one better
#Figure 3A
GR_merge_dt$dMML_relative_more_less=GR_merge_dt$dMML_relative
GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dMML_relative_more_less=GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dMML_relative*sign(GR_merge_dt[GR_merge_dt$CpGdiff!=0]$CpGdiff)
t.test(GR_merge_dt[CpGdiff!=0&dMML_pval<=pval_cutoff]$dMML_relative_more_less,alternative="less")
pdf('../downstream/output/graphs/Figure3/Figure3A_CpG_number_MML.pdf',width=7,height=7)
ggplot(GR_merge_dt[dMML_pval<=pval_cutoff],aes(y=dMML_relative_more_less,x=CpG_stat,fill=CpG_stat))+
  geom_violin()+xlab("")+
  theme_glob+ylab('relative dMML')+theme(legend.position = "none")

dev.off()
pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dMML_ratio.pdf',width=7,height=7)
ggplot(GR_merge_dt[dMML_pval<=pval_cutoff&density_diff!=1],aes(x=as.factor(round(density_diff,digits = 2)),y=dMML_relative))+geom_violin(fill='light blue')+
  xlab("CpG density ratio")+ylab("relative dMML")+stat_summary(fun=median, geom="point")+theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('../downstream/output/graphs/Figure3/FigureS4_CpG_density_dMML_quantile.pdf',width=7,height=7)
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
GR_merge_dt_S4=GR_merge_dt[dMML_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_S4$density_diff=findInterval(GR_merge_dt_S4$density_diff,quantile(GR_merge_dt_S4$density_diff,prob=seq(0,0.9,0.1)))
GR_merge_dt_S4$density_diff=factor(paste0(GR_merge_dt_S4$density_diff*10,'%'),levels=paste0( seq(0,100,10),"%"))
ggplot(GR_merge_dt_S4,aes(x=density_diff,y=dMML_relative))+geom_violin(fill='light blue')+
  xlab("CpG density difference quantile")+ylab("relative dMML")+stat_summary(fun=median, geom="point")+theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# ggplot(GR_merge_dNME[dNME_pval<=pval_cutoff&CpGdiff!=0],aes(x=dNME_relative_more_less))+geom_density()+xlab("relative dNME")+theme_glob
# ggplot(GR_merge_dNME[dNME_pval<=pval_cutoff&CpGdiff!=0],aes(y=dNME_relative_more_less,x=as.factor(density_diff)))+
#   geom_boxplot()+xlab("CpG difference")+
#   theme_glob+ylab('relative dNME')

# # check NME vs NME from allele specific way -------------------------------
GR_merge_tb=readRDS('../downstream/output/GR_merge_ASM_comp.rds')
GR_merge_tb=GR_merge_tb[!is.na(MML_ASM)&!is.na(MML)&!is.na(NME)]#3284912/52263042,total ASM regions=3332744
pdf('../downstream/output/graphs/Figure3/Fig-S5A-dNME_NME_all_pt.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb,aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb$NME_ASM,GR_merge_tb$NME)
pdf('../downstream/output/graphs/Figure3/Fig-S5B-dNME_NME_all_dMML_ASM_pt.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb[dMML_pval<=pval_cutoff],aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb[dMML_pval<=pval_cutoff]$NME_ASM,GR_merge_tb[dMML_pval<=pval_cutoff]$NME)
pdf('../downstream/output/graphs/Figure3/Fig-S5C-dNME_NME_non_dMML.pdf',width=3.5,height=3.5)
ggplot(GR_merge_tb[dMML_pval>pval_cutoff],aes(x=NME_ASM,y=NME))+geom_bin2d(bins=200)+xlab('mean allelic NME')+ylab('allele-agnostic NME')+
  geom_abline(slope=1,size=1,color='red')+ xlim(0,1)+ylim(0,1)+theme_glob+theme(legend.position = 'bottom')+  
  scale_fill_gradient2(high = "darkblue",low="lightblue")
dev.off()
cor.test(GR_merge_tb[dMML_pval>pval_cutoff]$NME_ASM,GR_merge_tb[dMML_pval>pval_cutoff]$NME)


# Figure 3B: allele-agnotic density ---------------------------------------

NME_in=readRDS(NME_agnostic_file)
CpG_hg19=readRDS('../downstream/input/human_analysis/CpG_hg19.rds')
NME_in$CG_hg19=countOverlaps(NME_in,CpG_hg19)
NME_in_gr=unique(granges(NME_in))
gr_seq=getSeq(Hsapiens,NME_in_gr,as.character=T)
NME_in_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
NME_in_olap=findOverlaps(NME_in,NME_in_gr,type='equal')
NME_in$CGcont_exp[queryHits(NME_in_olap)]=NME_in_gr$CGcont_exp[subjectHits(NME_in_olap)]
NME_in$density=NME_in$CG_hg19/NME_in$CGcont_exp
NME_in=readRDS('../downstream/input/human_analysis/NME_allele_agnostic_merge_20k_homogeneous_excluding_dMML2_CG.rds')
NME_in=NME_in[seqnames(NME_in)%in%paste0("chr",1:22)]
cor.test(NME_in$density,NME_in$NME,method='pearson')
#Make boxplot of this
#quantile(NME_in$density,prob=seq(0,1,0.2),na.rm=T)
NME_in$density_quant=findInterval(NME_in$density,seq(0,1,0.1))
#NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
NME_in$density_quant=factor(quant_conv[NME_in$density_quant],levels=quant_conv)
pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot_CG_exp.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in)),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()   

png('../downstream/output/graphs/Figure3/FigureS_CpG_density_NME_line.png',width=1080,height=1080)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in[sample(1:length(NME_in),round(length(NME_in)/5))])),aes(x=density, y=NME))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
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

pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot_CG_exp_non_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in[-queryHits(olap_islands)])),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()   

pdf('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME_boxplot_CG_exp_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(NME_in[queryHits(olap_islands)])),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off() 
png('../downstream/output/graphs/FigureS6/NME_density_feature_box.png')
# my_comparisons <- list( c("islands", "shores"), c("islands", "shelf"), c("islands", "open sea"),
#                         c("shores", "shelf"),c("shores", "open sea"),c("shelf", "open sea"))
ggplot(CpG_density_NME,aes(x=feature,y=NME))+
  geom_boxplot()+
  #geom_violin()+ 
  theme_glob+xlab("genomic features")
#stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()

# #same analysis on dMML-ASM
MML_in=readRDS(MML_agnostic_file)
MML_in$CG_hg19=countOverlaps(MML_in,CpG_hg19)
MML_in_gr=unique(granges(MML_in))
gr_seq=getSeq(Hsapiens,MML_in_gr,as.character=T)
MML_in_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
#Do it for MML
MML_in_olap=findOverlaps(MML_in,MML_in_gr,type='equal')
MML_in$CGcont_exp[queryHits(MML_in_olap)]=MML_in_gr$CGcont_exp[subjectHits(MML_in_olap)]
MML_in$density=MML_in$CG_hg19/MML_in$CGcont_exp
MML_in=readRDS('../downstream/input/human_analysis/MML_allele_agnostic_merge_20k_homogeneous_excluding_dMML2_CG.rds')
# MML_in$density=MML_in$CG_hg19/width(MML_in)
# MML_in=MML_in[seqnames(MML_in)%in%paste0("chr",1:22)]
# MML_in$density_quant=findInterval(MML_in$density,seq(0,0.1,0.01))
# MML_in$density_quant=factor(quant_conv[MML_in$density_quant],levels=quant_conv)
MML_in$density_quant=findInterval(MML_in$density,seq(0,1,0.1))
#NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
MML_in$density_quant=factor(quant_conv[MML_in$density_quant],levels=quant_conv)
olap=findOverlaps(MML_in,genomic_features$`CpG island`)
cor.test(MML_in$density,MML_in$MML,method='pearson')#-0.346
cor.test(MML_in[queryHits(olap)]$density,MML_in[queryHits(olap)]$MML,method='pearson')#-0.395
cor.test(MML_in[-queryHits(olap)]$density,MML_in[-queryHits(olap)]$MML,method='pearson')#0.0059
#Inverse correlation between MML an d CpG content: 19325872 [-0.46]
pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in)),aes(x=density_quant, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in[queryHits(olap)])),aes(x=density_quant, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot_non_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in[-queryHits(olap)])),aes(x=density_quant, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_line.png',width=1080,height=1080)#Totally having 69530406 points
ggplot(as.data.frame(mcols(MML_in[MML_in$density<=1][sample(1:length(MML_in[MML_in$density<=1]),round(length(MML_in)/3))])),aes(x=density, y=MML))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Figure S5
olap_islands=findOverlaps(MML_in,genomic_features$`CpG island`)
olap_shores=findOverlaps(MML_in,genomic_features$`CpG shore`)
olap_shelf=findOverlaps(MML_in,genomic_features$`CpG shelf`)
olap_open_sea=findOverlaps(MML_in,genomic_features$`CpG open sea`)

CpG_density_MML=rbind(data.table(MML=MML_in$MML[queryHits(olap_islands)],feature='islands'),
                      data.table(MML=MML_in$MML[queryHits(olap_shores)],feature='shores'),
                      data.table(MML=MML_in$MML[queryHits(olap_shelf)],feature='shelf'),
                      data.table(MML=MML_in$MML[queryHits(olap_open_sea)],feature='open sea'))



pdf('../downstream/output/graphs/FigureS6/MML_density_feature_box.pdf')
# my_comparisons <- list( c("islands", "shores"), c("islands", "shelf"), c("islands", "open sea"),
#                         c("shores", "shelf"),c("shores", "open sea"),c("shelf", "open sea"))
ggplot(CpG_density_MML,aes(x=feature,y=MML))+
  geom_boxplot(outlier.shape = NA)+
  #geom_violin()+ 
  theme_glob+xlab("genomic features")
#stat_compare_means(comparisons=my_comparisons,method = "wilcox.test")
dev.off()

# In mouse context --------------------------------------------------------
mml=readRDS('../downstream/input/mouse_analysis/MML_matrix_mouse_all_dedup_N2.rds')
nme=readRDS('../downstream/input/mouse_analysis/NME_matrix_mouse_all_dedup_N2.rds')
mm10_CpG=getCpgSitesmm10()
nme$CG_mm10=countOverlaps(nme,mm10_CpG)
mml$CG_mm10=countOverlaps(mml,mm10_CpG)
nme_gr=unique(c(granges(nme),granges(mml)))
gr_seq=getSeq(Mmusculus,nme_gr,as.character=T)
nme_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
saveRDS(nme_gr,'../downstream/output/mouse_analysis/region_CG_exp.rds')
nme_gr=readRDS('../downstream/output/mouse_analysis/region_CG_exp.rds')
nme_olap=findOverlaps(nme,nme_gr,type='equal')
nme$CGcont_exp[queryHits(nme_olap)]=nme_gr$CGcont_exp[subjectHits(nme_olap)]
nme$density=nme$CG_mm10/nme$CGcont_exp
mml_olap=findOverlaps(mml,nme_gr,type='equal')
mml$CGcont_exp[queryHits(mml_olap)]=nme_gr$CGcont_exp[subjectHits(mml_olap)]
mml$density=mml$CG_mm10/mml$CGcont_exp
nme=nme[seqnames(nme) %in% paste0("chr",1:19)]
mml=mml[seqnames(mml) %in% paste0("chr",1:19)]
nme$density_quant=findInterval(nme$density,seq(0,1,0.1))
mml$density_quant=findInterval(mml$density,seq(0,1,0.1))
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
nme$density_quant=factor(quant_conv[nme$density_quant],levels=quant_conv)

mml$density_quant=factor(quant_conv[mml$density_quant],levels=quant_conv)
enhancer=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')
tss=get_mm10_tss()
pdf('../downstream/output/graphs/Figure3/mouse_MML_density_enhancer_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(mml,enhancer),stat_name="MML")
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_MML_density_promoter_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(mml,tss,maxgap = 2000),stat_name="MML")
dev.off()

pdf('../downstream/output/graphs/Figure3/mouse_NME_density_enhancer_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(nme,enhancer),stat_name="NME")
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_NME_density_promoter_boxplot.pdf',width=3.5,height=3.5)#Totally having 69530406 points
density_mouse_calc(subsetByOverlaps(nme,tss,maxgap = 2000),stat_name="NME")
dev.off()

#Looking only looking at ones in cluster result
cluster_out=readRDS('../downstream/input/mouse_analysis/clustering/tissue_specific/kmeans_10run/uc_11.rds')

nme_dt=convert_GR(nme,direction="DT")
nme_dt=melt.data.table(nme_dt,id.vars = c("CG_mm10","CGcont_exp","density","region","density_quant"),variable.name = "Sample",value.name="stat_in")
nme_dt=nme_dt[!is.na(stat_in)]
mml_dt=convert_GR(mml,direction="DT")
mml_dt=melt.data.table(mml_dt,id.vars = c("CG_mm10","CGcont_exp","density","region","density_quant"),variable.name = "Sample",value.name="stat_in")
mml_dt=mml_dt[!is.na(stat_in)]
nme_dt_clu=data.table()
mml_dt_clu=data.table()
for(ts in names(cluster_out)){
  nme_dt_clu=rbind(nme_dt_clu,nme_dt[grepl(ts,Sample) & (region %in%names(cluster_out[[ts]]))])
  mml_dt_clu=rbind(mml_dt_clu,mml_dt[grepl(ts,Sample) & (region %in%names(cluster_out[[ts]]))])
  
}


pdf('../downstream/output/graphs/Figure3/mouse_MML_density_enhancer_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
enhancer_olap=findOverlaps(convert_GR(mml_dt_clu$region),enhancer)
print(ggplot(mml_dt_clu[queryHits(enhancer_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_MML_density_promoter_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
promoter_olap=findOverlaps(convert_GR(mml_dt_clu$region),tss,maxgap=2000)
print(ggplot(mml_dt_clu[queryHits(promoter_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

pdf('../downstream/output/graphs/Figure3/mouse_NME_density_enhancer_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
enhancer_olap=findOverlaps(convert_GR(nme_dt_clu$region),enhancer)
print(ggplot(nme_dt_clu[queryHits(enhancer_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
pdf('../downstream/output/graphs/Figure3/mouse_NME_density_promoter_boxplot_clu.pdf',width=3.5,height=3.5)#Totally having 69530406 points
promoter_olap=findOverlaps(convert_GR(nme_dt_clu$region),tss,maxgap=2000)
print(ggplot(nme_dt_clu[queryHits(promoter_olap)],aes(x=density_quant, y=stat_in))+
        ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
        ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

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

# pdf('../downstream/output/graphs/Figure3/FigureS_CpG_density_MML_boxplot_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
# ggplot(as.data.frame(mcols(MML_in[queryHits(olap)])),aes(x=density_quant, y=MML))+
#   ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
#   ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()  
# pdf('../downstream/output/graphs/Figure3/FigureS6_CpG_density_MML_boxplot_non_island.pdf',width=3.5,height=3.5)#Totally having 69530406 points
# 
# ggplot(as.data.frame(mcols(MML_in[-queryHits(olap)])),aes(x=density_quant, y=MML))+
#   ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
#   ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()  


# Plot density vs NME with points -----------------------------------------
# NME_in_dt=data.table(NME=NME_in$NME,density=NME_in$density)
# digits_round=4
# NME_in_dt_agg=NME_in_dt[, list(NME=median(NME),Bottom25=quantile(NME,probs=0.25),
#                                    top25=quantile(NME,probs=0.75)), 
#                             by = list(density = round(density,digits=digits_round))]
# NME_in_dt_agg$Bottom25= predict(loess(Bottom25~density,NME_in_dt_agg),newdata=NME_in_dt_agg$density)
# NME_in_dt_agg$top25= predict(loess(top25~density,NME_in_dt_agg),newdata=NME_in_dt_agg$density)
# png('../downstream/output/graphs/Figure3/Figure3D_CpG_density_NME.png',width=1080,height=1080)#Totally having 69530406 points
# ggplot(NME_in_dt_agg,aes(x=density, y=NME))+
#   ylim(c(0,1))+geom_smooth(method="loess",se=FALSE)+theme_glob+xlab("CpG density")+
#   ylab("NME")+geom_ribbon(aes(ymin = Bottom25, ymax = top25), fill = "blue", alpha = .1)+
#   theme(axis.title.x=element_text(hjust=0.5,size=48,face="bold"),
#         axis.title.y=element_text(hjust=0.5,size=48,face="bold"),
#         axis.text.x=element_text(size=46),
#         axis.text.y=element_text(size=46))+
#   geom_point(data=NME_in_dt,size=0.1,aes(x=density,y=NME),alpha=0.1)
# dev.off()              


# #LiftOver
# GR_merge=readRDS(GR_merge_file)
# ch = import.chain('../downstream/data/hg19ToMm10.over.chain')
# cur=granges(unique(GR_merge[GR_merge$dMML_pval<=pval_cutoff]))
# seqlevelsStyle(cur) = "UCSC"  # necessary
# cur19 = unlist(liftOver(cur, ch))
# overlap=list()
# overlap_dat=data.table()
# for(fn in dir('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/',pattern='all.csv')){
#   csv_in=fread(paste0('../downstream/output/mm10_result/chromHMM_enhancer/all_gene_list/',fn))
#   region_in=convert_GR(unique(csv_in$region))
#   ts=gsub('_all.csv','',fn)
#   region_olap=subsetByOverlaps(region_in,cur19)
#   region_olap=paste0(seqnames(region_olap),':',start(region_olap),'-',end(region_olap))
#   overlap[[ts]]=csv_in[region%in%region_olap]
#   overlap_dat=rbind(overlap_dat,data.table(ts=ts,total_regions=length(region_in),overlap=length(region_olap),
#                                            olap_GO=sum(region_olap_GO%in%csv_in[GO_result!=""]$region)))
# }

#Find which one lose CG increase Entropy and motif prefer higher ent
motif_gene <- readRDS(motif_gene_file)

olap=findOverlaps(variant_HetCpG_meta,motif_gene)
variant_HetCpG_meta_motif=variant_HetCpG_meta[queryHits(olap)]
variant_HetCpG_meta_motif$ref_score=motif_gene$scoreRef[subjectHits(olap)]
variant_HetCpG_meta_motif$alt_score=motif_gene$scoreRef[queryHits(olap)]
variant_HetCpG_meta_lose_CG=variant_HetCpG_meta_motif[grepl('CG',variant_HetCpG_meta_motif$REF_tri)&!(grepl('CG',variant_HetCpG_meta_motif$ALT_tri))]

binom.test(sum(variant_HetCpG_meta_lose_CG$ref_score<variant_HetCpG_meta_lose_CG$alt_score),length(variant_HetCpG_meta_lose_CG),
          p=sum(variant_HetCpG_meta_motif$ref_score<variant_HetCpG_meta_motif$alt_score)/length(variant_HetCpG_meta_motif))
#MML and NME curve for 
MML_in=readRDS('../downstream/input/MML_allele_agnostic_merge_20k_homogeneous_excluding_dMML2_CG.rds')
NME_in=readRDS('../downstream/input/NME_allele_agnostic_merge_20k_homogeneous_excluding_dMML2_CG.rds')
NME_in$stat_type="NME"
MML_in$stat_type="MML"
plot_dt=as.data.table(rbind(mcols(NME_in)[,c("score","density","stat_type")],mcols(MML_in)[,c("score","density","stat_type")]))
saveRDS(plot_dt,'../downstream/output/plot_dt_NME_MML.rds')

rm(MML_in)
rm(NME_in)
plot_dt=readRDS('../downstream/output/plot_dt_NME_MML.rds')
pdf('../downstream/output/FigureS_CpG_density_MML_NME_line_density_small1.pdf',width=7,height=7)#Totally having 69530406 points
ggplot(plot_dt[density<=1],aes(x=density, y=score,color=stat_type))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="bottom")
dev.off()
pdf('../downstream/output/FigureS_CpG_density_MML_NME_line_density_all.pdf',width=7,height=7)#Totally having 69530406 points
ggplot(plot_dt,aes(x=density, y=score,color=stat_type))+
  ylim(c(0,1))+geom_smooth()+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="bottom")
dev.off()