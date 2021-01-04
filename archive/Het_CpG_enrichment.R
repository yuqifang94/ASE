
##CpG open sea##
subjects=unique(NME_all$Subject)
open_sea_OR=data.frame(subject=subjects,OR=0,CpG_type='All regions')
for(subj in subjects){
  open_sea_OR$OR[open_sea_OR$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Subject==subj],genomic_features$`CpG open sea`,'dNME')$estimate
}  
#With het CpG
open_sea_OR_het=data.frame(subject=subjects,OR=0,CpG_type='Het CpG')
for(subj in subjects){
  open_sea_OR_het$OR[open_sea_OR_het$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Subject==subj & NME_all$CpGdiff!=0],genomic_features$`CpG open sea`,'dNME')$estimate
}  
#Without het CpG
open_sea_OR_nonhet=data.frame(subject=subjects,OR=0,CpG_type='Non Het CpG')
for(subj in subjects){
  open_sea_OR_nonhet$OR[open_sea_OR_nonhet$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Subject==subj & NME_all$CpGdiff==0],genomic_features$`CpG open sea`,'dNME')$estimate
} 
open_sea=rbind(open_sea_OR,open_sea_OR_het,open_sea_OR_nonhet)
ggplot(open_sea,aes(x=subject,y=OR,fill=CpG_type)) + ylim(0,1.5)+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle('Subject OR for open sea')

##CpG islands##
subjects=unique(NME_all$Subject)
islands_OR=data.frame(subject=subjects,OR=0,CpG_type='All regions')
for(subj in subjects){
  islands_OR$OR[islands_OR$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Subject==subj],genomic_features$`CpG island`,'dNME')$estimate
}  
#With het CpG
islands_OR_het=data.frame(subject=subjects,OR=0,CpG_type='Het CpG')
for(subj in subjects){
  islands_OR_het$OR[islands_OR_het$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Subject==subj & NME_all$CpGdiff!=0],genomic_features$`CpG island`,'dNME')$estimate
}  
#Without het CpG
islands_OR_nonhet=data.frame(subject=subjects,OR=0,CpG_type='Non Het CpG')
for(subj in subjects){
  islands_OR_nonhet$OR[islands_OR_nonhet$subject==subj]=testEnrichmentFeature_stat(NME_all[NME_all$Subject==subj & NME_all$CpGdiff==0],genomic_features$`CpG island`,'dNME')$estimate
} 
islands=rbind(islands_OR,islands_OR_het,islands_OR_nonhet)
ggplot(islands,aes(x=subject,y=OR,fill=CpG_type)) + ylim(0,1.5)+
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle('Subject OR for CpG islands')

##CpG open sea##
subjects=unique(MML_all$Subject)
open_sea_OR=data.frame(subject=subjects,OR=0,CpG_type='All regions')
for(subj in subjects){
  open_sea_OR$OR[open_sea_OR$subject==subj]=testEnrichmentFeature_stat(MML_all[MML_all$Subject==subj],genomic_features$`CpG open sea`,'dMML')$estimate
}  
#With het CpG
open_sea_OR_het=data.frame(subject=subjects,OR=0,CpG_type='Het CpG')
for(subj in subjects){
  open_sea_OR_het$OR[open_sea_OR_het$subject==subj]=testEnrichmentFeature_stat(MML_all[MML_all$Subject==subj & MML_all$CpGdiff!=0],genomic_features$`CpG open sea`,'dMML')$estimate
}  
#Without het CpG
open_sea_OR_nonhet=data.frame(subject=subjects,OR=0,CpG_type='Non Het CpG')
for(subj in subjects){
  open_sea_OR_nonhet$OR[open_sea_OR_nonhet$subject==subj]=testEnrichmentFeature_stat(MML_all[MML_all$Subject==subj & MML_all$CpGdiff==0],genomic_features$`CpG open sea`,'dMML')$estimate
} 
open_sea=rbind(open_sea_OR,open_sea_OR_het,open_sea_OR_nonhet)
ggplot(open_sea,aes(x=subject,y=OR,fill=CpG_type)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle('Subject OR for open sea')

##CpG islands##
subjects=unique(MML_all$Subject)
islands_OR=data.frame(subject=subjects,OR=0,CpG_type='All regions')
for(subj in subjects){
  islands_OR$OR[islands_OR$subject==subj]=testEnrichmentFeature_stat(MML_all[MML_all$Subject==subj],genomic_features$`CpG island`,'dMML')$estimate
}  
#With het CpG
islands_OR_het=data.frame(subject=subjects,OR=0,CpG_type='Het CpG')
for(subj in subjects){
  islands_OR_het$OR[islands_OR_het$subject==subj]=testEnrichmentFeature_stat(MML_all[MML_all$Subject==subj & MML_all$CpGdiff!=0],genomic_features$`CpG island`,'dMML')$estimate
}  
#Without het CpG
islands_OR_nonhet=data.frame(subject=subjects,OR=0,CpG_type='Non Het CpG')
for(subj in subjects){
  islands_OR_nonhet$OR[islands_OR_nonhet$subject==subj]=testEnrichmentFeature_stat(MML_all[MML_all$Subject==subj & MML_all$CpGdiff==0],genomic_features$`CpG island`,'dMML')$estimate
} 
islands=rbind(islands_OR,islands_OR_het,islands_OR_nonhet)
ggplot(islands,aes(x=subject,y=OR,fill=CpG_type)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+ggtitle('Subject OR for CpG islands')