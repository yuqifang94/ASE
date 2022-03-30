#cgmaptools
#/ibox/afeinbe2/yfang/trim_bam/analyzed_reads
#has bug: fix it in R https://github.com/guoweilong/cgmaptools/blob/master/src/mCTanghulu.pl
#Do a agnostic version of https://github.com/markrobinsonuzh/DAMEfinder/blob/master/R/methyl_circle_plot.R
# cgmaptools tanghulu -t CG -r mm10.fa -b mm10_heart_day14_5_all.sub.trimmed.bam -l chr13:43131840-43132089 -o Edn1_mm10_GATA4_E14_5.pdf -s /software/apps/samtools/1.9/gcc/5.5/bin/samtools
#Need to download R file so it can run under R/3.6.3
source('mainFunctions_sub.R')
source('methyl_circle_plot_agnostic.R')
# methyl_circle_plotCpG <- function(cpgsite = cpgsite, bamFile = bamFile, 
#                                   pointSize = 3, refFile = refFile, dame = NULL, order = FALSE, 
#                                   sampleName = NULL, sampleReads = FALSE, numReads = 20)


CG_mm10=getCpgSitesmm10()
genome <- BSgenome.Mmusculus.UCSC.mm10

#Bam file location:
bam_dir="J:/mouse_trimmed/"
plot_two_samples<-function(region_plot,genome,gene_in,stage_1,stage_2,tissue,bam_dir="J:/mouse_trimmed/",
                           folder_out='../downstream/output/mouse_analysis/examples/motif_region_example/',CG_mm10=CG_mm10,
                           pointSize=pointSize,lwd=lwd){
  chr_name=gsub(':.*','',region_plot)
  region_plot=convert_GR(region_plot)
  
 
  dna <- DNAStringSet(genome[[chr_name]], use.names = TRUE)
  names(dna) <- chr_name
  CpG_site=subsetByOverlaps(CG_mm10,region_plot)
  
  pdf(paste0(folder_out,tissue,'_',stage_1,'_',gene_in,'.pdf'),height=0.5,width=2)
  print(methyl_circle_plotCpG(cpgsite=CpG_site,bamFile = paste0(bam_dir,"mm10_",tissue,"_",gsub("E","day",gsub("\\.5","_5",stage_1)),"_all.sub.trimmed.bam"),dame=region_plot,sampleName = paste0(tissue," ",stage_1),
                        refFile=dna,pointSize=pointSize,lwd=lwd))
  dev.off()
  
  pdf(paste0(folder_out,tissue,'_',stage_2,'_',gene_in,'.pdf'),height=0.5,width=2)
  print(methyl_circle_plotCpG(cpgsite=CpG_site,bamFile = paste0(bam_dir,"mm10_",tissue,"_",gsub("E","day",gsub("\\.5","_5",stage_2)),"_all.sub.trimmed.bam"),dame=region_plot,sampleName = paste0(tissue," ",stage_2),
                        refFile=dna,pointSize=pointSize,lwd=lwd))
  dev.off()
  
  
}

# dNME examples -----------------------------------------------------------
#heart
#Region: chr13:43131840-43132089, edn1
plot_two_samples("chr13:43131840-43132089",genome,"Edn1","E14.5","E15.5","heart",CG_mm10=CG_mm10)
#chr12:80208821-80209070: Actn1, 14.5-15.5
plot_two_samples("chr12:80208821-80209070",genome,"Actn1","E14.5","E15.5","heart",CG_mm10=CG_mm10)
#chr11:65303752-65304001, Myocd E14.5-E15.5
plot_two_samples("chr11:65303752-65304001",genome,"Myocd","E14.5","E15.5","heart",CG_mm10=CG_mm10)
#chr13:115005436-115005685, igta1 E14.5-E15.5
plot_two_samples("chr13:115005436-115005685",genome,"Igta1","E14.5","E15.5","heart",CG_mm10=CG_mm10)
#chr4:100226672-100226921, Ror1 ,E13.5-E14.5
plot_two_samples("chr4:100226672-100226921",genome,"Ror1","E13.5","E14.5","heart",CG_mm10=CG_mm10)
#chr5:107794007-107794256,Tgfbr3,E14.5-E15.5
plot_two_samples("chr5:107794007-107794256",genome,"Tgfbr3","E14.5","E15.5","heart",CG_mm10=CG_mm10)
#chr3:86687916-86688168
plot_two_samples("chr3:86687916-86688168",genome,"Mab21l2","E14.5","E15.5","heart",CG_mm10=CG_mm10)
#chr7:81134783-81135032, Alpk3, E14.5, E15.5
plot_two_samples("chr7:81134783-81135032",genome,"Alpk3","E14.5","E15.5","heart",CG_mm10=CG_mm10,pointSize=0.35,lwd=0.075)
#chr2:180322867-180323133, GATA5, E11.5, E12.5
plot_two_samples("chr2:180322867-180323133",genome,"GATA5","E11.5","E12.5","heart",CG_mm10=CG_mm10,pointSize=0.35,lwd=0.075)

#forebrain
#chr14:28438995-28439244 forebrain       2  Wnt5a   dNME_only E10.5-E11.5     0.4116945           MEIS1
# chr3:108433651-108433900 forebrain       6 Celsr2   dNME_only E12.5-E13.5     0.3494379            KLF4

plot_two_samples("chr14:28438995-28439244",genome,"Wnt5a","E10.5","E11.5","forebrain",CG_mm10=CG_mm10,pointSize=0.15,lwd=0.075)
plot_two_samples("chr3:108433651-108433900",genome,"Celsr2","E12.5","E13.5","forebrain",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
#EFP
# chr16:44616342-44616591    EFP       2     Boc   dNME_only E10.5-E11.5     0.4416258            TBX5
# chr7:142516416-142516665    EFP       4    Igf2   dNME_only E11.5-E12.5     0.3919922             REL
# chr13:56142950-56143199    EFP       3   Tgfbi   dMML_only E11.5-E12.5     0.3269635            KLF4
# chr9:60961980-60962229    EFP       8    Tle3   dNME_only E12.5-E13.5     0.4685049             REL
#chr8:87822376-87822625    EFP       2  Zfp423   dNME_only E10.5-E11.5     0.4466273            KLF4
plot_two_samples("chr16:44616342-44616591",genome,"Boc","E10.5","E11.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr7:142516416-142516665",genome,"Igf2","E11.5","E12.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr13:56142950-56143199",genome,"Tgfbi","E11.5","E12.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr9:60961980-60962229",genome,"Tle3","E12.5","E13.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)#GD
plot_two_samples("chr8:87822376-87822625 ",genome,"Zfp423","E10.5","E11.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr11:95213671-95213920",genome,"Dlx3","E14.5","E15.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr8:83601830-83602079",genome,"Il15","E10.5","E11.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr18:69614251-69614500",genome,"Tcf4","E10.5","E11.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr4:140219094-140219343",genome,"Igsf21","E14.5","E15.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr1:180911257-180911506",genome,"Lefty1","E10.5","E11.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr8:25294235-25294484",genome,"Fgfr1","E12.5","E13.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
plot_two_samples("chr11:95213671-95213920",genome,"Nnat","E11.5","E12.5","EFP",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
#Limb

plot_two_samples("chr4:155283595-155283844",genome,"Ski","E14.5","E15.5","limb",CG_mm10=CG_mm10,pointSize=0.25,lwd=0.075)
#Forebrain old
#Region: chr1:75693286-75693535, Epha4
plot_two_samples("chr1:75693286-75693535",genome,"Epha4","E15.5","E16.5","forebrain")
#chr4:23636023-23636272,Pou3f2
plot_two_samples("chr4:23636023-23636272",genome,"Pou3f2","E11.5","E12.5","forebrain")
#chr16:22428600-22428849, Etv5, E13-E14
plot_two_samples("chr16:22428600-22428849",genome,"Etv5","E13.5","E14.5","forebrain")
#chr4:148854190-148854439,Mfn2, E13-E14
plot_two_samples("chr4:148854190-148854439",genome,"Mfn2","E13.5","E14.5","forebrain")
#chr13:25750988-25751237,SOx4,E14-E15
plot_two_samples("chr13:25750988-25751237",genome,"Sox4","E14.5","E15.5","forebrain")
#chr4:82464491-82464740,Nfib,E12-E13
plot_two_samples("chr4:82464491-82464740",genome,"Nfib","E12.5","E13.5","forebrain")

plot_two_samples("chr14:63207424-63207673",genome,"Gata4","E11.5","E12.5","forebrain",CG_mm10=CG_mm10)
plot_two_samples("chr14:63212488-63212737",genome,"Gata4","E15.5","E16.5","forebrain",CG_mm10=CG_mm10)


# chr1:75443504-75443755,Asic4, E13.5-E14.5
plot_two_samples("chr1:75443504-75443755",genome,"Asic4","E13.5","E14.5","forebrain",CG_mm10=CG_mm10)
#Limb
#
#Liver
#chr18:38278051-38278300
plot_two_samples("chr18:38278051-38278300",genome,"Fgf1","E11.5","E12.5","liver")

#Human examples
plot_two_samples_human<-function(region_plot,genome,gene_in,tissue_in,sample_in,bam_dir="/dcs04/feinberg/data/personal/yfang/allele_specific_roadmap_CEPL/work_archive/CpelAsm/data/",
                                 folder_out='../downstream/output/mouse_analysis/examples/motif_region_example/',CG_hg19=CG_hg19,pointSize=pointSize,lwd=lwd){
  bam_dir=paste0(bam_dir,tissue_in,'/bam/')
  chr_name=gsub(':.*','',region_plot)
  region_plot=convert_GR(region_plot)
  
  
  dna <- DNAStringSet(genome[[chr_name]], use.names = TRUE)
  names(dna) <- chr_name
  CpG_site=subsetByOverlaps(CG_hg19,region_plot)
  
  pdf(paste0(folder_out,tissue_in,'_',sample_in,'_genome1_',gene_in,'.pdf'),height=0.5,width=2, useDingbats=FALSE)
  print(methyl_circle_plotCpG(cpgsite=CpG_site,bamFile = paste0(bam_dir,sample_in,'.sort.genome1.bam'),dame=region_plot,sampleName = paste0(sample_in,"_genome1"),
                              refFile=dna,pointSize=pointSize,lwd=lwd,softclip=0))
  dev.off()
  
  pdf(paste0(folder_out,tissue_in,'_',sample_in,'_genome2_',gene_in,'.pdf'),height=0.5,width=2, useDingbats=FALSE)
  print(methyl_circle_plotCpG(cpgsite=CpG_site,bamFile = paste0(bam_dir,sample_in,'.sort.genome2.bam'),dame=region_plot,sampleName = paste0(sample_in,"_genome2"),
                              refFile=dna,pointSize=pointSize,lwd=lwd,softclip=0))
  dev.off()
  
  
}
bam_dir='/dcs04/feinberg/data/personal/yfang/allele_specific_roadmap_CEPL/work_archive/CpelAsm/data/'
CG_hg19=getCpgSiteshg19()
genome <- BSgenome.Hsapiens.UCSC.hg19
seqlevels(genome)=gsub('chr','',seqlevels(genome))
seqlevels(CG_hg19)=gsub('chr','',seqlevels(CG_hg19))
plot_two_samples_human("6:144329572-144329772",genome,'PLAG1','150','150_Brain_Hippocampus_middle_paired_phased',CG_hg19=CG_hg19,pointSize=0.35,lwd=0.075)
plot_two_samples_human("1:18993859-18994059",genome,'PAX7','HUES64','HUES64_endoerm_27_paired_phased',CG_hg19=CG_hg19,pointSize=3)
plot_two_samples_human("9:4023395-4023595",genome,'GLIS3','HUES64','HUES64_stem_27_undifferentiated_paired_phased',CG_hg19=CG_hg19,pointSize=3)
plot_two_samples_human("11:2720817-2721033",genome,'KCNQ1','HUES64','HUES64_stem_27_undifferentiated_paired_phased',CG_hg19=CG_hg19,pointSize=0.25,lwd=0.075)
##chr4:7,308,694-7,309,119
plot_two_samples_human("4:7308694-7309119",genome,'SORCS2','HUES64','HUES64_stem_27_undifferentiated_paired_phased',CG_hg19=CG_hg19,pointSize=0.35,lwd=0.075)
plot_two_samples_human("14:104566670-104566870",genome,'ASPG','HUES64','HUES64_stem_27_undifferentiated_paired_phased',CG_hg19=CG_hg19,pointSize=0.5,lwd=0.075)
plot_two_samples_human("1:175568947-175569147",genome,'TNR','HUES64','HUES64_stem_27_undifferentiated_paired_phased',CG_hg19=CG_hg19,pointSize=0.35,lwd=0.075)
plot_two_samples_human("19:54668387-54668587",genome,'TMC4','HUES64','HUES64_stem_27_undifferentiated_paired_phased',CG_hg19=CG_hg19,pointSize=0.35,lwd=0.075)

#Select human sample
bam_dir='/dcs04/feinberg/data/personal/yfang/allele_specific_roadmap_CEPL/work_archive/CpelAsm/data/'
CG_hg19=getCpgSiteshg19()
genome <- BSgenome.Hsapiens.UCSC.hg19
seqlevels(genome)=gsub('chr','',seqlevels(genome))
seqlevels(CG_hg19)=gsub('chr','',seqlevels(CG_hg19))
outDir='../downstream/output/human_analysis/example/'
bam_dir='/dcs04/feinberg/data/personal/yfang/allele_specific_roadmap_CEPL/work_archive/CpelAsm/data/'
GR_merge=readRDS(GR_merge_file)
GR_merge=convert_GR(GR_merge,direction="DT")
GR_merge_sub= GR_merge[dMML<0.1&dNME_pval<0.1&N_hg19>=5&(!is.na(genes_promoter)|!is.na(genes_body))][order(dNME,decreasing=T),
    list(dNME,dMML,N,region,Sample,genes_promoter,genes_body,MML1,MML2,NME1,NME2)]
 source('mouse_example_selection/methyl_circle_plot_agnostic.R')
 GR_merge_sub[grepl("ESC",Sample)]$Sample=gsub("ESC","ESC_paired",GR_merge_sub[grepl("ESC",Sample)]$Sample)
for (i in 1:nrow(GR_merge_sub) ){
  region_select=as.data.frame(GR_merge_sub)[i,]
  gene = na.omit(c(unlist(region_select$genes_promoter),unlist(region_select$genes_body)))[1]
  tissue=gsub('.*- ','',region_select$Sample)
    sample=paste0(tissue,'_',gsub(' -.*','',region_select$Sample),'_phased')
  cat("Processing:",sample,"\n")
  if(tissue != "H1"&!sample %in% c('STL003_Spleen_single_phased')){

  plot_two_samples_human(gsub("chr","",region_select$region),genome,gene,tissue,sample,CG_hg19=CG_hg19,pointSize=0.35,lwd=0.075,folder_out=outDir,
  bam_dir=bam_dir)
  }
}
