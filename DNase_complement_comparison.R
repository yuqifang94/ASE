source('mainFunctions_sub.R')
DNAase=readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp.rds')
control=c(readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp_control.rds'),
          readRDS('../downstream/input/mouse_analysis/mm10_PRC.rds'),
          readRDS('../downstream/output/mouse_analysis/mm10_allele_agnostic_analysis_compliment.rds'))
NME_old=readRDS('../downstream/input/mouse_analysis/NME_agnostic_mouse_all_merged.rds')
NME_old=NME_old[grepl('heart',NME_old$Sample)]
NME_new=readRDS('../downstream/input/mouse_analysis/NME_agnostic_mouse_all_merged_complement_heart.rds')
NME_heart=c(NME_old[NME_old$N>=2],NME_new[NME_new$N>=2])
NME_heart$stat_type=NULL
NME_heart$score=NULL
rm(NME_old)
rm(NME_new)
DNase_olap=findOverlaps(NME_heart,DNAase,type='equal')
NME_heart_DNase=NME_heart[queryHits(DNase_olap)]
NME_heart_non_DNase=NME_heart[-queryHits(DNase_olap)]
NME_heart_dt=rbind(data.table(NME=NME_heart_DNase$NME,Sample=NME_heart_DNase$Sample,region_type="DNase"),
                   data.table(NME=NME_heart_non_DNase$NME,Sample=NME_heart_non_DNase$Sample,region_type="non-DNase"))
pdf('../downstream/output/mouse_analysis/DNase_non_DNase/all_samples_DNase_non_DNase.pdf')
ggplot(NME_heart_dt,aes(x=NME,color=region_type))+geom_density()
dev.off()
pdf('../downstream/output/mouse_analysis/DNase_non_DNase/each_samples_DNase_non_DNase.pdf')
for(sp in unique(NME_heart_dt$Sample)){
  print(ggplot(NME_heart_dt[Sample==sp],aes(x=NME,color=region_type))+geom_density()+ggtitle(sp))
}
dev.off()
#MML
MML_old=readRDS('../downstream/input/mouse_analysis/MML_agnostic_mouse_all_merged.rds')
MML_old=MML_old[grepl('heart',MML_old$Sample)]
MML_new=readRDS('../downstream/input/mouse_analysis/MML_agnostic_mouse_all_merged_complement_heart.rds')
MML_heart=c(MML_old[MML_old$N>=2],MML_new[MML_new$N>=2])
MML_heart$stat_type=NULL
MML_heart$score=NULL
rm(MML_old)
rm(MML_new)
DNase_olap=findOverlaps(MML_heart,DNAase,type='equal')
MML_heart_DNase=MML_heart[queryHits(DNase_olap)]
MML_heart_non_DNase=MML_heart[-queryHits(DNase_olap)]
MML_heart_dt=rbind(data.table(MML=MML_heart_DNase$MML,Sample=MML_heart_DNase$Sample,region_type="DNase"),
                   data.table(MML=MML_heart_non_DNase$MML,Sample=MML_heart_non_DNase$Sample,region_type="non-DNase"))
pdf('../downstream/output/mouse_analysis/DNase_non_DNase/all_samples_DNase_non_DNase.pdf')
ggplot(MML_heart_dt,aes(x=MML,color=region_type))+geom_density()
dev.off()
pdf('../downstream/output/mouse_analysis/DNase_non_DNase/each_samples_DNase_non_DNase.pdf')
for(sp in unique(MML_heart_dt$Sample)){
  print(ggplot(MML_heart_dt[Sample==sp],aes(x=MML,color=region_type))+geom_density()+ggtitle(sp))
}
dev.off()

#UC
UC_old=readRDS('../downstream/input/mouse_analysis/UC_agnostic_mouse_matrix_dedup_N2_all_merged_ls_fix.rds')
UC_old=UC_old$heart
UC_new=readRDS('../downstream/input/mouse_analysis/UC_agnostic_mouse_matrix_dedup_N2_heart_complement.rds')
colnames(mcols(UC_new$heart))=sub('.5-heart-E','.5-E',colnames(mcols(UC_new$heart)))
UC_heart=c(UC_new$heart,UC_old)
DNase_olap=findOverlaps(UC_heart,DNAase,type='equal')
UC_heart_DNase=UC_heart[queryHits(DNase_olap)]
UC_heart_non_DNase=UC_heart[-queryHits(DNase_olap)]
UC_heart_DNase_dt=melt.data.table(convert_GR(UC_heart_DNase,direction="DT"),value.name="UC",variable.name = "Sample",id.vars = 'region')
UC_heart_DNase_dt$region_type="DNase"
UC_heart_non_DNase_dt=melt.data.table(convert_GR(UC_heart_non_DNase,direction="DT"),value.name="UC",variable.name = "Sample",id.vars = 'region')
UC_heart_non_DNase_dt$region_type="non-DNase"
UC_heart_dt=rbind(UC_heart_DNase_dt[!grepl("P0",Sample)],UC_heart_non_DNase_dt[!grepl("P0",Sample)])[!is.na(UC)]
pdf('../downstream/output/mouse_analysis/DNase_non_DNase/all_samples_DNase_non_DNase_UC_log.pdf')
ggplot(UC_heart_dt,aes(x=log(UC) ,color=region_type))+geom_density()
dev.off()

pdf('../downstream/output/mouse_analysis/DNase_non_DNase/all_samples_DNase_non_DNase_UC_log_max.pdf')
ggplot(UC_heart_dt[,list(UC=max(UC)),by=list(region,region_type)],aes(x=log(UC) ,color=region_type))+geom_density()
dev.off()

#Find number of extra enhancer genes
enhancer=readRDS('../downstream/output/mouse_analysis/enhancers/bin_enhancer.rds')
UC_heart$DNase=FALSE
UC_heart$DNase[queryHits(DNase_olap)]=TRUE
enhancer_olap=findOverlaps(UC_heart,enhancer)
UC_heart$gene=NA
UC_heart[queryHits(enhancer_olap)]$gene=enhancer$`Target Gene`[subjectHits(enhancer_olap)]
UC_heart_01_sel=UC_heart[rowMaxs(as.matrix(mcols(UC_heart)[,grepl("heart",colnames(mcols(UC_heart)))]),na.rm = T)>=0.1]
UC_heart_dt=melt.data.table(convert_GR(UC_heart,direction = "DT"),id.vars = c('region','gene','DNase'),variable.name = 'Sample',value.name = 'UC')
UC_heart_dt_gene=UC_heart_dt[,list(UC=median(UC)),by=list(Sample,gene,DNase)]
UC_heart_dt_gene=UC_heart_dt_gene[!is.na(UC)]
UC_heart_dt=UC_heart_dt[!is.na(UC)]
UC_heart_dt=UC_heart_dt[!grepl("P0",Sample)]
UC_heart_dt=UC_heart_dt[!grepl("chrX|chrY",region)]
out=data.table()
for(g in unique(UC_heart_dt$gene)) {

  if(sum(UC_heart_dt$gene==g,na.rm = T)>6){
   gene_dt=dcast.data.table(UC_heart_dt[gene==g],Sample~DNase,value.var = 'UC',fun.aggregate=median)
   if(ncol(gene_dt)==3){
    colnames(gene_dt)=c("Sample","Non_DNase","DNase")
    out=rbind(out,data.table(cor=cor(gene_dt$Non_DNase,gene_dt$DNase),gene=g))
   }
  }
  }

#5609 all,4398 non-DNase, 5024 total, 4109 DNase gene  
write(unique(UC_heart_01_sel[UC_heart_01_sel$DNase==FALSE]$gene),'../downstream/output/mouse_analysis/DNase_non_DNase/non_DNase_gene_heart_01.txt')
write(unique(UC_heart_01_sel[UC_heart_01_sel$DNase==TRUE]$gene),'../downstream/output/mouse_analysis/DNase_non_DNase/DNase_gene_heart_01.txt')

#Human
GR_merge=readRDS(GR_merge_file)
GR_merge$NME_mean=(GR_merge$NME1+GR_merge$NME2)/2
DNase=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp_control.rds')
DNase_olap=findOverlaps(GR_merge,DNase)
GR_merge$DNase="Non_DNase"
GR_merge$DNase[queryHits(DNase_olap)]="DNase"
pdf("DNase_human_density.pdf")
ggplot(as.data.table(mcols(GR_merge)),aes(color=DNase,x=NME_mean))+geom_density()
dev.off()
