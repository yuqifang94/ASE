source('mainFunctions_sub.R')
informME_dir='../downstream/data/JSD_heart_mouse_informME/'
JSD_in=GRanges()
for(JSD_fn in dir(informME_dir,pattern='JSD')){
    JSD_in_sample=import.bedGraph(paste0(informME_dir,JSD_fn))
    sample_in=paste(unlist(strsplit(gsub('JSD-|heart_|_all|.bed|mm10_','',JSD_fn),'-VS-')),collapse='-')
    sample_in=paste0('heart-',gsub('day','E',gsub('_','-',gsub('_5','.5',sample_in))),'-all')
    JSD_in_sample$Sample=sample_in
    JSD_in=c(JSD_in,JSD_in_sample)
}
JSD_in_dt=convert_GR(JSD_in,dir='DT')
JSD_in_dt=dcast.data.table(data=JSD_in_dt,region~Sample,value.var='score')
saveRDS(JSD_in_dt,'../downstream/output/mouse_analysis/JSD_informME_heart.rds')
JSD_in_dt=readRDS('../downstream/output/mouse_analysis/JSD_informME_heart.rds')
UC=readRDS(UC_in_matrix_ls_file)
UC=UC$heart
olap=findOverlaps(UC,convert_GR(JSD_in_dt$region),minoverlap=100)
UC_dt=convert_GR(UC[queryHits(olap)],dir='DT')
colnames(UC_dt)[22]='region_UC'
UC_dt$region_JSD=JSD_in_dt[subjectHits(olap)]$region
JSD_in_dt=JSD_in_dt[subjectHits(olap)]
colnames(JSD_in_dt)[1]='region_JSD'
JSD_in_dt$region_UC=UC_dt$region_UC
JSD_in_dt=melt.data.table(id.vars=c('region_UC','region_JSD'),data=JSD_in_dt,
    value.name='score',variable.name='Sample')
JSD_in_dt$stat_type='JSD'
colnames(UC_dt)=gsub('-5','\\.5',gsub('\\.','-',colnames(UC_dt)))
UC_dt=melt.data.table(UC_dt,id.vars=c('region_UC','region_JSD'),
    value.name='score',variable.name='Sample')
UC_dt$stat_type='UC'

UC_JSD_comp=rbind(UC_dt,JSD_in_dt)
UC_JSD_comp_dc=dcast.data.table(UC_JSD_comp,region_UC+region_JSD+Sample~stat_type,value.var='score',
fun.aggregate=mean)
UC_JSD_comp_dc=UC_JSD_comp_dc[!is.na(JSD)]
saveRDS(UC_JSD_comp_dc,'../downstream/output/mouse_analysis/UC_CPEL_JSD_informME_heart.rds')
UC_JSD_comp_dc[UC>0.1,list(cor=cor.test(UC,JSD,method='pearson')$estimate),by=Sample]
#                   Sample          cor
# 1: heart-E12.5-E13.5-all -0.055417879
# 2: heart-E13.5-E14.5-all -0.066865054
# 3: heart-E11.5-E12.5-all -0.052692067
# 4: heart-E14.5-E15.5-all -0.051065299
# 5: heart-E15.5-E16.5-all -0.044829555
# 6: heart-E10.5-E11.5-all -0.004559908

UC_JSD_comp_dc[,list(cor=cor.test(UC,JSD,method='pearson')$estimate),by=Sample]

#                   Sample         cor
# 1: heart-E10.5-E11.5-all -0.03964284
# 2: heart-E11.5-E12.5-all -0.05235849
# 3: heart-E12.5-E13.5-all -0.04315038
# 4: heart-E13.5-E14.5-all -0.03953295
# 5: heart-E14.5-E15.5-all -0.05063979
# 6: heart-E15.5-E16.5-all -0.05054658

NME_in=data.table()
for(NME_fn in dir(informME_dir,pattern='^NME')){
    NME_in_sample=import.bedGraph(paste0(informME_dir,NME_fn))
    NME_in_sample=convert_GR(NME_in_sample,dir='DT')
    sample_in=gsub('NME-|heart_|_all|.bed|mm10_','',NME_fn)
    sample_in=paste0('heart-',gsub('day','E',gsub('_','-',gsub('_5','.5',sample_in))),'-all')
    NME_in_sample$Sample=sample_in
    NME_in=rbind(NME_in,NME_in_sample)
}
NME_in=dcast.data.table(data=NME_in,region~Sample,value.var='score')
saveRDS(NME_in,'../downstream/output/mouse_analysis/NME_informME_heart.rds')
NME_in=readRDS('../downstream/output/mouse_analysis/NME_informME_heart.rds')
NME_CPEL=readRDS(NME_matrix_file)
NME_CPEL=convert_GR(NME_CPEL,dir='DT')
colnames(NME_CPEL)=gsub('-5','\\.5',gsub('\\.','-',colnames(NME_CPEL)))
NME_CPEL=NME_CPEL[,c(48,which(grepl('heart',colnames(NME_CPEL)))),with=F]
olap=findOverlaps(convert_GR(NME_in$region,dir='GR'),convert_GR(NME_CPEL$region,dir='GR'),minoverlap=100)
NME_in=NME_in[queryHits(olap)]
NME_CPEL=NME_CPEL[subjectHits(olap)]
colnames(NME_in)[1]='region_informME'
NME_in$region_CPEL=NME_CPEL$region
colnames(NME_CPEL)[1]='region_CPEL'
NME_CPEL$region_informME=NME_in$region_informME
NME_CPEL=melt.data.table(NME_CPEL,id.vars=c('region_informME','region_CPEL'),value.name='score',variable.name='Sample')
NME_in=melt.data.table(NME_in,id.vars=c('region_informME','region_CPEL'),value.name='score',variable.name='Sample')
NME_in$stat_type='informME'
NME_CPEL$stat_type='CPEL'
NME_comp=rbind(NME_in,NME_CPEL)
NME_comp$Sample=as.character(NME_comp$Sample)
NME_comp_dc=dcast.data.table(NME_comp,region_informME+region_CPEL+Sample~stat_type,value.var='score',
fun.aggregate=mean)
saveRDS(NME_comp_dc,'../downstream/output/mouse_analysis/NME_CPEL_informME_heart.rds')
NME_comp_dc[,list(cor=cor.test(CPEL,informME,method='spearman')$estimate),by=Sample]
#            Sample       cor
#1: heart-E10.5-all 0.6322133
#2: heart-E11.5-all 0.6010391
#3: heart-E12.5-all 0.6097146
#4: heart-E13.5-all 0.5946103
#5: heart-E14.5-all 0.5881557
#6: heart-E15.5-all 0.6092266
#7: heart-E16.5-all 0.6005338

NME_comp_dc[,list(cor=cor.test(CPEL,informME,method='pearson')$estimate),by=Sample]
#            Sample       cor
#1: heart-E10.5-all 0.6744848
#2: heart-E11.5-all 0.6680603
#3: heart-E12.5-all 0.6760570
#4: heart-E13.5-all 0.6667906
#5: heart-E14.5-all 0.6440991
#6: heart-E15.5-all 0.6628684
#7: heart-E16.5-all 0.6543901