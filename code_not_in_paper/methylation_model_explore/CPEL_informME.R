source('mainFunctions_sub.R')
informME_dir='../downstream/data/JSD_heart_mouse_informME/'
informME_QC_fn='../downstream/output/mouse_analysis/model_QC/JSD_informME_heart.rds'
informME_in=GRanges()
for(informME_fn in dir(informME_dir,pattern='JSD')){
    informME_in_sample=import.bedGraph(paste0(informME_dir,informME_fn))
    sample_in=paste(unlist(strsplit(gsub('JSD-|heart_|_all|.bed|mm10_','',informME_fn),'-VS-')),collapse='-')
    sample_in=paste0('heart-',gsub('day','E',gsub('_','-',gsub('_5','.5',sample_in))),'-all')
    informME_in_sample$Sample=sample_in
    informME_in=c(informME_in,informME_in_sample)
}
informME_in_dt=convert_GR(informME_in,dir='DT')
informME_in_dt=dcast.data.table(data=informME_in_dt,region~Sample,value.var='score')
mm10_CpG=getCpgSitesmm10()
informME_in_dt$N_informME=countOverlaps(convert_GR(informME_in_dt$region,direction="GR"),mm10_CpG)
saveRDS(informME_in_dt,informME_QC_fn)
informME_in_dt=readRDS(informME_QC_fn)
CPEL=readRDS(JSD_in_fn)
CPEL=CPEL[,grepl('heart|region|N',colnames(CPEL)),with=F]
colnames(CPEL)[which(colnames(CPEL)=="N")]="N_CPEL"
olap=findOverlaps(convert_GR(CPEL$region,direction="GR"),convert_GR(informME_in_dt$region,direction="GR"),minoverlap=100)
CPEL=CPEL[queryHits(olap)]
informME_in_dt=informME_in_dt[subjectHits(olap)]
CPEL$region_informME=informME_in_dt$region
CPEL$N_informME=informME_in_dt$N_informME
colnames(informME_in_dt)[which(colnames(informME_in_dt)=="region")]='region_informME'
informME_in_dt$region_CPEL=CPEL$region
informME_in_dt$N_CPEL=CPEL$N_CPEL
informME_in_dt=melt.data.table(id.vars=c('region_CPEL','region_informME',"N_CPEL","N_informME"),data=informME_in_dt,
    value.name='score',variable.name='Sample')
informME_in_dt$stat_type='informME'
# colnames(CPEL)=gsub('-5','\\.5',gsub('\\.','-',colnames(CPEL)))
colnames(CPEL)[which(colnames(CPEL)=="region")]='region_CPEL'
CPEL=melt.data.table(CPEL,id.vars=c('region_CPEL','region_informME',"N_CPEL","N_informME"),
    value.name='score',variable.name='Sample')
CPEL$stat_type='CPEL'
# CPEL$Sample=gsub('-heart|_all','',CPEL$Sample)
CPEL=CPEL[Sample %in% informME_in_dt$Sample]
CPEL_informME_comp=rbind(CPEL,informME_in_dt)
CPEL_informME_comp_dc=dcast.data.table(CPEL_informME_comp,region_CPEL+region_informME+N_CPEL+N_informME+Sample~stat_type,value.var='score',
fun.aggregate=mean)
CPEL_informME_comp_dc=CPEL_informME_comp_dc[!is.na(informME)&!is.na(CPEL)]
saveRDS(CPEL_informME_comp_dc,'../downstream/output/mouse_analysis/model_QC/CPEL_informME_JSD_heart_N.rds')


CPEL_informME_comp_dc[N_informME>1&N_CPEL>1][,list(cor=cor.test(CPEL[CPEL>quantile(CPEL,prob=0.99)&informME>quantile(informME,prob=0.99)],
                                                                    informME[CPEL>quantile(CPEL,prob=0.99)&informME>quantile(informME,prob=0.99)],method='pearson')$estimate),by=Sample]

#                   Sample        cor
# 1: heart-E10.5-E11.5-all 0.15584101
# 2: heart-E11.5-E12.5-all 0.15074483
# 3: heart-E12.5-E13.5-all 0.15884542
# 4: heart-E13.5-E14.5-all 0.06526965
# 5: heart-E14.5-E15.5-all 0.24524360
# 6: heart-E15.5-E16.5-all 0.13212972         


#Mean: 0.15
#Preprocessing NME and MML
informME_CPEL_NME_MML<-function(statIn,CPEL_file){
    tt1=proc.time()[[3]]
    cat("Start loading data\n")
    informME_dir='../downstream/data/JSD_heart_mouse_informME/'
    dat_in=data.table()
    for(dat_fn in dir(informME_dir,pattern=paste0('^',statIn))){
        dat_in_sample=import.bedGraph(paste0(informME_dir,dat_fn))
        dat_in_sample=convert_GR(dat_in_sample,dir='DT')
        sample_in=gsub('dat-|heart_|_all|.bed|mm10_','',dat_fn)
        sample_in=paste0('heart-',gsub('day','E',gsub('_','-',gsub('_5','.5',sample_in))),'-all')
        sample_in=gsub('-NME|-MML','',sample_in)
        dat_in_sample$Sample=sample_in
        dat_in=rbind(dat_in,dat_in_sample)
    }
    dat_in=dcast.data.table(data=dat_in,region~Sample,value.var='score')
   
    saveRDS(dat_in,paste0('../downstream/output/mouse_analysis/model_QC/',statIn,'_informME_heart.rds'))
    tt2=proc.time()[[3]]
    cat("Finish loading data in:", tt2-tt1,"\n")
    cat("Start reshaping samples\n")
    dat_CPEL=readRDS(CPEL_file)
    dat_CPEL=convert_GR(dat_CPEL,dir='DT')
    colnames(dat_CPEL)=gsub('-5','\\.5',gsub('\\.','-',colnames(dat_CPEL)))
    dat_CPEL=dat_CPEL[,c(which(grepl('region',colnames(dat_CPEL))),which(grepl('heart',colnames(dat_CPEL)))),with=F]
    dat_CPEL=dat_CPEL[rowSums(is.na(dat_CPEL))==0]
    olap=findOverlaps(convert_GR(dat_in$region,dir='GR'),convert_GR(dat_CPEL$region,dir='GR'),minoverlap=100)
    dat_in=dat_in[queryHits(olap)]
    dat_CPEL=dat_CPEL[subjectHits(olap)]
    colnames(dat_in)[1]='region_informME'
    dat_in$region_CPEL=dat_CPEL$region
    colnames(dat_CPEL)[1]='region_CPEL'
    dat_CPEL$region_informME=dat_in$region_informME
    dat_CPEL=melt.data.table(dat_CPEL,id.vars=c('region_informME','region_CPEL'),value.name='score',variable.name='Sample')
    dat_in=melt.data.table(dat_in,id.vars=c('region_informME','region_CPEL'),value.name='score',variable.name='Sample')
    dat_in$stat_type='informME'
    dat_CPEL$stat_type='CPEL'
    dat_comp=rbind(dat_in,dat_CPEL)
    dat_comp$Sample=as.character(dat_comp$Sample)
    dat_comp_dc=dcast.data.table(dat_comp,region_informME+region_CPEL+Sample~stat_type,value.var='score',
    fun.aggregate=mean)
    saveRDS(dat_comp_dc,paste0('../downstream/output/mouse_analysis/model_QC/',statIn,'_CPEL_informME_heart.rds'))
    tt3=proc.time()[[3]]
    cat("Finish reshaping samples in: ",tt3-tt2,"\n")
    return(list(spearman_cor=dat_comp_dc[,list(cor=cor.test(CPEL,informME,method='spearman')$estimate),by=Sample],
                pearson_cor=dat_comp_dc[,list(cor=cor.test(CPEL,informME,method='pearson')$estimate),by=Sample]))

}
NME_cor=informME_CPEL_NME_MML("NME",NME_matrix_file)
saveRDS(NME_cor,'../downstream/output/mouse_analysis/model_QC/NME_cor_informME_CPEL.rds')
MML_cor=informME_CPEL_NME_MML("MML",MML_matrix_file)
saveRDS(MML_cor,'../downstream/output/mouse_analysis/model_QC/MML_cor_informME_CPEL.rds')
theme_glob=theme_classic()+
        theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14),
                                 legend.position="bottom")
                                 
informME_CPEL_NME_MML_plot<-function(stat_comp_dc){
    ggplot(stat_comp_dc,aes(x=CPEL,y=informME))+geom_bin2d(bins=100)+
        geom_smooth(colour="white")+
        geom_abline(intercept=0,slope=1,colour="red")+
        scale_fill_continuous(type = "viridis")

}

NME_comp_dc=readRDS('../downstream/output/mouse_analysis/model_QC/NME_CPEL_informME_heart.rds')
pdf('../downstream/output/mouse_analysis/model_QC/NME_CPEL_informME_heart_cor.pdf')
    print(informME_CPEL_NME_MML_plot(NME_comp_dc)+ggtitle("NME")+theme_glob)
dev.off()

MML_comp_dc=readRDS('../downstream/output/mouse_analysis/model_QC/MML_CPEL_informME_heart.rds')
pdf('../downstream/output/mouse_analysis/model_QC/MML_CPEL_informME_heart_cor.pdf')
    print(informME_CPEL_NME_MML_plot(MML_comp_dc)+ggtitle("MML")+theme_glob)
dev.off()

CPEL_informME_comp_dc=readRDS('../downstream/output/mouse_analysis/model_QC/CPEL_informME_JSD_heart_N.rds')


pdf('../downstream/output/mouse_analysis/model_QC/JSD_CPEL_informME_heart_cor.pdf')
    print(informME_CPEL_NME_MML_plot(CPEL_informME_comp_dc[N_informME>1&N_CPEL>1])+xlab("CPEL")+ylab("informME")+theme_glob)
dev.off()

pdf('../downstream/output/mouse_analysis/model_QC/JSD_CPEL_informME_heart_cor_top.pdf')
    print(informME_CPEL_NME_MML_plot(CPEL_informME_comp_dc[N_informME>1&N_CPEL>1&CPEL>quantile(CPEL,prob=0.99)&informME>quantile(informME,prob=0.99)])+xlab("CPEL")+ylab("informME")+theme_glob)
dev.off()

