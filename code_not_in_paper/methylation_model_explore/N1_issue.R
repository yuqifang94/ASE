source('mainFunctions_sub.R')
theme_glob=theme_classic()+
        theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14),
                                 legend.position="bottom")
#CPEL vs informME in N  CpG
informME_QC_fn='../downstream/output/mouse_analysis/model_QC/JSD_informME_heart.rds'
JSD_in_dt=readRDS(informME_QC_fn)
JSD_in_dt$N=JSD_in_dt$N_informME
JSD_in_dt$N_informME=NULL
CPEL_JSD=readRDS(JSD_in_fn)
CPEL_JSD=CPEL_JSD[,grepl("heart|N|region",colnames(CPEL_JSD)),with=F]
CPEL_JSD=CPEL_JSD[,!grepl("NT-",colnames(CPEL_JSD)),with=F]
CPEL_UC=readRDS(UC_in_QC_fn)
CPEL_UC=CPEL_UC[,grepl("heart|region",colnames(CPEL_UC)),with=F]
mm10_CpG=getCpgSitesmm10()
CPEL_UC$N=countOverlaps(convert_GR(CPEL_UC$region,direction="GR"),mm10_CpG,minoverlap=2)


plot_boxPlot_N_CpG<-function(datIn,quantile=F,statName){
    datIn=melt.data.table(datIn,id.vars=c("N","region"),variable.name="Sample",value.name ="statIn")
    datIn=datIn[N<=17&!is.na(statIn)]
    datIn$N_f = factor(datIn$N,levels=sort(unique(datIn$N)))
    datIn[as.numeric(N_f)>=10]$N_f=">=10"
    if(quantile){
        quant=ecdf(sample(datIn$statIn,1000000))
        datIn$statIn=quant(datIn$statIn)

    }
 
    ggplot(datIn,aes(x=N_f,y=statIn))+geom_boxplot(outlier.shape=NA)+xlab("N")+ylab(statName)+theme_glob+ylim(c(0,quantile(datIn$statIn,prob=0.99,na.rm=T)))+
        ggtitle(paste0("Correlation\nWith N=1 included:",round(cor.test(datIn$N,datIn$statIn)$estimate,digits=2),"\n","Without N=1 included:",
                    round(cor.test(datIn[N>1]$N,datIn[N>1]$statIn)$estimate,digits=2)))
}


JSD_informME=plot_boxPlot_N_CpG(JSD_in_dt[N>0],statName="JSD informME")
JSD_CPEL=plot_boxPlot_N_CpG(CPEL_JSD[N>0],statName="JSD CPEL")
UC_CPEL=plot_boxPlot_N_CpG(CPEL_UC[N>0],statName="UC CPEL")
pdf("../downstream/output/mouse_analysis/model_QC/JSD_UC_informME_N.pdf",height=20,width=7)
    grid.arrange(UC_CPEL,JSD_CPEL,JSD_informME,nrow=3)
dev.off()

JSD_informME_q=plot_boxPlot_N_CpG(JSD_in_dt[N>0],quantile=T,statName="JSD informME")
JSD_CPEL_q=plot_boxPlot_N_CpG(CPEL_JSD[N>0],quantile=T,statName="JSD CPEL")
UC_CPEL_q=plot_boxPlot_N_CpG(CPEL_UC[N>0],quantile=T,statName="UC CPEL")
pdf("../downstream/output/mouse_analysis/model_QC/JSD_UC_informME_N_quant.pdf",height=20,width=7)
    grid.arrange(UC_CPEL_q,JSD_CPEL_q,JSD_informME_q,nrow=3)
dev.off()
