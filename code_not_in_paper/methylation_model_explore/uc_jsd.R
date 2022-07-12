source('mainFunctions_sub.R')
library(Gmisc)
# UC_in_dir='../downstream/data/compliment_UC_non_MDS_mouse/'
# UC_in_dir_analyzed='../downstream/data/DNase_control_PRC_non_MDS_mouse/'
# gff_in_DNase=import.gff3(mouse_DNase_control_gff_file)
# gff_in_DNase=paste0(seqnames(gff_in_DNase),':',start(gff_in_DNase),'-',end(gff_in_DNase))
# UC_in_analyzed_MDS=data.table(region=gff_in_DNase)
# UC_in_analyzed_MDS_UC=fastDoCall('cbind',
#                                  mclapply(dir(UC_in_dir,pattern = '.*uc.bedGraph'),function(x){
#                                    read.agnostic.mouse.uc(paste(UC_in_dir_analyzed,x,sep=''),matrix=T,fileter_N=1,gff_in=gff_in_DNase)},mc.cores=10))
# UC_in_analyzed_MDS=cbind(UC_in_analyzed_MDS,UC_in_analyzed_MDS_UC)
# saveRDS(UC_in_analyzed_MDS,UC_in_QC_fn)

cat("Loading dat\n")
tt1=proc.time()[[3]]
reformat_dt<-function(datIn){
  rn=datIn$region
  datIn=datIn[,-which(colnames(jsd)=="region"|colnames(jsd)=="N"),with=F]
  datIn=apply(as.matrix(datIn),2,as.numeric)
  rownames(datIn)=rn
  return(datIn)

}
jsd=readRDS(JSD_in_fn)
jsd=reformat_dt(jsd)
uc=readRDS(UC_in_QC_fn)
colnames(uc)=gsub(paste(paste0('-',unique(gsub('-.*','',colnames(jsd)))),collapse='|'),'',gsub('_','-',colnames(uc)))
uc=reformat_dt(uc)

region_inter=intersect(rownames(jsd),rownames(uc))
sample_inter=intersect(colnames(jsd),colnames(uc))

jsd=jsd[region_inter,sample_inter]
uc=uc[region_inter,sample_inter]
tt2=proc.time()[[3]]
cat("Finishing loading and reshaping data in: ",tt2-tt1,"\n")
JSD_UC=data.table(JSD=as.vector(jsd),UC=as.vector(uc))[!is.na(JSD)&!is.na(UC)]
saveRDS(JSD_UC,'../downstream/output/mouse_analysis/QC/JSD_UC_CPEL.rds')
cat("Start calculating correlation\n")
cut_num=2000
percent_cut=quantile(1:cut_num,prob=seq(0.1,1,0.1))
split_data=cut(1:length(region_inter),cut_num,label=FALSE)
tt2=proc.time()[[3]]
cor_output=fastDoCall('rbind',lapply(1:cut_num,function(x){
  if(x %in% percent_cut){
    cat("Finished:",x/cut_num*100,'%','in ',proc.time()[[3]]-tt2,'\n')
    }
  return(data.table(region=region_inter[split_data==x],
                    cor=diag(cor(t(jsd[region_inter[split_data==x],]),t(uc[region_inter[split_data==x],]),method="spearman"))))
}))

saveRDS(cor_output,'../downstream/output/mouse_analysis/QC/cor_JSD_UC_CPEL.rds')
tt3=proc.time()[[3]]
cat("Finish calculating correlation in: ",tt3-tt2,'\n')
theme_glob=theme_classic()+
        theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.text.x=element_text(size=14),
                                 axis.text.y=element_text(size=14),
                                 legend.position="bottom")
                                 
#Plotting
cor_output=readRDS('../downstream/output/mouse_analysis/QC/cor_JSD_UC_CPEL.rds')
cor_output=cor_output[!is.na(cor)]
pdf('../downstream/output/mouse_analysis/QC/cor_JSD_UC_CPEL.pdf')
cor_all=ggplot(cor_output,aes(x=cor))+geom_histogram(fill='white',color='blue')+xlab("Correlation")+theme_glob+ggtitle("All correlation")
cor_high=ggplot(cor_output[cor>=0.95],aes(x=cor))+geom_histogram(fill='white',color='blue')+xlab("Correlation")+theme_glob+ggtitle("Correlation >=0.95")
grid.arrange(cor_all,cor_high,nrow=2)
dev.off()
print(mean(cor_output$cor,na.rm=T))
#This takes extremely long time
JSD_UC=readRDS('../downstream/output/mouse_analysis/QC/JSD_UC_CPEL.rds')
JSD_UC=fastDoCall('rbind',JSD_UC
pdf('../downstream/output/mouse_analysis/QC/JSD_UC_CPEL.pdf')
      ggplot(JSD_UC)[sample(1:nrow(JSD_UC),100000000)],aes(x=UC,y=JSD))+geom_bin2d(bins=100)+
        geom_smooth(colour="white",method='lm')+
        #geom_abline(intercept=0,slope=1,colour="red")+
        scale_fill_continuous(type = "viridis")+theme_glob
      #Need large mem here
      # lapply(names(JSD_UC),function(x){
      # print(ggplot(JSD_UC[[x]],aes(x=UC,y=JSD))+geom_bin2d(bins=100)+
      #   geom_smooth(colour="white",method='lm')+
      #   #geom_abline(intercept=0,slope=1,colour="red")+
      #   scale_fill_continuous(type = "viridis")+theme_glob+ggtitle(x)
      # )

      # })

dev.off()
#Check /users/yfang/yfang_dcs04/allele_agnostic_mouse_all/CPEL/CPEL_agnostic/cpelasm/cpelasm/MARCC_archive_20210830/jsd_output_CPEL
