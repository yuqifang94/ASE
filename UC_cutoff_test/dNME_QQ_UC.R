#QQ plot

source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
UC_in_max_loc=readRDS(UC_merge_max_loc_file)
dNME_in=mclapply(UC_in,function(x) x[,grepl('dNME',colnames(x))],mc.cores=7)
dMML_in=mclapply(UC_in,function(x) x[,grepl('dMML',colnames(x))],mc.cores=7)
UC_in_only=mclapply(UC_in,function(x) x[,grepl('UC',colnames(x))],mc.cores=7)
rm(UC_in)
for (ts in names(UC_in_max_loc)){
    max_ts=UC_in_max_loc[[ts]][,c('UC_max_pair','dNME_max_UC_pair','dMML_max_UC_pair')]
    max_ts=cbind(data.table(regions=rownames(max_ts)),as.data.table(max_ts))
    #type='cairo' is needed for JHPCE
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dNME_max_UC_pair,x=UC_max_pair))+geom_point(alpha=0.01,color='blue')+geom_vline(xintercept=0.1,color='red')+
    xlab('UC')+ylab('dNME'))
    dev.off()
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dMML_max_UC_pair,x=UC_max_pair))+geom_point(alpha=0.01,color='blue')+geom_vline(xintercept=0.1,color='red')+
    xlab('UC')+ylab('dMML'))
    dev.off()
    max_ts$dNME_quant=ecdf(max_ts$dNME_max_UC_pair)(max_ts$dNME_max_UC_pair)
    max_ts$dMML_quant=ecdf(max_ts$dMML_max_UC_pair)(max_ts$dMML_max_UC_pair)
    max_ts$UC_quant=ecdf(max_ts$UC_max_pair)(max_ts$UC_max_pair)
    #ecdf(max_ts$UC_max_pair)(0.1)=0.8060806
    UC_cutoff_quant=ecdf(max_ts$UC_max_pair)(0.1)
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_quant.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dNME_quant*100,x=UC_quant*100))+geom_point(alpha=0.001,color='blue')+
    geom_vline(xintercept=UC_cutoff_quant,color='red')+xlab('UC quantile')+ylab('dNME quantile'))
    dev.off()
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_quant.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts,aes(y=dMML_quant*100,x=UC_quant*100))+geom_point(alpha=0.001,color='blue')+
    geom_vline(xintercept=UC_cutoff_quant,color='red')+xlab('UC quantile')+ylab('dMML quantile'))
    dev.off()
    max_ts$high_UC="smaller than 0.1"
    max_ts[UC_max_pair>0.1]$high_UC="greater than 0.1"
      png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_density_04.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts[dMML_max_UC_pair>0.4])+geom_density(aes(x=dMML_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0.4,1)
    dev.off()
        png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_density_all.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts)+geom_density(aes(x=dMML_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0,1)
    dev.off()
         png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_density_04.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts[dNME_max_UC_pair>0.4])+geom_density(aes(x=dNME_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0.4,1)
    dev.off()
        png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_density_all.png'),type='cairo',width=1000,height=1000)
    print(ggplot(max_ts)+geom_density(aes(x=dNME_max_UC_pair,color=high_UC))+
        xlab('dNME'))+xlim(0,1)
    dev.off()
    png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dNME_QQ.png'),type='cairo',width=1000,height=1000)
    qqplot(x=max_ts$UC_max_pair, y=max_ts$dNME_max_UC_pair,xlab='UC',ylab='dNME')
    abline(0,1,col='black')
    abline(v=0.1,col='red')
    dev.off()
      png(paste0('../downstream/output/mouse_analysis/UC_dNME_olap/',ts,'_UC_vs_dMML_QQ.png'),type='cairo',width=1000,height=1000)
    qqplot(x=max_ts$UC_max_pair, y=max_ts$dMML_max_UC_pair,xlab='UC',ylab='dMML')
    abline(0,1,col='black')
    abline(v=0.1,col='red')
    dev.off()
}