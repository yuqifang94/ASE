af <- list.files('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/MOCA/data/saver')
for (f in af) {
  expr <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/MOCA/data/saver/',f))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.frame(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  saveRDS(res,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/MOCA/res/',f))
}

#Plot mean vs var and MAV
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=12,face="bold"),
                                 axis.text.x=element_text(size=10),
                                 axis.text.y=element_text(size=10))
HESC_scRNA=readRDS(paste0(scRNA_dir,"HESC_1.rds"))
HESC_scRNA=as.data.frame(HESC_scRNA)
pdf("../downstream/output/human_analysis/QC/MAV_mean.pdf",height=3.5,width=7)
MAV_mean=ggplot(HESC_scRNA,aes(x=log(mean),y=hypervar_logvar))+geom_point(alpha=0.1,size=0.5)+geom_smooth()+xlab("log(mean)")+ylab("MAV")+theme_glob
var_mean=ggplot(HESC_scRNA,aes(x=log(mean),y=log(var)))+geom_point(alpha=0.1,size=0.5)+geom_smooth()+xlab("log(mean)")+ylab("log(var)")+theme_glob
print(grid.arrange(var_mean,MAV_mean,nrow=1))
dev.off()