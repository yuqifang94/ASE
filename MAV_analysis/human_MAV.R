#This code is copied from Dr.Jason Ji only changing file path
#Load gene matrix
library(Matrix)
library(data.table)
m <- readMM('../downstream/data/human_expression_MOCA/gene_count.txt.gz')
g <- fread('../downstream/data/human_expression_MOCA/gene_annotate.csv',data.table=F)
c <- fread('../downstream/data/human_expression_MOCA/cell_annotate.csv',data.table=F)
id <- which(!duplicated(g[,3]))
m <- m[id,]
rownames(m) <- g[id,3]
colnames(m) <- paste0(gsub(' ','',c$Main_cell_type),':',c$development_stage,'_cell',1:nrow(c))
saveRDS(m,file='../downstream/input/human_analysis/NME_expression_var/proc/count.rds')
#Runing SAVER
library(SAVER)
m <- readRDS('../downstream/input/human_analysis/NME_expression_var/proc/count.rds')
p <- sub(':.*','',colnames(m))
m <- m[,p %in% c('Cardiacmusclelineages','Limbmesenchyme','NeuralTube')]
m <- m[,colSums(m > 0) >= 500]
m <- m[rowSums(m) > 0,]
rc <- colSums(m)
rc <- rc/median(rc)
m <- t(t(m)/rc)
p <- sub(paste0('_.*'),'',colnames(m))
print(sort(table(p)))
for (tp in unique(p)) {
  d <- saver(m[,p==tp],ncores=20,size.factor=1)$estimate  
  saveRDS(d,file=paste0('../downstream/input/human_analysis/NME_expression_var/saver/',tp,'.rds'))
}

af <- list.files('../downstream/input/human_analysis/NME_expression_var/saver')
for (f in af) {
  expr <- readRDS(paste0('../downstream/input/human_analysis/NME_expression_var/saver/',f))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.frame(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  saveRDS(res,file=paste0('../downstream/input/human_analysis/NME_expression_var/scRNA/',f))
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