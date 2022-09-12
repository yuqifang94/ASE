source('mainFunctions_sub.R')
suppressMessages(library(GenomicRanges))
d <- readRDS(UC_in_MDS_all_file)
rn <- d[,1]
d <- as.matrix(d[,-1])
rownames(d) <- rn
d <- d[,!grepl('P0',colnames(d))]

mm <- c("EFP","forebrain","heart","hindbrain","limb","liver","midbrain")
library(RColorBrewer)
library(viridis)
library(gridExtra)
ct <- sub('_.*','',do.call(rbind,strsplit(sub('-all','',colnames(d)),'-')))
d <- d[,rowMeans(cbind(ct[,1] %in% mm,ct[,2] %in% mm))==1]
d <- d[rowMeans(is.na(d)) == 0,]
d <- d[rowSums(d > 0.1) >= 1,]
m <- colMeans(d)
names(m) <- sub('-all','',names(m))
st <- do.call(rbind,strsplit(names(m),'-'))
us <- unique(as.vector(st))
ma <- matrix(NA,nrow=length(us),ncol=length(us),dimnames = list(us,us))
for (i in 1:nrow(st)) ma[st[i,1],st[i,2]] <- m[i]
for (i in 1:nrow(st)) ma[st[i,2],st[i,1]] <- m[i]
diag(ma) <- 0
p <- cmdscale(as.dist(ma))
library(ggplot2)
pd <- data.frame(MDS1=p[,1],MDS2=p[,2],tissue=sub('_.*','',rownames(p)),time=sub('.*_','',rownames(p)),stringsAsFactors = F)
m <- aggregate(pd[,1:2],list(pd$tissue),mean)
colnames(m)[1] <- 'tissue'
cv <- brewer.pal(7,'Set1')
cv[6] <- 'goldenrod1'
names(cv) <- setdiff(unique(pd$tissue),'Center')
ut <- seq(0.7,1,length.out = length(unique(pd$time)))
names(ut) <- unique(pd$time)

pd <- rbind(pd,data.frame(MDS1=-0.015,MDS2=0,tissue='Center',time='Center'))

contcv <- rev(viridis_pal()(length(unique(pd$time))-1))
contcv <- c(contcv,'red')
names(contcv) <- c(setdiff(unique(pd$time),'Center'),'Center')
saveRDS(pd,"../downstream/output/mouse_analysis/MDS/pd.rds")
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                 axis.title.x=element_text(hjust=0.5,size=7,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=7,face="bold"),
                                 axis.text.x=element_text(size=7),
                                 axis.text.y=element_text(size=7))
pdf('../downstream/output/mouse_analysis/MDS/UC_tissue.pdf',width=4,height=3)
ggplot() + geom_point(data=pd[pd$tissue!='Center',],aes(x=MDS1,y=MDS2,col=tissue)) + theme_glob + theme(legend.title = element_blank())
dev.off()
pdf('../downstream/output/mouse_analysis/MDS/UC_time.pdf',width=4,height=3)
ggplot() + geom_point(data=pd,aes(x=MDS1,y=MDS2,col=time)) + theme_glob+ scale_color_manual(values=contcv) + theme(legend.title = element_blank())
dev.off()


