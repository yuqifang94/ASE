suppressMessages(library(GenomicRanges))
d <- readRDS('/home-4/zji4@jhu.edu/scratch/tmp/UC_in_all.rds')
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
rn <- rownames(d)

type <- commandArgs(trailingOnly = T)
type="MML"
d <- readRDS(paste0('../downstream/output/mouse_analysis/CPEL_outputs/',type,'_matrix_mouse_all_dedup_N2_all_regions.rds'))
dname <- paste0(as.character(seqnames(d)),':',start(d),'-',end(d))
d <- as.matrix(mcols(d))
rownames(d) <- dname
d <- d[rn,]
d <- d[,!grepl('P0',colnames(d))]
d<-d[rowSums(is.na(d))==0,]
us <- colnames(d)
ma <- matrix(NA,nrow=length(us),ncol=length(us),dimnames = list(us,us))
for (i in 1:(length(us)-1)) {
  for (j in (i+1):length(us)) {
    ma[us[i],us[j]] <- ma[us[j],us[i]] <- mean(abs(d[,us[i]]-d[,us[j]]))
  }
}
diag(ma) <- 0
set.seed(12345)
p <- cmdscale(as.dist(ma))
library(ggplot2)
library(scales)
pd <- data.frame(MDS1=p[,1],MDS2=p[,2],tissue=sub('-.*','',rownames(p)),time=sub('.*-','',sub('-all','',rownames(p))),stringsAsFactors = F)
m <- aggregate(pd[,1:2],list(pd$tissue),mean)
colnames(m)[1] <- 'tissue'
cv <- brewer.pal(8,'Set1')
cv[6] <- 'goldenrod1'
ut <- seq(0.7,1,length.out = length(unique(pd$time)))
names(ut) <- unique(pd$time)

pdf(paste0('../downstream/output/mouse_analysis/MDS/',type,'_tissue.pdf'),width=5.5,height=4.5)
ggplot() + geom_point(data=pd,aes(x=MDS1,y=MDS2,col=tissue)) + theme_classic() + scale_color_manual(values=cv) + theme(legend.title = element_blank())
dev.off()
pdf(paste0('../downstream/output/mouse_analysis/MDS/',type,'_time.pdf'),width=5.5,height=4.5)
ggplot() + geom_point(data=pd,aes(x=MDS1,y=MDS2,col=time)) + theme_classic() + scale_color_manual(values=rev(viridis_pal()(length(unique(pd$time))))) + theme(legend.title = element_blank())
dev.off()

