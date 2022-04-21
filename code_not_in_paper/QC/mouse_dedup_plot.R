library(ggplot2)
library(data.table)
mouseDupIn=fread('../downstream/output/mouse_analysis/QC/duplication_rate_fn.csv')
ggplot(mouseDupIn,aes(x=Sample,y=duplication*100))+geom_point()+geom_text(label=mouseDupIn$Sample,angle=270)+
  theme(axis.text.x=element_blank(),axis.text.y= element_text(size=14),axis.title.y=element_text(size=20))+
  ylab("duplication(%)")+xlab('')+ylim(c(0,65))

