fhat.exp.ggplot=data.table()
for(i in 1:length(fhat.exp$x)){
  for(j in 1:length(fhat.exp$y)){
    fhat.exp.ggplot=rbind(fhat.exp.ggplot,
                          data.table(x=fhat.exp$x[i],y=fhat.exp$y[j],z_trans=fhat.exp$z_trans[i,j]))
    
    
  }
  
}
print(ggplot()+xlim(c(-1,1))+ylim(c(-1,1))+  geom_tile(data=fhat.exp.ggplot,aes(x=x,y=y,fill=z_trans))+
        scale_fill_gradient(low = "yellow", high = "dark red"))