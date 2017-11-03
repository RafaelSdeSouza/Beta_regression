require(reshape2)
require(ggplot2)
require(ggthemes)
K=read.table("GAM.txt",header=T)
plot(K$x,K$y,col="Green",type="l",lwd=5)
lines(K$x,K$hat,lwd=5,col="red")
for( i in 1:10) lines(K$x, K[,paste0("B.",i)],lwd=2,lty=2)


g_data <- melt(K,id.vars = "x")

ggplot(g_data,aes(x = x,y = value,group = variable,linetype = variable,color=variable)) +
  geom_line(size=1.25) + scale_linetype_manual(name="",values=c("solid","dashed",rep("dotted",10))) +
  scale_color_manual(name="",values=c("green","red",rep("gray",10))) +
  theme_hc() + ylab("h(x)") + xlab("x")
  