require(reshape2)
require(ggplot2)
require(ggthemes)
K=read.table("GAM.txt",header=T)
plot(K$x,K$y,col="Green",type="l",lwd=5)
lines(K$x,K$hat,lwd=5,col="red")
for( i in 1:10) lines(K$x, K[,paste0("B.",i)],lwd=2,lty=2)


g_data <- melt(K,id.vars = "x")

pdf("hx.pdf",width = 6,height = 4)
ggplot(g_data,aes(x = x,y = value,group = variable,linetype = variable,color=variable)) +
  geom_line(size=1.25) + scale_linetype_manual(name="",values=c("solid","dashed",rep("dotted",10))) +
  scale_color_manual(name="",values=c("#3698BF","#A63921",rep("#D9D384",10))) +
  scale_size_manual(name="",values=c(1.5,1.25,rep(1,10))) +
  theme_economist_white() + ylab("h(x)") + xlab("x") + 
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
#        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
#        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 17.5))
dev.off()