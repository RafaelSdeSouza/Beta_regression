require(ggplot2)
require(ggthemes)
rank <- c(8838.27,6494.3,1385.06,1065.71,520.97,67.07)
names <- c("Q[HI]","f[b]","C","M[200]","sSFR","lambda") 
grank <- data.frame(rank,names)
#levels(grank$names) <- names

grank$names <- factor(grank$names, levels = rev(names))

grank$names <- parse(grank$names)
ggplot(grank,aes(x=names,y=rank)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels=rev(c(expression(Q[HI]),expression(f[b]),
                  "C",expression(M[200]),"sSFR",expression(lambda)))) +
  theme_economist_white()+
  coord_flip() +
  ylab("Test statistic")+
  xlab("")+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))