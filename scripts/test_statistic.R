require(ggplot2)
require(ggthemes)
require(dplyr)


#names <- c("Q[HI]","f[b]","C","M[200]","M['*']","sSFR","lambda") 


dat <- read.csv("importance.csv") %>% melt(id="X")
dat$X <- factor(dat$X, levels = rev(c("QHI","baryon_fraction","C","Mvir","Mstar","ssfr_stars","spin")))



pdf("test_stat_hurdle.pdf",width = 12,height = 5)
ggplot(dat,aes(x=X,y=value,group=variable,fill=variable)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels=rev(c(expression(Q[HI]),expression(f[b]),
                                "C",expression(M[200]),expression(M['*']),"sSFR",expression(lambda)))) +
  theme_economist_white()+
  scale_fill_manual(values=c("#477187","#D4B86A","#802D15"))+
  coord_flip() +
  ylab("Test statistic")+
  xlab("")+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif")) +
  facet_wrap(.~variable)
dev.off()


rank <- c(8838.27,6494.3,1385.06,1065.71, 720.38,520.97,67.07)




names <- c("Q[HI]","f[b]","C","M[200]","M['*']","sSFR","lambda") 
grank <- data.frame(rank,names)


#levels(grank$names) <- names

grank$names <- factor(grank$names, levels = rev(names))



pdf("test_stat_hurdle.pdf",width = 6,height = 5)
ggplot(grank,aes(x=names,y=rank)) +
  geom_bar(stat="identity",fill="#3698BF") +
  scale_x_discrete(labels=rev(c(expression(Q[HI]),expression(f[b]),
                  "C",expression(M[200]),expression(M['*']),"sSFR",expression(lambda)))) +
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
dev.off()