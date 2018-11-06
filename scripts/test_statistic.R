require(ggplot2)
require(ggthemes)
require(dplyr)



#names <- c("Q[HI]","f[b]","C","M[200]","M['*']","sSFR","lambda") 


dat <- read.csv("importance.csv") %>% melt(id="X") %>%
 mutate(alpha = (value/max(value))^{1/10})

dat$X <- factor(dat$X, levels = rev(c("QHI","baryon_fraction","C","Mvir","Mstar","sfr_stars","spin")))



pdf("test_stat_hurdle.pdf",width = 12,height = 5)
ggplot(dat,aes(x=X,y=value,group=variable,fill=X)) +
  geom_bar(stat="identity") +
  scale_x_discrete(labels=rev(c(expression(Q[HI]),expression(f[b]),
                                "C",expression(M[200]),expression(M['*']),"sSFR",expression(lambda)))) +
  theme_economist_white()+
  scale_fill_manual(values=rev(c('#a6761d','#66a61e','#666666','#d95f02','#1b9e77',
                             '#7570b3','#e6ab02')))+


#  scale_fill_manual(values=c("#477187","#D4B86A","#802D15"))+
#  scale_y_log10() +
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




