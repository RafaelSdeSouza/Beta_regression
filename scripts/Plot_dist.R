# Plot Normal, Bernoulli and Beta distributions for paper
require(ggplot2)
require(ggthemes)

x <- seq(0, 1, length = 100)

a <- c(0.5,1,5,1,2,2)
b <- c(0.5,1,1,3,2,5)
  

y <- list()
for (i in 1:6){  
y[[i]] <- dbeta(x, a[i], b[i])
}


b1 <- data.frame(x,y=y[[1]])
b2 <- data.frame(x,y=y[[2]])
b3 <- data.frame(x,y=y[[3]])
b4 <- data.frame(x,y=y[[4]])
b5 <- data.frame(x,y=y[[5]])
b6 <- data.frame(x,y=y[[6]])


b_all <- rbind(b1,b2,b3,b4,b5,b6)

b_all$case <- rep(c("a = 0.5 b = 0.5", "a = 1 b = 1",  "a = 5 b = 1","a = 1 b = 3","a = 2 b = 2","a = 2 b = 5"),each=100)

pdf("beta_dist.pdf",width = 10,height = 8)
ggplot(b_all,aes(x=x,y=y,group=case,color=case,linetype=case)) +
  geom_line(name="",size=1.25) + scale_color_stata(name="") +
  scale_linetype_stata(name="") +
  theme_economist_white() +
  theme(legend.key.width = unit(2.5, "cm"),legend.text=element_text(size=12.5),legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position=c(0.5,0.875),
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif")) +
    ylab("Beta Probability Density Function") +
  guides(col = guide_legend(nrow=3))
dev.off()

