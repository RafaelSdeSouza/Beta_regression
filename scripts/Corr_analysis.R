# Correlation Analysis to sub-select predictors
rm(list=ls(all=TRUE))
library(caret);
library(corrplot);library(reshape);
library(corrplot);library(caret);
library(Hmisc);
require(PerformanceAnalytics);
require(psych);library(corrr);
require(superheat);require(RColorBrewer)


Data=read.csv("..//data/FiBY.csv",header=T)

index <- sample(seq_len(nrow(Data)),replace=F, size = 40000)
data.2 = Data[,]


## fEsc is the variable of interest
Data <- as.data.frame(data.2[,c("Mstar","Mvir","ssfr_stars","baryon_fraction","spin","QHI","C")])
X    <- as.matrix(Data)
trans       <- preProcess(X,method = c("YeoJohnson", "center", "scale","spatialSign"))                                                                                       # distribute things around.  
Xtrans      <- predict(trans,X) # The "new" Transformed X                 


Xf <- cbind(X,fEsc=data.2$fEsc)
cxf<-cor(Xf)[-8, "fEsc"]


## Examine relationships between predictors before and after:
## It is advisable to analyze (scientifically and statistically) the relationships between predictors apart from the response variable
corrplot(cor(X,method="spearman"), order = "hclust",diag=F)
corrplot(cor(Xtrans,method="spearman"), order = "hclust",diag=F)

## Labels
names <- c("Mstar","M200", "sSFR",    
           "fb", "spin","QHI","C")   


require(tabplot)
X2 <- as.data.frame(cbind(redshift= data.2$redshift, fEsc=data.2$fEsc,X))
colnames(X2) <- c("z","fesc","Mstar","M200", "sSFR",    
                           "fb", "spin","QHI","C") 
tab<-tableplot(X2,fontsize = 10, legend.lines = 8)
tableSave(tab, filename="super_table.pdf",width=12,height = 6)



cx <- cor(X)
row.names(cx) <- names
colnames(cx) <- names


superheat(cx, 
          # color pallete
          heat.pal = brewer.pal(5, "PuBuGn"),
          # order rows/cols based on heirarchical clustering
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          # place dendrograms on columns and rows 
          row.dendrogram = F, 
          col.dendrogram = F,
          # make gridlines white for enhanced prettiness
          grid.hline.col = "white",
          grid.vline.col = "white",
          yt = cxf,
          yt.plot.type = "bar",
          yt.plot.size = 0.6,
          yt.axis.name = "Correlation with fesc",
          # rotate bottom label text
          bottom.label.text.angle = 90,
          legend.width=2.75)
quartz.save("superheat.pdf",type = "pdf",height=9,width=7.5)


## Cluster Analysis:
plot(varclus(X, similarity="spear")) 
plot(varclus(Xtrans, similarity="spear")) 


## Check boxplots
ggplot(data=melt(as.data.frame(scale(X))), aes(variable, value)) +
  geom_boxplot(fill="#33a02c",outlier.size = 0.5,outlier.colour = "grey")+
  theme_bw()+xlab("")+
  scale_x_discrete(labels=c(expression(M[star]),expression(M[gas]),
                            expression(M[200]),expression(SFR[gas]),expression(SSFR[gas]),expression(SFR[star]),expression(SSFR[star]),expression(f[b]),expression(tau),expression(lambda),expression(N[H]),expression(Q[HI]),"C")) +
  ylab("Untransformed  values")+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"),axis.text.x = element_text(angle = 90, hjust = 1))
# Print pdf

quartz.save(type = 'pdf', file = '../figures/box_raw.pdf',width = 9, height = 6)

ggplot(data=melt(as.data.frame(Xtrans)), aes(variable, value)) + geom_boxplot(fill="#33a02c",outlier.size = 0.5,outlier.colour = "grey")+
  theme_bw()+xlab("")+
  scale_x_discrete(labels=c(expression(M[star]),expression(M[gas]),
                            expression(M[200]),expression(SFR[gas]),expression(SSFR[gas]),expression(SFR[star]),expression(SSFR[star]),expression(f[b]),expression(tau),expression(lambda),expression(N[H]),expression(Q[HI]),"C")) +
  ylab("Transformed  values")+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = 0),
        text = element_text(size = 20,family="serif"),axis.text.x = element_text(angle = 90, hjust = 1))

quartz.save(type = 'pdf', file = '../figures/box_transf.pdf',width = 9, height = 6)



