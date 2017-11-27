# Correlation Analysis to sub-select predictors
rm(list=ls(all=TRUE))
library(caret);
library(corrplot);library(reshape2);
library(Hmisc);
require(PerformanceAnalytics);
require(psych);library(corrr);
require(superheat);require(RColorBrewer)
require(scales);require(GGally);require(ggthemes)


Data=read.csv("..//data/FiBY.csv",header=T)


data.1 = Data[Data$redshift < 25,]

# Log modulus transformation
L_M <-function(x){sign(x)*log10(abs(x) + 1)}

## fEsc is the variable of interest
data.2 <- as.data.frame(data.1[,c("Mstar","Mvir","ssfr_stars","baryon_fraction","spin","QHI","C")])
data.2$Mstar <- log10(data.2$Mstar)
data.2$Mvir <- log10(data.2$Mvir)
#data.2$ssfr_gas <- L_M(data.2$ssfr_gas)
data.2$ssfr_stars <- L_M(data.2$ssfr_stars)
data.2$spin <- log10(data.2$spin)
data.2$QHI  <- log10(data.2$QHI)
data.2$C  <- log10(data.2$C)


#index <- sample(seq_len(nrow(Data)),replace=F, size = 50000)


# Histogram of f_esc

fEsc <- data.frame(fEsc=data.1$fEsc)
fEsc[fEsc < 10^-3] = 0


ggplot(fEsc, aes(x=fEsc)) + geom_histogram(aes(y=..count../sum(..count..)),size=1.5,breaks = c(0,seq(0.0001,1,by=0.1)),fill="#3698BF",colour="#D9D384") + theme_classic() +
xlab(expression(f[esc])) + ylab("Fraction of Galaxies")  +
# scale_y_continuous(trans = 'log10', breaks = trans_breaks('log10', function(x) 10^x),
#labels = trans_format('log10', math_format(10^.x))) +
  theme(text = element_text(size = 20,family="serif")) 
#+coord_cartesian(ylim=c(1e0,5e4))




quartz.save(type = 'pdf', file = '../figures/hist_fesc.pdf',width = 9, height = 6)


## fEsc is the variable of interest
X    <- data.frame(data.2,fEsc = fEsc)

my_fn <- function(data, mapping, ...){
 ggplot(data = data, mapping = mapping) + 
    geom_density2d(...) 
}

my_bin <- function(data, mapping, ..., low = "#3698BF", high = "#D97C2B") {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    scale_fill_gradient(low = low, high = high,trans = log10_trans()) +
    theme_bw()
}


my_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(bins = 10,fill="#3698BF",colour="#D9D384",...) +
    theme_bw()
}


my_custom_cor_color <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  
  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ct <- cor.test(x,y)
  
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]
  tt <- as.character(rt)
  
  # plot the cor value
  p <- ggally_text(
    label = tt, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 6,
    color=color,
    ...
  ) 
  corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdYlGn")[2:6]
  
  if (r <= -0.8) {
    corCol <- corColors[1]
  } else if (r <= -0.6) {
    corCol <- corColors[2]
  } else if (r < 0.6) {
    corCol <- corColors[3]
  } else if (r < 0.8) {
    corCol <- corColors[4]
  } else {
    corCol <- corColors[5]
  }
  p <- p + theme(
    panel.background = element_rect(fill= corCol)
  )
  
  p
}





pm <- ggpairs(
  X,
  upper = list(continuous = my_custom_cor_color ),
  diag=list(continuous= my_hist),
  lower = list(continuous = my_bin)
)


p2 <- pm +  
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 10,family="serif")) 

p2







Xf <- data.frame(X,fEsc=fEsc)
cxf<-cor(Xf)[-8, "x"]


## Examine relationships between predictors before and after:
## It is advisable to analyze (scientifically and statistically) the relationships between predictors apart from the response variable
corrplot(cor(X,method="spearman"), order = "hclust",diag=F)
corrplot(cor(Xtrans,method="spearman"), order = "hclust",diag=F)

## Labels
names <- c("Mstar","M200", "sSFR",    
           "fb", "spin","QHI","C")   


require(tabplot)
X2 <- as.data.frame(cbind(redshift= data.2$redshift, fEsc=fEsc,X))
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



