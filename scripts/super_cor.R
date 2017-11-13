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





# Histogram of f_esc

fEsc <- data.frame(fEsc=data.1$fEsc)
fEsc[fEsc < 10^-3] = 0


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

  ct <- cor.test(x,y,method = "spearman" )

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
  corColors <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9",
                 "#74add1","#4575b4","#313695")

  if (r <= -0.8) {
    corCol <- corColors[1]
  } else if (r > -0.8 & r <= -0.6  ) {
    corCol <- corColors[2]
  } else if (r > -0.6 & r <= -0.4) {
    corCol <- corColors[3]
  } else if (r > -0.4 & r <= -0.2) {
    corCol <- corColors[4]
  } else if (r > -0.2 & r <= 0) {
    corCol <- corColors[5]
  } else if (r > 0 & r <= 0.2) {
    corCol <- corColors[6]
  } else if (r > 0.2 & r <= 0.4) {
    corCol <- corColors[7]
  } else if (r > 0.4 & r <= 0.6) {
    corCol <- corColors[8]
  } else if (r > 0.6 & r <= 0.8) {
    corCol <- corColors[9]
    } else {
    corCol <- corColors[10]
  }
  p <- p +
    theme(panel.background = element_rect(fill= corCol),
                                          panel.grid.major = element_line(),
                                          panel.grid.major.x = element_blank(),
                                          panel.grid.major.y = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          panel.grid.minor.x = element_blank(),
                                          panel.grid.minor.y = element_blank()
  )

  p
}



## Labels
names <- c("M[star]","M[200]", "sSFR",
           "f[b]", "lambda","Q[HI]","C","f[esc]")



pm <- ggpairs(
  X,
  upper = list(continuous = my_custom_cor_color ),
  diag=list(continuous= my_hist),
  lower = list(continuous = my_bin),
  columnLabels = names,
  labeller = "label_parsed"
)


p2 <- pm +
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 18.5,family="serif"),
        axis.text.x = element_text(face="bold", color="gray15",
                                   size=12, angle=45),
        axis.text.y = element_text(face="bold", color="gray15",
                                   size=12))


png("../figures/super_cor.png",width = 1.75*480, height = 1.75*480)
p2
dev.off()





