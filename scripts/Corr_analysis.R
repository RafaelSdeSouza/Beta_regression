# Correlation Analysis to sub-select predictors
rm(list=ls(all=TRUE))
library(corrplot);library(reshape2);
library(Hmisc);
require(PerformanceAnalytics);
require(psych);
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



ggplot(fEsc, aes(x=fEsc)) + geom_histogram(aes(y=..count../sum(..count..)),size=1.5,breaks = c(0,seq(0.0001,1,by=0.05)),fill="#4271AE",colour = "gray80") + theme_classic() +
#ggplot(fEsc, aes(x=fEsc)) + geom_histogram(aes(y=..count../sum(..count..)),size=1.5,breaks = c(0,seq(0.0001,1,by=0.05)),fill="#3698BF",colour="#D9D384") + theme_classic() +
#ggplot(fEsc, aes(x=fEsc)) + geom_histogram(aes(y=..count../sum(..count..)),size=1.5,breaks = c(0,seq(0.0001,1,by=0.1)),fill="#4271AE",colour = "gray80") + theme_classic() +
xlab(expression(f[esc])) + ylab("Fraction of Galaxies")  +
# scale_y_continuous(trans = 'log10', breaks = trans_breaks('log10', function(x) 10^x),
#labels = trans_format('log10', math_format(10^.x))) +
  theme(text = element_text(size = 20,family="serif")) 
#+coord_cartesian(ylim=c(1e0,5e4))

require(kdensity)

kde = kdensity(fEsc[[1]], kernel = "beta")
plot(kde, lwd = 2, col = "blue")

rug(diff(LakeHuron))



quartz.save(type = 'pdf', file = '../figures/hist_fesc.pdf',width = 9, height = 6)


## fEsc is the variable of interest
X    <- data.frame(data.2,fEsc = fEsc)


equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    round(seq(min(x), max(x), length=n),1)
  }
}



my_fn <- function(data, mapping, ...){
 ggplot(data = data, mapping = mapping) + 
    geom_density2d(...) 
}

my_bin <- function(data, mapping, ..., low = "#3698BF", high = "#D97C2B") {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    scale_fill_gradient(low = low, high = high,trans = log10_trans()) +
    theme_bw() +
    scale_y_continuous(round(seq(min(data), max(data), length.out = 3),1))
}


my_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(bins = 10,fill="#4271AE",colour="#1F3552",...) +
    theme_void() + theme( panel.grid.minor=element_blank(),
                        panel.grid.major=element_blank())
}


my_custom_cor_color <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  
  # get the x and y data to use the other code
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  
  ct <- cor.test(x, y, method = "spearm")
  
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]
  tt <- as.character(rt)
  
  
  # plot the cor value
  p <- ggally_text(
    label = tt, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 6,
    color = color,
    ...
  ) 
  corColors <- c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",
                 "#d1e5f0","#92c5de","#4393c3","#2166ac","#053061")
  
  if (r <= -0.8) {
    corCol <- corColors[1]
   
 } else if (r > -0.8  && r <= -0.6) {
    corCol <- corColors[2]
  } 
  else if (r > -0.6  && r <= -0.4) {
    corCol <- corColors[3]
  } 
  else if (r > -0.4  && r <= -0.2) {
    corCol <- corColors[4]
  } 
  else if (r > -0.2  && r <= 0) {
    corCol <- corColors[5]
  } 
  else if (r > 0  && r <= 0.2) {
    corCol <- corColors[6]
  } 
  else if (r > 0.2  && r <= 0.4) {
    corCol <- corColors[7]
  } 
  else if (r > 0.4  && r <= 0.6) {
    corCol <- corColors[8]
  }
  else if (r > 0.6  && r <= 0.8) {
    corCol <- corColors[9]
  } 
  else {
    corCol <- corColors[10]
  }
  p <- p + theme(
    panel.background = element_rect(fill = corCol),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) 
  
  p
}






pm <- ggpairs(
  X, columnLabels = c("M[star]","M[200]", "sSFR","f[b]", "lambda","Q[HI]","C","f[esc]"), 
  labeller = "label_parsed",
  upper = list(continuous = my_custom_cor_color ),
  diag = list(continuous = my_hist),
  lower = list(continuous = my_bin
  )
)



p2 <- pm +  
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 16,family="serif")) 

png("../figures/super_cor.png", width = 12, height = 10.5, units = 'in',res = 1200)
p2
dev.off()


quartz.save(type = 'pdf', file = '../figures/super_cor.pdf',width = 11, height = 10.5)




