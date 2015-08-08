# Exploration of the dataset
require(AMADA)
require(caret)
require(reshape)
require(ggthemes)
library(scales)
data.1= read.table(file="FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10")


trainIndex <- createDataPartition(data.1$redshift, p = .3,
                                  list = FALSE,
                                  times = 1)
#data.2<-data.1[data.1$redshift==8.86815,]
data.2<-data.1[trainIndex,]
#data.2<-data.1[data.1$redshift==8.86815,]
#data.2<-data.1
N<-nrow(data.2)


data.2$Y<-(data.2$fEsc*(N-1)+0.5)/N
data.2$Y[data.2$Y>=0.1]<-1
data.2$Y[data.2$Y<0.1]<-0
data.2$Y<-as.factor(data.2$Y)
ggdata<-melt(data.2[,-2],id=("Y"))
g1<-ggplot(ggdata,aes(x=Y,y=value),group=Y,colour=Y,fill=Y)

g2<-g1+geom_boxplot(data=ggdata,aes(group=Y,fill=Y),
                    outlier.colour="gray60",colour="gray60",fatten=2.5, outlier.size=1.5,notch=F)+
  scale_y_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  scale_fill_stata(guide="none")+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  facet_wrap(~variable,scales="free",ncol=4)+
  xlab("")+theme(plot.title = element_text(hjust=0.5),axis.title=element_text(vjust=-0.25),text = element_text(size=20))+
  ylab("")

pdf("box.pdf",height = 12,width = 12)
g2
dev.off()

cor1<-Corr_MIC(data.2,method="pearson")
plotgraph(cor1)
