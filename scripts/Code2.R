rm(list=ls(all=TRUE))
library(caret);library(visreg);library(mgcv);library(ggplot2);library(corrplot);library(reshape);require(ggthemes)
Data=read.table("..//data/FiBY_escape_data_all.dat",header=F)
data.1 = Data
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10","C")
data.2 = data.1[data.1$redshift<=10,]
## fEsc is the variable of interest
Data <- as.data.frame(data.2[,c("Mstar","Mgas","Mvir","sfr_gas","ssfr_gas","sfr_stars","ssfr_stars","baryon_fraction","age_star_mean","spin","NH_10","QHI","C")])
X    <- as.matrix(Data)
y    <- data.2$fEsc; 
y[ y < 10^-3] = 0

# Transform the columns of the design matrix (X)  to lighten the skewness, reduce the effect of outliers and reduce pairwise correlations if possible
trans       <- preProcess(X,method = c("YeoJohnson", "center", "scale","spatialSign")) # Yeo-Johnson followed by centering and scaling. "spatialSign" is a bit complicated but it seems useful here
                                                                                       # to reduce outliers. Also I noticed that many of the variables are highly concetrated at one point. "spatailSign" will
                                                                                       # distribute things around.  
Xtrans      <- predict(trans,X) # The "new" Transformed X                 
## Check skewness Before and After transformation
apply(X,2,function(x) skewness(x)) # All are exteremly skewed except perhaps age_star_mean
apply(Xtrans,2,function(x) skewness(x)) # Much Improved. 
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

#quartz.save(type = 'pdf', file = '../figures/box_raw.pdf',width = 9, height = 6)

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

#quartz.save(type = 'pdf', file = '../figures/box_transf.pdf',width = 9, height = 6)

### Two Models: 1) Model the probability that y > 0. 2) Model the Average of Y if y > 0 
### Some graphics:
### Plot the data (transformed) with a smoother:
### Good practice to go through all predictors (notice non-linearity):
non.zero <- ifelse(y > 0, 1, 0)
d        <- data.frame(Xtrans,y=y,non.zero)
ggplot(d, aes(x=Mstar, y=y, colour = as.factor(non.zero))) + geom_point(size=3,alpha = .5,pch=20) + geom_smooth(lwd=1.5,col="Blue")+ theme_bw()
ggplot(d, aes(x=sfr_gas, y=y, colour = as.factor(non.zero))) + geom_point(size=3,alpha = .5,pch=20) + geom_smooth(lwd=1.5,col="Blue")+ theme_bw()
ggplot(d, aes(x=NH_10, y=y, colour = as.factor(non.zero))) + geom_point(size=3,alpha = .5,pch=20) + geom_smooth(lwd=1.5,col="Blue")+ theme_bw() #  weird !

####
#### Modelling:
#### Using Non-parametric regression function we will adopt more general regression functions without restricting to specific functional assumption 
#### 1) Model Prob(y>0)
mod_linear_nzero   <- gam(non.zero ~ Mstar + Mgas + Mvir + sfr_gas + ssfr_gas + sfr_stars+ ssfr_stars+  baryon_fraction + age_star_mean + spin + NH_10 + QHI +C, data = d, family = binomial(link = logit))
r                  <- 30
M_non.zero         <- gam(non.zero ~ s(Mstar,bs="cr",k=r)    + s(Mgas,bs="cr",k=r) + s(Mvir,bs="cr",k=r) + s(sfr_gas,bs="cr",k=r) + s(baryon_fraction,bs="cr",k=r) +
                    s(ssfr_gas,bs="cr",k=r) + s(age_star_mean,bs="cr",k=r) + s(spin,bs="cr",k=r) + s(NH_10,bs="cr",k=r) + s(QHI,bs="cr",k=r),
                    data=d,family=binomial(link="logit"))



anova.gam(mod_linear_nzero,M_non.zero,test="Chisq") # Test the simple model against the more complicated one
summary(M_non.zero)
plot(M_non.zero,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) # Non-linearity does not seem obvious for the some of the variables 
gam.check(M_non.zero) # Residual analysis


# Plot using visreg
visreg(M_non.zero,"Mstar",ylab = expression(f[esc] > 0),line=list(col="#33a02c"), points=list(cex=0.25, pch=2,col="grey80"),
       fill.par=list(col=c('blue')),scale = "response",rug = 2,partial = TRUE)

### 2) Model Average y when y > 0. 
## Competing models     1) log normal (simple but can produce predicted values greater than 1)
##                      2) beta    (correct scale but more complicated) 
##                      3) Gamma models (can handle skewness naturally)
# Log normal
M_L_y    <- gam(log(y) ~ s(Mstar,bs="cr",k=100)    + s(Mgas,bs="cr",k=100) + s(Mvir,bs="cr",k=100) + s(sfr_gas,bs="cr",k=100) + s(baryon_fraction,bs="cr",k=100) +
                s(ssfr_gas,bs="cr",k=100) + s(age_star_mean,bs="cr",k=100) + s(spin,bs="cr",k=100) + s(NH_10,bs="cr",k=100) + s(QHI,bs="cr",k=100),
                subset=y>0,data=d,gamma=1.4)

summary(M_L_y)
plot(M_L_y,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) 
gam.check(M_L_y) # Residual analysis


# Beta  
M_Beta_y <- gam(y ~ s(Mstar,bs="cr",k=100)    + s(Mgas,bs="cr",k=100) + s(Mvir,bs="cr",k=100) + s(sfr_gas,bs="cr",k=100) + s(baryon_fraction,bs="cr",k=100) +
                    s(ssfr_gas,bs="cr",k=100) + s(age_star_mean,bs="cr",k=100) + s(spin,bs="cr",k=100) + s(NH_10,bs="cr",k=100) + s(QHI,bs="cr",k=100),
                     subset=y>0,data=d,family=betar(link="logit"),gamma=1.4)

summary(M_Beta_y)
plot(M_Beta_y,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) 
gam.check(M_Beta_y) # Residual analysis

visreg(M_Beta_y,"QHI",scale = "response",rug = 2,ylab = expression(f[esc])) # Plot using visreg

# Gamma 
M_Gamma_y <- gam(y ~ s(Mstar,bs="cr",k=100)    + s(Mgas,bs="cr",k=100) + s(Mvir,bs="cr",k=100) + s(sfr_gas,bs="cr",k=100) + s(baryon_fraction,bs="cr",k=100) +
                     s(ssfr_gas,bs="cr",k=100) + s(age_star_mean,bs="cr",k=100) + s(spin,bs="cr",k=100) + s(NH_10,bs="cr",k=100) + s(QHI,bs="cr",k=100),
                     subset=y>0,data=d,family=Gamma(link="log"),gamma=1.4)

summary(M_Gamma_y)
plot(M_Gamma_y,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) 
gam.check(M_Gamma_y) # Residual analysis


#### 
#### Prediction :
#### Say at xp = X[c(100,363,987),] 
xp= X[c(100,363,987),]; xp=as.data.frame(xp);colnames(xp) = colnames(X)
xp.trans <- predict(trans,xp) # The "new" Transformed xp      
# Find pr(y = 0) at xp 
predict(M_non.zero,newdata=xp.trans,type="response")
# Predict the average of y at xp: 
predict(M_non.zero,newdata=xp.trans,type="response")*predict(M_Beta_y,newdata=xp.trans,type="response")



## Fitted values:
fit = predict(M_non.zero,newdata=as.data.frame(Xtrans),type="response")*predict(M_Beta_y,newdata=as.data.frame(Xtrans),type="response")
cor(fit,y)
plot(fit,y,pch=20,cex=0.5)
abline(0,1,lty=2,col="Black",lwd=3)







