#rm(list=ls(all=TRUE))
library(caret);library(visreg);
library(mgcv);library(ggplot2);
library(corrplot);library(reshape);
require(ggthemes);library(e1071);
library(scales);library(MASS);library(Hmisc)
library(corrplot)
Data=read.csv("..//data/FiBY.csv",header=T)


#index <- sample(seq_len(nrow(Data)),replace=F, size = 40000)
data.1 = Data[Data$redshift < 15,]


## fEsc is the variable of interest
data.2 <- as.data.frame(data.1[,c("Mstar","Mvir","ssfr_gas","ssfr_stars","baryon_fraction","spin","NH_10","QHI","C")])
X    <- as.matrix(data.2)
y    <- data.1$fEsc; 
y[ y < 10^-3] = 0

# Transform the columns of the design matrix (X)  to lighten the skewness, reduce the effect of outliers and reduce pairwise correlations if possible
#trans       <- preProcess(X,method = c("YeoJohnson", "center", "scale","spatialSign")) # Yeo-Johnson followed by centering and scaling. "spatialSign" is a bit complicated but it seems useful here
                                                                                       # to reduce outliers. Also I noticed that many of the variables are highly concetrated at one point. "spatailSign" will

trans       <- preProcess(X,method = c("YeoJohnson", "center", "scale","spatialSign"))                                                                                       # distribute things around.  
Xtrans      <- predict(trans,X) # The "new" Transformed X                 

### Two Models: 1) Model the probability that y > 0.  2) Model the Average of Y if y > 0 
### Some graphics:
### Plot the data (transformed) with a smoother:
### Good practice to go through all predictors (notice non-linearity):
non.zero <- ifelse(y > 0, 1, 0);fesc = as.factor(non.zero);levels(fesc) = c(expression("fesc = 0"), expression(fesc>0))
Dat.trans        <- data.frame(Xtrans,y=y,non.zero)
ggplot(Dat.trans, aes(x=Mstar, y=y, colour = fesc))   + geom_point(size=3,alpha = .5,pch=20) + geom_smooth(lwd=1.5,col="Blue")+ theme_bw()
ggplot(Dat.trans, aes(x=ssfr_gas, y=y, colour = fesc)) + geom_point(size=3,alpha = .5,pch=20) + geom_smooth(lwd=1.5,col="Blue")+ theme_bw()
ggplot(Dat.trans, aes(x=NH_10, y=y, colour = fesc))   + geom_point(size=3,alpha = .5,pch=20) + geom_smooth(lwd=1.5,col="Blue")+ theme_bw() #  weird !

####
#### Modelling using non-parametric hurdle regression model:
#### Non-parametric regression function we will adopt more general regression functions without restricting to specific functional assumption 
#### 1) Model Prob(y>0) using nonparametric logistic regression

simple_nzero       <- bam(non.zero ~ Mstar +  Mvir +  ssfr_gas + ssfr_stars+  baryon_fraction +  spin + NH_10 + QHI +C, data = Dat.trans, family = binomial(link = logit))
r                  <- 30
npar_nzero         <- bam(non.zero ~ s(Mstar,bs="cr",k=r) + s(Mvir,bs="cr",k=r)  + s(ssfr_gas,bs="cr",k=r)  + 
                             s(ssfr_stars,bs="cr",k=r)   + s(baryon_fraction,bs="cr",k=r) +
                            s(spin,bs="cr",k=r)      + s(NH_10,bs="cr",k=r)   + s(QHI,bs="cr",k=r) + s(C,bs="cr",k=r),
                    data=Dat.trans,family=binomial(link="logit"),gamma=1.4)  # Please use bam with cautious as it might diverege. For example
                                                                             # for the beta regression below it did not converge although it did not directly mention that.
                                                                             # For this reason I abandoned it and retruned to gam as it is more thorough and better built 
                                                                             # In any case we can still use bam if it works but we just need to be careful.

anova.gam(simple_nzero,npar_nzero,test="Chisq") # Test the simple model against the more complicated one
summary(npar_nzero)
plot(npar_nzero,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) # Non-linearity does not seem obvious for the some of the variables 
gam.check(npar_nzero) # Residual analysis

######################################################################################
#### Produce the previous plots on the original scale of the predictors
### Mvir is an example
nn = 10^4+1; XMvir = matrix(apply(X,2,median),nrow=1); XMvir=XMvir%x% rep(1,nn);colnames(XMvir) = colnames(X); XMvir = as.data.frame(XMvir)    
XMvir$Mvir       = seq(min(X[,"Mvir"]),max(X[,"Mvir"]),length=nn)
XMvir.trans      = predict(trans,XMvir) # Convert to the transformed scale
# Predict and Produce confidence intervals:
Preds_nzero <- predict(npar_nzero,newdata = XMvir.trans,type="link",se=T,unconditional=T) 
fit.link    <- Preds_nzero$fit 
se          <- Preds_nzero$se  # Standard errors 
CI.L        <- fit.link-qnorm(0.975)*se 
CI.R        <- fit.link+qnorm(0.975)*se 
CI          <- cbind(fit.link,CI.L,CI.R) 
CI          <- exp(CI)/(1+exp(CI)) # The first column correponds to the estimated probability of being non-zero.
colnames(CI) <- c("Predictions","CI_L","CI_R")
### One can ago ahead and plot Mvir against CI 
# Using ggplot2

# Data for ggplot2
gg_mvir <- as.data.frame(cbind(CI,Mvir=XMvir$Mvir))
gg_original <- data.frame(x=data.2$Mvir,y=non.zero)

# Plot  via ggplot2
ggplot(gg_mvir,aes(x=Mvir,y=Predictions))+
#  geom_hex(data=gg_original,bins = 100,aes(x=x,y=y),alpha=0.3,fill=c("#000000"))+
  geom_point(data = gg_original,aes(x = x,y = y),size = 1,alpha = 0.1,col = "gray75",position = position_jitter(height  = 0.025))+
  geom_ribbon(aes(x = Mvir,y=Predictions,ymin=CI_L, ymax=CI_R),fill = c("#33a02c")) +
  geom_line(col="white",size=1.5)+
  theme_bw()+
  ylab(expression(paste(f[esc] > 0.1,"%",sep="")))+
  xlab(expression(M[200]))+
  scale_x_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))
        
        
### On the transformed scale:
###
gg_mvir_trans     <- as.data.frame(cbind(CI,Mvir=XMvir.trans$Mvir))
gg_trans <- data.frame(x=Dat.trans$Mvir,y=y)
# Plot  via ggplot2
ggplot(gg_mvir_trans,aes(x=Mvir,y=Predictions))+
  geom_point(data=gg_trans,aes(x=x,y=y),size=1,alpha=0.2,col="orange2",position = position_jitter (h = 0.025))+
  geom_ribbon(aes(x=Mvir,y=Predictions,ymin=CI_L, ymax=CI_R),fill=c("#33a02c")) +
  geom_line(col="white",size=1.5)+
  theme_bw()+
  ylab(expression(f[esc]))+
  xlab(expression(M[200]))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))


# Or
### If you still to use visreg. do the following
#Plot = visreg(npar_nzero,"Mvir",scale = "response",nn=nn,partial=T)
# replace the last three columns of Plot$fit by CI and replace Mvir as well 
#D  = dim(Plot$fit)[2]
#Plot$fit[,(D-2):D] = CI ;  Plot$fit$Mvir = log10(XMvir$Mvir) # Transformation is recommneded so the scale can make some sense instead from 0-10^24
                                                             #  
                                                             # Also, we can change the ticks of the x-axis to 10^x  
#plot(Plot,ylab = expression(paste(f[esc] > 0.1,"%",sep="")),line=list(col="white"), points=list(cex=0.05, pch=3,col="orange"),
#       fill.par=list(col=c('#33a02c')),rug = 2,xlab=expression(M[200]),partial=T)


### 2) Model Average y when y > 0 using non parametric beta regression model 
simple_Beta_y   <- gam(y ~ Mstar +  Mvir  + ssfr_gas +  ssfr_stars+  baryon_fraction +  spin + NH_10 + QHI +C,subset=y>0,data=Dat.trans,family=betar(link="logit"),gamma=1.4)
r <- 30
npar_Beta_y <- gam(y ~ s(Mstar,bs="cr",k=r)         +  s(Mvir,bs="cr",k=r)       + 
                       s(ssfr_gas,bs="cr",k=r)      +  s(ssfr_stars,bs="cr",k=r) + s(baryon_fraction,bs="cr",k=r) +
                        s(spin,bs="cr",k=r)       + s(NH_10,bs="cr",k=r)      + s(QHI,bs="cr",k=r) + s(C,bs="cr",k=r),
                       subset=y>0,data=Dat.trans,family=betar(link="logit"),gamma=1.4)

anova.gam(simple_Beta_y,npar_Beta_y,test="Chisq") # Test the simple model against the more complicated one
summary(npar_Beta_y)
plot(npar_Beta_y,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) 
gam.check(npar_Beta_y) # Residual analysis

# Predict and Produce confidence intervals:
Preds_y     <- predict(npar_Beta_y,newdata = XMvir.trans,type="link",se=T,unconditional=T) 
fit.link    <- Preds_nzero$fit 
se          <- Preds_nzero$se  # Standard errors 
CI.L        <- fit.link-qnorm(0.975)*se 
CI.R        <- fit.link+qnorm(0.975)*se 
CI          <- cbind(fit.link,CI.L,CI.R) 
CI           <- exp(CI)/(1+exp(CI)) # The first column correponds to the estimated probability of being non-zero.
colnames(CI) <- c("Predictions","CI_L","CI_R")
## On the transformed scale:
gg_mvir_trans     <- as.data.frame(cbind(CI,Mvir=XMvir.trans$Mvir))
gg_trans <- data.frame(x=Dat.trans$Mvir,y=y)

## For the original scale: replace the above two lines with
##gg_mvir <- as.data.frame(cbind(CI,Mvir=XMvir$Mvir))
##gg_original <- data.frame(x=Data$Mvir,y=y)

# Plot  via ggplot2
ggplot(gg_mvir_trans,aes(x=Mvir,y=Predictions))+
  geom_point(data=gg_trans,aes(x=x,y=y),size=1,alpha=0.2,col="orange2",position = position_jitter (h = 0.025))+
  geom_ribbon(aes(x=Mvir,y=Predictions,ymin=CI_L, ymax=CI_R),fill=c("#33a02c")) +
  geom_line(col="white",size=1.5)+
  theme_bw()+
  ylab(expression(f[esc]))+
  xlab(expression(M[200]))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))




#visreg(npar_Beta_y,"Mvir",scale = "response",rug = 2,ylab = expression(f[esc]),line=list(col="white"),
#       points=list(cex=0.05, pch=3,col="orange"), fill.par=list(col=c('#33a02c')),partial=T) # Plot using visreg

#### Finally: 
#### Prediction :
#### Say at xp = X[c(100,363,987),] 
xp= X[c(100,363,987),]; xp=as.data.frame(xp);colnames(xp) = colnames(X)
xp.trans <- predict(trans,xp) # The "new" Transformed xp      
# Find pr(y > 0) at xp 
predict(npar_nzero,newdata=xp.trans,type="response")
# Predict the average of y at xp (Important) : 
predict(npar_nzero,newdata=xp.trans,type="response")*predict(npar_Beta_y,newdata=xp.trans,type="response")

## Fitted values:
fit = predict(npar_nzero,newdata=as.data.frame(Xtrans),type="response")*predict(npar_Beta_y,newdata=as.data.frame(Xtrans),type="response")
cor(fit,y)
plot(fit,y,pch=20,cex=0.5)
abline(0,1,lty=2,col="Black",lwd=3)


#### To monitor the changes in the average of y as one variables is changing keeping all other predictors fixed
#### These are important plots as one can notice the importance of each predictor to y
### Mvir is an example:
Preds_nzero <- predict(npar_nzero ,newdata = XMvir.trans,type="response",se=T,unconditional=T)
Preds_y     <- predict(npar_Beta_y,newdata = XMvir.trans,type="response",se=T,unconditional=T)
mu1         <- Preds_nzero$fit;mu2 <- Preds_y$fit;
se1         <- Preds_nzero$se ;se2 <- Preds_y$se;
Predictions <- mu1*mu2
se          <- sqrt(mu2^2*se1^2 + mu1^2*se2^2 + se1^2*se2^2)
CI.L        <- Predictions-qnorm(0.975)*se 
CI.R        <- Predictions+qnorm(0.975)*se 
CI          <- cbind(Predictions,CI.L,CI.R) 
colnames(CI) <- c("Predictions","CI_L","CI_R")

gg_mvir <- as.data.frame(cbind(CI,Mvir=XMvir$Mvir))
gg_original <- data.frame(x=data.1$Mvir,y=y)


# Plot  via ggplot2
ggplot(gg_mvir,aes(x=Mvir,y=Predictions))+
 # geom_point(data=gg_original,aes(x=x,y=y),size=1,alpha=0.2,col="orange2",position = position_jitter (h = 0.025))+
  geom_hex(data=gg_original,bins = 50,aes(x=x,y=y))+
  scale_fill_continuous(low = "gray90", high = "gray10", trans = log10_trans())+
  geom_ribbon(aes(x=Mvir,y=Predictions,ymin=CI_L, ymax=CI_R),fill=c("#33a02c")) +
  geom_line(col="white",size=1.5)+
  theme_bw()+
  ylab(expression(f[esc]))+
  xlab(expression(M[200]))+
  scale_x_continuous(trans = 'log10',
                     breaks=trans_breaks("log10",function(x) 10^x),
                     labels=trans_format("log10",math_format(10^.x)))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))

### On the transformed scale:

###
gg_mvir_trans     <- as.data.frame(cbind(CI,Mvir=XMvir.trans$Mvir))
gg_trans <- data.frame(x=Dat.trans$Mvir,y=y)
# Plot  via ggplot2
ggplot(gg_mvir_trans,aes(x=Mvir,y=Predictions))+
  geom_point(data=gg_trans,aes(x=x,y=y),size=1,alpha=0.2,col="orange2",position = position_jitter (h = 0.025))+
  geom_ribbon(aes(x=Mvir,y=Predictions,ymin=CI_L, ymax=CI_R),fill=c("#33a02c")) +
  geom_line(col="white",size=1.5)+
  theme_bw()+
  ylab(expression(f[esc]))+
  xlab(expression(M[200]))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))




