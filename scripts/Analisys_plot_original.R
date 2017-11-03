#rm(list=ls(all=TRUE))
library(caret);library(visreg);
library(mgcv);library(ggplot2);
library(corrplot);library(reshape2);
require(ggthemes);library(e1071);
library(scales);library(MASS);library(Hmisc)
library(corrplot);library(gridExtra)

# Read data
Data=read.csv("..//data/FiBY.csv",header=T)

#subsampling 
#index <- sample(seq_len(nrow(Data)),replace=F, size = 40000)

# Cut in redshift
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

#data.2$spin <- log10(data.2$spin)
X    <- as.matrix(data.2)
y    <- data.1$fEsc; 
y[ y < 10^-3] = 0 # transform to zero everything below 1e-2
#y_b <- (y*(nrow(data.2)-1)+0.5)/nrow(data.2)
  


# Transform the columns of the design matrix (X)  to lighten the skewness, reduce the effect of outliers and reduce pairwise correlations if possible
#trans       <- preProcess(X,method = c("YeoJohnson", "center", "scale","spatialSign")) # Yeo-Johnson followed by centering and scaling. "spatialSign" is a bit complicated but it seems useful here
                                                                                       # to reduce outliers. Also I noticed that many of the variables are highly concetrated at one point. "spatailSign" will

trans       <- preProcess(X,method = c("YeoJohnson","center", "scale"))                                                                                       # distribute things around.  
Xtrans      <- predict(trans,X) # The "new" Transformed X                 


box1 = melt(as.data.frame(scale(X)));box1$case <- rep("Original  scale",nrow(box1))
box2 = melt(as.data.frame(scale(Xtrans)));box2$case <- rep("Transformed  scale",nrow(box2))
boxjoint <- rbind(box1,box2)
boxjoint$case <- as.factor(boxjoint$case)

## Check boxplots
ggplot(data=boxjoint, aes(variable, value)) +
  geom_boxplot(fill="#ba122b",outlier.size = 0.5,outlier.colour = "grey80",colour="#CCCC99")+
  theme_bw()+xlab("")+
  scale_x_discrete(labels=c(expression(M[star]),expression(M[200]),expression(sSFR[stars]),expression(f[b]),expression(log(lambda)),expression(Q[HI]),"C")) +
  ylab("")+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"),axis.text.x = element_text(angle = 90, hjust = 1))+
        facet_wrap(~case,ncol=1,scale="free_y")
# Print pdf

quartz.save(type = 'pdf', file = '../figures/box_raw.pdf',width = 9, height = 8)

ggplot(data=melt(as.data.frame(Xtrans)), aes(variable, value)) + 
  geom_boxplot(fill="#ba122b",outlier.size = 0.5,outlier.colour = "grey80",colour="#CCCC99")+
  theme_bw()+xlab("")+
  scale_x_discrete(labels=c(expression(M[star]),expression(M[200]),
                            expression(sSFR[gas]),expression(sSFR[stars]),expression(f[b]),expression(log(lambda)),expression(Q[HI]),"C")) +
  ylab("Transformed  scale")+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="top",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = 0),
        text = element_text(size = 20,family="serif"),axis.text.x = element_text(angle = 90, hjust = 1))

quartz.save(type = 'pdf', file = '../figures/box_transf.pdf',width = 9, height = 6)



### Two Models: 1) Model the probability that y > 0.  2) Model the Average of Y if y > 0 
### Some graphics:
### Plot the data (transformed) with a smoother:
### Good practice to go through all predictors (notice non-linearity):
non.zero <- ifelse(y > 0, 1, 0);fesc = as.factor(non.zero);levels(fesc) = c(expression("fesc = 0"), expression(fesc>0))
Dat.trans        <- data.frame(Xtrans,y=y,non.zero)



####
#### Modelling using non-parametric hurdle regression model:
#### Non-parametric regression function we will adopt more general regression functions without restricting to specific functional assumption 
#### 1) Model Prob(y>0) using nonparametric logistic regression

simple_nzero       <- gam(non.zero ~ Mstar +  Mvir + ssfr_stars+  baryon_fraction +  spin +  QHI +C , data = Dat.trans, family = binomial(link = logit))
r                  <- 30
npar_nzero         <- gam(non.zero ~ s(Mstar,bs="cr",k=r) + s(Mvir,bs="cr",k=r)  + 
                             s(ssfr_stars,bs="cr",k=r)   + s(baryon_fraction,bs="cr",k=r) +
                            s(spin,bs="cr",k=r)        + s(QHI,bs="cr",k=r) + s(C,bs="cr",k=r),
                    
                    data=Dat.trans,family=binomial(link="logit"),gamma=1.4)  # Please use bam with cautious as it might diverge. For example
                                                                             # for the beta regression below it did not converge although it did not directly mention that.
                                                                             # For this reason I abandoned it and retruned to gam as it is more thorough and better built 
                                                                             # In any case we can still use bam if it works but we just need to be careful.

#anova.gam(simple_nzero,npar_nzero,test="Chisq") # Test the simple model against the more complicated one
#summary(npar_nzero)
#plot(npar_nzero,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) # Non-linearity does not seem obvious for the some of the variables 
#gam.check(npar_nzero) # Residual analysis


######################################################################################
#### Produce the previous plots on the original scale of the predictors
### loop over all 
#names for plot 
#names <- c("log~(M[star]/M[sun])","log~(M[200]/M[sun])", "logM~(sSFR[stars]/Gyrs^-1)",    
#"f[b]", "log~(lambda)","log~Q[HI]/s^-1","C")    

names <- c("M[star]","M[200]", "sSFR",    
"f[b]", "lambda","Q[HI]","C")    


gg<-list()
gg_x <- list()
gg_original <-list()
for(i in 1:ncol(X)){
nn = 10^4+1; XX = matrix(apply(X,2,mean),nrow=1); 
XX=XX%x% rep(1,nn);colnames(XX) = colnames(X); 
XX = as.data.frame(XX)
XX[,i] <- seq(min(X[,i]),max(X[,i]),length=nn)
XX.trans      = predict(trans,XX) # Convert to the transformed scale

# Predict and Produce confidence intervals:
Preds_nzero <- predict(npar_nzero,newdata = XX.trans,type="link",se=T,unconditional=T) 
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
gg_x[[i]] <- data.frame(cbind(CI,x=XX[,i]),var = rep(names[i],length(XX[,i])))
gg_original[[i]] <- data.frame(x=data.2[,i],y=non.zero,var = rep(names[i],length(data.2[,i])))
}        

# put altogether for facets
ggg_x<-c()
for(i in 1:7){
ggg_x <- rbind(ggg_x,gg_x[[i]])
}  
ggg_original <- c()
for(i in 1:7){
  ggg_original <- rbind(ggg_original,gg_original[[i]])
}  

# Plot  via ggplot2
pdf("logit_part_3.pdf",width = 16,height = 8)
ggplot(ggg_x,aes(x=x,y=Predictions))+
  geom_hex(data=ggg_original,bins = 75,aes(x=x,y=y))+
  scale_fill_continuous(low = "#D9D3B4", high = "#441D0D", trans = log10_trans())+
  geom_ribbon(aes(ymin=CI_L, ymax=CI_R),fill = c("#3698BF"),alpha=0.75) +
  geom_line(col="#D97C2B",size=0.75)+
  theme_economist_white()+
  ylab(expression(paste("Probability of ",~f[esc] > 0,sep="")))+
  xlab("")+
  theme(legend.background = element_rect(fill="white"),
       legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))+
  facet_wrap(~var,scales = "free_x",ncol=4,labeller = label_parsed)
dev.off()






### 2) Model Average y when y > 0 using non parametric beta regression model 
simple_Beta_y   <- gam(y ~ Mstar +  Mvir  +   ssfr_stars+  baryon_fraction +  spin +  QHI +C,subset=y>0,data=Dat.trans,family=betar(link="logit"),gamma=1.4)
r <- 30
npar_Beta_y <- gam(y ~ s(Mstar,bs="cr",k=r)         +  s(Mvir,bs="cr",k=r)       + 
                       s(ssfr_stars,bs="cr",k=r) + s(baryon_fraction,bs="cr",k=r) +
                        s(spin,bs="cr",k=r)          + s(QHI,bs="cr",k=r) + s(C,bs="cr",k=r),
                       subset=y>0,data=Dat.trans,family=betar(link="logit"),gamma=1.4)

anova.gam(simple_Beta_y,npar_Beta_y,test="Chisq") # Test the simple model against the more complicated one
summary(npar_Beta_y)
plot(npar_Beta_y,pages=1,residuals=F,scheme=1,rug=FALSE,lwd=3,shade=TRUE,seWithMean=TRUE) 
gam.check(npar_Beta_y) # Residual analysis





#### Produce the previous plots on the original scale of the predictors
### loop over all 
gg2<-list()
pp_x <- list()
pp_original <- list()
for(i in 1:ncol(X)){
  nn = 10^4+1; XX = matrix(apply(X,2,mean),nrow=1); 
  XX=XX%x% rep(1,nn);colnames(XX) = colnames(X); 
  XX = as.data.frame(XX)
  XX[,i] <- seq(min(X[,i]),max(X[,i]),length=nn)
  XX.trans      = predict(trans,XX) # Convert to the transformed scale
# Predict and Produce confidence intervals:
Preds_y     <- predict(npar_Beta_y,newdata = XX.trans,type="link",se=T,unconditional=T) 
fit.link    <- Preds_y$fit 
se          <- Preds_y$se  # Standard errors 
CI.L        <- fit.link-qnorm(0.975)*se 
CI.R        <- fit.link+qnorm(0.975)*se 
CI          <- cbind(fit.link,CI.L,CI.R) 
CI           <- exp(CI)/(1+exp(CI)) # The first column correponds to the estimated probability of being non-zero.
colnames(CI) <- c("Predictions","CI_L","CI_R")
## On the transformed scale:
#gg_mvir_trans     <- as.data.frame(cbind(CI,Mvir=XMvir.trans$Mvir))
#gg_trans <- data.frame(x=Dat.trans$Mvir,y=y)


## For the original scale: replace the above two lines with
pp_x[[i]] <- data.frame(cbind(CI,x=XX[,i]),var = rep(names[i],length(XX[,i])))
pp_original[[i]] <- data.frame(x=data.2[,i],y=y,var = rep(names[i],length(data.2[,i])))

}

# put altogether for facets
ppp_x<-c()
for(i in 1:7){
  ppp_x <- rbind(ppp_x,pp_x[[i]])
}  
ppp_original <- c()
for(i in 1:7){
  ppp_original <- rbind(ppp_original,pp_original[[i]])
} 
# Plot  via ggplot2
pdf("beta_part_3.pdf",width = 16,height = 8)
ggplot(ppp_x,aes(x=x,y=Predictions))+
  geom_hex(data=ppp_original,bins = 50,aes(x=x,y=y),alpha=0.25)+
  scale_fill_continuous(low = "#D9D3B4", high = "#441D0D",trans = log_trans())+
  #  geom_point(data=gg_original,aes(x=x,y=y),size=1,alpha=0.2,col="orange2",position = position_jitter (h = 0.025))+
  geom_ribbon(aes(ymin=CI_L, ymax=CI_R),fill= c("#3698BF"),alpha=0.75) +
  geom_line(col="#D97C2B",size=0.75)+
  theme_economist_white()+
  ylab(expression(f[esc]))+
  xlab("")+
  #  scale_x_continuous(trans = 'log10',
  #                     breaks=trans_breaks("log10",function(x) 10^x),
  #                     labels=trans_format("log10",math_format(10^.x)))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))+
  facet_wrap(~var,scales = "free_x",ncol=4,labeller = label_parsed)
dev.off()







#### To monitor the changes in the average of y as one variables is changing keeping all other predictors fixed
#### These are important plots as one can notice the importance of each predictor to y
gg3<-list()
vv_x <- list()
vv_original <-list()
for(i in 1:ncol(X)){
  nn = 10^4+1; XX = matrix(apply(X,2,mean),nrow=1); 
  XX=XX%x% rep(1,nn);colnames(XX) = colnames(X); 
  XX = as.data.frame(XX)
  XX[,i] <- seq(min(X[,i]),max(X[,i]),length=nn)
  XX.trans      = predict(trans,XX) # Convert to the transformed scale
##
Preds_nzero <- predict(npar_nzero ,newdata = XX.trans,type="response",se=T,unconditional=T)
Preds_y     <- predict(npar_Beta_y,newdata = XX.trans,type="response",se=T,unconditional=T)
mu1         <- Preds_nzero$fit;mu2 <- Preds_y$fit;
se1         <- Preds_nzero$se ;se2 <- Preds_y$se;
Predictions <- mu1*mu2
se          <- sqrt(mu2^2*se1^2 + mu1^2*se2^2 + se1^2*se2^2)
CI.L        <- Predictions-qnorm(0.975)*se 
CI.R        <- Predictions+qnorm(0.975)*se 
CI          <- cbind(Predictions,CI.L,CI.R) 
colnames(CI) <- c("Predictions","CI_L","CI_R")

vv_x[[i]] <- data.frame(cbind(CI,x=XX[,i]),var = rep(names[i],length(XX[,i])))
vv_original[[i]] <- data.frame(x=data.2[,i],y=y,var = rep(names[i],length(data.2[,i])))

}

# put altogether for facets
vvv_x<-c()
for(i in 1:7){
  vvv_x <- rbind(vvv_x,vv_x[[i]])
}  
vvv_original <- c()
for(i in 1:7){
  vvv_original <- rbind(vvv_original,vv_original[[i]])
} 

# Plot  via ggplot2
pdf("hurdle_3.pdf",width = 16,height = 8)
ggplot(vvv_x,aes(x=x,y=Predictions))+
  # geom_point(data=gg_original,aes(x=x,y=y),size=1,alpha=0.2,col="orange2",position = position_jitter (h = 0.025))+
  geom_hex(data=vvv_original,bins = 50,aes(x=x,y=y),alpha=0.25)+
  scale_fill_continuous(low = "#D9D3B4", high = "#441D0D",trans = log10_trans())+
  geom_ribbon(aes(ymin=CI_L, ymax=CI_R),fill=c("#3698BF"),alpha=0.75) +
  geom_line(col="#D97C2B",size=0.75)+
  theme_economist_white()+
  ylab(expression(f[esc]))+
  xlab("")+
  #    scale_y_continuous(trans = 'log10',
  #                       breaks=trans_breaks("log10",function(x) 10^x),
  #                       labels=trans_format("log10",math_format(10^.x)))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))+
  facet_wrap(~var,scales = "free_x",ncol=4,labeller = label_parsed)
dev.off()




#### Finally: 
#### Prediction :
#### Say at xp = X[c(100,363,987),] 
#xp= X[c(100,363,987),];
xp= X; xp=as.data.frame(xp);colnames(xp) = colnames(X)
xp.trans <- predict(trans,xp) # The "new" Transformed xp      
# Find pr(y > 0) at xp 
predict(npar_nzero,newdata=xp.trans,type="response")
# Predict the average of y at xp (Important) : 
predict(npar_nzero,newdata=xp.trans,type="response")*predict(npar_Beta_y,newdata=xp.trans,type="response")

## Fitted values:
fit = predict(npar_nzero,newdata=as.data.frame(Xtrans),type="response",se=T)$fit*predict(npar_Beta_y,newdata=as.data.frame(Xtrans),type="response",se=T)$fit
 



cor(fit,y)
plot(fit,y,pch=20,cex=0.5)
abline(0,1,lty=2,col="Black",lwd=3)







