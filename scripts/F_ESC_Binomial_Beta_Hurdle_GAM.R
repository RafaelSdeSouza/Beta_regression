rm(list=ls(all=TRUE))
library(mgcv);library(ggplot2);library(reshape2);library(ggthemes);library(MASS);library(hexbin);library(scales)
# Read data
Data=read.csv("..//data/FiBY.csv")
Data=subset(Data,redshift < 25)
##
## Log modulus transformation
L_M <-function(x){sign(x)*log10(abs(x) + 1)}

## fEsc is the variable of interest
## x is a vector of covariates
x         = c("Mstar","Mvir","ssfr_stars","baryon_fraction","spin","QHI","C") # with variable names 
var.names <- c("M[star]/M[sun]","M[200]/M[sun]", "sSFR/Gyrs^-1","f[b]", "lambda","Q[HI]/s^-1","C")   
Data      <-  Data[,c("fEsc",x)]
##
## Log Transform each variable except C
Data$Mstar       <- log10(Data$Mstar)
Data$Mvir        <- log10(Data$Mvir)
Data$ssfr_stars  <- L_M(Data$ssfr_stars)
Data$spin        <- log10(Data$spin)
Data$QHI         <- log10(Data$QHI)
Data$C           <- log10(Data$C)
# Transform to zero everything below 1e-3 
Data$fEsc[Data$fEsc < 10^-3] = 0
colnames(Data)[1] = "f_esc"
# Add 0/1 variable indicating f_esc = 0 vs f_esc > 0 
Data$non.zero.f_esc <- ifelse(Data$f_esc> 0, 1, 0)
n         = nrow(Data) # Number of Observations
n0        = sum(Data$f_esc == 0) # Number of zeros
p         = length(var.names) # Number of covariates
###
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
####
#### Modelling using Hurdle Bionmial_Beta_GAM 
#### Two stages
#### 1) Model Prob(f_esc>0) through Bionmial_GAM with logistic link 
## Binomial_GAM

Binomial_GAM         <- gam(non.zero.f_esc ~ s(Mstar,bs="cr",k=12) +  s(Mvir,bs="cr",k=12)  +  s(ssfr_stars,bs="cr",k=12) + s(baryon_fraction,bs="cr",k=25) +
                                             s(spin,bs="cr")  +  s(QHI,bs="cr",k=25)   +  s(C,bs="cr",k=20),                    
                                             data=Data,family= binomial(link="logit"),method="REML") 

summary(Binomial_GAM) 
######################################################################################
##    
##
gg          <-list()
gg_x        <- list()
gg_original <-list()
for(i in 1:p){
nn     = 3*10^4; 
R      = matrix(apply(Data[,x],2,median),nrow=1); 
R      =  R%x% rep(1,nn);colnames(R) = x; 
R      = as.data.frame(R)
a      = quantile(Data[,x[i]],0.001); b= quantile(Data[,x[i]],0.999); I = Data[,x[i]] > a & Data[,x[i]] < b
R[,i]  = seq(a,b,length=nn)    
#
# Predict and Produce confidence intervals:
#
Preds_nzero <- predict(Binomial_GAM,newdata = R,type="link",se=T,unconditional=T) 
fit.link    <- Preds_nzero$fit 
se          <- Preds_nzero$se  
CI.L        <- fit.link-2*se 
CI.R        <- fit.link+2*se 
CI          <- cbind(fit.link,CI.L,CI.R) 
CI          <- exp(CI)/(1+exp(CI)) # The first column correponds to the estimated probability of being non-zero.
colnames(CI) <- c("Predictions","CI_L","CI_R")
##
gg_x[[i]] <- data.frame(cbind(CI,x=R[,i]),var = rep(var.names[i],nn))
gg_original[[i]] <- data.frame(x=Data[I,x[i]],y=Data[I,"non.zero.f_esc"],var = rep(var.names[i],sum(I)))
}        

# put altogether for facets
ggg_x<-c()
for(i in 1:p){
ggg_x <- rbind(ggg_x,gg_x[[i]])
}  
ggg_original <- c()
for(i in 1:p){
  ggg_original <- rbind(ggg_original,gg_original[[i]])
}  

# Plot  via ggplot2
pdf("Binomial_GAM.pdf",width = 16,height = 8)
ggplot(ggg_x,aes(x=x,y=Predictions))+
  geom_hex(data=ggg_original,alpha=0.25,bins = 75,aes(x=x,y=y))+
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
####
####

#
#######################################################################################
#######################################################################################
## The Second Stage
### 2) Model  f_esc when f_esc > 0 throug Beta_GAM with logistic link
###
r <- 35
Beta_GAM    <- gam(f_esc ~ s(Mstar,bs="cr",k=r) +  s(Mvir,bs="cr",k=r)  +  s(ssfr_stars,bs="cr",k=r) + s(baryon_fraction,bs="cr",k=r) +
                           s(spin,bs="cr",k=r)  +  s(QHI,bs="cr",k=r)   +  s(C,bs="cr",k=r),
                           subset=f_esc>0,data=Data,family=betar(link="logit"),method="REML")

summary(Beta_GAM)
###########################
## 
##
gg          <-list()
gg_x        <- list()
gg_original <-list()
for(i in 1:p){
nn     = 3*10^4+1; 
R      = matrix(apply(Data[,x],2,median),nrow=1); 
R      =  R%x% rep(1,nn);colnames(R) = x; 
R      = as.data.frame(R)
a      = quantile(Data[,x[i]],0.001); b= quantile(Data[,x[i]],0.999); I = Data[,x[i]] > a & Data[,x[i]] < b
R[,i]  = seq(a,b,length=nn)
#
# Predict and Produce confidence intervals:
#
Preds_fesc  <- predict(Beta_GAM,newdata = R,type="link",se=T,unconditional=T) 
fit.link    <- Preds_fesc$fit 
se          <- Preds_fesc$se 
CI.L        <- fit.link-2*se 
CI.R        <- fit.link+2*se 
CI          <- cbind(fit.link,CI.L,CI.R) 
CI          <- exp(CI)/(1+exp(CI)) # The first column correponds to the estimated average of f_esc when f_esc > 0
colnames(CI) <- c("Predictions","CI_L","CI_R")
##
gg_x[[i]] <- data.frame(cbind(CI,x=R[,i]),var = rep(var.names[i],nn))
gg_original[[i]] <- data.frame(x=subset(Data,f_esc >0&I)[,x[i]],y=subset(Data,f_esc >0&I)[,"f_esc"],var = rep(var.names[i],sum(Data$f_esc>0&I)))
}        

# put altogether for facets
ggg_x<-c()
for(i in 1:p){
ggg_x <- rbind(ggg_x,gg_x[[i]])
}  
ggg_original <- c()
for(i in 1:p){
  ggg_original <- rbind(ggg_original,gg_original[[i]])
}  

# Plot  via ggplot2
pdf("Beta_GAM.pdf",width = 16,height = 8)
ggplot(ggg_x,aes(x=x,y=Predictions))+
  geom_hex(data=ggg_original,alpha=0.25,bins = 75,aes(x=x,y=y))+
  scale_fill_continuous(low = "#D9D3B4", high = "#441D0D", trans = log10_trans())+
  geom_ribbon(aes(ymin=CI_L, ymax=CI_R),fill = c("#3698BF"),alpha=0.75) +
  geom_line(col="#D97C2B",size=0.75)+
  theme_economist_white()+
  ylab(expression(paste("Average of ",~f[esc]," given that ", ~f[esc] > 0 ,sep="")))+
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


#######################################################################################
#######################################################################################
### Hurdle model:Putting things together
###
### Fitted values:
####
#### Produce hurdle plots, using delta method
gg          <-list()
gg_x        <- list()
gg_original <-list()
for(i in 1:p){
nn     = 3*10^4+1; 
R      = matrix(apply(Data[,x],2,median),nrow=1); 
R      = R%x% rep(1,nn);colnames(R) = x; 
R      = as.data.frame(R)
a      = quantile(Data[,x[i]],0.001); b= quantile(Data[,x[i]],0.999); I = Data[,x[i]] > a & Data[,x[i]] < b
R[,i]  = seq(a,b,length=nn)
#
fit_Binomial = predict(Binomial_GAM,newdata=R,type="response",se=T,unconditional=T) 
fit_Beta     = predict(Beta_GAM,newdata=R,type="response",se=T,unconditional=T)     
mu_Binomial  = fit_Binomial$fit
mu_Beta      = fit_Beta$fit
se_Binomial  = fit_Binomial$se 
se_Beta      = fit_Beta$se
##
mu_Hurdle    = mu_Binomial*mu_Beta
se_Hurdle    = sqrt(se_Binomial^2*mu_Beta^2  + mu_Binomial^2*se_Beta^2 + se_Binomial^2*se_Beta^2)
## 
phi <- Beta_GAM$family$getTheta(TRUE)
sd.y       <-  sqrt(mu_Binomial*(mu_Beta*(1-mu_Beta)/(1+phi) + (1-mu_Binomial)*mu_Beta^2))
#
CI.L        <- mu_Hurdle-2*se_Hurdle 
CI.R        <- mu_Hurdle+2*se_Hurdle 
CI          <- cbind(mu_Hurdle,CI.L,CI.R,sd.y) 
colnames(CI) <- c("Predictions","CI_L","CI_R","SD")
##
gg_x[[i]] <- data.frame(cbind(CI,x=R[,i]),var = rep(var.names[i],nn))
gg_original[[i]] <- data.frame(x=Data[I,x[i]],y=Data[I,"f_esc"],var = rep(var.names[i],sum(I)))
}        

# put altogether for facets
ggg_x<-c()
for(i in 1:p){
ggg_x <- rbind(ggg_x,gg_x[[i]])
}  
ggg_original <- c()
for(i in 1:p){
  ggg_original <- rbind(ggg_original,gg_original[[i]])
}  

# Plot  via ggplot2
pdf("Hurdle_GAM.pdf",width = 16,height = 8)
ggplot(ggg_x,aes(x=x,y=Predictions))+
  geom_hex(data=ggg_original,alpha=0.25,bins = 75,aes(x=x,y=y))+
  geom_line(aes(x=x,y=SD),size=0.75,linetype="dashed")+
  scale_fill_continuous(low = "#D9D3B4", high = "#441D0D", trans = log10_trans())+
  geom_ribbon(aes(ymin=CI_L, ymax=CI_R),fill = c("#3698BF"),alpha=0.75) +
  geom_line(col="#D97C2B",size=0.75)+
  theme_economist_white()+
  ylab(expression(paste(~f[esc],sep="")))+
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













