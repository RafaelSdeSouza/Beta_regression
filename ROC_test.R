#  Required libraries
library(rjags)
library(ggmcmc)
library(ggplot2)
library(ggthemes)
library(pander)
library(Cairo)
library(plyr)
library(MASS)
library(scales)
library(plyr)
require(gdata)
require(runjags)
require(gdata)
require(caret)
require(pROC)
require(plyr)
require(LOGIT)
require(usdm)

#Read the  dataset

data.1= read.table(file="FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10")


#trainIndex <- createDataPartition(data.1$redshift, p = .25,
#                                  list = FALSE,
#                                  times = 1)
data.2<-data.1[data.1$redshift<=10,]
data.2<-data.1[trainIndex,]
#data.2<-data.1[data.1$redshift==8.86815,]
#data.2<-data.1
N<-nrow(data.2)



data.2$Y<-(data.2$fEsc*(N-1)+0.5)/N
data.2$Y[data.2$Y>=0.1]<-1
data.2$Y[data.2$Y<0.1]<-0
#data.2$Y<-as.factor(data.2$Y)


# Prepare data for JAGS
data.2$Mstar<-(data.2$Mstar-mean(data.2$Mstar))/sd(data.2$Mstar)
data.2$Mgas<-(data.2$Mgas-mean(data.2$Mgas))/sd(data.2$Mgas)
data.2$Mvir<-(data.2$Mvir-mean(data.2$Mvir))/sd(data.2$Mvir)
data.2$sfr_gas<-(data.2$sfr_gas-mean(data.2$sfr_gas))/sd(data.2$sfr_gas)
data.2$baryon_fraction<-(data.2$baryon_fraction-mean(data.2$baryon_fraction))/sd(data.2$baryon_fraction)
data.2$QHI<-(data.2$QHI-mean(data.2$QHI))/sd(data.2$QHI)
data.2$ssfr_gas<-(data.2$ssfr_gas-mean(data.2$ssfr_gas))/sd(data.2$ssfr_gas)
data.2$age_star_mean<-(data.2$age_star_mean-mean(data.2$age_star_mean))/sd(data.2$age_star_mean)
data.2$spin<-(data.2$spin-mean(data.2$spin))/sd(data.2$spin)
data.2$NH_10<-(data.2$NH_10-mean(data.2$NH_10))/sd(data.2$NH_10)
data.2$sfr_stars<-(data.2$sfr_stars-mean(data.2$sfr_stars))/sd(data.2$sfr_stars)

x<-as.matrix(data.2[,c("Mvir","sfr_gas","baryon_fraction","ssfr_gas","age_star_mean","spin","NH_10")])
vifcor(x,th=0.75)
#fit=glm(Y~QHI+baryon_fraction+redshift,data=data.2,family=binomial("probit"))

fit<-glmnet(x,y=data.2$Y,alpha=1,family="binomial")


plot(fit,xvar="lambda",label = TRUE)
plot(fit, xvar = "dev", label = TRUE)
cv.glmmod <- cv.glmnet(x,y=data.2$Y,alpha=1,family="binomial",type.measure = "class")
plot(cv.glmmod)
best_lambda <- cv.glmmod$lambda.min

coef.min = coef(cv.glmmod, s = "lambda.min")
active.min = coef.min[which(abs(coef.min) > 0.1)]
index.min = coef.min[active.min]



fit=glm(Y~Mstar+Mvir,data=data.2,family=binomial("probit"))
ROCtest(fit,10,"ROC")


fit2=glm(Y~sfr_gas,data=data.2,family=binomial("logit"))


ROCtest(fit2,10,"ROC")


formula1 = Y~QHI+baryon_fraction+f(redshift,model="ar1")
mod.1 = inla(formula1,data=data.2,family="binomial")

formula2 = Y~QHI+baryon_fraction
mod.2 = inla(formula2,data=data.2,family="binomial")


plot(mod.surg)
hyper = inla.hyperpar(mod.surg)
summary(hyper)
plot(hyper)


mydata <- as.data.frame(mod.surg$marginals.fixed)

ggplot(mydata) +   
  geom_line(aes(ID, `0.5quant`)) +   
  geom_line(aes(ID, `0.025quant`), linetype="dashed") +   
  geom_line(aes(ID, `0.975quant`), linetype="dashed")

model.mcmc = MCMClogit(Y~QHI+baryon_fraction,data=data.2,, mcmc=5000)
