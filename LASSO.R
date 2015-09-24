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
library(lme4)
library(nlme)
library(arm)
require(gam)
require(glmnet)
require(AMADA)
require(mgcv)



#Read the  dataset
data.1= read.table(file="FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10")

data.1$Mvir<-log(data.1$Mvir,10)
data.1$Mstar<-log(data.1$Mstar,10)
data.1$Mgas<-log(data.1$Mgas,10)

#trainIndex <- createDataPartition(data.1$redshift, p = .25,
#                                  list = FALSE,
#                                  times = 1)
data.2<-data.1[data.1$redshift<=20,]


# Exploratory plot AMADA correlations 


#cor1<-Corr_MIC(data.2[,-2],method="pearson")
#plotdendrogram(cor1,type="p")

# remove collinearity 

VIF<-vifcor(data.2[,c(-1,-2)],th=0.7)
remain<-levels(VIF@results$Variables)

data.3<-as.data.frame(scale(data.2[,remain]))

N<-nrow(data.2)
data.3$Y<-(data.2$fEsc*(N-1)+0.5)/N

# Reionization cutoff

cut<-0.5

data.3$Y[data.3$Y>=cut]<-1
data.3$Y[data.3$Y<cut]<-0
#data.2$Y<-as.factor(data.2$Y)


# Fit LASSO binomial regression 

x2<-as.matrix(data.3[,remain])

fit<-glmnet(x2,y=data.3$Y,alpha=1,family="binomial")
plot(fit,xvar="lambda",label = TRUE)
plot(fit, xvar = "dev", label = TRUE)
cv.glmmod <- cv.glmnet(x2,y=data.3$Y,alpha=1,family="binomial",type.measure = "auc")

# Plot AUC vs Lambda
plot(cv.glmmod,xlab=expression(log~lambda),family="serif")


best_lambda <- cv.glmmod$lambda.min
coef.min = coef(cv.glmmod, s = "lambda.1se")
coef.min



# Fit GLM with selected features



fit=glm(Y~Mgas +ssfr_stars+age_star_min ,data=data.2,family=binomial("probit"))
ROCtest(fit,10,"ROC")
fit2=gam(Y~te(Mgas,ssfr_stars),data=data.2,family=binomial("logit"))
plot(fit2)
vis.gam(fit2,type="response",plot.type = "persp",color="topo", border=NA, n.grid=500,theta=-60,phi=30)
ROCtest(fit2,10,"ROC")


formula1 = Y~QHI+baryon_fraction+f(redshift,model="ar1")
mod.1 = inla(formula1,data=data.2,family="binomial")

formula2 = Y~Mstar+Mgas+ssfr_gas+ssfr_stars+baryon_fraction+spin+age_star_max+age_star_min+NH_10+f(redshift,model="ar1")
mod.2 = inla(formula2,data=data.2,family="binomial")


plot(mod.2)
hyper = inla.hyperpar(mod.2)
summary(hyper)
plot(hyper)


mydata <- as.data.frame(mod.surg$marginals.fixed)

ggplot(mydata) +   
  geom_line(aes(ID, `0.5quant`)) +   
  geom_line(aes(ID, `0.025quant`), linetype="dashed") +   
  geom_line(aes(ID, `0.975quant`), linetype="dashed")

model.mcmc = MCMClogit(Y~QHI+baryon_fraction,data=data.2,, mcmc=5000)





# Prepare data 
#data.2$Mstar<-(data.2$Mstar-mean(data.2$Mstar))/sd(data.2$Mstar)
#data.2$Mgas<-(data.2$Mgas-mean(data.2$Mgas))/sd(data.2$Mgas)
#data.2$Mvir<-(data.2$Mvir-mean(data.2$Mvir))/sd(data.2$Mvir)
#data.2$sfr_gas<-(data.2$sfr_gas-mean(data.2$sfr_gas))/sd(data.2$sfr_gas)
#data.2$ssfr_stars<-(data.2$ssfr_stars-mean(data.2$ssfr_stars))/sd(data.2$ssfr_stars)
#data.2$baryon_fraction<-(data.2$baryon_fraction-mean(data.2$baryon_fraction))/sd(data.2$baryon_fraction)
#data.2$QHI<-(data.2$QHI-mean(data.2$QHI))/sd(data.2$QHI)
#data.2$ssfr_gas<-(data.2$ssfr_gas-mean(data.2$ssfr_gas))/sd(data.2$ssfr_gas)
#data.2$age_star_mean<-(data.2$age_star_mean-mean(data.2$age_star_mean))/sd(data.2$age_star_mean)
#data.2$spin<-(data.2$spin-mean(data.2$spin))/sd(data.2$spin)
#data.2$NH_10<-(data.2$NH_10-mean(data.2$NH_10))/sd(data.2$NH_10)
#data.2$sfr_stars<-(data.2$sfr_stars-mean(data.2$sfr_stars))/sd(data.2$sfr_stars)
#data.2$age_star_max<-(data.2$age_star_max-mean(data.2$age_star_max))/sd(data.2$age_star_max)
#data.2$age_star_min<-(data.2$age_star_min-mean(data.2$age_star_min))/sd(data.2$age_star_min)
