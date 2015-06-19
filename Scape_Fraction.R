# Bayesian  Logistic Regression using JAGS

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
#Read the already clean dataset

data.1= read.table(file="FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min")

data.2<-data.1[data.1$redshift==15.20530,]
data.3<-data.2[data.2$fEsc>0,]

jags.data <- list(Y= data.3$fEsc,
                  N = nrow(data.3),
                  x=(data.3$Mgas-mean(data.3$Mgas))/sd(data.3$Mgas)
#                  x=log(data.3$Mvir,10)
)


model<-"model{
#1. Priors
b0~dnorm(0,1e-5)
b1~dnorm(0,1e-5)
#c0~dnorm(0,1e-5)
#c1~dnorm(0,1e-5)
phi~dgamma(0.01,0.01)
#2. Likelihood
for(i in 1:N){
logit(mu[i]) <- b0 + b1*x[i]
#   sigma[i]<-exp(c0+c1*x[i])

#  p[i] <-max(0.000001,((1-mu[i])*mu[i]*mu[i]-mu[i]*sigma[i]^2)/sigma[i])
#  q[i] <- max(0.00001,(1-mu[i])*(mu[i]-mu[i]*mu[i]-sigma[i])/sigma[i])

p[i]<-mu[i]*phi
q[i]<-(1-mu[i])*phi
  Y[i] ~ dbeta(p[i],q[i])

#3. Prediction
prediction[i]~dbeta(p[i],q[i])
}
}"

params <- c("b0","b1","prediction")

inits1=list(b0=rnorm(1,0,1),b1=rnorm(1,0,1),c0=rnorm(1,0,1),c1=rnorm(1,0,1))
inits2=list(b0=rnorm(1,0,1),b1=rnorm(1,0,1),c0=rnorm(1,0,1),c1=rnorm(1,0,1))
inits3=list(b0=rnorm(1,0,1),b1=rnorm(1,0,1),c0=rnorm(1,0,1),c1=rnorm(1,0,1))

library(parallel)
cl <- makeCluster(3)
jags.logit <- run.jags(method="rjparallel", 
                       data = jags.data, 
                       inits = list(inits1,inits2,inits3),
                       model=model,
                       n.chains = 3,
                       adapt=1500,
                       monitor=c(params),
                       burnin=1000,
                       sample=5000,
                       summarise=FALSE,
                       plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)
summary<-extend.jags(jags.logit,drop.monitor=c("prediction"), summarise=TRUE)
print(summary)

# Diagnostics Confusion Matrix (very unlikely to be good)
predscape<-summary(as.mcmc.list(jags.logit, vars="prediction"))
predscape<-predscape$quantiles


plot(x,predscape[,3])
plot(x,data.3$fEsc)
