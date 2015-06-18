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




jags.data <- list(Y= data.1$fEsc,
                  N = nrow(data.1),
                  x=data.1$Mvir
)


model<-"model{
#1. Priors

tau.R<-pow(sdBeta,-1)
sdBeta ~ dgamma(0.001,0.001)
#alpha1~dnorm(0,0.00001)
#tau.z~dgamma(0.001,0.001)

# Random intercept 
for (j in 1:Ntype){
ranef[j]~ddexp(0,tau.R)
#ranef[j]~dnorm(alpha1,tau.z)
}

beta.0~dnorm(0,0.001)
beta.1~dnorm(0,0.001)
beta.2~dnorm(0,0.001)

#2. Likelihood
for (i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- beta.0+beta.1*mag_g[i]+beta.2*bar[i]+ranef[galtype[i]]
#eta[i] <- beta.0+beta.1*mag_g[i]+ranef[galtype[i]]
#3. Prediction
prediction[i]~dbern(pi[i])
}
}"
