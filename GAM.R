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

covariate<-round(runif(250,20,70))
response<-rnorm(250,13.5,1)

knot0<-c()
for(i in 1:15){
  knot0[i]<-i/(15+1)
}

nos<-quantile(covariate,knot0)
jags.data <-list(n=250,
                 nknots=15,
                 degree=2,
                 knot=as.numeric(nos))

model<-"model{ 
#Begin model

#Likelihood of the model
for (i in 1:n){
response[i]~dnorm(m[i],taueps)
m[i]<-inprod(beta[],X[i,])+inprod(b[],Z[i,])}

#Prior distributions of the random effects parameters
for (k in 1:nknots){b[k]~dnorm(0,taub)}

#Prior distribution of the fixed effects parameters
for (l in 1:degree+1){beta[l]~dnorm(0,1.0E-6)}

#Prior distributions of the precision parameters
taueps~dgamma(1.0E-3,1.0E-3) 
taub~dgamma(1.0E-3,1.0E-3)

#Construct the design matrix of fixed effects
for (i in 1:n){
for (l in 1:degree+1){
X[i,l]<-pow(covariate[i],l-1)
                     }
              }

#Construct the design matrix of random effects
for (i in 1:n){
for (k in 1:nknots){
u[i,k]<-(covariate[i]-knot[k])*step(covariate[i]-knot[k])
Z[i,k]<-pow(u[i,k],degree)
                   }
              }

#Deterministic transformations. Obtain the standard deviations and
#the smoothing parameter
sigmaeps<-1/sqrt(taueps);sigmab<-1/sqrt(taub)
lambda<-pow(sigmab,2)/pow(sigmaeps,2)

#Predicting new observations
for (i in 1:n){
epsilonstar[i]~dnorm(0,taueps)
ystar[i]<-m[i]+epsilonstar[i]
              }
}"

params<-c("beta","ystar","b")

inits0  <- function () {
  list(beta  = rnorm(1, 0, 0.01),
       b = rnorm(1, 0, 0.01)
       
  )  }

inits1=inits0()
inits2=inits0()
inits3=inits0()

inits = list(inits1,inits2,inits3)
bugs(data = jags.data,inits,params,model,n.chains = 3,n.iter=5000)

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