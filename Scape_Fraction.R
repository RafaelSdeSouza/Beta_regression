#  JAGS script  Scape_Fraction.R
#  Copyright (C) 2015  Rafael S. de Souza
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
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
#Read the  dataset

data.1= read.table(file="FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10")

data.2<-data.1[data.1$redshift==8.86815,]
N<-nrow(data.2)


data.2$Y<-(data.2$fEsc*(N-1)+0.5)/N


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



X<-model.matrix(~Mvir+baryon_fraction+age_star_mean+ssfr_gas+NH_10,data=data.2)
# Scale

K<-ncol(X)


jags.data <- list(Y= data.2$Y,
                  N = nrow(data.2),
                  X=X,
                  b0 = rep(0,K),
                  B0=diag(1e-4,K),
                  Npred = K
)


model<-"model{
#1. Priors 
#beta~dmnorm(b0[],B0[,]) # Normal Priors
# Jefreys priors for sparseness 
for(j in 1:Npred)   {
      lnTau[j] ~ dunif(-50, 50)   
      TauM[j] <- exp(lnTau[j])
      beta[j] ~ dnorm(0, TauM[j]) 
      Ind[j] <- step(abs(beta[j]) - 0.05)
}

theta~dgamma(0.01,0.01)
#2. Likelihood

for(i in 1:N){

Y[i] ~ dbeta(shape1[i],shape2[i])
shape1[i]<-theta*pi[i]
shape2[i]<-theta*(1-pi[i])
logit(pi[i]) <- eta[i]
eta[i]<-inprod(beta[],X[i,])
ExpY[i]<-pi[i]
VarY[i]<-pi[i]*(1-pi[i])/(theta+1)
PRes[i]<-(Y[i]-ExpY[i])/sqrt(VarY[i])


#3. Discrepancy measures

newY[i]~dbeta(shape1[i],shape2[i])
PResNew[i]<-(newY[i]-ExpY[i])/sqrt(VarY[i])
D[i]<-pow(PRes[i],2)
DNew[i]<-pow(PResNew[i],2)

}
Fit<-sum(D[1:N])
newFit<-sum(DNew[1:N])
}"

params <- c("beta","theta","PRes","Fit","newFit","newY","Ind")

inits0  <- function () {
  list(beta  = rnorm(K, 0, 0.01)  #Regression parameters
      
  )  }

inits1=inits0()
inits2=inits0()
inits3=inits0()

library(parallel)
cl <- makeCluster(3)
jags.logit <- run.jags(method="rjparallel", 
                       data = jags.data, 
                       inits = list(inits1,inits2,inits3),
                       model=model,
                       n.chains = 3,
                       adapt=1000,
                       monitor=c(params),
                       burnin=1000,
                       sample=5000,
                       summarise=FALSE,
                       plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)
summary<-extend.jags(jags.logit,drop.monitor=c("PRes","Fit","newFit","newY"), summarise=TRUE)
print(summary)

L.factors <- data.frame(
  Parameter=paste("beta[", seq(1:6), "]", sep=""),
  Label=c("(Intercept)","Mvir","baryon_fraction","age_star_mean","ssfr_gas","NH_10"))
beta_post<-ggs(jagssamples,par_labels=L.factors,family=c("beta"))

pdf("beta_scape.pdf")
ggs_caterpillar(beta_post)+theme_hc()+ylab("")+geom_vline(xintercept=c(0,0), linetype="dotted")+
theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))
dev.off()



require(scales)
Pres<-summary(as.mcmc.list(jags.logit, vars="PRes"),quantiles=0.5)$quantiles
plot(data.2$baryon_fraction,Pres)

Dispersion = sum(Pres^2)/(N-6)# beta.0, beta.1 and k, 3 parameters




# Diagnostics Confusion Matrix (very unlikely to be good)
predscape<-summary(as.mcmc.list(jags.logit, vars="newY"))
predscape<-predscape$quantiles


plot(data.2$fEsc,predscape[,3])
plot(data.2$baryon_fraction,predscape[,3])
plot(log(data.2$baryon_fraction,10),data.2$fEsc)



