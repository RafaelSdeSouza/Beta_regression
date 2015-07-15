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
require(caret)
#Read the  dataset


data.1= read.table(file="FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10")


trainIndex <- createDataPartition(data.1$redshift, p = .1,
                                  list = FALSE,
                                  times = 1)
#data.2<-data.1[data.1$redshift==8.86815,]
data.2<-data.1[trainIndex,]
N<-nrow(data.2)


data.2$Y<-(data.2$fEsc*(N-1)+0.5)/N


# Prepare data for JAGS
data.2$Mstar<-(data.2$Mstar-mean(data.2$Mstar))/sd(data.2$Mstar)
data.2$Mgas<-(data.2$Mgas-mean(data.2$Mgas))/sd(data.2$Mgas)
data.2$Mvir<-(data.2$Mvir-mean(data.2$Mvir))/sd(data.2$Mvir)
data.2$sfr_gas<-(data.2$sfr_gas-mean(data.2$sfr_gas))/sd(data.2$sfr_gas)
#data.2$baryon_fraction<-(data.2$baryon_fraction-mean(data.2$baryon_fraction))/sd(data.2$baryon_fraction)
data.2$QHI<-(data.2$QHI-mean(data.2$QHI))/sd(data.2$QHI)
data.2$ssfr_gas<-(data.2$ssfr_gas-mean(data.2$ssfr_gas))/sd(data.2$ssfr_gas)
data.2$age_star_mean<-(data.2$age_star_mean-mean(data.2$age_star_mean))/sd(data.2$age_star_mean)
data.2$spin<-(data.2$spin-mean(data.2$spin))/sd(data.2$spin)
data.2$NH_10<-(data.2$NH_10-mean(data.2$NH_10))/sd(data.2$NH_10)

# Random effects with redshift

re<-as.numeric(as.factor(data.2$redshift))
Nredshift<-length(unique(data.2$redshift))



#inla.m<-inla(Y~baryon_fraction+f(redshift,model="ar1"),data=data.2,family="beta",control.predictor = list(compute = TRUE),control.compute = list(config=TRUE))
#inla.m$summary.fitted.values$mean


#inla.m$summary.fitted.values$mean
# GAM 
library(splines)
df<-3
X.bs <- bs(data.2$baryon_fraction, df = df,
           intercept = FALSE)

X<-model.matrix(~baryon_fraction,data=data.2)
# Scale

K<-ncol(X)


jags.data <- list(Y= data.2$Y,
                  N = nrow(data.2),
                  X=X,
                  re=re,
                  X.bs=X.bs,
                  b0 = rep(0,K),
                  B0=diag(1e-4,K),
                  a0=rep(0,df),
                  A0=diag(1e-4,df),
#                  c0=rep(0,Nredshift),
#                  C0=diag(1,Nredshift),
                  Npred = K,
                  Nredshift = Nredshift
)


model<-"model{
#1. Priors 
beta~dmnorm(b0[],B0[,]) # Normal Priors
b~dmnorm(a0[],A0[,]) # Normal Priors
#c~dmnorm(c0[],tau*C0[,])
#tau<-1/(sigma*sigma)
sigma~dgamma(0.01,0.01)

tau.R<-pow(sdBeta,-1)
sdBeta ~ dgamma(0.01,0.01)


# Jefreys priors for sparseness 
#for(j in 1:Npred)   {
#lnTau[j] ~ dunif(-50, 50)   
#TauM[j] <- exp(lnTau[j])
#beta[j] ~ dnorm(0, TauM[j]) 
#Ind[j] <- step(abs(beta[j]) - 0.05)
#}

theta~dgamma(0.01,0.01)

# AR(1) model
for (j in 1:Nredshift){

ranef[j]~ddexp(0,tau.R)
}

#2. Likelihood


for(i in 1:N){

Y[i] ~ dbeta(shape1[i],shape2[i])
shape1[i]<-theta*pi[i]
shape2[i]<-theta*(1-pi[i])
logit(pi[i]) <- eta[i]
eta[i]<-inprod(beta[],X[i,])+inprod(b[],X.bs[i,])+ranef[re[i]]
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

params <- c("beta","theta","PRes","Fit","newFit","newY","Ind","c","b")

inits0  <- function () {
  list(beta  = rnorm(K, 0, 0.01)  #Regression parameters
       
  )  }

inits1=inits0()
inits2=inits0()
inits3=inits0()

library(parallel)
cl <- makeCluster(6)
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
                       plots=TRUE
)

jagssamples <- as.mcmc.list(jags.logit)
summary<-extend.jags(jags.logit,drop.monitor=c("PRes","Fit","newFit","newY"), summarise=TRUE)
print(summary)

L.factors <- data.frame(
  Parameter=paste("beta[", seq(1:5), "]", sep=""),
  Label=c("(Intercept)","Mstar","Mgas","sfr_gas","baryon_fraction"))
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

Dispersion = sum(Pres^2)/(N-7)# beta.0, beta.1 and k, 3 parameters




# Diagnostics Confusion Matrix (very unlikely to be good)
predscape<-summary(as.mcmc.list(jags.logit, vars="newY"),quantiles=c(0.005,0.025,0.25,0.5,0.75,0.975, 0.995))
predscape<-predscape$quantiles


plot(data.2$fEsc,predscape[,3])
plot(data.2$baryon_fraction,predscape[,5])
plot(data.2$baryon_fraction,data.2$fEsc)

pred2<-data.frame(f_gas=data.2$baryon_fraction,mean=predscape[,4],lwr1=predscape[,1],lwr2=predscape[,2],
                  lwr3=predscape[,3],upr1=predscape[,5],upr2=predscape[,6],upr3=predscape[,7])

CairoPDF("f_scape.pdf",height=8,width=9)
ggplot(data.2,aes(x=baryon_fraction,y=fEsc))+
  geom_ribbon(data=pred2,aes(x=f_gas,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill="gray") +
  geom_ribbon(data=pred2,aes(x=f_gas,y=mean,ymin=lwr2, ymax=upr2), alpha=0.35, fill="gray") +
  geom_ribbon(data=pred2,aes(x=f_gas,y=mean,ymin=lwr3, ymax=upr3), alpha=0.25, fill="gray") +
  geom_point(size=0.75,alpha=0.7)+
  geom_line(data=pred2,aes(x=f_gas,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_colour_gdocs()+
  theme_hc()+
  ylab(expression(f[scape]))+
  xlab(expression(f[gas]))+theme(axis.title.y=element_text(vjust=0.75),
                                                axis.title.x=element_text(vjust=-0.25),
                                                text = element_text(size=25))
dev.off()


inla.m$summary.fitted.values$`0.975quant`
pred2<-data.frame(f_gas=data.2$baryon_fraction,mean=inla.m$summary.fitted.values$mean,lwr1=inla.m$summary.fitted.values$`0.025quant`,upr1=inla.m$summary.fitted.values$`0.975quant`)

CairoPDF("f_scape_all.pdf",height=8,width=9)
ggplot(data.2,aes(x=baryon_fraction,y=fEsc))+
  geom_ribbon(data=pred2,aes(x=f_gas,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill="gray") +
  geom_point(size=0.75,alpha=0.7)+
  geom_line(data=pred2,aes(x=f_gas,y=mean),colour="gray25",linetype="dashed",size=1.2)+
  scale_colour_gdocs()+
  theme_hc()+
  ylab(expression(f[scape]))+
  xlab(expression(f[gas]))+theme(axis.title.y=element_text(vjust=0.75),
                                 axis.title.x=element_text(vjust=-0.25),
                                 text = element_text(size=25))
dev.off()
