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
require(R2jags)
library(ggmcmc)
library(ggplot2)
library(ggthemes)
#Read the  dataset

data.1= read.table(file="..//data/FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10")

data.2<-data.1[data.1$redshift<=10,]
N<-nrow(data.2)


data.2$Y[data.2$fEsc>=0.1]<-1
data.2$Y[data.2$fEsc<0.1]<-0


# Prepare data for JAGS
x.scale<-as.data.frame(scale(data.2[,c("Mstar","Mgas","Mvir","sfr_gas","baryon_fraction","ssfr_gas","age_star_mean",
                        "spin","NH_10","QHI")]))

X<-model.matrix(~baryon_fraction+QHI,data=x.scale)

# Scale
K<-ncol(X)


jags.data <- list(Y= data.2$Y,
                  N = nrow(data.2),
                  X=X,
                  b0 = rep(0,K),
                  B0=diag(1e-4,K),
                  Npred = K
)


LOGIT<-"model{
#1. Priors 
beta~dmnorm(b0[],B0[,]) # Normal Priors

#2. Likelihood

for(i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- inprod(beta[], X[i,])


#3. Prediction
NewPred[i]~dbern(pi[i])
}

}"

params <- c("beta","pi","NewPred")

inits  <- function () {
  list(beta = rnorm(K, 0, 0.1))}

jags.logit <- jags( 
                    data = jags.data, 
                    inits      = inits,
                    parameters = params,
                    model      = textConnection(LOGIT),
                    n.chains   = 3,
                    n.iter     = 500,
                    n.thin     = 1,
                    n.burnin   = 250)



print(jags.logit,justify = "left", digits=2)
