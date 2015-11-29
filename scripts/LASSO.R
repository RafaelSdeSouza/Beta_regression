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
data.1= read.table(file="..//data/FiBY_escape_data_all.dat",header=FALSE)
colnames(data.1)<-c("redshift","fEsc","Mvir","Mstar","Mgas","QHI","sfr_gas",
                    "sfr_stars","ssfr_gas","ssfr_stars","baryon_fraction",
                    "spin","age_star_mean","age_star_max","age_star_min","NH_10","clumping_factor")

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



# Labels
lbs_fun <- function(fit, ...) {
  L <- length(fit$lambda)
  x <- log(fit$lambda[L])+0.25
  y <- fit$beta[, L]
  labs <- names(y)
  text(x, y,pos = 1, offset = 0.5, labels=labs, ...)
}
plot(fit,xvar="lambda",label = TRUE)

lbs_fun(fit)

plot(fit, xvar = "dev", label = TRUE)
cv.glmmod <- cv.glmnet(x2,y=data.3$Y,alpha=1,family="binomial",type.measure = "auc")

# Plot AUC vs Lambda
CairoPDF("AUCvsLambda.pdf")
plot(cv.glmmod,xlab=expression(log~lambda),family="serif")
dev.off()

best_lambda <- cv.glmmod$lambda.1se

coef.min = coef(cv.glmmod, s = "lambda.1se")
coef.min

# Set only 2 most importants by hand 

coef.2 = coef(cv.glmmod, s = 0.0110300)
coef.2


# Fit GLM with selected features



fit=glm(Y~Mgas+ssfr_stars,data=data.3,family=binomial("logit"))


ROCtest(fit,10,"ROC")


# Plot 3D
# Plot 

x <-range(data.2$ssfr_stars)
x <- seq(x[1], 0.15*x[2], length.out=50)
y <- range(data.2$Mgas)
y <- seq(0.5*y[1], 1.25*y[2], length.out=50)

z <- outer(x,y,
           function(ssfr_stars,Mgas)
             predict(fit, data.frame(ssfr_stars,Mgas),type = 'response'))
library(rsm)
library(lattice)
YlOrBr <- c("#00A3DB")
#p<-persp(x,y,z, theta=150, phi=20,
#         expand = 0.5,shade = 0.1,
#         xlab="Z", ylab=expression(NII.Ha), zlab=expression(log10.EW.Ha),ticktype='detailed',
#         col = YlOrBr,border=NA,xlog=T,ylog=T)
cairo_pdf("logit3D.pdf",width = 9,height = 8)
trellis.par.set("axis.line",list(axis.text=list(cex=20),col=NA,lty=1,lwd=2))
par(mar=c(1,1,1,1))
wireframe(z~x+y,data=data.frame(x=x, y=rep(y, each=length(x)), z=z),
          par.settings = list(regions=list(alpha=0.6)),
          col.regions =YlOrBr,drape=T,light.source = c(5,5,5),colorkey = FALSE,
          ylab=list(label=expression(log~M[gas]/M['\u0298']),cex=1.25),
          xlab=list(label=expression(sSFR[stars]/Gyr),cex=1.25),
          zlab=list(rot=90,label=expression(P[f[scape]]>=0.1),cex=1.25,dist=-1,rot=0),
          scale=list(tck=0.75,arrows=FALSE,distance =c(0.75, 0.75, 0.75)))

dev.off()










fit2=gam(Y~te(Mgas,ssfr_stars),data=data.3,family=binomial("logit"))
plot(fit2)
vis.gam(fit2,type="response",plot.type = "persp",color="topo", border=NA, n.grid=500,theta=-60,phi=30)
ROCtest(fit2,10,"ROC")


formula1 = Y~QHI+baryon_fraction+f(redshift,model="ar1")
mod.1 = inla(formula1,data=data.2,family="binomial")

formula2 = Y~Mgas+ssfr_stars+f(redshift,model="ar1")
mod.2 = inla(formula2,data=data.3,family="binomial")


plot(mod.2)
hyper = inla.hyperpar(mod.2)
summary(hyper)
plot(hyper)


mydata <- as.data.frame(mod.2$marginals.fixed)

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
