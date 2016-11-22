#rm(list=ls(all=TRUE))
library(caret);library(visreg);
library(mgcv);library(ggplot2);
library(corrplot);library(reshape);
require(ggthemes);library(e1071);
library(scales);library(MASS);library(Hmisc)
library(corrplot);library(gridExtra)

# Read data
Data=read.csv("..//data/FiBY.csv",header=T)

#subsampling 
#index <- sample(seq_len(nrow(Data)),replace=F, size = 40000)

# Cut in redshift
data.1 = Data[Data$redshift < 30,]

# Log modulus transformation
L_M <-function(x){sign(x)*log10(abs(x) + 1)}

## fEsc is the variable of interest
data.2 <- as.data.frame(data.1[,c("Mvir","QHI")])
data.2$Mvir <- log10(data.2$Mvir)
data.2$QHI <- log10(data.2$QHI)

#data.2$spin <- log10(data.2$spin)
X    <- as.matrix(data.2)
y    <- data.1$fEsc; 
y[ y < 10^-3] = 0 # transform to zero everything below 1e-2
#y_b <- (y*(nrow(data.2)-1)+0.5)/nrow(data.2)



# Transform the columns of the design matrix (X)  to lighten the skewness, reduce the effect of outliers and reduce pairwise correlations if possible
#trans       <- preProcess(X,method = c("YeoJohnson", "center", "scale","spatialSign")) # Yeo-Johnson followed by centering and scaling. "spatialSign" is a bit complicated but it seems useful here
# to reduce outliers. Also I noticed that many of the variables are highly concetrated at one point. "spatailSign" will

trans       <- preProcess(X,method = c("YeoJohnson", "center", "scale","spatialSign"))                                                                                       # distribute things around.  
Xtrans      <- predict(trans,X) # The "new" Transformed X                 



### Two Models: 1) Model the probability that y > 0.  2) Model the Average of Y if y > 0 
### Some graphics:
### Plot the data (transformed) with a smoother:
### Good practice to go through all predictors (notice non-linearity):
non.zero <- ifelse(y > 0, 1, 0);fesc = as.factor(non.zero);levels(fesc) = c(expression("fesc = 0"), expression(fesc>0))
Dat.trans        <- data.frame(Xtrans,y=y,non.zero)



####
#### Modelling using non-parametric hurdle regression model:
#### Non-parametric regression function we will adopt more general regression functions without restricting to specific functional assumption 
#### 1) Model Prob(y>0) using nonparametric logistic regression
r                  <- 30
npar_nzero         <- gam(non.zero ~  s(Mvir,bs="cr",k=r)  + s(QHI,bs="cr",k=r),
                          data=Dat.trans,family=binomial(link="logit"),gamma=1.4)  # Please use bam with cautious as it might diverge. For example

######################################################################################
#### Produce the previous plots on the original scale of the predictors
### loop over all 
#names for plot 
names <- c("log~(M[200]/M[sun])",    
           "QHI")    

### 2) Model Average y when y > 0 using non parametric beta regression model 
r <- 30
npar_Beta_y <- gam(y ~     s(Mvir,bs="cr",k=r)       + 
                     s(QHI,bs="cr",k=r),
                   subset=y>0,data=Dat.trans,family=betar(link="logit"),gamma=1.4)


x1 <- Xtrans[,1]
x2 <- Xtrans[,2]
grid.lines = 50                # predict values on regular xy grid
x1.pred <- seq(1.01*min(x1), 0.99*max(x1), length.out = grid.lines)
x2.pred <- seq(1.01*min(x2), 0.99*max(x2), length.out = grid.lines)
x1x2 <- expand.grid(Mvir = x1.pred, QHI = x2.pred)

y.pred <- matrix(predict(npar_nzero, newdata = x1x2,type="response")*predict(npar_Beta_y, newdata = x1x2,type="response"), 
                 nrow = grid.lines, ncol = grid.lines)


x1_o <- X[,1]
x2_o <- X[,2]
grid.lines = 50                # predict values on regular xy grid
x1.pred_o <- seq(1.01*min(x1_o), 0.99*max(x1_o), length.out = grid.lines)
x2.pred_o <- seq(1.01*min(x2_o), 0.99*max(x2_o), length.out = grid.lines)
x1x2_o <- expand.grid(Mvir = x1.pred_o, QHI = x2.pred_o)



fit = predict(npar_nzero,newdata=as.data.frame(Xtrans),type="response")*predict(npar_Beta_y,newdata=as.data.frame(Xtrans),type="response")
cor(fit,y)


require("plot3D")


scatter3D_fancy(x, y, z,colvar = as.integer(CLUST$classification),col = c("#D46A6A","#D4B16A","#764B8E"),
                colkey=F,
                box = T,ticktype = "detailed",theta=40,phi=20,
                zlab = "LogOIII_Hb",ylab="LogNII_Ha", d=20,
                xlab="EWHa",bty = "u",col.panel = "gray95",col.grid = "gray35",contour = T)

# scatter plot with regression plane
pdf("3D.pdf",width = 10,height = 10)
scatter3D(x1_o, x2_o, y, pch = 19,                    
          cex = 0, cex.lab=1.5,  
          box = T,
          theta = 70, phi = 25, ticktype = "detailed",
          col="red2",bty = "u",col.panel = "gray95",
          col.grid = "gray35",
          xlab="log Mvir",
          ylab="log QHI",
          zlab="f_esc", 
          zlim=c(0,0.2),
          surf = list(col="red3",x = x1.pred_o, y = x2.pred_o, z = y.pred,  
                      facets = NA,lwd=1.5,lty=2),colkey = FALSE)
dev.off()


#### To monitor the changes in the average of y as one variables is changing keeping all other predictors fixed
#### These are important plots as one can notice the importance of each predictor to y
gg3<-list()
vv_x <- list()
vv_original <-list()
for(i in 1:ncol(X)){
  nn = 10^4+1; XX = matrix(apply(X,2,mean),nrow=1); 
  XX=XX%x% rep(1,nn);colnames(XX) = colnames(X); 
  XX = as.data.frame(XX)
  XX[,i] <- seq(min(X[,i]),max(X[,i]),length=nn)
  XX.trans      = predict(trans,XX) # Convert to the transformed scale
  ##
  Preds_nzero <- predict(npar_nzero ,newdata = XX.trans,type="response",se=T,unconditional=T)
  Preds_y     <- predict(npar_Beta_y,newdata = XX.trans,type="response",se=T,unconditional=T)
  mu1         <- Preds_nzero$fit;mu2 <- Preds_y$fit;
  se1         <- Preds_nzero$se ;se2 <- Preds_y$se;
  Predictions <- mu1*mu2
  se          <- sqrt(mu2^2*se1^2 + mu1^2*se2^2 + se1^2*se2^2)
  CI.L        <- Predictions-qnorm(0.975)*se 
  CI.R        <- Predictions+qnorm(0.975)*se 
  CI          <- cbind(Predictions,CI.L,CI.R) 
  colnames(CI) <- c("Predictions","CI_L","CI_R")
  
  vv_x[[i]] <- data.frame(cbind(CI,x=XX[,i]),var = rep(names[i],length(XX[,i])))
  vv_original[[i]] <- data.frame(x=data.2[,i],y=y,var = rep(names[i],length(data.2[,i])))
  
}

# put altogether for facets
vvv_x<-c()
for(i in 1:8){
  vvv_x <- rbind(vvv_x,vv_x[[i]])
}  
vvv_original <- c()
for(i in 1:8){
  vvv_original <- rbind(vvv_original,vv_original[[i]])
} 

# Plot  via ggplot2
pdf("hurdle_3.pdf",width = 16,height = 8)
ggplot(vvv_x,aes(x=x,y=Predictions))+
  # geom_point(data=gg_original,aes(x=x,y=y),size=1,alpha=0.2,col="orange2",position = position_jitter (h = 0.025))+
  geom_hex(data=vvv_original,bins = 50,aes(x=x,y=y))+
  scale_fill_continuous(low = "gray80", high = "gray10", trans = log10_trans())+
  geom_ribbon(aes(ymin=CI_L, ymax=CI_R),fill=c("#ba122b")) +
  geom_line(col="#CCCC99",size=1.5)+
  theme_bw()+
  ylab(expression(f[esc]))+
  xlab("")+
  #    scale_y_continuous(trans = 'log10',
  #                       breaks=trans_breaks("log10",function(x) 10^x),
  #                       labels=trans_format("log10",math_format(10^.x)))+
  theme(legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white",color = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(vjust = 0.1,margin=margin(0,10,0,0)),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 20,family="serif"))+
  facet_wrap(~var,scales = "free_x",ncol=4,labeller = label_parsed)
dev.off()




#### Finally: 
#### Prediction :
#### Say at xp = X[c(100,363,987),] 
xp= X[c(100,363,987),]; xp=as.data.frame(xp);colnames(xp) = colnames(X)
xp.trans <- predict(trans,xp) # The "new" Transformed xp      
# Find pr(y > 0) at xp 
predict(npar_nzero,newdata=xp.trans,type="response")
# Predict the average of y at xp (Important) : 
predict(npar_nzero,newdata=xp.trans,type="response")*predict(npar_Beta_y,newdata=xp.trans,type="response")

## Fitted values:
fit = predict(npar_nzero,newdata=as.data.frame(Xtrans),type="response",se=T)$fit*predict(npar_Beta_y,newdata=as.data.frame(Xtrans),type="response",se=T)$fit




cor(fit,y)
plot(fit,y,pch=20,cex=0.5)
abline(0,1,lty=2,col="Black",lwd=3)







