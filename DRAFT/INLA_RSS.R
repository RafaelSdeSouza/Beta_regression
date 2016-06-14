require(INLA)
# Read data
dat1<-read.csv("..//data/N_body.csv",header=T)
dat1$bin<-dat1$SFR
dat1$bin[dat1$bin > 0]<-1
dat2<-dat1[dat1$redshift < 10,]


################################
##INLA
################################
##Generate binary data
y <- log(dat2$Xmol,10)
x <- log(1e10*dat2$M_dm,10)
x0<-as.data.frame(dat2$x)
y0<-as.data.frame(dat2$y)
xy<-data.frame(x0,y0)
coords <-as.matrix(xy)


mesh2 <- inla.mesh.2d(loc=coords,max.edge=c(20, 40),cutoff=0.1)
plot(mesh2)
points(coords, pch=19, bg=1, col="blue", cex=0.7)


#Define the n × m projector matrix to project the process at the mesh nodes to locations
A.est <- inla.spde.make.A(mesh=mesh2, loc=coords)
dim(A.est)

#Build the SPDE model on the mesh. Exponential correlation function, α = 2
spde <- inla.spde2.matern(mesh=mesh2, alpha=2)



#Create a stack data 
stk.e <- inla.stack(tag='pred',  data=list(y=y), ## response
                    A=list(A.est, 1), ## two projector matrix 
                    effects=list(## two elements:
                    s=1:spde$n.spde, ## RF index 
                    data.frame(b0=1,x=x)))
                    
                    
# Fit 

formula <- y ~ 0 + b0 + x +f(s, model=spde)

res <- inla(formula, 
                data = inla.stack.data(stk.e),
                control.predictor=list(A=inla.stack.A(stk.e),compute=TRUE))


round(res$summary.fixed, 4)

post.s2e <-inla.tmarginal( function(x) 1/x, res$marginals.hyperpar$'Theta1 for s')

plot(post.s2e, type='l', ylab='Density',
     xlab=expression('parameter'))

rf <- inla.spde.result(
  inla=res, ## the inla() output
  name='s', ## name of RF index set
  spde=spde, ## SPDE model object
  do.transf=TRUE) 


plot(rf$marginals.var[[1]], ty = "l", xlab = expression('parameter1'), yla = "Density")

plot(rf$marginals.kap[[1]], type = "l", xlab = expression('parameter2'), ylab = "Density")

plot(rf$marginals.range[[1]], type = "l", xlab = "range nominal", ylab = "Density")

gproj <- inla.mesh.projector(mesh2, xlim = min(x0):max(x0), ylim = min(y0):max(y0), dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, res$summary.random$s$mean)
g.sd <- inla.mesh.project(gproj, res$summary.random$s$sd)



library(lattice)
library(gridExtra)
trellis.par.set(regions=list(col=terrain.colors(16)))
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean'),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd'), nrow=1)

#tcoo <- rbind(c(0.3, 0.3), c(0.5, 0.5), c(0.7, 0.7))
#dim(Ap <- inla.spde.make.A(mesh = mesh2, loc = tcoo))

#x0 <- c(0.5, 0.5, 0.5)

stk.pred <- inla.stack(tag='pred', A=list(Ap, 1), data=list(y=NA), effects=list(s=1:spde$n.spde, data.frame(x=x0, b0=1)))

stk.full <- inla.stack(stk.e, stk.pred)
p.res <- inla(formula, data=inla.stack.data(stk.full), ## full stack
              control.predictor=list(compute=TRUE, ## compute the predictor
                                     A=inla.stack.A(stk.full)))

pred.ind <- inla.stack.index(stk.e, tag = "pred")$data
ypost <- res$marginals.fitted.values[pred.ind]

names(ypost) <- paste('y', seq_along(ypost), sep='_')

library(plyr)

xyplot(y~x | .id, ldply(ypost), panel='llines', xlab='y', ylab='Density')
