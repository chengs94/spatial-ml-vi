# PLS + TPRS covariance

library(SpatioTemporal)
library(matrixStats)
library(pls)

source("spatTreeModified.R")
source("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/spatialRF/aux_functions.R")

dat = readRDS("dat.rda")
set.seed(2333)

n <- nrow(dat)
Y.pls <- Y.spatrf.cv <- Y.spatrf.pl <- rep(NA, n)

numfold <- 5
grid <- list(lambert_x = dat$lambert_x,
             lambert_y = dat$lambert_y)
smpl <- sample(n)
ind.stop <- round(n/numfold*1:numfold)
ind.strt <- c(1,ind.stop[1:numfold-1]+1)

X = dat[,-(1:21)]

logLikeTPRS <- function(theta,y,train.coords,x){
  
  #spatBas <- .makeSpatBas(pars=list(psill=theta,nugget=1-theta),
  #                        coords.train=train.coords,
  #                        coords.test = test.coords, out.Z = TRUE,
  #                        cov.type="TPRS",cov.opts=list(k=0,m=2))
  spatBas <- .makeSpatBas(pars=list(psill=1,nugget=1),
                          coords.train=train.coords, out.Z = TRUE,
                          cov.type="TPRS",cov.opts=list(k=0,m=2))
  R <- spatBas$Zb %*% t(spatBas$Zb)
  n = dim(R)[1]
  V <- chol(theta * R +  diag(rep(1-theta, n)))
  i.V = chol2inv(V)
  Om <- i.V - i.V %*% x %*% .fastSolve( t(x) %*% i.V %*% x) %*% t(x) %*% i.V
  
  cur.lklhd = 2 * sum(log(diag(V))) + n*log( t(y) %*% Om %*% y/ n)
  
  #cur.lklhd = n*log( t(y) %*% Om %*% y/ n) - 2 * sum ( log ( diag (i.Sig) ) )
  
  cur.lklhd
  
  ##n <- length(y)
  ##obj <- .detInv(c(sll,nggt),exp(-dist.mat/range))
  #obj = .makeSpatBas(pars=list(psill=sll,
  #                             nugget=nggt),
  #                   coords.train=train.coords,
  #                   coords.test = NULL,
  #                   cov.type="TPRS",cov.opts = list(k=0,m=2))
  
  ##Om <- obj$i.Sig
  #Om <- obj$i.sig
  #if(!is.null(x)){
  #  p <- dim(x)[2]
  #  Om <- Om - Om %*% x %*% .fastSolve( t(x) %*% Om %*% x) %*% t(x) %*% Om
  #}
  
  #detInvCov = base::det(Om)
  ##if(length(theta)==1){
  ##  out <- n*log( t(y) %*% Om %*% y/ n) + obj$det.Sig
  ##} else{
  ##out <- obj$det.Sig + t(y) %*% Om %*% y
  #out <- -detInvCov + t(y) %*% Om %*% y
  ##}
  
  #out
  
}

run = function(fold, annavg){
  cv.ind.i <- c(smpl[(ind.strt[fold]:ind.stop[fold])])
  n.test <- length(cv.ind.i)
  n.train <- n-n.test
  Y.train <- annavg[-cv.ind.i]
  Y.test <- annavg[cv.ind.i]
  
  ## Standardize Coordinates
  train.coords <- cbind(grid$lambert_x[-cv.ind.i],grid$lambert_y[-cv.ind.i])
  scale.params <- cbind(apply(train.coords,2,mean),apply(train.coords,2,sd))
  x1train <- train.coords[,1] <- (train.coords[,1] - scale.params[1,1])/scale.params[1,2]
  x2train <-train.coords[,2] <- (train.coords[,2] - scale.params[2,1])/scale.params[2,2]
  test.coords <- cbind(grid$lambert_x[cv.ind.i],grid$lambert_y[cv.ind.i])
  x1test <- test.coords[,1] <- (test.coords[,1] - scale.params[1,1])/scale.params[1,2]
  x2test <- test.coords[,2] <- (test.coords[,2] - scale.params[2,1])/scale.params[2,2]
  
  X.train = as.matrix(X[-cv.ind.i,])
  X.test = as.matrix(X[cv.ind.i,])
  
  dmat <- as.matrix(SpatioTemporal::crossDist(train.coords))
  dtestmat <- as.matrix(SpatioTemporal::crossDist(test.coords,train.coords))
  
  # UK - PLS
  strt <- proc.time()
  
  X.train.pc = X.train
  X.test.pc = X.test

  num.pcs = 3
  pc.obj <- get.pcs(X.train.pc,"pls",num.pcs,Y.train)
  #theta = optim(1,logLikeTPRS,method="L-BFGS-B",
  #              lower=-30, upper =15,
  #              y=Y.train, train.coords=train.coords, x=pc.obj$X)$par
  new.x.train <- pc.obj$X
  means<-apply(X.train.pc,2,mean)
  mean.matrix<-matrix(rep(means,n.test),nrow=n.test,byrow=T)
  sds <- apply(X.train.pc,2,sd)
  new.x.test <- (X.test.pc - mean.matrix) %*% diag(1/sds) %*% pc.obj$proj
  #new.x.train = as.data.frame(new.x.train[])
  #new.x.test = as.data.frame(new.x.test)
  
  thetas = seq(.05, .95, .05)
  lklhds = sapply(thetas, logLikeTPRS, 
                  y=Y.train, train.coords=train.coords, 
                  x=new.x.train)
  
  theta = thetas[which.min(lklhds)]
  spatBas <- .makeSpatBas(pars=list(psill=1,nugget=1),
                          coords.train=train.coords, out.Z = TRUE,
                          cov.type="TPRS",cov.opts=list(k=0,m=2))
  R <- spatBas$Zb %*% t(spatBas$Zb)
  n = dim(R)[1]
  V <- chol(theta * R +  diag(rep(1-theta, n)))
  i.V = chol2inv(V)
  
  fullspatbas <- .makeSpatBas(pars=list(psill=theta,
                                        nugget=1-theta),
                              coords.train=train.coords,
                              coords.test = test.coords,
                              cov.type="TPRS",cov.opts = list(k=0,m=2))
  beta.hat.pcr <- solve(t(new.x.train) %*% i.V %*% new.x.train,
                        t(new.x.train) %*% i.V %*% Y.train)
  beta.hat <- solve(t(train.coords) %*% fullspatbas$i.sig %*% train.coords,
                    t(train.coords) %*% fullspatbas$i.sig %*% 
                      (Y.train - new.x.train%*%beta.hat.pcr ) )
  Y.pls[cv.ind.i] <- predict.cor( new.x.test%*%beta.hat.pcr + 
                                    test.coords%*% beta.hat, 
                                  Y.train - train.coords %*% beta.hat - 
                                    new.x.train%*%beta.hat.pcr, 
                                  theta = theta, 
                                  R.test=fullspatbas$sig.test, 
                                  i.sig = fullspatbas$i.sig )
  
  #beta.hat.pcr = lklhds[[which(thetas == theta)]]$beta.hat.pcr
  #Y.pls[cv.ind.i] = lklhds[[which(thetas == theta)]]$Y.pred
  
  pls.time <- proc.time()-strt
  
  return(list(# pls-uk
    Y.pls = Y.pls, beta.hat.pls = beta.hat.pcr,
    pc.obj = pc.obj, theta.pls = theta, pls.time = pls.time))
}

plsTPRS.ufp.uw <- lapply(1:numfold, FUN = run, annavg = dat$ufp_uw)
plsTPRS.bc.uw <- lapply(1:numfold, FUN = run, annavg = dat$bc_uw)
plsTPRS.bc.primary <- lapply(1:numfold, FUN = run, annavg = dat$bc_primary)
plsTPRS.ufp.primary <- lapply(1:numfold, FUN = run, annavg = dat$ufp_primary)