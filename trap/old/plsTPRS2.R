# PLS + TPRS covariance

library(SpatioTemporal)
library(matrixStats)
library(pls)
library(glmnet)

source("spatTreeModified.R")
source("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/spatialRF/aux_functions.R")

dat = readRDS("dat.rda")
set.seed(2333)

n <- nrow(dat)
Y.pls <- Y.pls.mean <- Y.spatrf.cv <- Y.spatrf.pl <- rep(NA, n)

numfold <- 5
grid <- list(lambert_x = dat$lambert_x,
             lambert_y = dat$lambert_y)
smpl <- sample(n)
ind.stop <- round(n/numfold*1:numfold)
ind.strt <- c(1,ind.stop[1:numfold-1]+1)

X = dat[,-(1:21)]

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
  #pc.pars <- matrix(NA,ncol=length(num.pcs),nrow=3)
  pc.obj <- get.pcs(X.train.pc,"pls",num.pcs,Y.train)
  strt.vl <- optim(c(1,-1,1),logLikeExp,method="L-BFGS-B",
                   lower=c(-30,-20,-10), upper =c(15,15,30),
                   y=Y.train, dist.mat=dmat, x=pc.obj$X)
  #theta = optim(1,logLikeTPRS,method="L-BFGS-B",
  #              lower=-30, upper =15,
  #              y=Y.train, train.coords=train.coords, x=pc.obj$X)$par
  new.x.train <- pc.obj$X
  means<-apply(X.train.pc,2,mean)
  mean.matrix<-matrix(rep(means,n.test),nrow=n.test,byrow=T)
  sds <- apply(X.train.pc,2,sd)
  new.x.test <- (X.test.pc - mean.matrix) %*% diag(1/sds) %*% pc.obj$proj
  
  #beta.hat = solve(t(new.x.train)%*%new.x.train, t(new.x.train)%*%Y.train)
  #Y.pls.mean[cv.ind.i] = new.x.test %*% beta.hat
  
  pls.err = Inf
  pred.candidate = rep(NA, length(n))
  ncomps = seq(1,15)
  best.ncomps = NA
  mod = plsr(Y.train ~ X.train, ncomp = max(ncomps), scale = TRUE)
  for (j in ncomps){
    #new.x.trn = as.matrix(mod$scores)
    #proj = mod$projection
    #new.x.tst = (X.test.pc - mean.matrix) %*% diag(1/sds) %*% proj
    #mod.linear = lm(Y.train ~ new.x.trn)
    #bt.pls = solve(t(new.x.trn)%*%new.x.trn, t(new.x.trn)%*%Y.train)
    pred.candidate[cv.ind.i] = predict(mod, X.test, ncomp = j)
    #pred.candidate[cv.ind.i] = new.x.tst %*% bt.pls
    if (mean((pred.candidate[cv.ind.i] - Y.test)^2) < pls.err){
      pls.err = mean((pred.candidate[cv.ind.i] - Y.test)^2)
      Y.pls.mean[cv.ind.i] = pred.candidate[cv.ind.i]
      best.ncomps = j
    }
  }
  
  #theta.pls <- strt.vl$par
  
  #sig.pls <- exp(theta.pls[1])
  #nug.pls <- exp(theta.pls[2])
  #range.pls <- exp(theta.pls[3])
  
  #i.V <- .fastSolve(sig.pls * exp(-dmat/range.pls) + 
  #                    nug.pls * diag(length(Y.train)))
  #beta.hat.pcr <- solve(t(new.x.train) %*% i.V %*% new.x.train,
  #                      t(new.x.train) %*% i.V %*% Y.train)
  #Y.pls.mean[cv.ind.i]<- c(new.x.test %*% beta.hat.pcr) + 
  #  sig.pls * exp(-dtestmat/range.pls) %*% i.V %*% 
  #  c(Y.train-c(new.x.train %*% beta.hat.pcr))
  
  #######
  #Z <- makeSigmaTPRS(c(1),train.coords,out.ZB=TRUE,
  #                   pen.form=FALSE,fit.tprs=TRUE,K=n.train)
  #Z.test <- makeSigmaTPRS(1,train.coords,test.coords,out.ZB=TRUE,
  #                        pen.form=FALSE,fit.tprs = FALSE,K=n.train)
  #d = dim(Z)[2]
  #mod = cv.glmnet(x = cbind(new.x.train, Z), y = Y.train, alpha = 0,
  #                penalty.factor = c(rep(0, dim(new.x.train)[2]), rep(1, d)))
  #lambda = mod$lambda[which.min(mod$cvm)]
  #mod.best = glmnet(x = cbind(new.x.train, Z), y = Y.train, alpha = 0, lambda = lambda,
  #                  penalty.factor = c(rep(0, dim(new.x.train)[2]), rep(1, d)))
  
  #Y.pls[cv.ind.i] <- predict(mod.best, newx = cbind(new.x.test, Z.test))
  trndat = as.data.frame(cbind(new.x.train, x1train, x2train))
  colnames(trndat) = c("pc1","pc2","pc3","x1","x2")
  tstdat = as.data.frame(cbind(new.x.test, x1test, x2test))
  colnames(tstdat) = c("pc1","pc2","pc3","x1","x2")
  mod = mgcv::gam(Y.train~pc1+pc2+pc3+s(x1,x2,bs="tp",k=floor(n.train/2),m=2), data = trndat)
  Y.pls[cv.ind.i] = predict(mod, tstdat)
  
  #beta.hat.pcr = lklhds[[which(thetas == theta)]]$beta.hat.pcr
  #Y.pls[cv.ind.i] = lklhds[[which(thetas == theta)]]$Y.pred
  
  pls.time <- proc.time()-strt
  
  return(list(# pls-uk
    Y.pls = Y.pls, Y.pls.mean = Y.pls.mean, best.ncomps = best.ncomps#, beta.hat.pls = mod.best$beta,
    #pc.obj = pc.obj, theta.pls = lambda, pls.time = pls.time
    ))
}

plsTPRS.ufp.uw <- lapply(1:numfold, FUN = run, annavg = dat$ufp_uw)
plsTPRS.bc.uw <- lapply(1:numfold, FUN = run, annavg = dat$bc_uw)
plsTPRS.bc.primary <- lapply(1:numfold, FUN = run, annavg = dat$bc_primary)
plsTPRS.ufp.primary <- lapply(1:numfold, FUN = run, annavg = dat$ufp_primary)