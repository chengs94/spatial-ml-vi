# decomposing the variance from different components

library(SpatioTemporal)
library(matrixStats)
library(pls)
library(inline)
library(parallel)
library(randomForest)

setwd("/home/students/chengsi/Desktop/exposurepred/more")

source("aux_functions.R")
source("spatTreeModified.R")

dat = readRDS("dat.rda")
set.seed(2333)

#hist(dat$ufp_uw) => log transformation?

n <- nrow(dat)
#Y.pls <- Y.spatrf.cv <- Y.spatrf.pl <- rep(NA, n)
Y.tprs <-  Y.tprs.rf <- Y.rf <- Y.rf.tprs <- rep(NA, n)

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
  
  #cleaned <- cleanGIS(as.matrix(X[-cv.ind.i,]),
  #                    as.matrix(X[cv.ind.i,]),
  #                    rm.outliers = FALSE)
  #X.train <- cleaned$X
  #X.test <- cleaned$X.test
  
  X.train = as.matrix(X[-cv.ind.i,])
  X.test = as.matrix(X[cv.ind.i,])
  
  dmat <- as.matrix(SpatioTemporal::crossDist(train.coords))
  dtestmat <- as.matrix(SpatioTemporal::crossDist(test.coords,train.coords))
  
  # TPRS
  
  #cleaned_for_pc <- cleanGIS(as.matrix(X[-cv.ind.i,]),
  #                           as.matrix(X[cv.ind.i,]),
  #                           rm.outliers = FALSE)
  #X.train.pc <- cleaned_for_pc$X
  #X.test.pc <- cleaned_for_pc$X.test
  
  X.train.pc = X.train
  X.test.pc = X.test
  
  strt <- proc.time()
  m <- 2
  mod <- mgcv::gam(Y.train~s(x1train,x2train,bs="tp",k=n.train,m=m))
  Y.tprs[cv.ind.i] <- predict(mod,data.frame(x1train=x1test,x2train=x2test))
  tprs.time <- proc.time() - strt
  
  # TPRS - RF
  strt <- proc.time()
  tprs.rf.mod <- randomForest(X.train,mod$residuals,xtest=X.test,
                              nodesize=5,importance = TRUE)
  Y.tprs.rf[cv.ind.i] <- Y.tprs[cv.ind.i] + tprs.rf.mod$test$predicted
  tprs.rf.time <- tprs.time + proc.time()-strt
  
  # RF; RF - TPRS
  strt <- proc.time()
  rf <- randomForest(X.train,Y.train,xtest=X.test,nodesize=5,importance = TRUE)
  mod <- mgcv::gam((Y.train-rf$predicted)~s(x1train,x2train,bs="tp",k=n.train,m=2,fx=FALSE))
  Y.rf.tprs[cv.ind.i]<- rf$test$predicted+ predict(mod,data.frame(x1train=x1test,x2train=x2test))
  Y.rf[cv.ind.i] <- rf$test$predicted
  rf.time <- proc.time()-strt
  
  return(list(
    Y.tprs = Y.tprs, Y.tprs.rf = Y.tprs.rf, Y.rf = Y.rf, Y.rf.tprs = Y.rf.tprs,
    tprs.time = tprs.time, tprs.rf.time = tprs.rf.time, rf.time = rf.time,
    rf.imp = rf$importance, rf.impSD = rf$importanceSD,
    tprs.rf.imp = tprs.rf.mod$importance, tprs.rf.impSD = tprs.rf.mod$importanceSD
    ))
}

numCores=detectCores()

ts.ufp.uw <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_uw,
                        mc.cores = numCores)
ts.bc.uw <- mclapply(1:numfold, FUN = run, annavg = dat$bc_uw,
                       mc.cores = numCores)
save.image("two_stage_RFimp.RData")

ts.bc.primary <- mclapply(1:numfold, FUN = run, annavg = dat$bc_primary,
                            mc.cores = numCores)
ts.ufp.primary <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_primary,
                             mc.cores = numCores)
save.image("two_stage_RFimp.RData")
