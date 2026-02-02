# exposure prediction models

#library(tidyverse)
library(SpatioTemporal)
library(matrixStats)
library(pls)
library(inline)
library(parallel)

setwd("/home/students/chengsi/Desktop/exposurepred/")

source("aux_functions.R")
source("spatTreeModified.R")

#source("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/spatialRF/aux_functions.R")
#source("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/spatialRF/spatTree.R")

dat = readRDS("dat.rda")
set.seed(2333)

#hist(dat$ufp_uw) => log transformation?

n <- nrow(dat)
Y.pls <- Y.spatrf.cv <- Y.spatrf.pl <- rep(NA, n)

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
  
  pca.obj = prcomp(as.matrix(X[-cv.ind.i,]), center = TRUE, scale = TRUE, retx = TRUE)
  imp = summary(pca.obj)$importance[2,]
  ncomp = varexpl = 0
  while (varexpl < .8){
    ncomp = ncomp + 1
    varexpl = varexpl + imp[ncomp]
  }
  newpca = prcomp(as.matrix(X[-cv.ind.i,]), center = TRUE, scale = TRUE, retx = TRUE,
                   rank. = ncomp)
  X.train = newpca$x
  X.test = predict(newpca, as.matrix(X[cv.ind.i, ]))
  
  #X.train = as.matrix(X[-cv.ind.i,])
  #X.test = as.matrix(X[cv.ind.i,])
  
  dmat <- as.matrix(SpatioTemporal::crossDist(train.coords))
  dtestmat <- as.matrix(SpatioTemporal::crossDist(test.coords,train.coords))
  
  # spatial RF
  strt <- proc.time()
  
  cv_obj <- cvSpatRF(Y=Y.train,X=X.train,coords=train.coords,
                     #cv.lambda = c(.6, .8),
                     #cv.lambda=seq(.88,.98,.02), #expit(seq(-2,4,1/3)),
                     #cv.r = seq(2, 10, 2), cv.m = floor(dim(X)[2]/seq(2.5, 4.5, .5)),
                     #cv.s = ceiling(length(Y.train)*seq(0.65, 0.9, 0.05)),
                     random.tune = TRUE, nums = 20,
                     range.r = c(2, 10), range.m = c(2, 5),
                     range.s = c(0.5, 0.9),
                     cov.type="TPRS",cov.opts=list(k=0,m=2),lklhd=TRUE,
                     Xtest = X.test, replace=TRUE,
                     coords.test = test.coords, t = 500, 
                     var.imp=TRUE,imp.msr = "perm")
  
  Y.spatrf.cv[cv.ind.i] <- cv_obj$spatRF$ftest + cv_obj$spatRF$ztest
  
  fullspatbas <- .makeSpatBas(pars=list(psill=cv_obj$lambda.lklhd,
                                        nugget=1-cv_obj$lambda.lklhd),
                              coords.train=train.coords,
                              coords.test = test.coords,
                              cov.type="TPRS",cov.opts = list(k=0,m=2))
  
  beta.hat <- solve(t(train.coords) %*% fullspatbas$i.sig %*% train.coords,
                    t(train.coords) %*% fullspatbas$i.sig %*% 
                      (Y.train - cv_obj$spatRF.lklhd$fpredicted ) )
  Y.spatrf.pl[cv.ind.i] <- predict.cor( cv_obj$spatRF.lklhd$ftest + 
                                          test.coords%*% beta.hat, 
                                        Y.train - train.coords %*% beta.hat - 
                                          cv_obj$spatRF.lklhd$fpredicted, 
                                        theta = cv_obj$lambda.lklhd, 
                                        R.test=fullspatbas$sig.test, 
                                        i.sig = fullspatbas$i.sig )
  
  spatrf.time <- proc.time()-strt
  
  return(list(# pls-uk
    #Y.pls = Y.pls, beta.hat.pls = beta.hat.pcr,
    #pc.obj = pc.obj, theta.pls = theta.pls, pls.time = pls.time,
    # spat RF
    spatrf.varnames = colnames(X.train),
    Y.spatrf.cv = Y.spatrf.cv, spatrf.cv.obj = cv_obj,
    Y.spatrf.pl = Y.spatrf.pl, beta.hat.spatrf = beta.hat, spatrf.time = spatrf.time))
}

numCores=detectCores()

rslt.ufp.uw <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_uw,
                        mc.cores = numCores)
rslt.bc.uw <- mclapply(1:numfold, FUN = run, annavg = dat$bc_uw,
                       mc.cores = numCores)
save.image("trap_cv_pca.RData")

rslt.bc.primary <- mclapply(1:numfold, FUN = run, annavg = dat$bc_primary,
                            mc.cores = numCores)
rslt.ufp.primary <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_primary,
                             mc.cores = numCores)
save.image("trap_cv_pca.RData")
#save.image("trap_cv_grid1.RData")