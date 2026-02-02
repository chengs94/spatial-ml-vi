# exposure prediction models

#library(tidyverse)
#library(sf)
#library(blockCV)
library(SpatioTemporal)
library(matrixStats)
library(pls)
library(randomForest)
library(inline)
library(parallel)

setwd("/home/users/chengsi/Desktop/exposurepred/")

source("aux_functions.R")
source("spatTreeModified.R")

#source("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/spatialRF/aux_functions.R")
#source("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/spatialRF/spatTree.R")

# dat = readRDS("dat.rda")
# dat = readRDS("dat_updated.rda")
dat = readRDS("dat_final.rda")
# load("folds.rda")
set.seed(2333)

#hist(dat$ufp_uw) => log transformation?

n <- nrow(dat)
Y.pls <- Y.spatrf.cv <- Y.spatrf.pl <- Y.tprs <- Y.tprs.rf <- Y.rf <- Y.rf.tprs <- rep(NA, n)

numfold <- 5
grid <- list(lambert_x = dat$lambert_x,
             lambert_y = dat$lambert_y)

#temp = st_as_sf(dat, coords = c("lambert_x", "lambert_y"))
#sb <- spatialBlock(speciesData = temp, # sf or SpatialPoints
                   #species = "Species", # the response column (binomial or multi-class)
                   #rasterLayer = myrasters, # a raster for backgoround (optional)
#                   theRange = 10000, # size of the blocks
#                   k = numfold, # the number of folds
#                   selection = "random",
#                   iteration = 100, # find evenly dispersed folds
#                   biomod2Format = TRUE)
#rm(temp)

smpl <- sample(n)
ind.stop <- round(n/numfold*1:numfold)
ind.strt <- c(1,ind.stop[1:numfold-1]+1)

# X = dat[,-(1:21)]
X = dat[,-(1:12)]

run = function(fold, annavg){
  cv.ind.i <- c(smpl[(ind.strt[fold]:ind.stop[fold])])
  # cv.ind.i <- folds[[fold]][[2]]
  # cv.ind.i <- which(folds == fold)
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
  
  # UK - PLS
  strt <- proc.time()
  
  #cleaned_for_pc <- cleanGIS(as.matrix(X[-cv.ind.i,]),
  #                           as.matrix(X[cv.ind.i,]),
  #                           rm.outliers = FALSE)
  #X.train.pc <- cleaned_for_pc$X
  #X.test.pc <- cleaned_for_pc$X.test
  
  X.train.pc = X.train
  X.test.pc = X.test
  
  num.pcs <- 1:5
  pc.pars <- matrix(NA,ncol=length(num.pcs),nrow=3)
  for(i in num.pcs){
    pc.obj <- get.pcs(X.train.pc,"pls",i,Y.train)
    strt.vl <- optim(c(1,-1,1),logLikeExp,method="L-BFGS-B",
                     lower=c(-30,-20,-10), upper =c(15,15,30),
                     y=Y.train, dist.mat=dmat, x=pc.obj$X)
    pc.pars[,which(num.pcs==i)] <- strt.vl$par
  }
  cv.pcs.pls <- cv.pcr.exp(X.train.pc, Y.train, dmat, cv.pcs = num.pcs, k = numfold, 
                           method="pls",pars.init = pc.pars)
  pc.obj <- get.pcs(X.train.pc,"pls",cv.pcs.pls$opt.pcs,Y.train)
  
  new.x.train <- pc.obj$X
  means<-apply(X.train.pc,2,mean)
  mean.matrix<-matrix(rep(means,n.test),nrow=n.test,byrow=T)
  sds <- apply(X.train.pc,2,sd)
  new.x.test <- (X.test.pc - mean.matrix) %*% diag(1/sds) %*% pc.obj$proj
  
  theta.pls <- pc.pars[,which(num.pcs==cv.pcs.pls$opt.pcs)]
  
  sig.pls <- exp(theta.pls[1])
  nug.pls <- exp(theta.pls[2])
  range.pls <- exp(theta.pls[3])
  
  i.V <- .fastSolve(sig.pls * exp(-dmat/range.pls) + 
                      nug.pls * diag(length(Y.train)))
  beta.hat.pcr <- solve(t(new.x.train) %*% i.V %*% new.x.train,
                        t(new.x.train) %*% i.V %*% Y.train)
  Y.pls[cv.ind.i]<- c(new.x.test %*% beta.hat.pcr) + 
    sig.pls * exp(-dtestmat/range.pls) %*% i.V %*% 
    c(Y.train-c(new.x.train %*% beta.hat.pcr))
  
  pls.time <- proc.time()-strt
  
  # spatial RF
  strt <- proc.time()
  
  cv_obj <- cvSpatRF(Y=Y.train,X=X.train,coords=train.coords,
                     #cv.lambda = c(.6, .8),
                     # cv.lambda=seq(.88,.98,.02), #expit(seq(-2,4,1/3)),
                     cv.lambda=seq(.8,.98,.02),
                     #cv.r = seq(2, 10, 2), cv.m = floor(dim(X)[2]/seq(2.5, 4.5, .5)),
                     #cv.s = ceiling(length(Y.train)*seq(0.65, 0.9, 0.05)),
                     random.tune = TRUE, nums = 20,
                     range.r = c(2, 10), range.m = c(2.5, 4.5),
                     range.s = c(0.65, 0.9),
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
  
  # TPRS
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
  
  return(list(# pls-uk
    Y.pls = Y.pls, beta.hat.pls = beta.hat.pcr,
    pc.obj = pc.obj, theta.pls = theta.pls, pls.time = pls.time,
    # spat RF
    spatrf.varnames = colnames(X.train),
    Y.spatrf.cv = Y.spatrf.cv, spatrf.cv.obj = cv_obj,
    Y.spatrf.pl = Y.spatrf.pl, beta.hat.spatrf = beta.hat, spatrf.time = spatrf.time,
    # two-stage
    Y.tprs = Y.tprs, Y.tprs.rf = Y.tprs.rf, Y.rf = Y.rf, Y.rf.tprs = Y.rf.tprs,
    tprs.time = tprs.time, tprs.rf.time = tprs.rf.time, rf.time = rf.time,
    rf.imp = rf$importance, rf.impSD = rf$importanceSD,
    tprs.rf.imp = tprs.rf.mod$importance, tprs.rf.impSD = tprs.rf.mod$importanceSD))
}

numCores=detectCores()

# rslt.bc.uw <- mclapply(1:numfold, FUN = run, annavg = dat$bc_uw,
#                            mc.cores = numCores)
# rslt.ufp.uw <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_uw,
#                             mc.cores = numCores)

rslt.ufp.uw <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_uw,
                        mc.cores = numCores)
rslt.bc.uw <- mclapply(1:numfold, FUN = run, annavg = dat$bc_uw,
                        mc.cores = numCores)
rslt.co2.uw <- mclapply(1:numfold, FUN = run, annavg = dat$co2_uw,
                        mc.cores = numCores)
rslt.no2.uw <- mclapply(1:numfold, FUN = run, annavg = dat$no2_uw,
                       mc.cores = numCores)
rslt.pm25.uw <- mclapply(1:numfold, FUN = run, annavg = dat$pm25_uw,
                        mc.cores = numCores)

save(rslt.ufp.uw, rslt.bc.uw, rslt.co2.uw,
     rslt.no2.uw, rslt.pm25.uw, file = "mm_wind_all.RData")
save.image("mm_wind_all_full.RData")
# save(rslt.ufp.uw, rslt.bc.uw, file = "mm_ufp_bc.RData")
# save.image("mm_ufp_bc_full.RData")

# folds = fold_block
# block.ufp.uw <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_uw,
#                         mc.cores = numCores)
# block.bc.uw <- mclapply(1:numfold, FUN = run, annavg = dat$bc_uw,
#                        mc.cores = numCores)
# save.image("trap_cv_block.RData")
# 
# folds = fold_strat
# strat.ufp.uw <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_uw,
#                         mc.cores = numCores)
# strat.bc.uw <- mclapply(1:numfold, FUN = run, annavg = dat$bc_uw,
#                        mc.cores = numCores)
# save.image("trap_cv_block.RData")

#rslt.bc.primary <- mclapply(1:numfold, FUN = run, annavg = dat$bc_primary,
#                            mc.cores = numCores)
#rslt.ufp.primary <- mclapply(1:numfold, FUN = run, annavg = dat$ufp_primary,
#                             mc.cores = numCores)
#save.image("trap_cv_morepars.RData")