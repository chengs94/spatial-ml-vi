# prediction on a grid / at cohort locations
# very similar to model.R, just assigning the grid covariates as test data

library(SpatioTemporal)
library(matrixStats)
library(pls)
library(inline)
library(randomForest)

# load CV result to get optimal tuning parameters
load("trap_cv_morepars_finer.RData")
load("trap_cv_morepars_finer2.RData")
#cov_grid = readRDS("cov_grid.rda")
cov_cohort = readRDS("cov_cohort.rda")
#source("spatTreeModifiedMean.R")
source("spatTreeModified.R")
#temp = read.csv("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/data/dr0311_grid_covars.csv")
#cov_grid = cbind(temp$longitude, temp$latitude, cov_grid)
#colnames(cov_grid)[1:2] = c("longitude", "latitude")
#saveRDS(cov_grid, "cov_grid.rda")

rslt.all = list(rslt.ufp.primary, rslt.ufp.uw, rslt.bc.primary, rslt.bc.uw)

bestNumComp = rep(0, length(rslt.all))
bestPLpar = matrix(0, ncol = 4, nrow = length(rslt.all))
bestCVpar = matrix(0, ncol = 4, nrow = length(rslt.all))
#colnames(bestSpatrfPar) = c("lambda", "r", "m", "s")

set.seed(2333)

for (j in 1:length(rslt.all)){
  for (k in 1:numfold){
    bestNumComp[j] = bestNumComp[j] + ncol(rslt.all[[j]][[k]]$pc.obj$loadings)
    bestCVpar[j,1] = bestCVpar[j,1] + rslt.all[[j]][[k]]$spatrf.cv.obj$lambda.min
    bestCVpar[j,2] = bestCVpar[j,2] + rslt.all[[j]][[k]]$spatrf.cv.obj$r.min
    bestCVpar[j,3] = bestCVpar[j,3] + rslt.all[[j]][[k]]$spatrf.cv.obj$m.min
    bestCVpar[j,4] = bestCVpar[j,4] + rslt.all[[j]][[k]]$spatrf.cv.obj$s.min
    bestPLpar[j,1] = bestPLpar[j,1] + rslt.all[[j]][[k]]$spatrf.cv.obj$lambda.lklhd
    bestPLpar[j,2] = bestPLpar[j,2] + rslt.all[[j]][[k]]$spatrf.cv.obj$r.lklhd
    bestPLpar[j,3] = bestPLpar[j,3] + rslt.all[[j]][[k]]$spatrf.cv.obj$m.lklhd
    bestPLpar[j,4] = bestPLpar[j,4] + rslt.all[[j]][[k]]$spatrf.cv.obj$s.lklhd
  }
}
bestNumComp = round(bestNumComp/numfold)
bestPLpar = bestPLpar/numfold
bestCVpar = bestCVpar/numfold
bestPLpar[,2:4] = round(bestPLpar[,2:4])
bestCVpar[,2:4] = round(bestCVpar[,2:4])

#X.test = cov_grid[,colnames(cov_grid) %in% colnames(X)]
X.test = cov_cohort[,colnames(cov_cohort) %in% colnames(X)]
X.test = X.test[names(X)]

grid_pred = function(annavg, numComp, PLpar, CVpar, X.train, X.test){
  Y.train = annavg
  n.test = nrow(X.test)
  n.train = nrow(X)
  #spatRF.par = rbind(PLpar, CVpar)
  
  train.coords <- cbind(dat$lambert_x, dat$lambert_y)
  scale.params <- cbind(apply(train.coords,2,mean),apply(train.coords,2,sd))
  x1train <- train.coords[,1] <- (train.coords[,1] - scale.params[1,1])/scale.params[1,2]
  x2train <-train.coords[,2] <- (train.coords[,2] - scale.params[2,1])/scale.params[2,2]
  # test.coords <- cbind(cov_grid$lambert_x, cov_grid$lambert_y)
  test.coords <- cbind(cov_cohort$lambert_x, cov_cohort$lambert_y)
  x1test <- test.coords[,1] <- (test.coords[,1] - scale.params[1,1])/scale.params[1,2]
  x2test <- test.coords[,2] <- (test.coords[,2] - scale.params[2,1])/scale.params[2,2]
  
  dmat <- as.matrix(SpatioTemporal::crossDist(train.coords))
  dtestmat <- as.matrix(SpatioTemporal::crossDist(test.coords,train.coords))
  
  # uk-pls
  pc.obj <- get.pcs(X.train,"pls",numComp,Y.train)
  
  new.x.train <- pc.obj$X
  means<-apply(X.train,2,mean)
  mean.matrix<-matrix(rep(means,n.test),nrow=n.test,byrow=T)
  sds <- apply(X.train,2,sd)
  new.x.test <- (X.test - mean.matrix) %*% diag(1/sds) %*% pc.obj$proj
  
  pc.obj <- get.pcs(X.train,"pls",numComp,Y.train)
  strt.vl <- optim(c(1,-1,1),logLikeExp,method="L-BFGS-B",
                   lower=c(-30,-20,-10), upper =c(15,15,30),
                   y=Y.train, dist.mat=dmat, x=pc.obj$X)
  pc.pars <- strt.vl$par
  theta.pls <- pc.pars
  
  sig.pls <- exp(theta.pls[1])
  nug.pls <- exp(theta.pls[2])
  range.pls <- exp(theta.pls[3])
  
  i.V <- .fastSolve(sig.pls * exp(-dmat/range.pls) + 
                      nug.pls * diag(length(Y.train)))
  beta.hat.pcr <- solve(t(new.x.train) %*% i.V %*% new.x.train,
                        t(new.x.train) %*% i.V %*% Y.train)
  Y.pls <- c(new.x.test %*% beta.hat.pcr) + 
    sig.pls * exp(-dtestmat/range.pls) %*% i.V %*% 
    c(Y.train-c(new.x.train %*% beta.hat.pcr))
  
  Y.pls.mean <- c(new.x.test %*% beta.hat.pcr)
  
  # spatRF
  spatRF.pl <- try(spatRF( pars=list(psill=PLpar[1],nugget=1-PLpar[1]), 
                           Y=Y.train, X=X.train, Xtest = X.test, 
                           coords.test = test.coords, coords=train.coords, 
                           cov.type = "TPRS", cov.opts = list(k=0,m=2), 
                           r = PLpar[2], m = PLpar[3], s = PLpar[4],
                           var.imp = TRUE, imp.msr = "imp",
                           t = 500, replace = TRUE))
  if (class(spatRF.pl) != "try-error"){
    fullspatbas <- .makeSpatBas(pars=list(psill=PLpar[1],
                                          nugget=1-PLpar[1]),
                                coords.train=train.coords,
                                coords.test = test.coords,
                                cov.type="TPRS",cov.opts = list(k=0,m=2))
    beta.hat <- solve(t(train.coords) %*% fullspatbas$i.sig %*% train.coords,
                      t(train.coords) %*% fullspatbas$i.sig %*% 
                        (Y.train - spatRF.pl$fpredicted ) )
    Y.spatrf.pl <- predict.cor( spatRF.pl$ftest + 
                                test.coords%*% beta.hat, 
                                Y.train - train.coords %*% beta.hat - 
                                spatRF.pl$fpredicted, 
                                #theta = 0, 
                                theta = PLpar[1], 
                                R.test=fullspatbas$sig.test, 
                                i.sig = fullspatbas$i.sig )
    #Y.spatrf.pl = spatRF.pl$ftest + test.coords%*% beta.hat
    #Y.spatrf.pl = spatRF.pl$fpredicted + spatRF.pl$zpredicted
    #cur.mse <- sum((Y-y.test)^2)/n
  } else{
    Y.spatrf.pl = rep(NA, n.test)
    #cur.mse = Inf
  }
  
  spatRF.cv <- try(spatRF( pars=list(psill=CVpar[1],nugget=1-CVpar[1]), 
                           Y=Y.train, X=X.train, Xtest = X.test, 
                           coords.test = test.coords, coords=train.coords, 
                           cov.type = "TPRS", cov.opts = list(k=0,m=2), 
                           r = CVpar[2], m = CVpar[3], s = CVpar[4],
                           var.imp = TRUE, imp.msr = "imp",
                           t = 500, replace = TRUE))
  
  if (class(spatRF.cv) != "try-error"){
    Y.spatrf.cv = spatRF.cv$ftest + spatRF.cv$ztest
  } else{
    Y.spatrf.cv = rep(NA, n.test)
    #cur.mse = Inf
  }
  
  # TPRS
  #strt <- proc.time()
  m <- 2
  mod <- mgcv::gam(Y.train~s(x1train,x2train,bs="tp",k=n.train,m=m))
  Y.tprs <- predict(mod,data.frame(x1train=x1test,x2train=x2test))
  #tprs.time <- proc.time() - strt
  
  # TPRS - RF
  #strt <- proc.time()
  tprs.rf.mod <- randomForest(X.train,mod$residuals,xtest=X.test,
                              nodesize=5,importance = TRUE)
  Y.tprs.rf <- Y.tprs + tprs.rf.mod$test$predicted
  #tprs.rf.time <- tprs.time + proc.time()-strt
  
  # RF; RF - TPRS
  #strt <- proc.time()
  rf <- randomForest(X.train,Y.train,xtest=X.test,nodesize=5,importance = TRUE)
  mod <- mgcv::gam((Y.train-rf$predicted)~s(x1train,x2train,bs="tp",k=n.train,m=2,fx=FALSE))
  Y.rf.tprs <- rf$test$predicted+ predict(mod,data.frame(x1train=x1test,x2train=x2test))
  Y.rf <- rf$test$predicted
  #rf.time <- proc.time()-strt
  
  out = list(Y.pls = Y.pls, Y.pls.mean = Y.pls.mean, 
             pc.obj = pc.obj, theta.pls = theta.pls,
             # spat RF
             spatrf.varnames = colnames(X.train),
             Y.spatrf.cv = Y.spatrf.cv, spatrf.cv = spatRF.cv,
             Y.spatrf.pl = Y.spatrf.pl, spatrf.pl = spatRF.pl, 
             # two-stage
             Y.tprs = Y.tprs, Y.tprs.rf = Y.tprs.rf, Y.rf = Y.rf, Y.rf.tprs = Y.rf.tprs,
             rf.imp = rf$importance, rf.impSD = rf$importanceSD,
             tprs.rf.imp = tprs.rf.mod$importance, tprs.rf.impSD = tprs.rf.mod$importanceSD)
  return(out)
}

grid.ufp.primary = grid_pred(annavg = dat$ufp_primary, numComp = bestNumComp[1],
                        PLpar = bestPLpar[1,], CVpar = bestCVpar[1,],
                        X.train = as.matrix(X), X.test = as.matrix(X.test))

grid.ufp.uw = grid_pred(annavg = dat$ufp_uw, numComp = bestNumComp[2],
                        PLpar = bestPLpar[2,], CVpar = bestCVpar[2,],
                        X.train = as.matrix(X), X.test = as.matrix(X.test))
#save.image("pred_grid_mean.RData")

grid.bc.primary = grid_pred(annavg = dat$bc_primary, numComp = bestNumComp[3],
                       PLpar = bestPLpar[3,], CVpar = bestCVpar[3,],
                       X.train = as.matrix(X), X.test = as.matrix(X.test))

grid.bc.uw = grid_pred(annavg = dat$bc_uw, numComp = bestNumComp[4],
                       PLpar = bestPLpar[4,], CVpar = bestCVpar[4,],
                       X.train = as.matrix(X), X.test = as.matrix(X.test))
#save.image("pred_grid_mean.RData")
#save.image("pred_grid_imp.RData")
#save.image("rslt_grid_uw2.RData")

save.image("pred_cohort_image.RData")

cov_cohort = cov_cohort[c("longitude", "latitude")]
save(grid.ufp.primary, grid.ufp.uw, grid.bc.primary, grid.bc.uw,
     cov_cohort, file = "pred_cohort.RData")
