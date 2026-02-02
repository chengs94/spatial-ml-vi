library(SpatioTemporal)
library(matrixStats)
library(pls)
library(inline)
library(parallel)
library(randomForest)

setwd("/home/users/chengsi/Desktop/exposurepred/national")

source("aux_functions.R")
source("spatTreeModified2.R")

job_no = 1
set.seed(job_no*5, kind="L'Ecuyer-CMRG")

pollutant <- c("EC","OC","S","Si")[job_no]
plltnt.mat <- read.csv(paste0("data/", pollutant,".txt"))
plltnt.mat <- plltnt.mat[-which(plltnt.mat$native_id=="PHOE5"),]

if(pollutant %in% c("S","Si")){
  annavg <- plltnt.mat[,2]
  X <- plltnt.mat[,c(14:25,27:613)]
} else{
  sum1 <- ifelse(is.na(plltnt.mat[,2]),0,plltnt.mat[,2]*plltnt.mat[,3])
  obs1 <-  ifelse(is.na(plltnt.mat[,3]),0,plltnt.mat[,3])
  sum2 <- ifelse(is.na(plltnt.mat[,5]),0,plltnt.mat[,5]*plltnt.mat[,6])
  obs2 <- ifelse(is.na(plltnt.mat[,6]),0,plltnt.mat[,6])
  annavg <- (sum1+sum2)/(obs1+obs2)
  X <- plltnt.mat[,c(18:29,31:617)]
}

# cleaned <- cleanGIS(as.matrix(X), NULL, rm.outliers = FALSE)
# X0 = X
# X = cleaned$X

annavg <- sqrt(annavg)
n <- length(annavg)
Y.spatrf.cv <- Y.spatrf.pl <- matrix(NA, nrow = n, ncol = ncol(X)*3)

numfold <- 3
grid <- list(lambert_x = plltnt.mat$lambert_x,
             lambert_y = plltnt.mat$lambert_y)

smpl <- sample(n)
ind.stop <- round(n/numfold*1:numfold)
ind.strt <- c(1,ind.stop[1:numfold-1]+1)

fold = 1

run = function(fold){
  
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
  
  cleaned <- cleanGIS(as.matrix(X[-cv.ind.i,]),
                      as.matrix(X[cv.ind.i,]),
                      rm.outliers = FALSE)
  X.train <- cleaned$X
  X.test <- cleaned$X.test
  # X.train = X[-cv.ind.i,]
  # X.test = X[cv.ind.i, ]
  
  # if (col_idx <= ncol(X.test)){
  #   X.test[,col_idx] = quantile(X.test[,col_idx], qt)
  # } else{
  #   return(NULL)
  # }
  
  dmat <- as.matrix(SpatioTemporal::crossDist(train.coords))
  dtestmat <- as.matrix(SpatioTemporal::crossDist(test.coords,train.coords))
  
  # spatial RF
  strt <- proc.time()
  
  spatRF.pl <- try(spatRF_trntst( pars=list(psill=.7,nugget=.3), 
                           Y=Y.train, X=X.train, Xtest = X.test, 
                           coords.test = test.coords, coords=train.coords, 
                           cov.type = "TPRS", cov.opts = list(k=0,m=2), 
                           var.imp = FALSE, imp.msr = "imp",
                           t = 500, replace = TRUE))
  
  # cv_obj <- cvSpatRF(Y=Y.train,X=X.train,coords=train.coords,
  #                    #cv.lambda=c(0.6,0.8),#expit(seq(-2,4,1/3)),
  #                    cv.lambda = .7,
  #                    cov.type="TPRS",cov.opts=list(k=0,m=2),lklhd=TRUE,
  #                    Xtest = X.test, replace=TRUE,
  #                    coords.test = test.coords, t = 500, 
  #                    var.imp=FALSE,imp.msr = "perm")
  
  Y.spatrf.cv[cv.ind.i, 1:ncol(spatRF.pl$ftest)] <- spatRF.pl$ftest + spatRF.pl$ztest
  
  fullspatbas <- .makeSpatBas(pars=list(psill=.7, nugget=.3),
                              coords.train=train.coords,
                              coords.test = test.coords,
                              cov.type="TPRS",cov.opts = list(k=0,m=2))
  
  beta.hat <- solve(t(train.coords) %*% fullspatbas$i.sig %*% train.coords,
                    t(train.coords) %*% fullspatbas$i.sig %*% 
                      (Y.train - spatRF.pl$fpredicted ) )
  for (j in 1:ncol(spatRF.pl$ftest)){
    Y.spatrf.pl[cv.ind.i, j] <- predict.cor( spatRF.pl$ftest[,j] + 
                                            test.coords%*% beta.hat, 
                                          Y.train - train.coords %*% beta.hat - 
                                            spatRF.pl$fpredicted, 
                                          theta = .7, 
                                          R.test=fullspatbas$sig.test, 
                                          i.sig = fullspatbas$i.sig )
  }
  
  spatrf.time <- proc.time()-strt

  
  return(list(
    # spat RF
    spatrf.varnames = colnames(X.train),
    Y.spatrf.cv = Y.spatrf.cv, spatrf.obj = spatRF.pl,
    Y.spatrf.pl = Y.spatrf.pl, beta.hat.spatrf = beta.hat, 
    spatrf.time = spatrf.time,
    X.train = X.train, X.test = X.test))
}

# numCores = detectCores() - 1

# rslt = mcmapply(run, mc.cores = 1, #mc.cores = numCores,
#                 col_idx = 1:ncol(X), qt = c(.25, .5, .75),
#                 MoreArgs = list(fold = 1))
# rslt <- lapply(1:numfold, FUN = run)
rslt = run(fold = fold)

save.image("national_ec_spatRF_full.RData")
save(X, annavg, plltnt.mat, rslt, file = "national_ec_spatRF.RData")
