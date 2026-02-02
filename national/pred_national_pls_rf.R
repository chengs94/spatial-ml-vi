library(SpatioTemporal)
library(matrixStats)
library(pls)
# library(inline)
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
Y.pls <- Y.rf <- Y.rf.tprs <- matrix(NA, nrow = n, ncol = ncol(X)*3)

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
  
  col.idx = 1:ncol(X.test)
  qts = c(.25, .5, .75)
  # n.test <- dim(X.test)[1]
  
  # UK - PLS
  strt <- proc.time()
  
  X.train.pc <- X.train
  X.test.pc <- X.test
  
  num.pcs <- 1:5
  pc.pars <- matrix(NA,ncol=length(num.pcs),nrow=3)
  for(i in num.pcs){
    pc.obj <- get.pcs(X.train.pc,"pls",i,Y.train)
    strt.vl <- optim(c(1,-1,1),logLikeExp,method="L-BFGS-B",
                     lower=c(-30,-20,-10), upper =c(15,15,30),
                     y=Y.train, dist.mat=dmat, x=pc.obj$X)
    pc.pars[,which(num.pcs==i)] <- strt.vl$par
  }
  cv.pcs.pls <- cv.pcr.exp(X.train.pc, Y.train, dmat, cv.pcs = num.pcs, k = 5, 
                           method="pls",pars.init = pc.pars)
  pc.obj <- get.pcs(X.train.pc,"pls",cv.pcs.pls$opt.pcs,Y.train)
  
  new.x.train <- pc.obj$X
  means<-apply(X.train.pc,2,mean)
  mean.matrix<-matrix(rep(means,n.test),nrow=n.test,byrow=T)
  sds <- apply(X.train.pc,2,sd)
  
  theta.pls <- pc.pars[,which(num.pcs==cv.pcs.pls$opt.pcs)]
  
  sig.pls <- exp(theta.pls[1])
  nug.pls <- exp(theta.pls[2])
  range.pls <- exp(theta.pls[3])
  
  i.V <- .fastSolve(sig.pls * exp(-dmat/range.pls) + 
                      nug.pls * diag(length(Y.train)))
  beta.hat.pcr <- solve(t(new.x.train) %*% i.V %*% new.x.train,
                        t(new.x.train) %*% i.V %*% Y.train)
  
  X.test.pc0 = X.test.pc
  k = 0
  
  for (l in col.idx){
    for (j in 1:length(qts)){
      X.test.pc = X.test.pc0
      X.test.pc[,l] = quantile(X.test.pc[,l], qts[j])
      new.x.test <- (X.test.pc - mean.matrix) %*% diag(1/sds) %*% pc.obj$proj
      k = k + 1
      
      Y.pls[cv.ind.i, k]<- c(new.x.test %*% beta.hat.pcr) + 
        sig.pls * exp(-dtestmat/range.pls) %*% i.V %*% 
        c(Y.train-c(new.x.train %*% beta.hat.pcr))
      
      if (k %% 100 == 0) message("PLS", k)
    }
  }

  pls.time <- proc.time()-strt
  
  # strt <- proc.time()
  # X.test0 = X.test
  # k = 0
  # 
  # for (l in col.idx){
  #   for (j in 1:length(qts)){
  #     X.test = X.test0
  #     X.test[,l] = quantile(X.test[,l], qts[j])
  #     k = k + 1
  #     
  #     rf <- randomForest(X.train,Y.train,xtest=X.test,nodesize=5,importance = TRUE)
  #     mod <- mgcv::gam((Y.train-rf$predicted)~s(x1train,x2train,bs="tp",k=n.train,m=2,fx=FALSE))
  #     
  #     Y.rf.tprs[cv.ind.i, k]<- rf$test$predicted+ predict(mod,data.frame(x1train=x1test,x2train=x2test))
  #     Y.rf[cv.ind.i, k] <- rf$test$predicted
  #     
  #     if (k %% 100 == 0) message(paste("rf", k))
  #   }
  # }
  # 
  # rf.time <- proc.time()-strt
  # 
  # message(paste("one round done, col_idx =", col_idx, "quantile = ", qt))
  
  return(list(# pls-uk
    Y.pls = Y.pls, beta.hat.pls = beta.hat.pcr,
    pc.obj = pc.obj, theta.pls = theta.pls, pls.time = pls.time,
    # RF and RF-TPRS
    # Y.rf = Y.rf, Y.rf.tprs = Y.rf.tprs, rf.time = rf.time,
    # others
    X.train = X.train, X.test = X.test))
}

numCores = detectCores() - 1

# pars = expand.grid(1:ncol(X), c(.25, .5, .75))
# 
# rslt = mcmapply(run, mc.cores = 1, #mc.cores = numCores,
#                 col_idx = pars[,1], qt = pars[,2],
#                 MoreArgs = list(fold = 1))

rslt <- run(fold = fold)

save.image("national_ec_pls_rf_full.RData")
save(X, annavg, plltnt.mat, rslt, file = "national_ec_pls_rf.RData")
