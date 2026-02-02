## Set paths to directory
path_to_code <- "~/Dropbox/UW/SpatRF/JCGS/JASA_code/"
path_to_data <- "~/Dropbox/UW/SpatRF/JCGS/JASA_code/"
job_no <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(job_no*5,kind="L'Ecuyer-CMRG")

source(paste0(path_to_code,"aux_functions.R"))
source(paste0(path_to_code,"spatTree.R"))
library(randomForest)
library(SpatioTemporal)
library(matrixStats)
library(pls)
library(inline)

pollutant <- c("EC","OC","S","Si")[job_no]

source(paste0(path_to_code,"aux_functions.R"))
source(paste0(path_to_code,"spatTree.R"))

plltnt.mat <- read.csv(paste0(path_to_data,pollutant,".txt"))

## Remove duplicate site
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

annavg <- sqrt(annavg)

n <- length(annavg)
Y.tprs <-  Y.tprs.rf <- Y.rf <- Y.rf.tprs <- Y.rf.withbas <- 
  Y.spatrf.cv <- Y.spatrf.pl <- Y.pls <- rep(NA,n)

numfold <- 10

grid <- list(lambert_x =plltnt.mat$lambert_x,
             lambert_y =plltnt.mat$lambert_y)

smpl <- sample(n)
ind.stop <- round(n/numfold*1:numfold)
ind.strt <- c(1,ind.stop[1:numfold-1]+1)
tprs <- TRUE

##############
# Begin Runs #
##############
for (j in 1:numfold){
  
  cv.ind.i <- c(smpl[(ind.strt[j]:ind.stop[j])])
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
  
  R <- makeSigmaTPRS(c(1),train.coords,out.ZB=FALSE,
                      pen.form=FALSE,fit.tprs=TRUE,K=n.train)
  R.test <- (makeSigmaTPRS(1,train.coords,test.coords,out.ZB=FALSE,
                             pen.form=FALSE,fit.tprs = FALSE))
  
  Z <- makeSigmaTPRS(c(1),train.coords,out.ZB=TRUE,
                     pen.form=FALSE,fit.tprs=TRUE,K=n.train)
  Z.test <- makeSigmaTPRS(1,train.coords,test.coords,out.ZB=TRUE,
                          pen.form=FALSE,fit.tprs = FALSE,K=n.train)
  
  #############################################################
  # Method 1: Thin Plate Penalized Regression Splines by MGCV #
  #############################################################
  ## TPRS by MGCV
  strt <- proc.time()
  m <- 2
  mod <- mgcv::gam(Y.train~s(x1train,x2train,bs="tp",k=n.train,m=m))
  Y.tprs[cv.ind.i] <- predict(mod,data.frame(x1train=x1test,x2train=x2test))
  tprs.time <- proc.time() - strt
  
  #####################################
  # Method 2: TPRS then Random Forest #
  #####################################
  
  strt <- proc.time()
  tprs.rf.mod <- randomForest(X.train,mod$residuals,xtest=X.test,
                              nodesize=5,importance = FALSE)
  Y.tprs.rf[cv.ind.i] <- Y.tprs[cv.ind.i] + tprs.rf.mod$test$predicted
  tprs.rf.time <- tprs.time + proc.time()-strt
  
  #######################################
  # Method 3-4: Random Forest - TPRS    #
  #######################################
  
  strt <- proc.time()
  rf <- randomForest(X.train,Y.train,xtest=X.test,nodesize=5,importance = TRUE)
  mod <- mgcv::gam((Y.train-rf$predicted)~s(x1train,x2train,bs="tp",k=n.train,m=2,fx=FALSE))
  Y.rf.tprs[cv.ind.i]<- rf$test$predicted+ predict(mod,data.frame(x1train=x1test,x2train=x2test))
  Y.rf[cv.ind.i] <- rf$test$predicted
  rf.time <- proc.time()-strt
  
  #############################################
  # Method 5: Random Forest with TPRS basis #
  #############################################
  
  strt <- proc.time()
  rf.withbas <- randomForest(cbind(X.train,Z),Y.train,xtest=cbind(X.test,Z.test))
  Y.rf.withbas[cv.ind.i] <- rf.withbas$test$predicted
  rf.time <- proc.time()-strt
  
  ###################################################
  # Methods 6-7: Spatial Random Forests           #
  # 2 Versions : Optimization by CV                 #
  #              Optimization via Pseudo-Likelihood #
  ###################################################
  
  strt <- proc.time()
  
  cv_obj <- cvSpatRF(Y=Y.train,X=X.train,coords=train.coords,
                     cv.lambda=c(0.6,0.8),#expit(seq(-2,4,1/3)),
                     cov.type="TPRS",cov.opts=list(k=0,m=2),lklhd=TRUE,
                     Xtest = X.test, replace=TRUE,
                     coords.test = test.coords)
  
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
  

  #######################
  # Method 15: UK - PLS #
  #######################
  
  strt <- proc.time()
  
  cleaned_for_pc <- cleanGIS(as.matrix(X[-cv.ind.i,]),
                             as.matrix(X[cv.ind.i,]),
                             rm.outliers = FALSE)
  X.train.pc <- cleaned_for_pc$X
  X.test.pc <- cleaned_for_pc$X.test
  
  num.pcs <- 1:10
  pc.pars <- matrix(NA,ncol=length(num.pcs),nrow=3)
  for(i in num.pcs){
    pc.obj <- get.pcs(X.train.pc,"pls",i,Y.train)
    strt.vl <- optim(c(1,-1,1),logLikeExp,method="L-BFGS-B",
                     lower=c(-30,-20,-10), upper =c(15,15,30),
                     y=Y.train, dist.mat=dmat, x=pc.obj$X)
    pc.pars[,which(num.pcs==i)] <- strt.vl$par
  }
  cv.pcs.pls <- cv.pcr.exp(X.train.pc, Y.train, dmat, cv.pcs = num.pcs, k = 10, 
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
  
  
}
####

out <- cbind(annavg,Y.tprs,Y.tprs.rf, Y.rf,Y.rf.tprs, Y.rf.withbas,
             Y.pls,Y.spatrf.cv,Y.spatrf.pl )

save(out,file=paste0(path_to_data,"national_",pollutant,".Rdata"))


