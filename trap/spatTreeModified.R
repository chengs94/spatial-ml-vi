############
# spatTree #
############
# Inputs: Y - Outcome
#         X - Covariates to build tree on
#         X.fix - Fixed Covariates to regress on in addition 
#                   to building the tree
#         Sig - Correlation matrix between observations
#         i.Sig - inverse correlation matrix, used to save time in loops
#         m - number of covariates to sample at each stage in the tree
#             only samples from X, not X.bas
#         r - minimum leaf size
#         k - max number of leaves to use
#         lambda - ridge penalty for the basis functions

spatTree <- function(Y, X, X.fix = NULL, Sig = NULL, i.Sig = NULL, 
                     m = ceiling(dim(X)[2]/3), r = 5, k = floor(length(Y)/r), 
                     lambda = 0,...){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  if(is.null(X.fix)){
    b <- NULL
  } else{
    p.fix <- dim(X.fix)[2]
    b <- rep(NA,p.fix)
  }
  
  if(is.null(i.Sig)) {
    if(is.null(Sig)){
      i.Sig <- diag(rep(1,n))
    } else{
      i.Sig <- .fastSolve(Sig)
    }
  }
  
  if(dim(X)[2] != 1){
    sortind <- t(t(matrix(order(rep(1:p,each=n),c(X)),n,p)) - c(0:(p-1)*n))
  } else{
    sortind <- t(t(matrix(order(rep(1:p,each=n),c(X)),n,p)))
  }
  covar.num <- split.val <- branch.loss <- leaf.split <- rep(NA,k)
  Ak <- Ck <- matrix(NA,n,k+1)
  Ak[,1] <- Ck[,1] <- 1
  a <- Ck.ind <- 1
  b <- NULL
  
  ## First, Start with the intercept estimate
  ## Calculate starting kernel estimate
  if(!is.null(X.fix)) {
      X.start <- cbind(X.fix,Ck[,1])
      pen <- c(rep(lambda,p.fix),0)
      Kern.C <- i.Sig - i.Sig %*% X.start %*% 
                  .fastSolve( diag(pen) + t(X.start) %*% i.Sig %*% X.start ) %*% 
                  t(X.start) %*% i.Sig
  } else {
    Kern.C <- i.Sig - rowSums(i.Sig) %*% t(colSums(i.Sig)) / c(sum(i.Sig))
  }
  
  loss.full <- t(Y) %*% Kern.C %*% Y
  omegaY <- as.numeric(Kern.C %*% Y)
  
  while( max(colSums2( Ak[,1:a,drop=F] )) > (2 * r - 1)  & a <= k ) {
    
    smpl <- sample(1:p,m)
    num.per.leaf <- colSums2(Ak[,1:a,drop=F])
    leaf.ind <- which(num.per.leaf > (2*r-1))
    n.k.full <- num.per.leaf[leaf.ind]
    Bk <- Ak[,leaf.ind,drop=FALSE]
    leaf.mem <- rowSums2(Bk %*% leaf.ind)
    leaf.mem.use <- match(leaf.mem,c(0,leaf.ind))-1
    leaf.use <- 1:length(leaf.ind)
    if(any(leaf.mem==0) ){
      n.k.full <- c(n-sum(n.k.full),n.k.full)
      leaf.use <- c(0,leaf.use)
    }
    
    split_C <- C_splits(x = as.numeric(X) ,numperleaf = as.integer(n.k.full),
                        splitleaf = as.integer(leaf.use) ,
                        leafmem = as.integer(leaf.mem.use),
                        smpl = as.integer(smpl),
                        r = as.integer(r),
                        sortind = as.integer(sortind),
                        kern = Kern.C,omegaY = omegaY
    )
    
    
    tot.splits <- split_C[[4]]
    if(tot.splits!=0){
      loss.den <- split_C[[6]][1:tot.splits]
      loss.num  <- (split_C[[5]][1:tot.splits])^2 
      loss <- loss.num/loss.den
      
      if(max(loss,na.rm=T) < 0){
        break
      } else{
        ind <- which.max( loss )
      }
      
      covar.num[a] <- split_C[[2]][ind]
      split.val[a] <- split_C[[1]][ind]
      leaf.split[a] <- leaf.ind[split_C[[3]][ind]]
      
      branch.loss[a] <- loss[ind]
      Ak[,a+1] <- Ck[,a+Ck.ind] <- 
        as.numeric(X[,covar.num[a]] <= split.val[a] & leaf.mem==leaf.split[a])
      
      Ak[,1:a] <- ifelse(Ak[,1:a,drop=F] - Ak[,a+1] < 0 , 
                          0, 
                          Ak[,1:a,drop=F] - Ak[,a+1]) 
      
      Kern.C <- Kern.C - 
        (Kern.C %*% Ck[,a+Ck.ind] %*% crossprod(Ck[,a+Ck.ind], Kern.C)) / 
        c( crossprod(Ck[,a+Ck.ind], Kern.C) %*% Ck[,a+Ck.ind] )
      omegaY <- as.numeric(Kern.C %*% Y)
      
      a <- a + 1

    } else{
      break
    }
  }

  Ck <- Ck[,1:a,drop=FALSE]
  
  if(is.null(X.fix)){
    muk <- solve(t(Ck) %*% i.Sig %*% Ck,t(Ck) %*% i.Sig %*% Y)
  } else{
    X.out <- cbind(X.fix,Ck)
    muk <- solve(t(X.out) %*% i.Sig %*% X.out + 
                   diag(c(rep(lambda,p.fix),rep(0,a))),
                 t(X.out) %*% i.Sig %*% Y)
    b <- muk[1:p.fix]
    muk <- muk[(p.fix+1):length(muk)]
  }
  
  ## Return the tree 
  return(list(covar=covar.num[1:(a-1)],split=split.val[1:(a-1)],
              leaf=leaf.split[1:(a-1)],loss=branch.loss[1:(a-1)],
              Ck=Ck,muk=muk,X=X.fix,b=b,max.loss = loss.full))
}

##########
# spatRF #
##########
# Inputs: pars - Parameters of the covariance function
#         Y - Observations at sites
#         X - Covariates to build Random Forest on
#         coords - Locations of observations
#         t - Number of Trees to aggregate over
#         cov.type - type of covariance function (see .makeSpatBas)
#         cov.opts - options for covariance function (see .makeSpatBas)
#         Xtest - Covariates for test observations
#         coords.test - locations of test sites
#         m - number of covariates to sample for each tree
#         k - max number of leaves to use
#         s - number of observations to use in each tree
#         replace - Sample with/without replacement

spatRF <- function( pars, Y, X, coords, t = 500, cov.type = "exp", 
                            cov.opts = NULL, Xtest = NULL, 
                            coords.test = NULL, m = floor(dim(X)[2]/3),  r = 5, 
                            k = floor(length(Y)/r), replace = FALSE,
                            s = ifelse (replace,length(Y), 
                                        ceiling(.632*length(Y))), 
                            var.imp=FALSE,imp.msr = c("perm","imp"),... ){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  times.selected <- times.trained <- rep(0,n)
  Yib <- fib <- zib <- rep(0,n)
  Yoob <- foob <- zoob <- rep(0,n)
  M <- matrix(0,n,n)
  
  if(cov.type=="TPRS" & !is.null(coords.test)){
    X.vldt <- as.matrix(coords.test)
  } else{
    X.vldt <- NULL
  }
  
  if(!is.null(Xtest)) {
    n.test <- dim(Xtest)[1]
    Ytest <- ftest <- ztest <- rep(0,n.test)
  }
  
  if(var.imp){
    var.imp.msr <- rep(0,p)
    imp.msr <- imp.msr[1]
    if(imp.msr=="perm"){
      perm.ind <- replicate(p,sample(1:n))
      perm.mat <- matrix(0,nrow=n,ncol=p)
    }
  }
  
  ## Build t Spatial trees
  for(i in 1:t){
    
    ## For each tree, subsample the data into training and test
    vldt.sample <- replace
    train <- sample(1:n,s,replace=replace)
    test <- (1:n)[-which(1:n %in% train)]
    Y.train <- Y[train]
    X.train <- X[train,,drop=FALSE]
    
    ## If TPRS, generate the fixed basis functions
    if(cov.type=="TPRS") {
      X.fix <- as.matrix(coords[train,])
      X.fix.test <- as.matrix(coords[test,])
    } else{
      X.fix <- NULL
      X.fix.test <- NULL
    }
    
    ## Generate the training and test covariance matrices
    spatBas <- try( .makeSpatBas(pars,coords.train=coords[train,],
                                coords.test = coords[test,],
                                cov.type=cov.type,cov.opts = cov.opts) )
    if(class(spatBas) != "try-error") vldt.sample <- FALSE
      
    while(vldt.sample){
      train <- sample(1:n,s,replace=replace)
      test <- (1:n)[-which(1:n %in% train)]
      Y.train <- Y[train]
      X.train <- X[train,,drop=FALSE]
      if(cov.type=="TPRS") {
        X.fix <- as.matrix(coords[train,])
        X.fix.test <- as.matrix(coords[test,])
      } else{
        X.fix <- NULL
        X.fix.test <- NULL
      }
        
      spatBas <- try( .makeSpatBas(pars,coords.train=coords[train,],
                                  coords.test = coords[test,],
                                  cov.type=cov.type,cov.opts = cov.opts) )
      if(class(spatBas) != "try-error") vldt.sample <- FALSE
    }
    
    i.Sig <- spatBas$i.sig
    
    ## Run the tree
    mod <- spatTree(Y=Y.train,X=X.train,X.fix = X.fix,
                    i.Sig =spatBas$i.sig,m=m,r=r,k=k)
    
    ## Compile predictions on held out samples
    C.test <- .makeCmat(mod,X[test,,drop=FALSE])
    times.selected[test] <- times.selected[test] + 1
    
    if(cov.type=="TPRS"){
      foob[test] <- foob[test] + C.test %*% mod$muk
      zoob[test] <- zoob[test] + X.fix.test %*% mod$b + 
                      spatBas$sig.test %*% spatBas$i.sig %*% 
                      (Y.train - cbind(mod$Ck,mod$X) %*% c(mod$muk,mod$b))
    } else{
      foob[test] <- foob[test] + C.test %*% mod$muk
      zoob[test] <- zoob[test] + spatBas$sig.test %*% spatBas$i.sig %*% 
                      (Y.train - mod$Ck %*% mod$muk )
    }
    
    ## Generate predictions at held out sites using the current tree
    if( !is.null(Xtest) ) {
      C.vldt <- .makeCmat(mod,Xtest)
      ftest <- ftest + C.vldt %*% mod$muk 
      globalSpatBas <- .makeSpatBas(pars,coords.train=coords[train,],
                                   coords.test = coords.test,
                                   cov.type=cov.type,cov.opts = cov.opts)
      if(cov.type=="TPRS"){
        ztest <- ztest + X.vldt %*% mod$b +
                  globalSpatBas$sig.test %*% globalSpatBas$i.sig %*% 
                  (Y.train - cbind(mod$Ck,mod$X) %*% c(mod$muk,mod$b))
      } else{
        ztest <- ztest + 
                  globalSpatBas$sig.test %*% globalSpatBas$i.sig %*% 
                  (Y.train - mod$Ck %*% mod$muk) 
      }
    }
    
    ## Generate variable importance measures if requested
    if(var.imp){
      if(imp.msr=="imp"){
        loss.explain <- rowsum(mod$loss/c(mod$max.loss),mod$covar)
        var.imp.msr[as.numeric(rownames(loss.explain))] <- 
          var.imp.msr[as.numeric(rownames(loss.explain))] + c(loss.explain)
      } else if(imp.msr=="perm"){
        for (w in unique(mod$covar)){
          X.perm <- X
          X.perm[,w] <- X[perm.ind[,w],w]
          C.perm <- .makeCmat(mod,X.perm[test,,drop=FALSE])
          perm.mat[test,w] <- perm.mat[test,w] + C.perm %*% mod$muk
        }
        notused <- (1:p)[-unique(mod$covar)]
        perm.mat[test,notused] <- perm.mat[test,notused] + c(C.test %*% mod$muk)
      }
    }
    
  }
  
  out <-list(fpredicted=foob/times.selected,zpredicted = zoob/times.selected)
  
  if(!is.null(Xtest)) {
    out$ftest <- ftest/t
    out$ztest <- ztest/t
  }
  
  if(var.imp){
    if(imp.msr=="perm"){
      for(w in 1:p){
        var.imp.msr[w] <- (sum((Y-(perm.mat[,w]+zoob)/times.selected)^2) - 
                             sum((Y-(foob + zoob)/times.selected)^2) )/n
      }
    }
    out$var.imp <- var.imp.msr
  }
  
  return(out)
}

############
# cvSpatRF #
############
# Inputs: Y - Observations at sites
#         X - Covariates to build Random Forest on
#         coords - Locations of observations
#         cv.lambda - "Penalty" parameters to optimize over
#         cov.type - type of covariance function (see .makeSpatBas)
#         cov.opts - options for covariance function (see .makeSpatBas)
#         lklhd - returns pseudo-likelihood estimates in addition to CV
#         ... - additional params to pass to  spatRF
cvSpatRF <- function(Y,X,coords,cov.type="exp", cov.opts = list(range=NULL),
                     cv.lambda = seq(0.05,0.95,0.05),
                     random.tune = FALSE, nums = 10, # random search for tuning
                     range.r = c(2, 12), range.m = c(1.5, 6),
                     range.s = c(0.3, 0.9),
                     cv.r = seq(2, 10, 2), cv.m = floor(dim(X)[2]/c(2, 3, 4)),
                     cv.s = ceiling(length(Y)*seq(0.4, 0.8, 0.1)),
                     lklhd=FALSE, ...) {
  
  n <- length(Y)
  p <- dim(X)[2]
  if (random.tune){
    cv.lambda = rep(cv.lambda, nums)
    cv.r = sample(range.r[1]:range.r[2], size = nums, replace = TRUE)
    cv.m = floor(dim(X)[2]/runif(nums, min = range.m[1], max = range.m[2]))
    cv.s = ceiling(length(Y)*runif(nums, min = range.s[1], max = range.s[2]))
    params = cbind(cv.lambda, cv.r, cv.m, cv.s)
  } else{
    params = expand.grid(cv.lambda, cv.r, cv.m, cv.s)
  }
  m = nrow(params)
  #m <- length(cv.lambda)
  MSE <- rep(NA, m)
  #names(MSE) <- cv.lambda
  y.full <- matrix(NA,n,m)
  #colnames(y.full) <- cv.lambda
  lklhd.vals <- MSE
  bestRF <- bestLklhd <- NULL
  
  if(lklhd) {
    spatBas <- .makeSpatBas(pars=list(psill=1,nugget=1),
                          coords.train=coords, out.Z = TRUE,
                          cov.type=cov.type,cov.opts = cov.opts)
    R <- spatBas$Zb %*% t(spatBas$Zb)
  }

  for(j in 1:m){
    lambda = params[j,1]
    r = params[j,2]
    m_ = params[j,3]
    s = params[j,4]
    spatRF.obj <- try(spatRF( pars=list(psill=lambda,nugget=1-lambda), Y=Y, X=X, 
                              coords=coords, cov.type = cov.type, 
                              cov.opts = cov.opts, 
                              r = r, m = m_, s = s, ...))
    if (class(spatRF.obj) != "try-error"){
      y.test <- spatRF.obj$fpredicted + spatRF.obj$zpredicted
      cur.mse <- sum((Y-y.test)^2)/n
    } else{
      y.test = rep(NA, nrow(y.full))
      cur.mse = Inf
    }
    
    if(is.null(bestRF) || all(is.na(MSE))){
      bestRF <- spatRF.obj
    } else{
      if( cur.mse < min(MSE,na.rm=T)){
        bestRF <- spatRF.obj
      }
    }
    
    #MSE[which(names(MSE)==lambda)] <- cur.mse
    #y.full[,which(colnames(y.full)==lambda)] <- y.test
    MSE[j] = cur.mse
    y.full[,j] = y.test
    
    if(lklhd){
      V <- chol(lambda * R +  diag(rep(1-lambda,n)))
      cur.lklhd <- 2 * sum(log(diag(V))) + 
        n * log( t(Y-spatRF.obj$fpredicted) %*% 
                   chol2inv(V) %*% (Y-spatRF.obj$fpredicted) / n ) 
      
      if(is.null(bestLklhd) || all(is.na(lklhd.vals))){
        bestLklhd <- spatRF.obj
      } else{
        if( cur.lklhd < min(MSE,na.rm=T)){
          bestLklhd <- spatRF.obj
        }
      }
      #klhd.vals[which(names(lklhd.vals)==lambda)] <-  cur.lklhd
      lklhd.vals[j] <-  cur.lklhd
    }
    
  }
    
  out <- list(mse=MSE, y.cv = y.full[,which.min(MSE)],
              #lambda.min=as.numeric(names(MSE)[which.min(MSE)]),
              lambda.min=params[which.min(MSE),1],
              r.min=params[which.min(MSE),2],
              m.min=params[which.min(MSE),3],
              s.min=params[which.min(MSE),4],
              spatRF = bestRF)
  if(lklhd) {
    out$lklhd <- lklhd.vals
    out$y.lklhd <- y.full[,which.min(lklhd.vals)]
    #out$lambda.lklhd <- as.numeric(names(lklhd.vals)[which.min(lklhd.vals)])
    out$lambda.lklhd=params[which.min(lklhd.vals),1]
    out$r.lklhd=params[which.min(lklhd.vals),2]
    out$m.lklhd=params[which.min(lklhd.vals),3]
    out$s.lklhd=params[which.min(lklhd.vals),4]
    out$spatRF.lklhd <- bestLklhd
  }
  
  return(out)
}

#######################
# Auxiliary Functions #
#######################
.fastSolve <- function(mat){
  out <- try(chol2inv(chol(mat)),silent=TRUE)
  if(class(out)=="try-error") {
    out <- try(solve(mat),silent=TRUE)
    if(class(out)=="try-error"){
      out <- svd(mat)
      out <- out$u %*% diag(1/out$d) %*% t(out$v)
    }
  }
  out
}

.makeCmat <- function(mod,X){
  num.branches <- length(mod$covar)
  n <- dim(X)[1]
  Ak <- Ck <- matrix(NA,n,num.branches+1)
  Ak[,1] <- Ck[,1] <- 1
  leaf.mem <- rep(1,n)
  for (a in 1:num.branches){
    
    Ak[,a+1] <- Ck[,a+1] <- as.numeric(X[,mod$covar[a]] <= mod$split[a] & 
                                         leaf.mem==mod$leaf[a])
    Ak[,1:a] <- ifelse( Ak[,1:a,drop=F] - Ak[,a+1] < 0 , 0, Ak[,1:a,drop=F] - Ak[,a+1] ) 
    leaf.mem <- Ak[,1:(a+1)] %*%  1:(a+1)
  }
  return(Ck)
}

.makeSpatBas <- function(pars,coords.train,coords.test=NULL, cov.type, 
                         cov.opts=list(NULL),out.Z = FALSE){
  n <- dim(coords.train)[1]
  n.test <- dim(coords.test)[1]
  if(pars$psill!=0){
    if(cov.type == "exp"){
      ## Full Rank Kriging with Specified Range
      train.dist.mat <- as.matrix(SpatioTemporal::crossDist(coords.train))
      Sigma <-  exp(-train.dist.mat/cov.opts$range) 
      i.sig <- .fastSolve(pars$psill * Sigma + pars$nugget * diag(n))
      if(!is.null(coords.test)){
        test.dist.mat <- as.matrix(SpatioTemporal::crossDist(coords.test,coords.train))
        Sig.test <- pars$psill * exp(-test.dist.mat/cov.opts$range)
      }
      if(out.Z){
        Zb <- .sqrtMat(Sigma)
      } 
    } else if ( cov.type == "LRK" ){
      ## Low Rank Kriging with Specified Range
      train.dist.mat <- as.matrix(SpatioTemporal::crossDist(coords.train,cov.opts$knots))
      Z1 <- exp(-train.dist.mat/cov.opts$range)
      Om <- exp(-SpatioTemporal::crossDist(cov.opts$knots)/cov.opts$range)
      i.sqrt.Om <- .iSqrtMat(Om,by.chol=TRUE)
      Zb <- Z1 %*% i.sqrt.Om
      i.sig <- .blockSolve(pars,Zb)
      if(!is.null(coords.test)){
        test.dist.mat <- as.matrix(SpatioTemporal::crossDist(coords.test,cov.opts$knots))
        Z2 <- exp(-test.dist.mat/cov.opts$range)
        Sig.test <- pars$psill*Z2 %*% i.sqrt.Om %*% t(Z1 %*% i.sqrt.Om)
      }
    } else if(cov.type == "TPRS"){
    
      if(cov.opts$k==0){
        cov.opts$k <- n
      }
      K <- cov.opts$k
      if(is.null(cov.opts$m)) cov.opts$m <- 2
      m <- cov.opts$m
      ## Thin Plate Regression Splines
      Xbas <- cbind(1, coords.train)
      if(m > 2){
        for(jj in 2:(m-1)){
          temp.T <- coords.train^jj
          Xbas <- cbind(Xbas, temp.T)
          rm(temp.T)
        }
      }
      Xbas <- as.matrix(Xbas)
      M <- choose(m+2-1, 2) #d=2
      E <- .etaFunc(r=crossDist(coords.train), d=2, m=m)
      Ek <- eigen(E)
      Uk <- Ek$vectors[, 1:K]
      Dk <- diag(Ek$values[1:K])
      tUkT <- t(Uk)%*%Xbas
      QR <- qr(tUkT)
      Zk <- qr.Q(QR, complete=TRUE)[, (M+1):K]
      O <- t(Zk)%*%Dk%*%Zk
      O <- svd(O)
      O$d <- sqrt(O$d)
      O <- solve(O$u%*%diag(O$d)%*%t(O$v))
      Zb <- Uk%*%Dk%*%Zk%*%O
      i.sig <- .blockSolve(pars,Zb)
      if(!is.null(coords.test)){
        E1 <- .etaFunc(r=SpatioTemporal::crossDist(coords.test, coords.train), d=2, m=m)
        E2 <- .etaFunc(r=SpatioTemporal::crossDist(coords.train, coords.train), d=2, m=m)
        Sig.test <- pars$psill * E1%*%Uk%*%Zk%*%O%*%O%*%t(Zk)%*%t(Uk)%*%t(E2)
      }
    }
  } else{
    i.sig <- diag(rep(1/pars$nugget,n))
    if(!is.null(coords.test)){
      Sig.test <- matrix(0,n.test,n)
    }
  }
  out <- list(i.sig = i.sig)
  if(!is.null(coords.test)){
    out$sig.test <- Sig.test
  }
  if(out.Z){
    out$Zb <- Zb
  }
  
  return(out)
}

.etaFunc <- function(m, d, r){
  if(d%%2 == 0){
    out <- (-1)^(m+1+d/2)/
      (2^(2*m-1)*pi^(d/2)*factorial(m-1)*factorial(m-d/2))*r^(2*m-d)*log(r)
  } else {
    out <- gamma(d/2-m)/(2^(2*m)*pi^(d/2)*factorial(m-1))*r^(2*m-d)
  }
  
  out[is.na(out)] <- 0
  return(out)
}

.iSqrtMat <- function(mat,by.chol=TRUE){
  if(by.chol){
    out <- try(chol(mat),silent=TRUE)
    if(class(out)=="try-error") {
      by.chol=FALSE
    } else{
      out <- (chol(chol2inv(out)))
    }
  }
  if(!by.chol){
    out <- svd(mat)
    out <- out$u %*% diag(1/sqrt(out$d)) %*% t(out$v)
  } 
  out
}

.sqrtMat <- function(mat,by.chol=TRUE){
  if(by.chol){
    out <- try(t(chol(mat)), silent = TRUE)
    if(class(out)=="try-error"){
      by.chol = FALSE
    }
  }
  if(!by.chol){
    out <- svd(mat)
    out <- out$u %*% diag(sqrt(out$d)) %*% t(out$v)
  }
  out
}

.blockSolve <- function(pars,Z){
  n <- dim(Z)[1]
  k <- dim(Z)[2]
  out <- try(diag(n) - Z %*% 
               .fastSolve (t(Z) %*% Z  + pars$nugget/pars$psill* diag(k)) %*% 
               t(Z)) / pars$nugget
  if(class(out) == "try-error"){
    out <- .fastSolve(pars$psill* Z %*% t(Z) + pars$nugget * diag(n))
  }
  out
}

############
# C_splits #
############
require(inline)
C_splits <- cfunction(c(x = "numeric",numperleaf="integer",
                        splitleaf = "integer",leafmem="integer",smpl="integer",
                        r="integer",sortind = "integer",kern = "numeric",
                        omegaY = "numeric"),"
                      int k = length(splitleaf);
                      int n = length(leafmem);
                      int s = length(smpl);
                      int numsplits=0;
                      int strtind[k];
                      int obsind,covarnum,leafnum,obsleaf,xnum,leafind,xnumnext,ynum,yind,yindold;
                      int nummaxsplit,yold,leaforder;
                      int count[k];
                      int strtj = 0;
                      int lmadjust = 1;
                      double runningnum, runningden;
                      
                      int *pnpl, *psl,*plm,*psmpl,*pr,*pcovar,*pleaf,*psortind,*pcli,*pns;
                      double *pctpt,*px,*plossnum,*plossden,*pkern,*poy ;
                      
                      pnpl = INTEGER(numperleaf);
                      psl = INTEGER(splitleaf);
                      plm = INTEGER(leafmem);
                      psmpl = INTEGER(smpl);
                      pr = INTEGER(r);
                      px = REAL(x);
                      psortind = INTEGER(sortind);
                      pkern = REAL(kern);
                      poy = REAL(omegaY);
                      
                      nummaxsplit = s* n; 
                      SEXP ctpt = PROTECT(allocVector(REALSXP, nummaxsplit));
                      SEXP covar = PROTECT(allocVector(INTSXP, nummaxsplit));
                      SEXP leaf = PROTECT(allocVector(INTSXP, nummaxsplit));
                      SEXP covleafind = PROTECT(allocVector(INTSXP, n));
                      SEXP vec = PROTECT(allocVector(VECSXP, 6));
                      SEXP ns = PROTECT(allocVector(INTSXP, 1));
                      SEXP lossnum = PROTECT(allocVector(REALSXP, nummaxsplit));
                      SEXP lossden = PROTECT(allocVector(REALSXP, nummaxsplit));
                      
                      pcli = INTEGER(covleafind);
                      pctpt = REAL(ctpt);
                      pcovar = INTEGER(covar);
                      pleaf = INTEGER(leaf);
                      pns = INTEGER(ns);
                      plossnum = REAL(lossnum);
                      plossden = REAL(lossden);
                      
                      // Set up counter indices for each leaf
                      strtind[0] = 0;
                      for(int t = 1; t < k; t++){
                      strtind[t] = strtind[t-1] + pnpl[t-1];
                      }
                      
                      // Check if there are leaves that shouldn't be split
                      if(psl[0] ==0){
                      strtj =1;
                      lmadjust = 0;
                      } 
                      
                      // For Each Covariate in the sample
                      for ( int i = 0; i < s; i++ ) {
                      covarnum = psmpl[i]-1;    
                      for(int d=0; d<k; d++){
                      count[d] = 0;
                      }
                      for(int check =0; check < n; check++){
                      pcli[check]=check;
                      }
                      
                      // Sort covariate first by leaf, then in order of value
                      for( int a=0; a < n; a++){
                      obsind = psortind[ a + n * covarnum ]-1;
                      obsleaf = plm[ obsind] - lmadjust;
                      leaforder = strtind[ obsleaf ] + count[ obsleaf ];
                      pcli[leaforder] = obsind;
                      count[obsleaf]+=1;
                      }
                      // covleafind now contains the order of observations, first sorted by leaf, then by covariate
                      
                      // Initialize the leaf counter
                      for( int j = strtj; j < k; j++ ){
                      leafnum = psl[j];
                      runningnum = 0;
                      runningden = 0;
                      
                      if( (*pr-1) > 0 ){
                      for(int ci = 0; ci < (*pr-1); ci++){
                      yind = pcli[ strtind[j] + ci ];
                      runningnum += poy[ yind ];
                      if( ci > 0){
                      for( int cii = 0; cii < ci; cii++){
                      yindold = pcli[ strtind[j] + cii ];
                      runningden += 2 * pkern[n * yind + yindold];
                      }
                      }
                      runningden += pkern[n * yind + yind];
                      
                      }
                      }
                      
                      // For each observation in the leaf
                      for( int b = *pr-1; b < (pnpl[j] - *pr); b++){
                      
                      // Get index of which observation is next lowest               
                      leafind = strtind[j] + b;
                      ynum = pcli[ leafind ];
                      xnum = n*covarnum + pcli[ leafind ];
                      xnumnext = n*covarnum + pcli[ leafind +1 ];
                      
                      // Update running numerator and denominator total  
                      runningnum += poy[ynum];
                      
                      if(b > 0){
                      for( int bi = 0; bi < b; bi++){
                      yold = pcli[ strtind[j] + bi ];
                      runningden += 2 * pkern[n * ynum + yold];
                      }
                      }
                      runningden +=  pkern[n * ynum + ynum];
                      
                      if(px[ xnum ] - px[ xnumnext] < 0){
                      pctpt[numsplits] = (px[ xnum ] + px[ xnumnext])/2;
                      pleaf[numsplits] = leafnum;
                      pcovar[numsplits] = covarnum+1;
                      plossnum[numsplits] = runningnum;
                      plossden[numsplits] = runningden;
                      numsplits += 1;
                      }
                      
                      // End loop for each observation in the leaf
                      } 
                      
                      // End loop for each leaf  
                      }
                      
                      // End loop for each covariate
                      } 
                      
                      pns[0] = numsplits;
                      SET_VECTOR_ELT(vec, 0, ctpt);
                      SET_VECTOR_ELT(vec, 1, covar);
                      SET_VECTOR_ELT(vec, 2, leaf);
                      SET_VECTOR_ELT(vec, 3, ns);    
                      SET_VECTOR_ELT(vec, 4, lossnum);
                      SET_VECTOR_ELT(vec, 5, lossden); 
                      UNPROTECT(8);
                      
                      return vec;
                      
                      ")

