# importance measure for PLS based on permutation

library(parallel)

setwd("/home/students/chengsi/Desktop/exposurepred/")
load("trap_cv_morepars_finer.RData")
load("trap_cv_morepars_finer2.RData")

Y.pls <- rep(NA, n)
set.seed(2333)

predPerm = function(fold, annavg, cv.rslt, id.x){
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
  X.test[,id.x] = sample(X.test[,id.x], n.test, replace = FALSE)
  
  X.train.pc = X.train
  X.test.pc = X.test
  
  dmat <- as.matrix(SpatioTemporal::crossDist(train.coords))
  dtestmat <- as.matrix(SpatioTemporal::crossDist(test.coords,train.coords))
  
  pc.obj = cv.rslt[[fold]]$pc.obj
  
  new.x.train <- pc.obj$X
  means<-apply(X.train.pc,2,mean)
  mean.matrix<-matrix(rep(means,n.test),nrow=n.test,byrow=T)
  sds <- apply(X.train.pc,2,sd)
  new.x.test <- (X.test.pc - mean.matrix) %*% diag(1/sds) %*% pc.obj$proj
  
  theta.pls = cv.rslt[[fold]]$theta.pls
  
  sig.pls <- exp(theta.pls[1])
  nug.pls <- exp(theta.pls[2])
  range.pls <- exp(theta.pls[3])
  
  i.V <- .fastSolve(sig.pls * exp(-dmat/range.pls) + 
                      nug.pls * diag(length(Y.train)))
  beta.hat.pcr = cv.rslt[[fold]]$beta.hat.pls
  #beta.hat.pcr <- solve(t(new.x.train) %*% i.V %*% new.x.train,
  #                      t(new.x.train) %*% i.V %*% Y.train)
  Y.pls[cv.ind.i]<- c(new.x.test %*% beta.hat.pcr) + 
    sig.pls * exp(-dtestmat/range.pls) %*% i.V %*% 
    c(Y.train-c(new.x.train %*% beta.hat.pcr))
  
  return(Y.pls)
}

numCores=detectCores()

pars = expand.grid(1:numfold, 1:ncol(X))

imp.ufp.uw = mcmapply(fold = pars[,1], id.x = pars[,2],
                      FUN = predPerm, mc.cores = numCores,
                      MoreArgs = list(annavg = dat$ufp_uw, cv.rslt = rslt.ufp.uw))
imp.bc.uw = mcmapply(fold = pars[,1], id.x = pars[,2],
                     FUN = predPerm, mc.cores = numCores,
                     MoreArgs = list(annavg = dat$bc_uw, cv.rslt = rslt.bc.uw))
save(imp.ufp.uw, imp.bc.uw, file = "pls_perm.RData")

imp.ufp.primary = mcmapply(fold = pars[,1], id.x = pars[,2],
                           FUN = predPerm, mc.cores = numCores,
                           MoreArgs = list(annavg = dat$ufp_primary, cv.rslt = rslt.ufp.primary))
imp.bc.primary = mcmapply(fold = pars[,1], id.x = pars[,2],
                          FUN = predPerm, mc.cores = numCores,
                          MoreArgs = list(annavg = dat$bc_primary, cv.rslt = rslt.bc.primary))
save(imp.ufp.uw, imp.bc.uw, imp.ufp.primary, imp.bc.primary,
     file = "pls_perm.RData")

## summarizing results to get importance
library(ggplot2)
library(gridExtra)

getImp = function(rslt.perm, rslt.cv, y.true){
  pars = data.frame(fold = pars[,1], id.x = pars[,2],
                    mse.cv = NA, mse.perm = NA, imp = NA)
  for (j in 1:ncol(rslt.perm)){
    cv.ind = which(!is.na(rslt.perm[,j]))
    pars$mse.cv[j] = mean((y.true[cv.ind] - rslt.cv[[pars$fold[j]]]$Y.pls[cv.ind])^2)
    pars$mse.perm[j] = mean((y.true[cv.ind] - rslt.perm[cv.ind, j])^2)
    pars$imp[j] = pars$mse.perm[j] - pars$mse.cv[j]
  }
  imp = aggregate(pars$imp, list(pars$id.x), mean)$x
  return(imp)
}

perm.ufp.primary = getImp(imp.ufp.primary, rslt.ufp.primary, dat$ufp_primary)
perm.ufp.uw = getImp(imp.ufp.uw, rslt.ufp.uw, dat$ufp_uw)
perm.bc.primary = getImp(imp.bc.primary, rslt.bc.primary, dat$bc_primary)
perm.bc.uw = getImp(imp.bc.uw, rslt.bc.uw, dat$bc_uw)

plotImp = function(impMat, pollutant){
  my.alpha = .3
  
  split_cov_name <- function(dt, cov) {
    dt <- suppressWarnings(dt %>%
                             rename(cov_full_name = cov) %>%
                             mutate(
                               buffer = substr(cov_full_name, nchar(cov_full_name)-4, nchar(cov_full_name)),
                               buffer = as.numeric(ifelse(!is.na(as.integer(buffer)), buffer, NA)),
                               cov = ifelse(is.na(buffer), cov_full_name, substr(cov_full_name, 1, nchar(cov_full_name)-5) )
                             )
    ) %>%
      select(contains("cov"), buffer, everything())
    
    # elevation
    dt$cov[grepl("^elev_.+_above$", dt$cov_full_name)] <- "elev_above"
    dt$cov[grepl("^elev_.+_below$", dt$cov_full_name)] <- "elev_below"
    dt$cov[grepl("^elev_.+_stdev$", dt$cov_full_name)] <- "elev_stdev"
    dt$cov[grepl("^elev_.+_at_elev$", dt$cov_full_name)] <- "elev_at_elev"
    
    dt$buffer[grepl("_1k_", dt$cov_full_name)] <- 1000
    dt$buffer[grepl("_5k_", dt$cov_full_name)] <- 5000
    
    return(dt)
  }
  
  impMat = as.matrix(impMat)
  rownames(impMat) = colnames(X)
  dtNew <- impMat %>%
    as.data.frame() %>%
    rownames_to_column(var = "cov") %>%
    # rename variables if buffers
    split_cov_name(cov = "cov")
  colnames(dtNew)[4] = "Importance"
  
  dtNew  %>%
    #buffered covariates
    drop_na(buffer) %>%
    ggplot(aes(x = Importance, y = cov)) +
    geom_point(aes(size=buffer),
               shape=1,
               alpha=my.alpha) +
    scale_size(breaks = c(min(dtNew$buffer, na.rm = T),
                          max(dtNew$buffer, na.rm = T)
    )) + #500, 5000, 10000,
    #non-buffered covariates
    geom_point(data = dtNew[is.na(dtNew$buffer),],
               alpha=my.alpha,
               aes(shape="")) +
    geom_vline(xintercept=0,
               linetype="solid",
               alpha=my.alpha) +
    labs(y = "Geocovariate",
         shape= "non-buffer", #"proximity,\nelevation",
         #title = paste("PLS Importance:", pollutant)) +
         title = pollutant) +
    theme(legend.position = "bottom", text = element_text(size = 16))
}

p.ufp.primary = plotImp(perm.ufp.primary, "UFP (Adjusted)")
p.ufp.uw = plotImp(perm.ufp.uw, "UFP (Unweighted)")
p.bc.primary = plotImp(perm.bc.primary, "BC (Adjusted)")
p.bc.uw = plotImp(perm.bc.uw, "BC (Unweighted)")

grid.arrange(arrangeGrob(p.ufp.primary, p.ufp.uw, p.bc.primary, p.bc.uw, nrow = 1), 
             top = textGrob("PLS importance based on permutation",
                            gp=gpar(fontsize=20,font=3)))
