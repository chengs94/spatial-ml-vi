# variable importance

library(SpatioTemporal)
library(matrixStats)
library(inline)
#library(parallel)

setwd("/home/students/chengsi/Desktop/exposurepred/")
#load("trap_cv_morepars.RData")
load("trap_cv_morepars_finer.RData")
load("trap_cv_morepars_finer2.RData")
source("aux_functions.R")
source("spatTreeModified.R")

rslt.all = list(rslt.ufp.primary, rslt.ufp.uw,
                rslt.bc.primary, rslt.bc.uw)

## summarize best parameters
bestNumComp = rep(0, length(rslt.all))
bestPLpar = matrix(0, ncol = 4, nrow = length(rslt.all))
bestCVpar = matrix(0, ncol = 4, nrow = length(rslt.all))
colnames(bestSpatrfPar) = c("lambda", "r", "m", "s")

#set.seed(2333)
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

set.seed(2333)

mod.ufp.primary = spatRF( pars=list(psill=bestPLpar[1,1],nugget=1-bestPLpar[1,1]), 
                          Y=dat$ufp_primary, X=as.matrix(X), 
                          Xtest = NULL, coords.test = NULL,
                          coords=cbind(grid$lambert_x, grid$lambert_y), 
                          cov.type = "TPRS", cov.opts=list(k=0,m=2),
                          r = bestPLpar[1,2], m = bestPLpar[1,3], s = bestPLpar[1,4],
                          replace=TRUE, t = 500, 
                          var.imp=TRUE,imp.msr = "imp")

mod.ufp.bc = spatRF( pars=list(psill=bestPLpar[2,1],nugget=1-bestPLpar[2,1]), 
                          Y=dat$ufp_uw, X=as.matrix(X), 
                          Xtest = NULL, coords.test = NULL,
                          coords=cbind(grid$lambert_x, grid$lambert_y), 
                          cov.type = "TPRS", cov.opts=list(k=0,m=2),
                          r = bestPLpar[2,2], m = bestPLpar[2,3], s = bestPLpar[2,4],
                          replace=TRUE, t = 500, 
                          var.imp=TRUE,imp.msr = "imp")

#names(rslt.ufp.primary[[1]]$spatrf.cv.obj)
getImp = function(num, name){
  # which pollutant; name in dataset
  #coords = cbind(grid$lambert_x, grid$lambert_y)
  loadings.pls = pls::plsr(dat[,which(colnames(dat) == name)] ~ as.matrix(X), 
                          ncomp = bestNumComp[num],scale=TRUE)$loadings[]
    #get.pcs(as.matrix(X), "pls", bestNumComp[num], 
    #                     dat[,which(colnames(dat) == name)])$loadings
  imp.spatrf = matrix(0, nrow = ncol(X), ncol = 2)
  colnames(imp.spatrf) = c("pl", "cv")
  for (j in 1:numfold){
    imp.spatrf[,2] = imp.spatrf[,2] + rslt.all[[num]][[j]]$spatrf.cv.obj$spatRF$var.imp
    imp.spatrf[,1] = imp.spatrf[,1] + rslt.all[[num]][[j]]$spatrf.cv.obj$spatRF.lklhd$var.imp
  }
  imp.spatrf = imp.spatrf/numfold
  return(list(loadings.pls = loadings.pls, imp.spatrf = imp.spatrf))
}

imp = mapply(FUN = getImp,
             num = 1:4, name = c("ufp_primary", "ufp_uw", "bc_primary", "bc_uw"))

## plot
my.alpha=0.3

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

plotImp = function(impMat, method, pollutant){
  if (method == "pls"){
    rownames(impMat) = colnames(X)
    dtNew <- impMat %>%
      as.data.frame() %>%
      rownames_to_column(var = "cov") %>%
      # rename variables if buffers
      split_cov_name(cov = "cov") %>%
      #make long format for faceting
      gather(key = "Component", value = "Loading", contains("Comp")) %>%
      mutate(Component = as.numeric(substr(Component, 6, nchar(Component)))) 
    
    dtNew  %>%
      #buffered covariates
      drop_na(buffer) %>%
      ggplot(aes(x = Loading, y = cov)) +
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
      facet_wrap(~Component,
                 labeller = "label_both",
                 nrow = 1) +
      labs(y = "Geocovariate",
           shape= "non-buffer", #"proximity,\nelevation",
           title = paste("PLS Geocovariate Component Loadings:", pollutant)) +
      theme(legend.position = "bottom")
    
  } else{
    colnames(impMat) = c("Method PL", "Method CV")
    rownames(impMat) = colnames(X)
    dtNew <- impMat %>%
      as.data.frame() %>%
      rownames_to_column(var = "cov") %>%
      # rename variables if buffers
      split_cov_name(cov = "cov") %>%
      #make long format for faceting
      gather(key = "Method", value = "Importance", contains("Method")) %>%
      mutate(Method = substr(Method, 8, nchar(Method))) 
    
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
      facet_wrap(~ Method,
                 labeller = "label_both",
                 nrow = 1) +
      labs(y = "Geocovariate",
           shape= "non-buffer", #"proximity,\nelevation",
           title = paste("Spatial RF Geocovariate Importance:", pollutant)) +
      theme(legend.position = "bottom")
  }
}

# UFP
plotImp(imp[1,1][[1]], "pls", "UFP (Primary)")
plotImp(imp[2,1][[1]], "spatrf", "UFP (Primary)")
plotImp(imp[1,2][[1]], "pls", "UFP (Unweighted)")
plotImp(imp[2,2][[1]], "spatrf", "UFP (Unweighted)")

# BC
plotImp(imp[1,3][[1]], "pls", "BC (Primary)")
plotImp(imp[2,3][[1]], "spatrf", "BC (Primary)")
plotImp(imp[1,4][[1]], "pls", "BC (Unweighted)")
plotImp(imp[2,4][[1]], "spatrf", "BC (Unweighted)")

#plotImp(impMat = impMat, method = "spatrf", pollutant = "UFP primary")

#bestModel = function(num, name){
#  coords = cbind(grid$lambert_x, grid$lambert_y)
#  pc.best = get.pcs(as.matrix(X), "pls", bestNumComp[num], dat[,which(colnames(dat) == name)])
#  lambda.pl = bestPLpar[num,1]
#  r.pl = bestPLpar[num,2]
#  m.pl = bestPLpar[num,3]
#  s.pl = bestPLpar[num,4]
#  spatrf.pl.best = cvSpatRF(Y=dat[,which(colnames(dat) == name)], 
#                            X=as.matrix(X),coords=coords,
#                            cv.lambda=lambda.pl, cv.r = r.pl,
#                            cv.s = s.pl, cv.m = m.pl,
#                            random.tune = FALSE,
#                            cov.type="TPRS",cov.opts=list(k=0,m=2),lklhd=FALSE,
#                            Xtest = as.matrix(X), replace=TRUE,
#                            coords.test = coords, t = 500, 
#                            var.imp=TRUE,imp.msr = "imp")
#  spatrf.pl.best = spatRF( pars=list(psill=lambda.pl, nugget=1-lambda.pl), 
#                           Y=dat[,which(colnames(dat) == name)], 
#                           X=as.matrix(X), coords = coords, t = 500, 
#                           cov.type = "TPRS", cov.opts=list(k=0,m=2), 
#                           Xtest = NULL, coords.test = NULL, 
#                           m = m.pl,  r = r.pl, 
#                           k = floor(nrow(X)/r.pl), replace = FALSE,
#                           s = s.pl, 
#                           var.imp=TRUE, imp.msr = "imp")
#}
