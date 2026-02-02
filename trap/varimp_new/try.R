library(SpatioTemporal)
library(matrixStats)
library(pls)
library(inline)
library(randomForest)
# library(tidyverse)
library(parallel)

# load CV result to get optimal tuning parameters
load("trap_cv_morepars_finer.RData")
load("trap_cv_morepars_finer2.RData")
#cov_grid = readRDS("cov_grid.rda")
cov_cohort = readRDS("cov_cohort.rda")
#source("spatTreeModifiedMean.R")
source("spatTreeModified.R")

numCores = detectCores() - 2

# project_crs <- 4326
# study_area_shp <- read_sf("oval_shp/oval_around_monitoring_area.shp") %>% 
#   st_transform(project_crs)
# 
# grid_covars_shp <- cov_grid %>%
#   st_as_sf(., coords=c("longitude","latitude"), remove = FALSE,
#            crs=project_crs)
# grid_covars_shp$in_study_area <- st_intersects(grid_covars_shp, study_area_shp,
#                                                sparse = FALSE) %>%
#   apply(., 1, any)

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
X.test = as.matrix(X.test[names(X)])

# X.test = X.test[grid_covars_shp$in_study_area,]
idx_subset = sample(1:nrow(X.test), 350, replace = FALSE)
X.test = X.test[idx_subset,]

col_subset = c(1, 19, 21, 33, 37, 82, 135, 144, 153, 164, 174, 165)
qts = c(.25, .5, .75)
k = 1
X_list = list()
for (i in 1:length(col_subset)){
  for (j in 1:length(qts)){
    X_list[[k]] = as.matrix(X.test)
    X_list[[k]][,col_subset[i]] = quantile(X_list[[k]][,col_subset[i]], qts[j])
    k = k + 1
  }
}

annavg = dat$ufp_uw
numComp = bestNumComp[2]
PLpar = bestPLpar[2,]
CVpar = bestCVpar[2,]
X.train = as.matrix(X)

Y.train = annavg
n.test = nrow(X.test)
n.train = nrow(X)
#spatRF.par = rbind(PLpar, CVpar)

train.coords <- cbind(dat$lambert_x, dat$lambert_y)
scale.params <- cbind(apply(train.coords,2,mean),apply(train.coords,2,sd))
x1train <- train.coords[,1] <- (train.coords[,1] - scale.params[1,1])/scale.params[1,2]
x2train <-train.coords[,2] <- (train.coords[,2] - scale.params[2,1])/scale.params[2,2]
# test.coords <- cbind(cov_grid$lambert_x, cov_grid$lambert_y)
test.coords <- cbind(cov_cohort[idx_subset,]$lambert_x, 
                     cov_cohort[idx_subset,]$lambert_y)
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
