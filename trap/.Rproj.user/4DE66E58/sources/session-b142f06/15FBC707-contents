# load("mm_ufp_bc.RData")
# load("mm_co2.RData")
# load("mm_no2.RData")
# load("mm_pm25.RData")
load("mm_wind_all.RData")

rslt.all = list(rslt.ufp.uw, rslt.bc.uw, rslt.co2.uw, rslt.no2.uw, rslt.pm25.uw)

bestNumComp = rep(0, length(rslt.all))
bestPLpar = matrix(0, ncol = 4, nrow = length(rslt.all))
bestCVpar = matrix(0, ncol = 4, nrow = length(rslt.all))

set.seed(2333)

numfold = length(rslt.ufp.uw)
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

temp = rep(NA, length(rslt.ufp.uw[[1]]$Y.spatrf.pl))
for (j in 1:length(rslt.ufp.uw)){
  .idx = which(!is.na(rslt.ufp.uw[[j]]$Y.spatrf.pl))
  temp[.idx] = rslt.ufp.uw[[j]]$Y.spatrf.pl[.idx]
}

save(bestNumComp, bestPLpar, bestCVpar, file = "mm_pars_wind.RData")