# refined variable importance plots

setwd("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/trap/")

dat = readRDS("dat.rda")
X = dat[,-(1:21)]

load("varimp_new/varimp_cohort_diff.RData")

cols = 1:ncol(X)
qs = c(.125, .25, .375, .5, .625, .75, .875)

avg.pls = avg.spatrf = avg.rf = avg.rf.tprs = matrix(NA, ncol = length(qs), nrow = ncol(X))
k = 1

for (i in 1:length(cols)){
  for (j in 1:length(qs)){
    avg.pls[i,j] = mean(pred.diff[[k]]$Y.pls)
    avg.spatrf[i,j] = mean(pred.diff[[k]]$Y.spatrf.pl)
    avg.rf[i,j] = mean(pred.diff[[k]]$Y.rf)
    avg.rf.tprs[i,j] = mean(pred.diff[[k]]$Y.rf.tprs)
    k = k + 1
  }
}

plotImp = function(impMat1, impMat2, impMat3, tt){
  # my.alpha = .3
  my.alpha = 1
  
  # standardize
  # impMat1.std = scale(impMat1) #%>% as.data.frame()
  # impMat2.std = scale(impMat2) #%>% as.data.frame()
  # impMat3.std = scale(impMat3) #%>% as.data.frame()
  # impMat3$diff = abs(impMat3[,1] - impMat3[,2])
  
  colnames(impMat1) = colnames(impMat2) = colnames(impMat3) =
    c("Method PLS", "Method SpatRF-PL", "Method RF", "Method RF-TPRS")
  impMat1 = impMat1[,c(1,2)]
  impMat2 = impMat2[,c(1,2)]
  impMat3 = impMat3[,c(1,2)]
  
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
  
  impMat = rbind(impMat1, impMat2, impMat3) %>% as.data.frame()
  # impMat.std = rbind(impMat1.std, impMat2.std, impMat3.std) %>% as.data.frame()
  impMat$Quantiles = rep(c("25% to 50%", "50% to 75%", "25% to 75%"), 
                         each = nrow(impMat1))
  impMat$Diff = abs(impMat[,1] - impMat[,2])
  # impMat$Diff = abs(impMat.std[,1] - impMat.std[,2])
  impMat$cov = rep(colnames(X), 3)
  
  # impMat_top = impMat %>% group_by(Quantiles) %>% 
  #   slice_max(order_by = Diff, n = 5)
  
  #rownames(impMat) = rep(colnames(X), 3)
  dtNew <- impMat %>%
    as.data.frame() %>%
    #rownames_to_column(var = "cov") %>%
    # rename variables if buffers
    split_cov_name(cov = "cov") %>%
    #make long format for faceting
    gather(key = "Method", value = "Importance", contains("Method")) %>%
    mutate(Method = substr(Method, 8, nchar(Method))) 
  
  set1 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "PLS") %>%
    slice_max(order_by = Diff, n = 5) %>% select(cov) %>% c() %>% unlist()
  set2 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "PLS") %>%
    slice_max(order_by = Diff, n = 5) %>% select(cov) %>% c() %>% unlist()
  set3 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "PLS") %>%
    slice_max(order_by = Diff, n = 5) %>% select(cov) %>% c() %>% unlist()
  subset = set1 %>% union(set2) %>% union(set3)
  dtNew = dtNew[dtNew$cov %in% subset, ]
  
  dtNew %>%
    #buffered covariates
    drop_na(buffer) %>%
    ggplot(aes(x = Importance, y = cov)) +
    geom_point(aes(size=buffer, color = Method), shape=1, alpha=my.alpha) +
    scale_size(breaks = c(min(dtNew$buffer, na.rm = T),
                          max(dtNew$buffer, na.rm = T))) + #500, 5000, 10000,
    #non-buffered covariates
    geom_point(data = dtNew[is.na(dtNew$buffer),],
               alpha=my.alpha, aes(shape="", color = Method)) +
    geom_vline(xintercept=0,
               linetype="solid",
               alpha=my.alpha) +
    facet_wrap(~ Quantiles, scales = "free_x",
               labeller = "label_both",
               nrow = 1) +
    labs(y = "Geocovariate",
         shape= "non-buffer", #"proximity,\nelevation",
         #title = paste("PLS Importance:", pollutant)) +
         title = tt) +
    theme(legend.position = "bottom", text = element_text(size = 16))
}

plotImp(cbind(exp(avg.pls[,4]) - exp(avg.pls[,2]), 
              exp(avg.spatrf[,4]) - exp(avg.spatrf[,2]),
              exp(avg.rf[,4]) - exp(avg.rf[,2]), 
              exp(avg.rf.tprs[,4]) - exp(avg.rf.tprs[,2])),
        cbind(exp(avg.pls[,6]) - exp(avg.pls[,4]), 
              exp(avg.spatrf[,6]) - exp(avg.spatrf[,4]),
              exp(avg.rf[,6]) - exp(avg.rf[,4]), 
              exp(avg.rf.tprs[,6]) - exp(avg.rf.tprs[,4])),
        cbind(exp(avg.pls[,6]) - exp(avg.pls[,2]), 
              exp(avg.spatrf[,6]) - exp(avg.spatrf[,2]),
              exp(avg.rf[,6]) - exp(avg.rf[,2]), 
              exp(avg.rf.tprs[,6]) - exp(avg.rf.tprs[,2])),
        "Mean difference in prediction, UFP")

load("varimp_new/bc_varimp_cohort_diff.RData")

avg.pls = avg.spatrf = avg.rf = avg.rf.tprs = matrix(NA, ncol = length(qs), nrow = ncol(X))
k = 1

for (i in 1:length(cols)){
  for (j in 1:length(qs)){
    avg.pls[i,j] = mean(bc.pred.diff[[k]]$Y.pls)
    avg.spatrf[i,j] = mean(bc.pred.diff[[k]]$Y.spatrf.pl)
    avg.rf[i,j] = mean(bc.pred.diff[[k]]$Y.rf)
    avg.rf.tprs[i,j] = mean(bc.pred.diff[[k]]$Y.rf.tprs)
    k = k + 1
  }
}

plotImp(cbind(exp(avg.pls[,4]) - exp(avg.pls[,2]), 
              exp(avg.spatrf[,4]) - exp(avg.spatrf[,2]),
              exp(avg.rf[,4]) - exp(avg.rf[,2]), 
              exp(avg.rf.tprs[,4]) - exp(avg.rf.tprs[,2])),
        cbind(exp(avg.pls[,6]) - exp(avg.pls[,4]), 
              exp(avg.spatrf[,6]) - exp(avg.spatrf[,4]),
              exp(avg.rf[,6]) - exp(avg.rf[,4]), 
              exp(avg.rf.tprs[,6]) - exp(avg.rf.tprs[,4])),
        cbind(exp(avg.pls[,6]) - exp(avg.pls[,2]), 
              exp(avg.spatrf[,6]) - exp(avg.spatrf[,2]),
              exp(avg.rf[,6]) - exp(avg.rf[,2]), 
              exp(avg.rf.tprs[,6]) - exp(avg.rf.tprs[,2])),
        "Mean difference in prediction, BC")
