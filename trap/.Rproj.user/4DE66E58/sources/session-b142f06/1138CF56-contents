var_labels_orig = c("pop_s" = "Population", 
                    "pop_dens" = "Population density", 
                    "ndvi_summer_a" = "NDVI (summer)",
                    "ndvi_q75_a" = "NDVI (annual 75%ile)",
                    "ndvi_q50_a" = "NDVI (annual median)",
                    "m_to_ry" = "Dist to railyard",
                    "m_to_rr" = "Dist to railroad",
                    "m_to_l_airp" = "Dist to large airport",
                    "m_to_airp" = "Dist to airport",
                    "m_to_a1" = "Dist to A1 road",
                    "lu_transport_p" = "Transportation",
                    "lu_resi_p" = "Residential",
                    "lu_mix_urban_p" = "Mixed urban",
                    "lu_industrial_p" = "Industrial",
                    "lu_comm_p" = "Commercial",
                    "lu_crop_p" = "Agricultural",
                    "lu_green_p" = "Forest land",
                    "lu_oth_urban_p" = "Other urban",
                    "ll_a3_s" = "A3 road length",
                    "tl_s" = "Truck route length",
                    "rlu_evergreen_p" = "Evergreen forest",
                    "rlu_water_p" = "Water",
                    "pop10_s" = "Population",
                    "ll_a1_s" = "A1 road length",
                    "ll_a23_s" = "A2,A3 road length",
                    "rlu_dev_hi_p" = "Dev high intensity",
                    "rlu_dev_med_p" = "Dev medium intensity",
                    "log_m_to_ry" = "Log(dist) to railyard",
                    "log_m_to_l_airp" = "Log(dist) to large airport",
                    "log_m_to_airp" = "Log(dist) to airport",
                    "log_m_to_m_port" = "Log(dist) to medium port",
                    "log_m_to_a1" = "Log(dist) to A1 road",
                    "log_m_to_a1_a1_intersect" = "Log(dist) to A1/A1 intersec",
                    "log_m_to_a2_a3_intersect" = "Log(dist) to A2/A3 intersec",
                    "log_m_to_a1_a3_intersect" = "Log(dist) to A1/A3 intersec",
                    "log_m_to_truck" = "Log(dist) to truck route",
                    "log_m_to_rr" = "Log(dist) to railroad",
                    "imp_a" = "Impervious surface",
                    "elev_elevation" = "Elevation",
                    "elev_below" = "High relative elevation",
                    "elev_at_elev" = "Intermediate relative elevation",
                    "elev_above" = "Low relative elevation")

var_labels = c("bus_s" = "sum of bus routes",
               "log_m_to_bus" = "log meters to bus route",
               "elev_above" = "low relative elevation",
               "elev_below" = "high relative elevation",
               "elev_at_elev" = "moderate relative elevation",
               "elev_stdev" = "std deviation of elevation",
               "elev_elevation" = "elevation above sea level",
               "imp_a" = "avg imperviousness",
               "intersect_a1_a3_s" = "A1-A3 intersections",
               "intersect_a2_a2_s" = "A2-A2 intersections",
               "intersect_a2_a3_s" = "A2-A3 intersections",
               "intersect_a3_a3_s" = "A3-A3 intersections",
               "ll_a1_s" = "length of A1 roads",
               "ll_a2_s" = "length of A2 roads",
               "ll_a3_s" = "length of A3 roads",
               "ll_a23_s" = "length of A2/A3 roads",
               "m_to_a1" = "meters to A1 road",
               "m_to_a1_a1_intersect" = "meters to A1-A1 intersection",
               "m_to_a1_a2_intersect"= "meters to A1-A2 intersection",
               "m_to_a1_a3_intersect" = "meters to A1-A3 intersection",
               "m_to_a123" = "meters to A1/A2/A3 road",
               "m_to_a2" = "meters to A2 road",
               "m_to_a2_a2_intersect" = "meters to A2-A2 intersection",
               "m_to_a2_a3_intersect"= "meters to A2-A3 intersection",
               "m_to_a23" = "meters to A2/A3 road",
               "m_to_a3" = "meters to A3 road",
               "m_to_a3_a3_intersect" = "meters to A3-A3 intersection",
               "m_to_airp" = "meters to airport",
               "m_to_coast" = "meters to coastline",
               "m_to_comm" = "meters to commercial area",
               "m_to_l_airp" = "meters to large airport",
               "m_to_l_port" = "meters to large port",
               "m_to_m_port" = "meters to medium port",
               "m_to_s_port" = "meters to small port",
               "m_to_rr" = "meters to railroad",
               "m_to_ry" = "meters to rail yard",
               "m_to_truck" = "meters to truck route",
               "m_to_waterway" = "meters to waterway",
               "ndvi_q25_a" = "NDVI (25th quantile)",
               "ndvi_q50_a" = "NDVI (50th quantile)",
               "ndvi_q75_a" = "NDVI (75th quantile)",
               "ndvi_summer_a" = "average summer time NDVI",
               "ndvi_winter_a" = "average winter time NDVI",
               "pop_s" = "population density within a buffer",
               "pop10_s" = "2010 population density",
               "pop90_s" = "1990 population density",
               "pop_dens" = "population density", 
               # "lu_barren_p" = "proportion of barren land",
               # "lu_decid_forest_p" = "proportion of deciduous forest",
               # "lu_dev_hi_p" = "proportion of highly developed land",
               # "lu_dev_lo_p" = "proportion of low developed land",
               # "lu_dev_med_p" = "proportion of medium developed land",
               # "lu_dev_open_p" = "proportion of developed open land",
               # "lu_evergreen_p" = "proportion of evergreen forest",
               # "lu_herb_wetland_p" = "proportion of herb wetland",
               # "lu_mix_forest_p" = "proportion of mixed forest",
               # "lu_water_p" = "proportion of water",
               # "lu_woody_wetland_p" = "proportion of woody wetland",
               "tl_s" = "length of truck routes",
               
               "emissions_s" = "sum of emissions (g)",
               
               "em_co_s" = "sum of CO stack emissions",
               "em_nox_s" = "sum of NOx stack emissions",
               "em_pm10_s" = "sum of PM10 stack emissions",
               "em_pm25_s" = "sum of PM2.5 stack emissions",
               "em_so2_s" = "sum of SO2 stack emissions",
               "no2_behr" = "avg columnar NO2",
               "no2_behr_2005" = "columnar NO2 in 2005",
               "no2_behr_2006" = "columnar NO2 in 2006",
               "no2_behr_2007" = "columnar NO2 in 2007",
               
               "lu_resi_p" = "proportion of land use: residential",      
               "lu_comm_p" = "proportion of land use: commercial",
               "lu_industrial_p" = "proportion of land use: industrial", 
               "lu_transport_p" = "proportion of land use: transportation",
               "lu_industcomm_p"  = "proportion of land use: industrial/commercial",
               "lu_mix_urban_p"   = "proportion of land use: mixed urban",
               "lu_oth_urban_p"   = "proportion of land use: other urban",
               "lu_crop_p" = "proportion of land use: cropland",
               "lu_grove_p" = "proportion of land use: groves", 
               "lu_herb_range_p" = "proportion of land use: herbaceous rangeland",
               "lu_shrub_p" = "proportion of land use: shrub",
               "lu_forest_p" = "proportion of land use: forest",
               "lu_mix_range_p" = "proportion of land use: mixed rangeland",
               "lu_mix_forest_p" = "proportion of land use: mixed forest",
               "lu_reservior_p" = "proportion of land use: reservior",
               "lu_bays_p" = "proportion of land use: bay",
               "lu_stream_p" = "proportion of land use: stream",
               "lu_green_p" = "proportion of land use: evergreen",
               "lu_wetland_p" = "proportion of land use: wetland",
               "lu_nf_wetland_p" = "proportion of land use: nonforested wetland",
               "lu_mine_p" = "proportion of land use: mine",
               "lu_transition_p" = "proportion of land use: transitional area",
               "lu_unspec_p" = "proportion of land use: unspecified",
               "lu_oth_agri_p" = "proportion of land use: other agricultural",
               "lu_lakes_p" = "proportion of land use: lake",
               "lu_rock_p" = "proportion of land use: rock",
               "lu_sandy_p" = "proportion of land use: sandy",
               "lu_beach_p" = "proportion of land use: beach",
               "lu_dry_salt_p" = "proportion of land use: dry salt",
               "lu_shrub_tun_p" = "proportion of land use: shrub tundra",
               "lu_feeding_p" = "proportion of land use: confined feeding",
               "lu_herb_tun_p" = "proportion of land use: herbaceous tundra",
               "lu_bare_tun_p" = "proportion of land use: bare tundra",
               "lu_mix_tun_p" = "proportion of land use: mixed tundra",
               "lu_wet_tun_p" = "proportion of land use: wet tundra",
               "lu_snowfield_p" = "proportion of land use: snowfield",
               "lu_glacier_p" = "proportion of land use: glacier",
               "lu_mix_barren_p" = "proportion of land use: mixed barren")

plotImp2 = function(impMat1, impMat2, impMat3, tt, n_top = 5, top_only = TRUE){
  # my.alpha = .3
  my.alpha = .8
  
  colnames(impMat1) = colnames(impMat2) = colnames(impMat3) =
    c("Method UK-PLS", "Method SpatRF", "Method Truth")
  
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
  impMat$cov = rep(colnames(X), 3)
  # impMat.std = rbind(impMat1.std, impMat2.std, impMat3.std) %>% as.data.frame()
  impMat$Quantiles = rep(c("25% to 50%", "50% to 75%", "25% to 75%"), 
                         each = nrow(impMat1))
  impMat$Diff = abs(impMat[,1] - impMat[,2])
  # impMat$Diff = abs(impMat.std[,1] - impMat.std[,2])
  # impMat$cov = rep(colnames(X), 3)
  
  # impMat_top = impMat %>% group_by(Quantiles) %>% 
  #   slice_max(order_by = Diff, n = 5)
  
  #rownames(impMat) = rep(colnames(X), 3)
  
  dtNew <- impMat %>% select(-`Method Truth`) %>%
    as.data.frame() %>%
    #rownames_to_column(var = "cov") %>%
    # rename variables if buffers
    split_cov_name(cov = "cov") %>%
    #make long format for faceting
    gather(key = "Method", value = "Importance", contains("Method")) %>%
    mutate(Method = substr(Method, 8, nchar(Method))) 

  if (top_only){
    set1 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set1.2 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2.2 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3.2 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set1.3 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
      # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2.3 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
      # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3.3 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
      # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    # subset = set1 %>% union(set2) %>% union(set3)
    subset = set1 %>% union(set2) %>% union(set3) %>% 
      union(set1.2) %>% union(set2.2) %>% union(set3.2) %>%
      union(set1.3) %>% union(set2.3) %>% union(set3.3)
    dtNew = dtNew[dtNew$cov %in% subset, ] %>% arrange(cov)
  }
  
  # dtNew$x1 = sign(dtNew$Importance) * log(abs(dtNew$Importance) + 1)
  
  dtNew %>%
    #buffered covariates
    drop_na(buffer) %>%
    ggplot(aes(x = Importance, y = cov)) + # x = Importance
    geom_point(aes(size=buffer, color = Method), shape=1, alpha=my.alpha) +
    scale_x_continuous(trans = ssqrt_trans) +
    scale_size(breaks = c(min(dtNew$buffer, na.rm = T),
                          max(dtNew$buffer, na.rm = T))) + #500, 5000, 10000,
    #non-buffered covariates
    geom_point(data = dtNew[is.na(dtNew$buffer),],
               alpha=my.alpha, aes(shape="", color = Method)) +
    geom_vline(xintercept=0,
               linetype="solid",
               alpha=.3) +
    scale_x_continuous(trans = ssqrt_trans) +
    scale_y_discrete(labels=var_labels) +
    facet_wrap(~ Quantiles, scales = "fixed", # scales = "free_x",
               labeller = "label_both",
               nrow = 1) +
    labs(y = "Geocovariate", x = "Predictor Contribution",
         shape= "non-buffer", #"proximity,\nelevation",
         #title = paste("PLS Importance:", pollutant)) +
         title = tt) +
    theme(legend.position = "bottom", text = element_text(size = 16))
}

plotmaxcorr = function(impMat1, impMat2, impMat3, tt, n_top = 5, top_only = TRUE,
                       corr_matrix){
  my.alpha = .5
  
  colnames(impMat1) = colnames(impMat2) = colnames(impMat3) =
    c("Method UK-PLS", "Method SpatRF", "Method Truth")
  
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
  impMat$Quantiles = rep(c("25% to 50%", "50% to 75%", "25% to 75%"), 
                         each = nrow(impMat1))
  impMat$Diff = abs(impMat[,1] - impMat[,2])
  impMat$cov = rep(colnames(X), 3)
  
  dtNew <- impMat %>%
    as.data.frame() %>%
    #rownames_to_column(var = "cov") %>%
    # rename variables if buffers
    split_cov_name(cov = "cov") %>%
    #make long format for faceting
    gather(key = "Method", value = "Importance", contains("Method")) %>%
    mutate(Method = substr(Method, 8, nchar(Method))) 
  
  if (top_only){
    set1 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set1.2 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2.2 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3.2 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set1.3 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
    # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2.3 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
    # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3.3 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
    # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    # subset = set1 %>% union(set2) %>% union(set3)
    subset = set1 %>% union(set2) %>% union(set3) %>% 
      union(set1.2) %>% union(set2.2) %>% union(set3.2) %>%
      union(set1.3) %>% union(set2.3) %>% union(set3.3)
    dtNew = dtNew[dtNew$cov %in% subset, ] %>% arrange(cov)
  }
  
  # dtNew$x1 = sign(dtNew$Importance) * log(abs(dtNew$Importance) + 1)
  
  impMat = apply(corr_matrix, 1, function(x) max(abs(x))) %>% as.data.frame()
  impMat$cov = colnames(X)
  colnames(impMat) = c("corr", "cov")
  
  dtNew <- impMat %>%
    as.data.frame() %>%
    split_cov_name(cov = "cov") # %>%
    #make long format for faceting
    # gather(key = "Method", value = "Importance", contains("Method")) %>%
    # mutate(Method = substr(Method, 8, nchar(Method))) 
  
  dtNew$Corr = "Max Abs Value"
  
  dtNew %>%
    #buffered covariates
    drop_na(buffer) %>%
    ggplot(aes(x = corr, y = cov)) + # x = Importance
    geom_point(aes(size=buffer), shape=1, alpha=my.alpha) +
    scale_x_continuous(trans = ssqrt_trans) +
    scale_size(breaks = c(min(dtNew$buffer, na.rm = T),
                          max(dtNew$buffer, na.rm = T))) + #500, 5000, 10000,
    #non-buffered covariates
    geom_point(data = dtNew[is.na(dtNew$buffer),],
               alpha=my.alpha, aes(shape="")) +
    # geom_vline(xintercept=0,
    #            linetype="solid",
    #            alpha=.3) +
    scale_x_continuous(trans = ssqrt_trans) +
    scale_y_discrete() +
    facet_wrap(~ Corr, scales = "fixed", # scales = "free_x",
               labeller = "label_both",
               nrow = 1) +
    theme(axis.text.y = element_blank()) +
    labs(y = " ", x = "Max Abs Corr",
         shape= "non-buffer", #"proximity,\nelevation",
         #title = paste("PLS Importance:", pollutant)) +
         title = tt) +
    theme(legend.position = "bottom", text = element_text(size = 16))
}

plotcorr = function(impMat1, impMat2, impMat3, tt, n_top = 5, top_only = TRUE,
                    corr_matrix){
  # my.alpha = .3
  my.alpha = .8
  
  colnames(impMat1) = colnames(impMat2) = colnames(impMat3) =
    c("Method UK-PLS", "Method SpatRF", "Method Truth")
  
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
  
  if (top_only){
    set1 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "UK-PLS") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set1.2 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2.2 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3.2 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "SpatRF") %>%
      slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set1.3 = dtNew %>% filter(Quantiles == "25% to 50%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
      # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set2.3 = dtNew %>% filter(Quantiles == "50% to 75%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
      # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    set3.3 = dtNew %>% filter(Quantiles == "25% to 75%" & Method == "Truth") %>%
      slice(which(Importance != 0)) %>% select(cov) %>% c() %>% unlist()
      # slice_max(order_by = abs(Importance), n = n_top) %>% select(cov) %>% c() %>% unlist()
    # subset = set1 %>% union(set2) %>% union(set3)
    subset = set1 %>% union(set2) %>% union(set3) %>% 
      union(set1.2) %>% union(set2.2) %>% union(set3.2) %>%
      union(set1.3) %>% union(set2.3) %>% union(set3.3)
    # dtNew = dtNew[dtNew$cov %in% subset, ] %>% arrange(cov)
  } else {
    subset = dtNew$cov
  }
  
  colnames(corr_matrix) =
    c("Var Distance to A1", "Var Population Density", "Var NDVI", "Var LU-Urban",
      "Var LU-Residential")
  
  impMat = corr_matrix %>% as.data.frame()
  impMat$Corr = "Active Predictors"
  impMat$cov = colnames(X)

  dtNew <- impMat %>%
    as.data.frame() %>%
    #rownames_to_column(var = "cov") %>%
    # rename variables if buffers
    split_cov_name(cov = "cov") %>%
    #make long format for faceting
    gather(key = "Var", value = "Correlation", contains("Var")) %>%
    mutate(Var = substr(Var, 5, nchar(Var))) %>%
    filter(cov %in% subset)
  
  dtNew %>%
    #buffered covariates
    drop_na(buffer) %>%
    ggplot(aes(x = Correlation, y = cov)) + # x = Importance
    geom_point(aes(size=buffer, color = Var), shape=1, alpha=my.alpha) +
    scale_x_continuous(trans = ssqrt_trans) +
    scale_size(breaks = c(min(dtNew$buffer, na.rm = T),
                          max(dtNew$buffer, na.rm = T))) + #500, 5000, 10000,
    #non-buffered covariates
    geom_point(data = dtNew[is.na(dtNew$buffer),],
               alpha=my.alpha, aes(shape="", color = Var)) +
    geom_vline(xintercept=0,
               linetype="solid",
               alpha=.3) +
    scale_x_continuous() + # trans = ssqrt_trans) +
    scale_y_discrete() +
    # scale_y_discrete(labels=c("pop_s" = "Population", 
    #                           "pop_dens" = "Population density", 
    #                           "ndvi_summer_a" = "NDVI (summer)",
    #                           "ndvi_q75_a" = "NDVI (annual 75%ile)",
    #                           "ndvi_q50_a" = "NDVI (annual median)",
    #                           "m_to_ry" = "Dist to railyard",
    #                           "m_to_rr" = "Dist to railroad",
    #                           "m_to_l_airp" = "Dist to large airport",
    #                           "m_to_airp" = "Dist to airport",
    #                           "m_to_a1" = "Dist to A1 road",
    #                           "lu_transport_p" = "Transportation",
    #                           "lu_resi_p" = "Residential",
    #                           "lu_mix_urban_p" = "Mixed urban",
    #                           "lu_industrial_p" = "Industrial",
    #                           "lu_comm_p" = "Commercial",
    #                           "lu_crop_p" = "Agricultural",
    #                           "lu_green_p" = "Forest land",
    #                           "lu_oth_urban_p" = "Other urban",
    #                           "ll_a3_s" = "A3 road length",
    #                           "tl_s" = "Truck route length",
    #                           "rlu_evergreen_p" = "Evergreen forest",
    #                           "rlu_water_p" = "Water",
    #                           "pop10_s" = "Population",
    #                           "ll_a1_s" = "A1 road length",
    #                           "ll_a23_s" = "A2,A3 road length",
    #                           "rlu_dev_hi_p" = "Dev high intensity",
    #                           "rlu_dev_med_p" = "Dev medium intensity",
    #                           "log_m_to_ry" = "Log(dist) to railyard",
    #                           "log_m_to_l_airp" = "Log(dist) to large airport",
    #                           "log_m_to_airp" = "Log(dist) to airport",
    #                           "log_m_to_m_port" = "Log(dist) to medium port",
    #                           "log_m_to_a1" = "Log(dist) to A1 road",
    #                           "log_m_to_a1_a1_intersect" = "Log(dist) to A1/A1 intersec",
    #                           "log_m_to_a2_a3_intersect" = "Log(dist) to A2/A3 intersec",
    #                           "log_m_to_a1_a3_intersect" = "Log(dist) to A1/A3 intersec",
    #                           "log_m_to_truck" = "Log(dist) to truck route",
    #                           "log_m_to_rr" = "Log(dist) to railroad",
    #                           "imp_a" = "Impervious surface",
    #                           "elev_elevation" = "Elevation",
    #                           "elev_below" = "High relative elevation",
    #                           "elev_at_elev" = "Intermediate relative elevation",
    #                           "elev_above" = "Low relative elevation")) +
    facet_wrap(~ Corr, scales = "fixed", # scales = "free_x",
               labeller = "label_both",
               nrow = 1) +
    theme(axis.text.y = element_blank()) +
    labs(y = "", x = "Corr w/ Active Predictors",
         shape= "non-buffer", #"proximity,\nelevation",
         #title = paste("PLS Importance:", pollutant)) +
         title = tt) +
    theme(legend.position = "bottom", text = element_text(size = 16))
}
