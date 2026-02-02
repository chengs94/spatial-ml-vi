# preprocessing for annual average pollution concentration

library(tidyverse)

source("preprocessing_function.R")

# file_path = "C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/"
file_path = "/Users/sicheng/Library/CloudStorage/OneDrive-UW/exposurepred/"
stops = read.csv(paste(file_path, "data/updated/stop_data.csv", sep = ""), header = TRUE)

annual0 = stops %>% group_by(location, variable) %>%
  mutate(med_trim = ifelse(median_value >= quantile(median_value, 0.05, na.rm = T) &
                           median_value <= quantile(median_value, 0.95, na.rm = T), median_value, NA),
         med_wind = ifelse(median_value == max(median_value), max(median_value[median_value!=max(median_value)]),
                           ifelse(median_value == min(median_value), min(median_value[median_value!=min(median_value)]),
                                           median_value))) %>%
  summarize(
    avg_trim_med = mean(med_trim, na.rm = TRUE),
    avg_wind_med = mean(med_wind)
  ) %>%
  gather("annual", "value", contains("med")) %>%
  ungroup()

# annual1 = annual0 %>% filter(annual == "avg_trim_med") %>% 
#   mutate(location = as.factor(location),
#          variable = as.factor(variable),
#          annual = as.factor(annual)) %>%
#   spread(key = "variable", value = "value")

annual1 = annual0 %>% filter(annual == "avg_wind_med") %>% 
  mutate(location = as.factor(location),
         variable = as.factor(variable),
         annual = as.factor(annual)) %>%
  spread(key = "variable", value = "value")

####### 
cov_mm0 <- read.csv(paste(file_path, "data/updated/dr0311_mobile_covars.csv", sep = ""), header = TRUE) %>%
  # fix labeling & duplicate issue w/ cohort pop covariates
  select(-contains(".y")) %>%
  rename_at(vars(contains(".x")), ~ gsub(".x", "", .))

#####
cov_mm0 = cov_mm0[!is.na(cov_mm0$pop10_s00500),]

cov_mm1 <- cov_mm0 %>%
  #drop columns if any NAs
  select_if(~!any(is.na(.))) %>%
  #don't need this
  select(-contains("region")) %>%
  #need this for fn below to work
  mutate(site_type = "FIXED")

#save mm locations
locations_mm <- cov_mm1 %>%
  select(native_id,
         latitude:lambert_y) %>%
  #duplicate column - for later
  mutate(site_id = native_id)

#only keep covariates
cov_mm1 <- cov_mm1 %>%
  ### --> check in future that this is correct
  select(-c(location_id:msa))

mm_names <- names(cov_mm1)

######################## covariate.preprocess() ########################
cov_preprocessed <- covariate.preprocess(covars.mon = cov_mm1, 
                                         covars.sub = NULL,
                                         #covars.sub = cov_act1,    
                                         region = "NW")

############ Drop additinoal variables I won't be using ####################
# * emissions data 
cov_mm <- cov_preprocessed$covars.mon %>%
  # drop AP predictions (these are based on other covariate models)
  select(-contains("log_em_"),  -contains("no2_")) 
final_cov <- names(cov_mm)

cov_mm <- cbind(locations_mm, cov_mm)

cov_mm <- cov_mm %>%
  # drop pop90, pop_ ( & ??NDVI)
  select(-contains(c("pop90_", "pop_"))) %>%
  #drop Roosevelt garage & stop w/ 1 obeservation that was replaced by MS0601
  filter(!site_id %in% c("MS0000", "MS0398"))

site_loc_vars <- cov_mm %>%
  select(native_id, site_id, latitude, longitude, lambert_x, lambert_y) %>%
  names()

#variable names in log and native scale
non_proximity_vars <- names(cov_mm)[!grepl("m_to_", names(cov_mm)) & !names(cov_mm) %in% site_loc_vars]
proximity_vars_log <- names(cov_mm)[grepl("m_to_", names(cov_mm))]
proximity_vars_native <- str_replace(string = proximity_vars_log, pattern = "log_", replacement = "")
cov_names_log <- append(non_proximity_vars, proximity_vars_log)
cov_names_native <- append(non_proximity_vars, proximity_vars_native)
proximity_vars_native <- cov_mm %>%
  select(site_id, proximity_vars_log) %>%
  mutate_at(proximity_vars_log, ~exp(.)) %>%
  rename_at(proximity_vars_log, ~sub(x = ., "log_", ""))

# annual = annual1 %>% mutate(ufp_uw = log(pnc_noscreen),
#                             bc_uw = log(ma200_ir_bc1),
#                             co2_uw = log(co2_umol_mol),
#                             no2_uw = log(no2),
#                             orig_pm25_uw = neph_bscat,
#                             pm25_uw = log(neph_bscat)) %>%
#   select(contains(c("location", "_uw"))) %>% 
#   left_join(cov_mm, by = c("location" = "native_id"))

# saveRDS(annual, "dat_updated.rda")

dat = readRDS("dat_updated.rda")
annual = annual1 %>% mutate(ufp_uw = log(pnc_noscreen),
                            bc_uw = log(ma200_ir_bc1),
                            co2_uw = log(co2_umol_mol),
                            no2_uw = log(no2),
                            orig_pm25_uw = neph_bscat,
                            pm25_uw = log(neph_bscat)) %>%
  select(contains(c("location", "_uw"))) %>% 
  left_join(dat %>% select(-(2:7)), by = "location")

saveRDS(annual, "dat_final.rda")
