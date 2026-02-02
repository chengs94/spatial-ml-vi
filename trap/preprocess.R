# replicating the TRAP analysis

library(tidyverse)

#source("https://raw.githubusercontent.com/magali17/TRAP/master/A2.0.3_functions_covariate_preprocessing_regional.R")
source("preprocessing_function.R")

#cov_mm0 <- read.csv("dr0311_mobile_covars.csv") %>%
#cov_mm0 <- read.csv("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/data/dr0311_grid_covars.csv") %>%
cov_mm0 <- read.csv("E:/dr0311_cohort_covars.csv") %>%
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

#saveRDS(cov_mm,file = "cov_mm.rda")

#### grid ######
# grid_covars0 <- read.csv("C:/Users/Si Cheng/Dropbox (Personal)/exposurepred/data/dr0311_grid_covars.csv")
# 
# cov_grid <- grid_covars0 %>%
#   combine_a123_m() %>%
#   combine_a23_m() %>%
#   combine_a23_ll() %>%
#   log_transform_distance() %>%
#   # duplicate column - for later
#   mutate(site_id = native_id)
# saveRDS(cov_mm,file = "cov_grid.rda")

############ annual averages ###############
annual0 = read.csv("annual_2020-10-02_loc_avg.csv")

#covariates
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

#annual <- annual0 %>%
  #convert to log
#  mutate_at(vars(contains(c("ufp", "bc"))), ~log(.)) %>%
#  mutate(
#    ufp_native_scale = exp(ufp_primary),
#    bc_native_scale = exp(bc_primary)) %>%
#  left_join(cov_mm) %>%
  # add native scale proximity variables
#  left_join(proximity_vars_native)

annual = merge(annual0, cov_mm, by = "site_id")

annual = annual[,-2]
annual[,5:16] = log(annual[,5:16])

saveRDS(annual, "dat.rda")

#cov_mm2 = cov_mm[,colnames(cov_mm) %in% colnames(dat)]
#saveRDS(cov_mm2, "cov_grid.rda")

#ufp_names <- names(annual)[grepl("ufp", names(annual))]
#bc_names <- names(annual)[grepl("bc", names(annual))]
#analysis_names <- c(ufp_names, bc_names)

##### cohort #####
cov_mm = cov_mm[,colnames(cov_mm) %in% colnames(cov_grid)]
cov_mm_ordered = cov_mm[colnames(cov_grid)]
saveRDS(cov_mm_ordered, "cov_cohort.rda")
