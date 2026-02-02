# exploratory analyses, plots, etc for the manuscript

library(dplyr)
# library(sf)
library(ggmap)
library(ggplot2)
library(ggpubr)
library(gstat)
library(ggspatial)

################ Seattle #################
dat = readRDS("dat_updated.rda")

# seaMap <- get_map(location=c(min(dat$longitude), min(dat$latitude),
#                              max(dat$longitude), max(dat$latitude)),
#                   source="stamen", maptype = "terrain", crop=FALSE)
seaMap <- suppressMessages(get_stamenmap(
  bbox=c(left = min(dat$longitude) - .05, bottom = min(dat$latitude) - .05,
             right = max(dat$longitude) + .01, top = max(dat$latitude) + .01),
  zoom = 11, 
  maptype = "toner-lite" #has airport symbol
  #maptype = "toner-background" #roads & water but no airport
))

map_x_labels <- c(seq(-122.5, -121.9, 0.2)) #0.2
map_y_labels <- c(seq(47.2, 48, 0.2))

p_ufp = ggmap(seaMap, darken = c(.5, "white")) + # labs(x = "longitude", y = "latitude") +
  geom_point(data = dat, aes(x = longitude, y = latitude, 
                                color = exp(ufp_uw)/1000, size = exp(ufp_uw)/1000, alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "log",
                        name = "1000pt/cm3", labels = scales::number_format(accuracy = 1)) +
  annotation_scale(location = "tr") +
  annotation_scale(location = "tr", unit_category ="imperial", pad_y = unit(0.55, "cm")) +
  annotation_north_arrow(location = "tr", pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  # theme_bw() +
  # theme(
  #   legend.justification=c(0,1),  
  #   legend.position=c(0,1),  
  #   legend.background =  element_blank()
  # ) +
  coord_sf(expand = F, crs = 4326) +
  # scale_x_continuous(breaks = map_x_labels, labels = map_x_labels ) +
  # scale_y_continuous(breaks = map_y_labels,
  #                    labels = format(map_y_labels,digits = 1, nsmall = 1)
  # ) +
  geom_text(aes(x=-Inf, y=-Inf, hjust=-0.01, vjust=-0.3, 
                label= "Map tiles by Stamen Design, under CC BY 3.0. \nData by OpenStreetMap, under ODbL."), size=2.5,
  ) +
  guides(size = "none", alpha = "none") + ggtitle("Annual Average UFP Concentration at Monitoring Locations") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16)) +
  labs(x = "Longitude", y = "Latitude")

p_bc = ggmap(seaMap, darken = c(.5, "white")) + labs(x = "longitude", y = "latitude") +
  geom_point(data = dat, aes(x = longitude, y = latitude, 
                             color = exp(bc_uw), size = exp(bc_uw), alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "log",
                        name = "ng/m3", labels = scales::number_format(accuracy = 1)) +
  guides(size = "none", alpha = "none") + ggtitle("BC") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16))+
  labs(x = "Longitude", y = "Latitude")

p_no2 = ggmap(seaMap, darken = c(.5, "white")) + labs(x = "longitude", y = "latitude") +
  geom_point(data = dat, aes(x = longitude, y = latitude, 
                             color = exp(no2_uw), size = exp(no2_uw), alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "log", name = "ppb", 
                        labels = scales::number_format(accuracy = 1)) +
  guides(size = "none", alpha = "none") + ggtitle("NO2") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16)) +
  labs(x = "Longitude", y = "Latitude")

p_co2 = ggmap(seaMap, darken = c(.5, "white")) + labs(x = "longitude", y = "latitude") +
  geom_point(data = dat, aes(x = longitude, y = latitude, 
                             color = exp(co2_uw), size = exp(co2_uw), alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "log",
                        name = "ppm", labels = scales::number_format(accuracy = 1)) +
  guides(size = "none", alpha = "none") + ggtitle("CO2") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16)) +
  labs(x = "Longitude", y = "Latitude")

p_pm25 = ggmap(seaMap, darken = c(.5, "white")) + labs(x = "longitude", y = "latitude") +
  geom_point(data = dat, aes(x = longitude, y = latitude, 
                             color = exp(pm25_uw), size = exp(pm25_uw), alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "log",
                        name = expression(mu * "g/m3"),  
                        labels = scales::number_format(accuracy = 1)) +
  guides(size = "none", alpha = "none") + ggtitle("PM2.5") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16)) +
  labs(x = "Longitude", y = "Latitude")

# ggarrange(p_ufp, p_bc, p_no2, p_co2, p_pm25, 
#           nrow = 1, ncol = 5, common.legend = FALSE, legend = "bottom") %>%
#   annotate_figure(top = text_grob("Annual Average Pollutant Concentration, Seattle", 
#                                   face = "bold", size = 16))

ggarrange(p_bc, p_no2, p_co2, p_pm25, 
          nrow = 1, ncol = 4, common.legend = FALSE, legend = "bottom") %>%
  annotate_figure(top = text_grob("Annual Average Pollutant Concentration, Seattle", 
                                  face = "bold", size = 16))

################ National #################
rm(list = ls())

dat.ec <- read.csv("/Users/sicheng/Library/CloudStorage/OneDrive-UW/exposurepred/spatialRF/data/EC.txt")
dat.s <- read.csv("/Users/sicheng/Library/CloudStorage/OneDrive-UW/exposurepred/spatialRF/data/S.txt")
dat.oc = read.csv("/Users/sicheng/Library/CloudStorage/OneDrive-UW/exposurepred/spatialRF/data/OC.txt")
dat.si <- read.csv("/Users/sicheng/Library/CloudStorage/OneDrive-UW/exposurepred/spatialRF/data/Si.txt")

sum1 <- ifelse(is.na(dat.ec[,2]),0,dat.ec[,2]*dat.ec[,3])
obs1 <-  ifelse(is.na(dat.ec[,3]),0,dat.ec[,3])
sum2 <- ifelse(is.na(dat.ec[,5]),0,dat.ec[,5]*dat.ec[,6])
obs2 <- ifelse(is.na(dat.ec[,6]),0,dat.ec[,6])
dat.ec$annavg <- (sum1+sum2)/(obs1+obs2)
dat.s$annavg <- dat.s[,2]
sum1 <- ifelse(is.na(dat.oc[,2]),0,dat.oc[,2]*dat.oc[,3])
obs1 <-  ifelse(is.na(dat.oc[,3]),0,dat.oc[,3])
sum2 <- ifelse(is.na(dat.oc[,5]),0,dat.oc[,5]*dat.oc[,6])
obs2 <- ifelse(is.na(dat.oc[,6]),0,dat.oc[,6])
dat.oc$annavg = (sum1+sum2)/(obs1+obs2)
dat.si$annavg = dat.si[,2]

# ecMap <- get_map(location=c(min(dat.ec$longitude) - 2, min(dat.ec$latitude) - 2,
#                             max(dat.ec$longitude) + 2, max(dat.ec$latitude) + 2) ,
#                  source="stamen", maptype = "terrain", crop=TRUE)
bbox = c(left = min(dat.ec$longitude) - 2, bottom = min(dat.ec$latitude) - 2,
         right = max(dat.ec$longitude) + 2, top = max(dat.ec$latitude) + 2)
ecMap = suppressMessages(get_stamenmap(
  bbox=bbox,
  zoom = 6,
  maptype = "toner-lite" #has airport symbol
  #maptype = "toner-background" #roads & water but no airport
))

p_ec = ggmap(ecMap, darken = c(.5, "white")) + labs(x = "Longitude", y = "Latitude") +
  geom_point(data = dat.ec, aes(x = longitude, y = latitude, 
                             color = annavg, size = annavg, alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "sqrt",
                        name = expression(mu * "g/m3"), 
                        labels = scales::number_format(accuracy = 0.01)) +
  guides(size = "none", alpha = "none") + ggtitle("EC") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16))

# sMap <- get_map(location=c(min(dat.s$longitude) - 2, min(dat.s$latitude) - 2,
#                            max(dat.s$longitude) + 2, max(dat.s$latitude) + 2) ,
#                 source="stamen", maptype = "terrain", crop=TRUE)
p_s = ggmap(ecMap, darken = c(.5, "white")) + labs(x = "Longitude", y = "Latitude") +
  geom_point(data = dat.s, aes(x = longitude, y = latitude, 
                                color = annavg, size = annavg, alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "sqrt",
                        name = expression(mu * "g/m3"), 
                        labels = scales::number_format(accuracy = 0.01)) +
  annotation_scale(location = "tr") +
  annotation_scale(location = "tr", unit_category ="imperial", pad_y = unit(0.55, "cm")) +
  annotation_north_arrow(location = "tr", pad_y = unit(0.5, "in"), pad_x = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  # theme_bw() +
  # theme(
  #   legend.justification=c(0,1),
  #   legend.position=c(0,1),
  #   legend.background =  element_blank()
  # ) +
  coord_sf(expand = F, crs = 4326) +
  geom_text(aes(x=-Inf, y=-Inf, hjust=-0.01, vjust=-0.3, 
                label= "Map tiles by Stamen Design, under CC BY 3.0. \nData by OpenStreetMap, under ODbL."), size=2.5,
  ) +
  guides(size = "none", alpha = "none") + ggtitle("Annual Average S Concentration at Monitoring Locations") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16))

# siMap <- get_map(location=c(min(dat.si$longitude) - 2, min(dat.si$latitude) - 2,
#                             max(dat.si$longitude) + 2, max(dat.si$latitude) + 2) ,
#                  source="stamen", maptype = "terrain", crop=TRUE)
p_si = ggmap(ecMap, darken = c(.5, "white")) + labs(x = "Longitude", y = "Latitude") +
  geom_point(data = dat.si, aes(x = longitude, y = latitude, 
                               color = annavg, size = annavg, alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "sqrt",
                        name = "ng/m3", 
                        labels = scales::number_format(accuracy = 0.01)) +
  guides(size = "none", alpha = "none") + ggtitle("Si") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16))

# ocMap <- get_map(location=c(min(dat.ec$longitude) - 2, min(dat.ec$latitude) - 2,
#                             max(dat.ec$longitude) + 2, max(dat.ec$latitude) + 2) ,
#                  source="stamen", maptype = "terrain", crop=TRUE)
p_oc = ggmap(ecMap, darken = c(.5, "white")) + labs(x = "Longitude", y = "Latitude") +
  geom_point(data = dat.oc, aes(x = longitude, y = latitude, 
                                color = annavg, size = annavg, alpha = .9)) +
  scale_color_viridis_c(option = "cividis", trans = "sqrt",
                        name = expression(mu * "g/m3"), 
                        labels = scales::number_format(accuracy = 0.01)) +
  guides(size = "none", alpha = "none") + ggtitle("OC") +
  theme(plot.margin=unit(rep(0,4), "pt"), plot.title = element_text(size=16))

ggarrange(p_ec, p_oc, p_si,
          nrow = 3, ncol = 1, common.legend = FALSE, legend = "right") %>%
  annotate_figure(top = text_grob("Annual Average Pollutant Concentration, National", 
                                  face = "bold", size = 16))
