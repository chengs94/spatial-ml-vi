library(tidyverse)
library(sf)
library(ggmap)
library(ggplot2)
library(gridExtra)
library(grid)
library(gstat)
library(ggspatial)

load("rslt_grid_uw.RData")

crs_m <- 32148
crs_deg <- 4326 #WGS84. in decimal degrees
## PROJ.4 string for crs_deg
crs_deg_proj4 <- "+proj=longlat +datum=WGS84 +no_defs"

monitoring_area_shp <- readRDS("/Users/sicheng/Library/CloudStorage/OneDrive-UW/exposurepred/data/monitoring_area_shp.rda") %>%
  #convert from 4269 
  st_transform(crs_deg)
monitoring_land_shp <- readRDS("/Users/sicheng/Library/CloudStorage/OneDrive-UW/exposurepred/data/monitoring_land_shp.rda")

grid_covars_shp <- cov_grid %>%
  st_as_sf(., coords=c("longitude","latitude"), remove = FALSE,
           crs=crs_deg)
grid_covars_shp$in_study_area <- st_intersects(grid_covars_shp, monitoring_land_shp,
                                               sparse = FALSE) %>% apply(., 1, any)
grid_covars_shp$Y.pls = exp(grid.ufp.uw$Y.pls)/1000
grid_covars_shp$Y.spatrf = exp(grid.ufp.uw$Y.spatrf.pl)/1000

grid_resolution <- 0.003
grid_idp <- 1.5 
grid_nmax <- 10

if (!file.exists("idw_df.RData")){
  finer_grid <- st_as_sf(grid_covars_shp, coords = c("longitude", "latitude"),  crs= crs_deg) %>%
    # make a rectangular box w/ evenly spaced points at ~500 m resolution
    st_bbox() %>% st_as_sfc() %>%
    st_make_grid(cellsize = grid_resolution, what = "centers") %>% 
    # view in df format
    st_as_sf()  
  
  idw.pls = idw(formula = Y.pls~1,  
                locations = st_as_sf(grid_covars_shp, coords = c("longitude", "latitude"),  
                                     remove = F, crs= crs_deg), 
                newdata = finer_grid,  
                #smaller inverse distnace powers produce more smoothing (i.e., give points furthewr away larger weights). default is the max, 2.0.
                idp=grid_idp, 
                # how many nearby points should be used to smooth
                nmax=grid_nmax #maxdist=0.1
  )
  idw.spatrf = idw(formula = Y.spatrf~1,  
                   locations = st_as_sf(grid_covars_shp, coords = c("longitude", "latitude"),  
                                        remove = F, crs= crs_deg), 
                   newdata = finer_grid,  
                   #smaller inverse distnace powers produce more smoothing (i.e., give points furthewr away larger weights). default is the max, 2.0.
                   idp=grid_idp, 
                   # how many nearby points should be used to smooth
                   nmax=grid_nmax #maxdist=0.1
  )
  idw_df = cbind(idw.pls, idw.spatrf$var1.pred)# [st_intersection(idw.pls, monitoring_land_shp),]
  idw_df = st_intersection(idw_df, monitoring_land_shp) %>% cbind(., st_coordinates(.))
  colnames(idw_df)[c(1,3)] = c("Y.pls", "Y.spatrf")
  
  saveRDS(idw_df, "idw_df.RData")
} else {
  idw_df = readRDS("idw_df.RData")
}

idw_df$Y.diff = idw_df$Y.spatrf - idw_df$Y.pls

bbox <- st_bbox(st_transform(st_buffer(st_transform(monitoring_area_shp, crs_m), 10000), crs_deg))
names(bbox) <- c("left", "bottom", "right", "top")
map0 <- suppressMessages(get_stamenmap(
  bbox = bbox, 
  zoom = 11, 
  maptype = "toner-lite" #has airport symbol
  #maptype = "toner-background" #roads & water but no airport
))
## map labels
map_x_labels <- c(seq(-122.5, -121.9, 0.2)) #0.2
map_y_labels <- c(seq(47.2, 48, 0.2))

idw_df_long = idw_df %>% gather(model, value, c("Y.pls", "Y.spatrf", "Y.diff")) %>%
  mutate(model = recode(model, Y.pls = "UK-PLS", Y.spatrf = "SpatRF (PL)", 
                        Y.diff = "Difference"))

p1 = ggmap(ggmap = map0, darken = c(.5, "white")) +
  #predictions
  geom_raster(data = idw_df_long %>% filter(model != "Difference"), 
              inherit.aes = F, 
              aes(fill= value, x=round(X, 6), y=round(Y, 6)),
              interpolate = T)  +
  geom_sf(data=monitoring_land_shp, inherit.aes = F, alpha=0, size=0.1) +
  # scale_fill_gradient(name = "1000pt/cm3", low = "yellow", high = "red") +
  scale_fill_viridis_c(option = "cividis", trans = "log",
                        name = "1000pt/cm3", labels = scales::number_format(accuracy = 1)) +
  facet_wrap(~model)  +
  # add scales & N arrow 
  annotation_scale(location = "tr") +
  annotation_scale(location = "tr", unit_category ="imperial", pad_y = unit(0.55, "cm")) +
  annotation_north_arrow(location = "tr", pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  theme_bw() +
  theme(
    legend.justification=c(0,1),  
    legend.position=c(0,1),  
    legend.background =  element_blank()
  ) +
  coord_sf(expand = F) +
  scale_x_continuous(breaks = map_x_labels, labels = map_x_labels ) +
  scale_y_continuous(breaks = map_y_labels,
                     labels = format(map_y_labels,digits = 1, nsmall = 1)
  ) +
  #add attribution/reference to bottom left
  geom_text(aes(x=-Inf, y=-Inf, hjust=-0.01, vjust=-0.3, 
                label= "Map tiles by Stamen Design, under CC BY 3.0. \nData by OpenStreetMap, under ODbL."), size=2.5,
  ) +
  labs(x = "Longitude", y = "Latitude") 

p2 = ggmap(ggmap = map0, darken = c(.5, "white")) +
  #predictions
  geom_raster(data = idw_df_long %>% filter(model == "Difference"), 
              inherit.aes = F, 
              aes(fill= value, x=round(X, 6), y=round(Y, 6)),
              interpolate = T)  +
  geom_sf(data=monitoring_land_shp, inherit.aes = F, alpha=0, size=0.1) +
  # scale_fill_gradient(name = "1000pt/cm3", low = "yellow", high = "red") +
  scale_fill_gradient2(mid = "grey",
                       # high = muted("blue"),
                       midpoint = 0, # trans = "log",
                       name = "1000pt/cm3") +
  facet_wrap(~model)  +
  # add scales & N arrow 
  annotation_scale(location = "tr") +
  annotation_scale(location = "tr", unit_category ="imperial", pad_y = unit(0.55, "cm")) +
  annotation_north_arrow(location = "tr", pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  theme_bw() +
  theme(
    legend.justification=c(0,1),  
    legend.position=c(0,1),  
    legend.background =  element_blank()
  ) +
  coord_sf(expand = F) +
  scale_x_continuous(breaks = map_x_labels, labels = map_x_labels ) +
  scale_y_continuous(breaks = map_y_labels,
                     labels = format(map_y_labels,digits = 1, nsmall = 1)
  ) +
  #add attribution/reference to bottom left
  geom_text(aes(x=-Inf, y=-Inf, hjust=-0.01, vjust=-0.3, 
                label= "Map tiles by Stamen Design, under CC BY 3.0. \nData by OpenStreetMap, under ODbL."), size=2.5,
  ) +
  labs(x = "Longitude", y = "Latitude") 

library(egg)
ggarrange(p1, p2, ncol = 2, widths = c(2,1))
