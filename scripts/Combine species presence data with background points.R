##################################################################################
### Combine species presence data with background points for habitat modelling ###
##################################################################################

# 2022-03-09
# updated 2023-04-05

# load libraries
require(tidyverse)
require(sf)
require(raster)

# define stereographic projection around Antarctica
prj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# load data
dat <- readRDS("data/combined_presence_data.rds")
dat <- dat %>% filter(lat < -40)

# append xy
dat[c("x", "y")] <- dat %>% st_as_sf(coords = c("lon", "lat")) %>% 
  sf::st_set_crs(4326) %>%
  st_transform(prj) %>%
  st_coordinates()

# generate 10000 background points for each species
source("R/pseudoAbs.R")
pabs <- pseudoAbs(x_min = -180,
                  x_max = 180,
                  y_min = -80,
                  y_max = -40,
                  projection = prj,
                  n = 10000 * length(unique(dat$species)))

# create dummy dataframe replicating structure of original data
# replace locations with pseudo absences coordinates
#pseudo <- dat
pseudo <- dat %>% group_by(group, species) %>% dplyr::select(group, species, date) %>% slice(rep(1, 10000)) %>% ungroup
pseudo$date <- NA
pseudo[c("x", "y")] <- pabs
pseudo[c("lon", "lat")] <- pseudo %>% st_as_sf(coords = c("x", "y")) %>% 
  sf::st_set_crs(prj) %>%
  st_transform(4326) %>%
  st_coordinates()
pseudo <- pseudo %>% dplyr::select("group", "species", "date", "lon", "lat", "x", "y")

# define 1s and 0s
dat$PresAbs <- 1
pseudo$PresAbs <- 0

# combine original data with pseudo absences
dat <- rbind(dat, pseudo)
dat <- dat %>% dplyr::select("group", "species", "date", "PresAbs", "lon", "lat", "x", "y")

# plot to check
ggplot() + 
  geom_point(aes(x = x, y = y, colour = factor(PresAbs)), data = dat)

# save outputted data frame
saveRDS(dat, "data/presence_absence_data_10k.rds")

# ends