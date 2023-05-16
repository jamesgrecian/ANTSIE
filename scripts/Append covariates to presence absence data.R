#####################################################################################
### Append contemporary covariates to the combined species presence/ absence data ###
#####################################################################################

# 2023-02-13
# split from the GAM evaluation script written at Tor 2022-12-24

# load libraries
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(raster)

# load contemporary environmental covariates
sst <- raster("data/contemporary covariates/temperature_025_aus_mean_wgs84")
sal <- raster("data/contemporary covariates/salinity_025_aus_mean_wgs84")
ssh <- raster("data/contemporary covariates/surfaceheight_025_aus_mean_wgs84")
mld <- raster("data/contemporary covariates/mixedlayerdepth_025_aus_mean_wgs84")
sic <- raster("data/contemporary covariates/seaice_025_aus_mean_wgs84") # sea ice is already polar projection
sst_grad <- terrain(sst, opt = "slope", unit = "radians")
ssh_grad <- terrain(ssh, opt = "slope", unit = "radians")

# mixed layer depth raster only extends to 78.5 S so extent to match
mld <- extend(mld, sst, value = NA)

# stack the covariates into single object - not sic as it's already projected
covs <- stack(sst, sst_grad, sal, ssh, ssh_grad, mld)
names(covs) <- c("sst", "sst_grad", "sal", "ssh", "ssh_grad", "mld")

# project to southern hemisphere polar projection
# need to clip to southern hemisphere first
ext <- extent(-180, 180, -90, -40)
covs <- crop(covs, ext)
template <- projectRaster(from = covs, to = sic, alignOnly = T) # make a template that aligns to the sic grid but wider for whole southern ocean
covs <- projectRaster(from = covs, to = template) 

sic <- extend(sic, covs, value = 0)
sic <- mask(sic, mask = subset(covs, 1), maskvalue = NA)
names(sic) <- "sic"
covs <- stack(covs, sic)

# download some bathymetry data as an example from marmap
# use 10 minute resolution as proxy of <0.25 degrees
# this can be downgraded to 25 km x 25 km
bathy <- marmap::getNOAA.bathy(lon1 = -180,
                               lon2 = 180,
                               lat1 = -90,
                               lat2 = -40,
                               resolution = 10)
bat <- marmap::as.raster(bathy)
bat <- projectRaster(from = bat, to = covs)
bat <- mask(bat, mask = subset(covs, 1), maskvalue = NA)
names(bat) <- "bat"
covs <- stack(covs, bat)
plot(covs)

# load species presence data
dat <- readRDS("data/presence_absence_data_10k.rds")
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

dat[c("x", "y")] <- dat %>% st_as_sf(coords = c("lon", "lat")) %>% 
  sf::st_set_crs(4326) %>%
  st_transform(prj) %>%
  st_coordinates()

# extract environmental covariates
pulled <- raster::extract(covs, as.matrix(dat[c("x", "y")])) %>% as_tibble()

# append to data frame
dat <- dat %>% bind_cols(pulled)

# drop NAs
dat <- dat %>% drop_na(sst, sst_grad, sal, ssh, ssh_grad, mld, sic, bat)

# output to file
saveRDS(dat, "data/presence_absence_data_10k_with_covariates_2023-04-05.rds")

# also output the covariates for use in spatial predictions
saveRDS(covs, "data/covariate_stack.rds")

# ends


# What to do with sea ice concentration
# Is it informative to just use concentration?
# This is summer concentration - what does that tell us?
# Is there a better proxy for what sea ice is doing?
# One option is distance to sea ice edge - i.e. distance to 15% contour...
# what about days of open water?
# what about timing of sea ice retreat?

#contour <- 15
#x <- sic
#x[x > contour] <- contour # reclassify ice to either 15% or open water
#x[x < contour] <- 0
#x[is.na(x)] <- 3 # replace NA with a different value for gridDistance function

# generate distance raster
# origin is where to calculate distance from
# omit is the land to ignore - set this as 3 as per
# https://www.r-bloggers.com/2020/02/three-ways-to-calculate-distances-in-r/
#d <- gridDistance(x,
#                  origin = contour,
#                  omit = 3)
#d <- d/1000 # metres to km

#sid <- d
#names(sid) <- "sid"

# ends


