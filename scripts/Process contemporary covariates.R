#########################################
### Processing contemporary data sets ###
#########################################

# 2022-11-01

# Characterise contemporary environmental conditions
# These will be used to drive species distribution models
# Use climatological means over Austral winter: October - March
# Find finest scale spatial resolution
# Find longest available time series

# load libraries
require(tidyverse)
require(raster)

# Temperature, salinity and mixed layer depth available from World Ocean Atlas 2018 database (Boyer et al. 2018).
# https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/
# https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DOC/woa18documentation.pdf
# t_an is objectively analyzed climatological mean
# level 1 is 0m depth
# download monthly files and then combine to create seasonal mean

# WOA temperature
fn <- list.files("/Volumes/Nifty 128/ANTSIE data/WOA", full.names = T, pattern = "t")
temp <- stack()
for (i in 1:length(fn)){
  foo <- brick(fn[i], varname = "t_an")
  foo <- subset(foo, 1)
  temp <- stack(temp, foo)
}
temp <- mean(temp)
plot(temp)
writeRaster(temp, "data/environmental covariates - present/temperature_025_aus_mean_wgs84", overwrite = TRUE)

# WOA salinity
fn <- list.files("/Volumes/Nifty 128/ANTSIE data/WOA", full.names = T, pattern = "s")
sal <- stack()
for (i in 1:length(fn)){
  foo <- brick(fn[i], varname = "s_an")
  foo <- subset(foo, 1)
  sal <- stack(sal, foo)
}
sal <- mean(sal)
plot(sal)
writeRaster(sal, "data/environmental covariates - present/salinity_025_aus_mean_wgs84", overwrite = TRUE)

# WOA mixed-layer depth
fn <- list.files("/Volumes/Nifty 128/ANTSIE data/WOA", full.names = T, pattern = ".csv")
mld <- stack()
for (i in 1:length(fn)){
  foo <- read_csv(fn[i], skip = 1)
  names(foo) <- c("lat", "lon", "mld")
  foo <- rasterFromXYZ(cbind(foo$lon, foo$lat, foo$mld))
  mld <- stack(mld, foo)
}
mld <- mean(mld)
plot(mld)
writeRaster(mld, "data/environmental covariates - present/mixedlayerdepth_025_aus_mean_wgs84", overwrite = TRUE)

# Sea surface height data were extracted from SSALTO/ DUACS and
# sea ice concentration data were extracted from NSIDC.
# From these data we calculated austral summer (October to March) climatologies.

#######################################################################
### Create sea ice climatology by averaging over years for each day ###
#######################################################################

# long-term average for each julian day

# output is a stack for a year - 1 layer per day
# each layer is the average ice conditions on that day over the whole period

# so load in all 1st Jan take average over 20 years
# then 2nd January etc
require(tidyverse)
require(stringr)
require(raster)

source("R/monthly_ice_average.R")
nsidc_01 <- monthly_ice_average(year_min = 1979, year_max = 2020, m = 1)
nsidc_02 <- monthly_ice_average(year_min = 1979, year_max = 2020, m = 2)
nsidc_03 <- monthly_ice_average(year_min = 1979, year_max = 2020, m = 3)
nsidc_10 <- monthly_ice_average(year_min = 1979, year_max = 2020, m = 10)
nsidc_11 <- monthly_ice_average(year_min = 1979, year_max = 2020, m = 11)
nsidc_12 <- monthly_ice_average(year_min = 1979, year_max = 2020, m = 12)
# January 1988
# December 1987
sic <- stack(nsidc_01, nsidc_02, nsidc_03, nsidc_10, nsidc_11, nsidc_12)
sic <- mean(sic)
writeRaster(sic, "data/environmental covariates - present/seaice_025_aus_mean_wgs84", overwrite = TRUE)

# AVISO Mean sea level height above geoid
# https://www.aviso.altimetry.fr/index.php?id=1526
fn <- list.files("/Volumes/Nifty 128/ANTSIE data/AVISO", full.names = T)
ssh <- stack(fn, varname = "adt")
ssh <- rotate(ssh)
plot(ssh)
ssh <- mean(ssh)
writeRaster(ssh, "data/environmental covariates - present/surfaceheight_025_aus_mean_wgs84", overwrite = TRUE)



  