#########################################################
### Generate covariate anomalies for each climatology ###
#########################################################

# Use 'change factor protocol' to base shift contemporary observed covariate data base on esgf data
# Calculate a climate anomaly as the difference between modelled historical period and modelled lgm period
# This anomaly is then added to the satellite observation data
# For more details see Freer et al. 2019 Diversity and Distributions supplement 1

# calculate anomaly as difference between historical and lgm for only models with both???
# this would reduce the number of models available...
# shist <- subset(shist, which(names(shist) %in% names(slgm)))

# instead calculate ensemble mean across all available models
# then calculate the anomaly from the ensemble mean
# this means more data available...

# load libraries
require(raster)

# load helper function
source("R/cmip_stack.r")

##############
## Salinity ##
##############

# create stack of historical period data
# one layer per model and named
shist <- cmip_stack(path = "data/CMIP6/",
                    var = "sos",
                    period = "historical")
shist_mean <- mean(shist)

# create stack of lgm period data
slgm <- cmip_stack(path = "data/CMIP6/",
                   var = "sos",
                   period = "lgm")
slgm_mean <- mean(slgm)

# create anomaly
sos_anom <- slgm_mean - shist_mean

# plot to check
par(mfrow = c(3, 1))
plot(shist_mean)
plot(slgm_mean)
plot(sos_anom)

#################
## Temperature ##
#################

# create stack of historical period data
# one layer per model and named
thist <- cmip_stack(path = "data/CMIP6/",
                    var = "tos",
                    period = "historical")
thist_mean <- mean(thist)

# create stack of lgm period data
tlgm <- cmip_stack(path = "data/CMIP6/",
                   var = "tos",
                   period = "lgm")
tlgm_mean <- mean(tlgm)

# create anomaly
tos_anom <- tlgm_mean - thist_mean

# plot to check
plot(thist_mean)
plot(tlgm_mean)
plot(tos_anom)

########################
## Sea surface height ##
########################

# create stack of historical period data
# one layer per model and named
zhist <- cmip_stack(path = "data/CMIP6/",
                    var = "zos",
                    period = "historical")
zhist_mean <- mean(zhist)

# create stack of lgm period data
zlgm <- cmip_stack(path = "data/CMIP6/",
                   var = "zos",
                   period = "lgm")
zlgm_mean <- mean(zlgm)

# create anomaly
zos_anom <- zlgm_mean - zhist_mean

# plot to check
plot(zhist_mean)
plot(zlgm_mean)
plot(zos_anom)

###########################
## Sea ice concentration ##
###########################

# create stack of historical period data
# one layer per model and named
ihist <- cmip_stack(path = "data/CMIP6/",
                    var = "siconc",
                    period = "historical")
ihist_mean <- mean(ihist)

# create stack of lgm period data
ilgm <- cmip_stack(path = "data/CMIP6/",
                   var = "siconc",
                   period = "lgm")
ilgm_mean <- mean(ilgm)

# create anomaly
siconc_anom <- ilgm_mean - ihist_mean

# plot to check
par(mfrow = c(3, 1))
plot(ihist_mean)
plot(ilgm_mean)
plot(siconc_anom)

#######################
## Mixed layer depth ##
#######################

# create stack of historical period data
# one layer per model and named
mhist <- cmip_stack(path = "data/CMIP6/",
                    var = "mlotst",
                    period = "historical")
mhist_mean <- mean(mhist)

# create stack of lgm period data
mlgm <- cmip_stack(path = "data/CMIP6/",
                   var = "mlotst",
                   period = "lgm")
mlgm_mean <- mean(mlgm)

# create anomaly
mlotst_anom <- mlgm_mean - mhist_mean

# plot to check
plot(mhist_mean)
plot(mlgm_mean)
plot(mlotst_anom)

####################
### Output files ###
####################

par(mfrow = c(2,3))
plot(sos_anom)
plot(tos_anom)
plot(zos_anom)
plot(siconc_anom)
plot(mlotst_anom)

lgm_anom <- stack(sos_anom, tos_anom, zos_anom, siconc_anom, mlotst_anom)
names(lgm_anom) <- c("sos", "tos", "zos", "siconc", "mlotst")
writeRaster(lgm_anom, "data/cmip_lgm_anom")

# ends
