#################################################
### Evaluate GAM models of habitat preference ###
#################################################

# 2022-12-24

# given a dataframe containing all the species location data and the environmental covariates
# run a GAM model and test the validity of it using the spatial k-folds cross validation proceedure...

# load libraries
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(raster)
require(mgcv)
require(gratia)

# load spatial cross-validation functions
source("R/spatial_kfolds.r")
source("R/kfoldCV.r")
source("R/rank_models.r")

# Process covariates - should put this somewhere else and then just load the raster stack...
# load environmental data
sst <- raster("data/environmental covariates - present/temperature_025_aus_mean_wgs84")
sal <- raster("data/environmental covariates - present/salinity_025_aus_mean_wgs84")
ssh <- raster("data/environmental covariates - present/surfaceheight_025_aus_mean_wgs84")
mld <- raster("data/environmental covariates - present/mixedlayerdepth_025_aus_mean_wgs84")
sic <- raster("data/environmental covariates - present/seaice_025_aus_mean_wgs84") # sea ice is already polar projection
sst_grad <- terrain(sst, opt = "slope", unit = "radians")
ssh_grad <- terrain(ssh, opt = "slope", unit = "radians")

# mixed layer depth raster only extends to 78.5 S
mld <- extend(mld, sst, value = NA)

# stack the covariates into single object - not sic
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

covs <- stack(covs, sic)

# download some bathymetry data as an example from marmap
# use 10 minute resolution as proxy of <0.25 degrees
# this can be downgraded to 25 km x 25 km
#bathy <- marmap::getNOAA.bathy(lon1 = -180,
#                               lon2 = 180,
#                               lat1 = -90,
#                               lat2 = -40,
#                               resolution = 10)
#bat <- marmap::as.raster(bathy)
bat <- raster("data/environmental covariates - present/depth_gr1")
bat <- projectRaster(from = bat, to = covs)
bat <- mask(bat, mask = subset(covs, 1), maskvalue = NA)
names(bat) <- "bat"
covs <- stack(covs, bat)
plot(covs)

# now can load in point data and run a habitat model or two
# have I already calculated offset covariates based on lgm - hist?

# load data
dat <- readRDS("data/combined_data_presence_absence.rds")
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

dat[c("x", "y")] <- dat %>% st_as_sf(coords = c("lon", "lat")) %>% 
  sf::st_set_crs(4326) %>%
  st_transform(prj) %>%
  st_coordinates()

pulled <- raster::extract(covs, as.matrix(dat[c("x", "y")])) %>% as_tibble()

dat <- dat %>% bind_cols(pulled)

dat <- dat %>% drop_na(sst, sst_grad, sal, ssh, ssh_grad, mld, sic, bat)

# check with VIF
# sst, ssh and sic or sid are all strongly correlated (with latitude)
faraway::vif(dat[c("sst", "sst_grad", "sal", "ssh", "ssh_grad", "mld", "bat", "sic")])
faraway::vif(dat[c("sst_grad", "sal", "ssh", "ssh_grad", "mld", "bat", "sic")])
faraway::vif(dat[c("sst", "sst_grad", "sal", "ssh_grad", "mld", "bat", "sic")])

cor(dat$sst_grad, dat$ssh_grad)


# need to go through all the possible combinations of 8 covariates
# check for correlations between pairs
# remove candidate models with correlated variables...
# ie sst and ssh_grad but exclude sst and ssh from same model etc...
variables <- c("sst", "sst_grad", "sal", "ssh", "ssh_grad", "mld", "sic", "bat")

# check pairwise correlations
test <- combn(variables, 2, simplify = T)
test <- t(test) %>% as_tibble()
# must be possible to append correlations to table then filter by cor > 0.7...?
test <- test %>% rowwise() %>% mutate(lala = cor(dat[c(V1)], dat[c(V2)]))
test %>% filter(lala > .7)

# can use all possible candidate models except those containing both sst and ssh

f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sst_grad,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(ssh_grad,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5) + s(sic,  bs = "cr", k = 5)

# how many terms does the formula have?
n.terms <- length(attr(terms(f), "term.labels"))
new_f <- f

for(i in 3:n.terms){
  t <- combn(attr(terms(f), "term.labels"), m = i, simplify = F)
  
  for(j in 1:length(t)){
    t2 <- t[[j]]
    temp_f <- reformulate(paste(t2, collapse = " + "), response = all.vars(attr(terms(f), "variables"))[1])
    new_f <- c(new_f, temp_f)
  }
}

sst_id <- grep("sst,", new_f)
ssh_id <- grep("ssh,", new_f)
ids <- sst_id[sst_id %in% ssh_id]
idx <- 1:length(new_f)
kp <- setdiff(idx, ids)
new_f <- new_f[kp]

# so come up with a way to generate all combinations and then filter them to exclude sst and ssh combos

# or eddie kenetic energy
# or geostrophic current anomaly

# quick demonstration of analysis approach
myct <- dat %>% filter(group == "myctophids")
#myct <- myct %>% filter(bat > -7000)
#myct <- myct %>% filter(sal > 32.5)

# fit maximum GAM to each species
# then do k-folds cross validation to select covariates
E_ant <- myct %>% filter(species == "Electrona antarctica")

# first run full model
m <- gam(new_f[[1]], E_ant, family = binomial, method = "REML", select = T)
foo <- kfoldCV(E_ant, new_f[[1]], k)
  
# create output tibble
  out <- tibble_row(model = 1,
                    formula = paste(substr(attr(terms(new_f[[1]]), "term.labels"), 3, 5), collapse = " + "),
                    AIC = AIC(m),
                    AUC = list(foo[[1]]$AUC),
                    TSS = list(foo[[1]]$TSS),
                    spatial_cv = list(foo[[2]]))
  
  # then loop through combination of terms
  for(i in 2:length(new_f)){
    m <- gam(new_f[[i]], data, family = binomial, method = "REML", select = T)
    foo <- kfoldCV(data, new_f[[i]], k)
    
    # append results to output dataframe
    out <- out %>%
      bind_rows(
        tibble_row(model = i+1,
                   formula = paste(substr(attr(terms(new_f[[i]]), "term.labels"), 3, 5), collapse = " + "),
                   AIC = AIC(m),
                   AUC = list(foo[[1]]$AUC),
                   TSS = list(foo[[1]]$TSS),
                   spatial_cv = list(foo[[2]])))
  }


  tri <- out %>% dplyr::select(model, formula, AIC, AUC, TSS)
  tri <- tri %>%
    unnest(c(AUC, TSS)) %>%
    group_by(model) %>%
    mutate(AUC = mean(AUC),
           TSS = mean(TSS)) %>%
    slice(1) %>%
    ungroup() %>% 
    arrange(desc(AUC), desc(TSS))
  
  









                
for (i in 1:nrow(test)){
  print(cor(dat[c(test$V1[i])], dat[c(test$V2[i])]))
}
  
  
cor(dat[test[[3]]])

combn(variables, 8, simplify = F)
combn(variables, 7, simplify = F)
combn(variables, 6, simplify = F)
combn(variables, 5, simplify = F)
combn(variables, 4, simplify = F)
combn(variables, 3, simplify = F)
combn(variables, 2, simplify = F)
combn(variables, 1, simplify = F)



# quick demonstration of analysis approach
myct <- dat %>% filter(group == "myctophids")
#myct <- myct %>% filter(bat > -7000)
#myct <- myct %>% filter(sal > 32.5)

# fit maximum GAM to each species
# then do k-folds cross validation to select covariates
E_ant <- myct %>% filter(species == "Electrona antarctica")
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5) + s(sic,  bs = "cr", k = 5)
output <- rank_combs(E_ant, f, k = 6)

# want to minimise AIC, maximise the AUC and TSS
trial4 <- output %>% dplyr::select(model, formula, AIC, AUC, TSS)
trial4 <- trial4 %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))


# generate plot for visualising the spatial cv
# define the bounding box
test <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = -20)) %>% 
  sf::st_as_sfc() %>%
  st_set_crs(4326)
# split the bounding box into k segments
test <- test %>% st_make_grid(n = c(k, 1))

# project the box into polar projection
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
test <- test %>% st_transform(prj)
test <- test %>% st_as_sf

ggplot() + 
  labs(title = "Coefficient of Variation",
       subtitle = "For each spatial k-fold") +
  theme_bw() +
  geom_sf(aes(), fill = NA, data = test) +
  geom_point(aes(x = x, y = y, colour = cv),
             data = output %>% unnest(spatial_cv)) +
  coord_sf(xlim = c(min(output %>% unnest(spatial_cv) %>% pull(x)),
                    max(output %>% unnest(spatial_cv) %>% pull(x))),
           ylim = c(min(output %>% unnest(spatial_cv) %>% pull(y)),
                    max(output %>% unnest(spatial_cv) %>% pull(y))),
           crs = prj) +
  scale_colour_viridis_c() +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ model)

###
### Second species...
###

E_car <- myct %>% filter(species == "Electrona carlsbergi")
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
output <- rank_models(E_car, f, k = 6)

E_car_output <- rank_combs(E_car, f, k = 6)

# want to minimise AIC, maximise the AUC and TSS
trial_2 <- E_car_output %>% dplyr::select(model, formula, AIC, AUC, TSS)
trial_2 <- trial_2 %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))



# all terms
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop bat
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop sid
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop mld
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop ssh
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop sal
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop sst
f <- PresAbs ~ s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)





