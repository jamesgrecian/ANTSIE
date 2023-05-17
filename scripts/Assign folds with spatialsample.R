#################################################
### Use spatialsample package to assign folds ###
#################################################

# 2023-04-10
# modified 2023-04-17 to make PresAbs a factor

# More info on package here: https://spatialsample.tidymodels.org/index.html
# spatialsample::spatial_block_cv function seems to work similarlily to block_CV
# however the function has fewer issues partitioning the data
# this includes being able to assign k = 10 folds rather than k = 6 folds
# literature suggests k = 10 folds is better...
# for now focus on 500 km block size - this could be varied later
# also question over whether block size should vary per species... adds complexity
# blockCV paper suggests matching block size to spatial autocorrelation...?

# spatal_block_CV will not partition data correctly on first go
# sometimes one block will contain no presences
# this will cause issues later down the line - can't get AUC if there are no presences to test against
# solution is to run while loop and only stop when each block has `some` presences in it...

# load libraries
require(tidyverse)
require(spatialsample)
require(sf)

# define projection
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# load data
dat <- readRDS("data/presence_absence_data_10k_with_covariates_2023-04-05.rds")

# make PresAbs a factor
dat <- dat |> mutate(PresAbs = factor(PresAbs))

# only consider cephalopod species included in the Xavier paper
drop_id <- c("Illex argentinus", "Lycoteuthis lorigera", "Teuthowenia pellucida",
             "Mastigoteuthis psychrophila", "Abraliopsis gilchristi", "Semirossia patagonica",
             "Batoteuthis skolops", "Chiroteuthis veranyi", "Galiteuthis suhmi",
             "Grimpoteuthis megaptera", "Histioteuthis miranda", "Moroteuthis knipovitchi",
             "Parateuthis tunicata", "Taningia danae")
dat <- dat |> filter(!species %in% drop_id)

# make sf object
dat_sf <- dat |> st_as_sf(coords = c("x", "y"), crs = prj)

# initialise list to store folds
folds <- NULL 

# when PresAbs is a factor values are stored in check as 1 2 not 0 1
# -1 to revert to 0 1 so that summing to 20 still works - not summing to 60!

for(k in c(1:19, 22:27)){ # issues with species 20 and 21
  print(k)
  check <- NULL
  while(sum(check-1) < 20){ # 10 folds and need presences in both analysis and assessment groups 10 * 2 = 20
  check <- NULL
  sub_df <- dat_sf %>% filter(species == unique(species)[k])
  block_folds <- spatial_block_cv(sub_df, v = 10, square = F, cellsize = 500000)
  for (j in 1:10){
    check <- c(check,
               block_folds$splits[[j]] %>% analysis %>% pull(PresAbs) %>% unique,
               block_folds$splits[[j]] %>% assessment %>% pull(PresAbs) %>% unique)
  }
}
folds[[k]] <- block_folds
}

# species 20 "Kondakovia longimana" will get assigned if given more time for while loop to run...
k <- 20
check <- NULL
while(sum(check-1) < 20){
  check <- NULL
  sub_df <- dat_sf %>% filter(species == unique(species)[k])
  block_folds <- spatial_block_cv(sub_df, v = 10, square = F, cellsize = 500000)
    for (j in 1:10){
      check <- c(check,
                 block_folds$splits[[j]] %>% analysis %>% pull(PresAbs) %>% unique,
                 block_folds$splits[[j]] %>% assessment %>% pull(PresAbs) %>% unique)
    }
  }

folds[[k]] <- block_folds

# species 21 "Loligo gahi" doesn't get assigned right - too clustered
# probably reason to exclude from further analysis...
#k <- 21
#check <- NULL
#while(sum(check) < 20){
#  check <- NULL
#  sub_df <- dat_sf %>% filter(species == unique(species)[k])
#  block_folds <- spatial_block_cv(sub_df, v = 10, square = F, cellsize = 500000)
#  for (j in 1:10){
#    check <- c(check,
#               block_folds$splits[[j]] %>% analysis %>% pull(PresAbs) %>% unique,
#               block_folds$splits[[j]] %>% assessment %>% pull(PresAbs) %>% unique)
#  }
#}
# folds[[k]] <- block_folds

# plot patagonian squid to illustrate clustering...
#p1 <- ggplot() + 
#  theme_bw() +
#  ggtitle(unique(dat_sf$species)[k]) +
#  geom_sf(aes(), data = world_shp) +
#  geom_sf(aes(colour = factor(PresAbs)), data = dat_sf %>% filter(species == unique(species)[k])) +
#  coord_sf(xlim = c(-6500000, 6500000), ylim = c(-6500000, 6500000)) +
#  facet_wrap(~PresAbs)

#quartz(width = 9, height = 5)
#print(p1)
#quartz.save(file = "test.jpg",
#            type = "jpeg",
#            dev = dev.cur(),
#            dpi = 500)
#dev.off()

# save folds object
saveRDS(folds, "data/folds.rds")

# ends