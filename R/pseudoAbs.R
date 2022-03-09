###############################################
### Function for generating pseudo-absences ###
###############################################

# WJ Grecian
# 2022-03-09

# inputs are bounding box lons and lats, projection and number of samples drawn
# function projects bounding box, removes land and then uses st_sample to generate n random points

pseudoAbs <- function(x_min = -180, x_max = 180, y_min = -80, y_max = -40, projection, n = 1000){
  
  # avoid sf manipulation errors
  sf::sf_use_s2(FALSE)
  # https://stackoverflow.com/questions/68478179/how-to-resolve-spherical-geometry-failures-when-joining-spatial-data
  
  # define region of interest to generate absences within
  # set this up as a rectangular polygon then project to form a donut around Antarctica
  # https://github.com/r-spatial/sf/issues/1078
  ext <- raster::extent(x_min, x_max, y_min, y_max)
  ext <- ext %>% sf::st_bbox(crs = 4326)
  ext <- ext %>% sf::st_as_sfc() 
  # remove the crs so we get planar ops, but restore immediate
  crs <- sf::st_crs(ext)
  ext <- sf::st_set_crs(ext, NA)
  ext <- ext %>% sf::st_segmentize(1) %>% sf::st_set_crs(crs)
  ext <- ext %>% sf::st_transform(crs = prj) 
  
  # cookie cut a land shapefile from this donut so absences aren't generated on land
  # Load in shapefile from rworldmap, clip to south and project
  world_shp <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
  
  # clip world shapefile to southern hemisphere
  CP <- sf::st_bbox(c(xmin = -180,
                      xmax = 180,
                      ymin = -90,
                      ymax = 0), crs = 4326) %>%
    sf::st_as_sfc()
  
  world_shp <- world_shp %>% sf::st_crop(CP)
  world_shp <- world_shp %>% sf::st_transform(prj) %>% sf::st_buffer(0)
  world_shp <- world_shp %>% sf::st_union() # merge to single feature to simplify intersection with donut
  
  # clip the ext donut by the southern hemisphere shapefile
  cookie <- sf::st_difference(ext, world_shp)
  
  # random sample of n points within cookie shape
  pabs <- cookie %>% sf::st_sample(n)
  pabs <- pabs %>% sf::st_coordinates()
  return(pabs)
}

# ends
