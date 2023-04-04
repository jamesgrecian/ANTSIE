########################################################################################
### Function to generate background points from at-sea buffer around input locations ###
########################################################################################

# 2023-03-31

# input sf dataframe of locations
# scaler is multiplier for maximum distance around buffer

buffered_points <- function(input_data, scaler = 2){
  
  # calculate distance for buffer
  range_distance <- input_data %>% st_distance() %>% as.numeric() # calculate distances between all points
  range_distance <- max(range_distance)/1000 # scale to km
  cut_off <- 2 * range_distance
  
  # create convex hull around range
  range_hull <- input_data %>% st_union() %>% st_convex_hull() # create convex hull around points
  
  # cookie cut out range_hull? so that it doesn't exclude Antarctica???
  
  # load in shapefile for background mask, clip and project
  world_shp <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
  CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = 4326) %>% sf::st_as_sfc()
  world_shp <- world_shp %>% sf::st_crop(CP)
  world_shp <- world_shp %>% sf::st_transform(projection(input_data))
  world_shp <- world_shp %>% st_union()
  
  # create raster of region of interest - make it a bit further north for contingency
  r <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 0, crs = 4326)
  r <- projectRaster(r, crs = projection(input_data), res = 50000)
  
  # calculate at-sea distances using Virginia's post on stack overflow
  # https://stackoverflow.com/questions/55815748/calculate-distance-raster-avoiding-land
  mask <- rasterize(sf:::as_Spatial(world_shp), r, field = -1) # create mask from land
  test <- rasterize(sf::as_Spatial(range_hull), mask, update = T) # append origin of interest
  x <- gridDistance(test, 1, omit = -1) # calculate distance from origin excluding land (-1)
  x <- x/1000 # convert from m to km
  
  x[mask == -1] <- NA # make sure that land is NA and hasn't been swamped by the convex hull
  
  x[x < cut_off] <- 9999 # identify all locations in x within cutoff distance
  out <- rasterToPolygons(x, fun = function(x){x == 9999}, dissolve = T)
  
  # sample 10*n random points within background mask shape
  pts <- out %>% st_as_sf() %>% sf::st_sample(10*nrow(input_data)) %>% st_coordinates() %>% as_tibble()
  
  # clip points to below latitude threshold
  pts[c("lon", "lat")] <- pts %>% 
    st_as_sf(coords = c("X", "Y")) %>% 
    sf::st_set_crs(projection(input_data)) %>%
    st_transform(4326) %>% 
    st_coordinates()
  
  pts <- pts %>% filter(lat <= -40) # remove points north of 40 degrees
  pts <- pts[1:nrow(input_data),1:2]
  return(pts)
}
