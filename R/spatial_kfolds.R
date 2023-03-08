#######################################################################################
### Function to generate indices for a spatial approach to k folds cross validation ###
#######################################################################################

# 2022-12-13

# given a data frame of points and k number of spatial folds

# define a bounding box for the southern ocean
# cut that shape into k slices
# project to correct projection for south pole
# query the location of each data point against k slices
# return the id and append which slice for each data point
# these numbers can then be used to segment the data for k-folds cross-validation

# 2023-02-27
# add option switch for blockCV instead

# 2023-03-06
# add set seed to blockCV function
# this prevents a new set of spatial CV blocks being selected at each run of the function
# seems to fix the issue of NA presences being predicted and causing dismo::evaluate to fail

spatial_kfolds <- function(dat, k, type = c("pie", "blockCV")){
  
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
  
  # convert data to polar projection
  dat_sf <- dat %>% st_as_sf(coords = c("x", "y"), crs = prj)
  
  if (type == "pie"){
  
    # find which segment each data point is within
    id_1 <- dat_sf %>% sf::st_within(test, sparse = F) # matrix of all answers
    id_1 <- which(id_1, arr.ind = T) # store TRUE answers
    id_1 <- id_1 %>% as_tibble() %>% arrange(row)
  
    # in very rare cases point might be on the line
    id_2 <- dat_sf %>% sf::st_touches(test, sparse = F)
    id_2 <- which(id_2, arr.ind = T)
    # st_touches will return the id of both sections - as points straddle both
    # pull the first one
    id_2 <- id_2 %>% as_tibble() %>% arrange(row)
    id_2 <- id_2 %>% filter(col == unique(col)[1])
  
    # append these to the first set of indices
    ids <- id_1 %>% bind_rows(id_2)
    ids <- ids %>% arrange(row)
    # return indices as per dismo kfold function
    id <- ids$col
    
  }
  
  if (type == "blockCV"){
    capture.output(sbl <- cv_spatial(x = dat_sf,
                      column = "PresAbs",
                      k = k,
                      selection = "random",
                      iteration = 50, # can up this to 250 to ensure equal points between each block - less likely to fail?
                      progress = F,
                      report = F,
                      plot = F,
                      seed = 1)) # has setting the seed fixed the issue?
    id <- sbl$folds_ids
  }
  return(id)
}
