#############################################
### Helper functions to process CMIP6 data ###
#############################################

# generate austral summer mean of product from PMIP lgm run
process_lgm <- function(path, project = T){
  # lgm files just need seasonal averages extracted
  # lgm files have data for 100 years * 12 months = 1200 layers
  # generate austral summer index - October to March
  index <- c(seq(1, 1200, 12),
             seq(2, 1200, 12),
             seq(3, 1200, 12),
             seq(10, 1200, 12),
             seq(11, 1200, 12),
             seq(12, 1200, 12))
  index <- sort(index)
  
  # process...
  # loop to retain dates for checking
  lgm <- raster::stack()
  for(i in 1:length(path)){
    lgm <- raster::stack(lgm, raster::stack(path[i]))
  }
  
  cat("object contains", raster::nlayers(lgm), "data layers\n")
  
  lgm <- raster::subset(lgm, index) # subset raster stack of 1200 months to only austral summer
  lgm <- raster::calc(lgm, mean) # calculate seasonal mean
  
  if (project == T){
    # crop and project
    lgm <- raster::crop(lgm, raster::extent(-180, 180, -90, -30)) # crop to southern hemisphere
    
    lgm <- raster::projectRaster(lgm,
                         res = 100,
                         method = "bilinear",
                         crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0")
    return(lgm)
  } else {
    return(lgm)
  }
}






# generate austral summer mean of product from CMIP historical run
process_hist <- function(path){
  
  # historical files have data from 1850 to 2014
  # 165 years * 12 months = 1980 layers
  # trim to 1950 to 2010
  # layer 1201 should be 1950...
  # layer 1932 is 2010.12.16 
  hist <- stack(path)
  hist <- subset(hist, 1201:1932)
  
  # then generate austral summer index - October to March for remaining 61 years
  index <- c(seq(1, 732, 12),
             seq(2, 732, 12),
             seq(3, 732, 12),
             seq(10, 732, 12),
             seq(11, 732, 12),
             seq(12, 732, 12))
  index <- sort(index)
  
  # process...
  hist <- subset(hist, index) # subset raster stack of 732 months to only austral summer
  hist <- calc(hist, mean) # calculate seasonal mean

  # rotate and crop
  hist <- rotate(hist) # convert from 0-360 to -180-180
  hist <- crop(hist, extent(-180, 180, -90, -30)) # crop to southern hemisphere and then plot
  
  # project
  hist <- projectRaster(hist,
                        res = 50,
                        method = "bilinear",
                        crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0")
  return(hist)
}