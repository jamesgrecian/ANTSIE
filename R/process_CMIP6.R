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
process_hist <- function(path, project = T){
  
  # historical files have data from 1850 to 2014
  # 165 years * 12 months = 1980 layers
  # trim to 1950 to 2010
  # layer 1201 should be 1950...
  # layer 1932 is 2010.12.16 
  hist <- stack(path)
  
  cat("object contains", raster::nlayers(hist), "data layers\n") # check it contains 1980 layers at start
  
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

  if (project == T){
    # crop and project
    hist <- raster::crop(hist, raster::extent(-180, 180, -90, -30)) # crop to southern hemisphere
    
    hist <- raster::projectRaster(hist,
                                  res = 100,
                                  method = "bilinear",
                                  crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0")
    return(hist)
  } else {
    return(hist)
  }
}

# write a function to automate the process_lgm function for each model
# saves having to write out the filename multiple times

# batch load data function
batch_lgm <- function(path, var){
  local_files <- list.files(path, recursive = T, full.names = T)
  local_files <- local_files[grepl("bil_1x1", local_files)]
  local_files <- local_files[grepl(var, local_files)]
  rast_wgs84 <- process_lgm(local_files, project = F)
  rast_stere <- process_lgm(local_files, project = T)
  
  # plot
  par(mfrow = c(2,1))
  plot(rast_wgs84, main = var)
  plot(rast_stere, main = var)
  
  # save outputs
  names <- tail(unlist(strsplit(local_files, split = "/")), 1)
  nm <- paste0(unlist(strsplit(names, split = "_"))[3:5], sep = "_", collapse = "")
  writeRaster(rast_wgs84, paste0("data/environmental covariates - past/", var, "_", nm, "bil_1x1_aus_mean_wgs84"), overwrite = TRUE)
  writeRaster(rast_stere, paste0("data/environmental covariates - past/", var, "_", nm, "bil_1x1_aus_mean_stere"), overwrite = TRUE)
}

# batch load data function
batch_hist <- function(path, var){
  local_files <- list.files(path, recursive = T, full.names = T)
  local_files <- local_files[grepl("bil_1x1", local_files)]
  local_files <- local_files[grepl(var, local_files)]
  rast_wgs84 <- process_hist(local_files, project = F)
  rast_stere <- process_hist(local_files, project = T)
  
  # plot
  par(mfrow = c(2,1))
  plot(rast_wgs84, main = var)
  plot(rast_stere, main = var)
  
  # save outputs
  names <- tail(unlist(strsplit(local_files, split = "/")), 1)
  nm <- paste0(unlist(strsplit(names, split = "_"))[3:5], sep = "_", collapse = "")
  writeRaster(rast_wgs84, paste0("data/environmental covariates - present/", var, "_", nm, "bil_1x1_aus_mean_wgs84"), overwrite = TRUE)
  writeRaster(rast_stere, paste0("data/environmental covariates - present/", var, "_", nm, "bil_1x1_aus_mean_stere"), overwrite = TRUE)
}



