################################################################
### Function to process NSIDC montly data into climatologies ###
################################################################

monthly_ice_average <- function(year_min, year_max, m){
  
  # Use RCurl library to query FTP server and list files
  url = "ftp://anonymous:wjg5@sidads.colorado.edu/DATASETS/NOAA/G02135/south/monthly/geotiff"
  
  # Time series is split into yearly and monthly folders
  months <- c("01_Jan", "02_Feb", "03_Mar", "04_Apr", "05_May", "06_Jun", "07_Jul", "08_Aug", "09_Sep", "10_Oct", "11_Nov", "12_Dec")
  years <- year_min:year_max
  
  # combine dates to construct file path
  # how to add zero infront of digit
  fn <- paste0(url, "/", months[m], "/", "S_", years, str_pad(m, 2, pad = "0"), "_concentration_v3.0.tif")
  
  writeLines("loading ice data")
  
  # set up progress bar
  pb1 <- txtProgressBar(min = 1, max = length(fn), style = 3)
  
  # Loop through years and load ice data into stack
  x <- raster::stack()
  
  for (i in 1:length(fn)){
    setTxtProgressBar(pb1, i) # update progress bar
    ice <- tryCatch(raster::raster(fn[i]),
                    error = function(err) NA)
    x <- tryCatch(raster::stack(x, ice),
                  error = function(err) NA)
  }
  # define projection
  raster::projection(x) = "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  # check number of layers processed
  writeLines(paste("\n", length(year_min:year_max), "layers requested"))
  writeLines(paste(nlayers(x), "layers processed"))
  
  writeLines("processing ice data...")
  
  # 0 is ocean; 2510 pole hole; 2530 coast line; 2540 land; 2550 missing
  # 0-1000 so divide by 10 to get percentage
  x[x == 2510] <- 1000 # make pole hole 100% ice cover
  x[x>1000] <- NA
  x <- x/10
  
  x <- raster::mean(x, na.rm = T)
  return(x)
}


