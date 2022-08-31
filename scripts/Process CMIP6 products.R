################################################################
### Process downloaded CMIP6 products ready for ENM analysis ###
################################################################

# CMIP6 products downloaded via ESGF are stored as .nc netcdf files
# These have been re-gridded from native grid (gn) to 1x1 degree using cdo bash script
# Format is monthly means for n years of data
# From these create seasonal averages for austral summer - October to March

# for lgm there are varying numbers of files:
# AWI has 10 files - 12 months and 10 years each
# MIROC has 1200 layers in one file
# MPI has 25 files each with 240 layers - 500 years of data!

# double check MPI file and delete all but last 100 years using bash script

# helper function will process CMIP6 data and save seasonal average
# based on indexing months for 1200 files
# ensure that 1200 files in month-year order are loaded

require(tidyverse)
require(raster)
require(ncdf4)
source("R/process_CMIP6.R")

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

# AWI-ESM-1-1-LR
batch_lgm(var = "sos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/AWI")
batch_lgm(var = "tos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/AWI")
batch_lgm(var = "zos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/AWI")
batch_lgm(var = "siconc", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/AWI")
batch_lgm(var = "mlotst", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/AWI")

# INM-CM4-8
batch_lgm(var = "tos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/INM")

# MIROC-ES2L
batch_lgm(var = "sos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MIROC")
batch_lgm(var = "tos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MIROC")
batch_lgm(var = "zos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MIROC")
batch_lgm(var = "siconc", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MIROC")

# MPI-ESM1-2-LR
batch_lgm(var = "sos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MPI-M")
batch_lgm(var = "tos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MPI-M")
batch_lgm(var = "zos", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MPI-M")
batch_lgm(var = "siconc", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MPI-M")
batch_lgm(var = "mlotst", path = "/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/MPI-M")

# ends

# alternative approach was to code the following for each variable/ model combo

## AWI surface temperature ##
local_files <- list.files("/Volumes/Nifty 128/ANTSIE data/CMIP6/PMIP/AWI", recursive = T, full.names = T)
local_files <- local_files[grepl("bil_1x1", local_files)]
local_files <- local_files[grepl("tos", local_files)]
tos_AWI_ESM_1_1_LR_lgm_r1i1p1f1_bil_1x1_aus_mean_wgs84 <- process_lgm(local_files, project = F)
tos_AWI_ESM_1_1_LR_lgm_r1i1p1f1_bil_1x1_aus_mean_stere <- process_lgm(local_files, project = T)

# plot
par(mfrow = c(2,1))
plot(tos_AWI_ESM_1_1_LR_lgm_r1i1p1f1_bil_1x1_aus_mean_wgs84)
plot(tos_AWI_ESM_1_1_LR_lgm_r1i1p1f1_bil_1x1_aus_mean_stere)

# save outputs
writeRaster(tos_AWI_ESM_1_1_LR_lgm_r1i1p1f1_bil_1x1_aus_mean_wgs84, "data/environmental covariates - past/tos_AWI-ESM-1-1-LR_lgm_r1i1p1f1_bil_1x1_aus_mean_wgs84", overwrite = TRUE)
writeRaster(tos_AWI_ESM_1_1_LR_lgm_r1i1p1f1_bil_1x1_aus_mean_stere, "data/environmental covariates - past/tos_AWI-ESM-1-1-LR_lgm_r1i1p1f1_bil_1x1_aus_mean_stere", overwrite = TRUE)

# even older code...

# need to process 4D files to single layer rasters for ENM analysis

# for PMIP products average over all pseudo-years for months of interest
# for CMIP products years are 1850 to 2014
# first extract years of interest (1950-2010?) before averaging over months 

# load libraries
require(tidyverse)
#remotes::install_github("eliocamp/rcmip6", force = T)
require(rcmip6)
require(ncdf4)
require(raster)

# check products stored locally
cmip_root_set("data")             # specify root directory
local_files <- cmip_available(root = cmip_root_get()) %>% as_tibble

# pull file paths
filenames <- local_files %>% pull(files) %>% unlist()
# pull only .nc files - drop the .log files
filenames <- filenames[grepl("[.]nc", filenames)]

# split filenames to historical, midHolocene and lgm as these will be processed differently
files_hist <- filenames[grepl("historical", filenames)]
files_lgm <- filenames[grepl("lgm", filenames)]
files_holo <- filenames[grepl("midHolocene", filenames)]

# helper function to process the CMIP6 data into seasonal average
source("R/process_CMIP6.R")

lgm_intpp  <- process_lgm(path = files_lgm[1])
lgm_o2os   <- process_lgm(path = files_lgm[2])
lgm_sos    <- process_lgm(path = files_lgm[3])
lgm_tos    <- process_lgm(path = files_lgm[4])
lgm_siconc <- process_lgm(path = files_lgm[5])

holo_intpp  <- process_lgm(path = files_holo[1])
holo_o2os   <- process_lgm(path = files_holo[2])
holo_sos    <- process_lgm(path = files_holo[3])
holo_tos    <- process_lgm(path = files_holo[4])
holo_siconc <- process_lgm(path = files_holo[5])

hist_intpp  <- process_hist(path = files_hist[1])
hist_o2os   <- process_hist(path = files_hist[2])
hist_sos    <- process_hist(path = files_hist[3])
hist_tos    <- process_hist(path = files_hist[4])
hist_siconc <- process_hist(path = files_hist[5])

# output these files to locally stored rasters
writeRaster(hist_intpp, "data/environmental covariates - present/intpp_Omon_MIROC-ES2L_historical_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(hist_o2os, "data/environmental covariates - present/o2os_Omon_MIROC-ES2L_historical_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(hist_sos, "data/environmental covariates - present/sos_Omon_MIROC-ES2L_historical_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(hist_tos, "data/environmental covariates - present/tos_Omon_MIROC-ES2L_historical_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(hist_siconc, "data/environmental covariates - present/siconc_SImon_MIROC-ES2L_historical_r1i1p1f2_gr1", overwrite = TRUE)

writeRaster(holo_intpp, "data/environmental covariates - past/intpp_Omon_MIROC-ES2L_midHolocene_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(holo_o2os, "data/environmental covariates - past/o2os_Omon_MIROC-ES2L_midHolocene_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(holo_sos, "data/environmental covariates - past/sos_Omon_MIROC-ES2L_midHolocene_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(holo_tos, "data/environmental covariates - past/tos_Omon_MIROC-ES2L_midHolocene_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(holo_siconc, "data/environmental covariates - past/siconc_SImon_MIROC-ES2L_midHolocene_r1i1p1f2_gr1", overwrite = TRUE)

writeRaster(lgm_intpp, "data/environmental covariates - past/intpp_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_o2os, "data/environmental covariates - past/o2os_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_sos, "data/environmental covariates - past/sos_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_tos, "data/environmental covariates - past/tos_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_siconc, "data/environmental covariates - past/siconc_SImon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)

# ends
