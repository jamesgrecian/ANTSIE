################################################################
### Process downloaded CMIP6 products ready for ENM analysis ###
################################################################

# CMIP6 products are stored as netcdf files
# need to process 4D files to single layer rasters for ENM analysis

# for PMIP products average over all pseudo-years for months of interest
# for CMIP products years are 1850 to 2014
# first extract years of interest (1950-2010?) before averaging over months 

# load libraries
require(tidyverse)
# remotes::install_github("eliocamp/rcmip6")
require(rcmip6)
require(ncdf4)
require(raster)

# check products stored locally
cmip_root_set("data/")             # specify root directory
local_files <- cmip_available(root = cmip_root_get()) %>% as_tibble

# pull file paths
filenames <- local_files %>% pull(files) %>% unlist()
# pull only .nc files - drop the .log files
filenames <- filenames[grepl("[.]nc", filenames)]
# split filenames to historical and lgm as these will be processed differently
files_hist <- filenames[grepl("historical", filenames)]
files_lgm <- filenames[grepl("lgm", filenames)]

# helper function to process the CMIP6 data into seasonal average
source("R/process_CMIP6.R")

lgm_intpp  <- process_lgm(path = files_lgm[1])
lgm_o2os   <- process_lgm(path = files_lgm[2])
lgm_sos    <- process_lgm(path = files_lgm[3])
lgm_tos    <- process_lgm(path = files_lgm[4])
lgm_siconc <- process_lgm(path = files_lgm[5])

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

writeRaster(lgm_intpp, "data/environmental covariates - past/intpp_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_o2os, "data/environmental covariates - past/o2os_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_sos, "data/environmental covariates - past/sos_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_tos, "data/environmental covariates - past/tos_Omon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)
writeRaster(lgm_siconc, "data/environmental covariates - past/siconc_SImon_MIROC-ES2L_lgm_r1i1p1f2_gr1", overwrite = TRUE)

# ends
