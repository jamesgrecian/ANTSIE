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

# Process data for the last glacial maximum

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


# Repeat the process the historical period data
# historical files have data from 1850 to 2014
# 165 years * 12 months = 1980 layers
# trim to 1950 to 2010
# layer 1201 should be 1950...
# layer 1932 is 2010.12.16 

# AWI-ESM-1-1-LR
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/AWI")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/AWI")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/AWI")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/AWI")
batch_hist(var = "mlotst", path = "/Users/home/CMIP6/CMIP/AWI")

# INM-CM4-8
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/INM")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/INM")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/INM")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/INM")

# MIROC
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/MIROC")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/MIROC")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/MIROC")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/MIROC")

# MPI-M
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/MPI-M")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/MPI-M")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/MPI-M")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/MPI-M")
batch_hist(var = "mlotst", path = "/Users/home/CMIP6/CMIP/MPI-M")

# IPSL
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/IPSL")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/IPSL")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/IPSL")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/IPSL")
batch_hist(var = "mlotst", path = "/Users/home/CMIP6/CMIP/IPSL")

# ACCESS
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/CSIRO")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/CSIRO")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/CSIRO")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/CSIRO")
batch_hist(var = "mlotst", path = "/Users/home/CMIP6/CMIP/CSIRO")

# MRI
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/MRI")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/MRI")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/MRI")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/MRI")
batch_hist(var = "mlotst", path = "/Users/home/CMIP6/CMIP/MRI")

# CESM2
batch_hist(var = "sos", path = "/Users/home/CMIP6/CMIP/NCAR")
batch_hist(var = "tos", path = "/Users/home/CMIP6/CMIP/NCAR")
batch_hist(var = "zos", path = "/Users/home/CMIP6/CMIP/NCAR")
batch_hist(var = "siconc", path = "/Users/home/CMIP6/CMIP/NCAR")
batch_hist(var = "mlotst", path = "/Users/home/CMIP6/CMIP/NCAR")

# next step is to create multi-model ensemble means from these products...

# ends
