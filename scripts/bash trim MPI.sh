###############################################
### Trim files downloaded for MPI-ESM1-2-LR ###
###############################################

# MPI models provide 500 years of data for lgm
# Trim files to only consider final 100 years of data
# list available files and delete first 80

# process netcdf files from AWI-ESM-1-1-LR for salinity
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/sos/gn/20190710/
find sos* | head -80 | xargs rm -rf

# process netcdf files from MPI-ESM1-2-LR for temperature
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/tos/gn/20190710/
find tos* | head -80 | xargs rm -rf

# process netcdf files from MPI-ESM1-2-LR for sea surface height
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/zos/gn/20190710/
find zos* | head -80 | xargs rm -rf

# process netcdf files from MPI-ESM1-2-LR for sea ice concentration
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/SImon/siconc/gn/20190710/
find siconc* | head -80 | xargs rm -rf

# process netcdf files from MPI-ESM1-2-LR for mixed layer depth
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/mlotst/gn/20190710/
find mlotst* | head -80 | xargs rm -rf




  