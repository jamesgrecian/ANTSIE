<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data

Environmental data for the LGM, midHolocene and historical periods were
download from the ESGF
[portal](https://esgf-node.llnl.gov/search/cmip6/) using the `rcmip6`
package in `R`.

Monthly mean sea surface salinity (sos), sea surface temperature (tos),
sea surface height (zos), sea ice concentration (siconc) and mixed layer
depth (mlotst) climatologies were downloaded on the native grid (gn) and
then regularised to a 1x1 degree grid using [Climate Data
Operators](https://code.mpimet.mpg.de/projects/cdo) bash scripts.

The outputed .nc files contain monthly global climatologies at a 1x1
degree spatial resolution. Historical period data run from 1850 to 2014
(12 months x 164 years = 1968 layers). LGM period data were for 100
years (12 months x 100 years = 1200 layers). MPI-ESM1-2-LR model outputs
for the LGM contained 500 replicate years. Processing was limited to the
final 100 years of data using bash script
[here](/scripts/bash%20trim%20MPI.sh). Due to the size of these files
they are not available in this repo.

From these data we generated Austral summer climatologies
(October-March) using a custom helper function
[process_CMIP6.R](/R/process_CMIP6.R). These files are stored in the
CMPI6 folder using [Git Large File Storage](https://git-lfs.com).
