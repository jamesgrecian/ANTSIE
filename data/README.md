<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data

## CMIP

Environmental data for the LGM, midHolocene and historical periods were
download from the ESGF
[portal](https://esgf-node.llnl.gov/search/cmip6/) with the `rcmip6`
package using the [Download CMIP6
script](scripts/Download%20CMIP6%20script.R).

Monthly mean sea surface salinity (sos), sea surface temperature (tos),
sea surface height (zos), sea ice concentration (siconc) and mixed layer
depth (mlotst) climatologies were downloaded on the native grid (gn) and
then regularised to a 1x1 degree grid using [Climate Data
Operators](https://code.mpimet.mpg.de/projects/cdo) bash
[scripts](scripts/).

The outputed .nc files contain monthly global climatologies at a 1x1
degree spatial resolution. Historical period data run from 1850 to 2014
(12 months x 164 years = 1968 layers). LGM period data were for 100
years (12 months x 100 years = 1200 layers). MPI-ESM1-2-LR model outputs
for the LGM contained 500 replicate years. Processing was limited to the
final 100 years of data using bash script
[here](/scripts/bash%20trim%20MPI.sh). Due to the size of these files
they are not available in this repo.

From these data we generated Austral summer climatologies
(October-March) using [Process CMIP6
products.R](scripts/Process%20CMIP6%20products.R). These files are
stored in the CMIP6 folder using [Git Large File
Storage](https://git-lfs.com).

## Contemporary covariates

Contemporary remotely sensed environmental data for temperature,
salinity and mixed layer depth are available from the [World Ocean Atlas
2018](https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/)
database. Data for sea surface height are available via
[AVISO](https://www.aviso.altimetry.fr/index.php?id=1526) and data for
sea ice concentration are available via
[NSIDC](https://nsidc.org/data/g02135/versions/3).

For each we calculated climatological means over the Austral summer
(October-March) for the longest time period available: for WOA
temperature, salinity and mixed layer depth this was 1981 to 2010; for
AVISO sea surface height this was 1993 to 2021; for NSIDC sea ice
concentration this was 1979 to 2020.

These data were processed using the [Process contemporary
covariates](/scripts/Process%20contemporary%20covariates.R) script and
outputed climatologies on 0.25 x 0.25 degree spatial resolution are
stored in the contemporary covariates folder.
