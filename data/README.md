<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data

Environmental data for the LGM, midHolocene and historical periods have
been download from the ESGF
[portal](https://esgf-node.llnl.gov/search/cmip6/) using the `rcmip6`
package in `R`.

These data have been downloaded on the native grid (gn) and then
processed using [Climate Data
Operators](https://code.mpimet.mpg.de/projects/cdo) bash scripts.

MPI-ESM1-2-LR model outputs contained more replicate years than
necessary. Processing was limited to the final 100 years of data using
bash script [here](/scripts/bash%20trim%20MPI.sh)

Processed data are available as monthly means for n years of data. From
these we generated Austral summer climatologies (October-March) for each
variable using a custom helper function
[process_CMIP6.R](/R/process_CMIP6.R)
