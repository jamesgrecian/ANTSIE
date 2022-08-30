########################################################################
### Download local copies of CMIP6 products using the rcmip6 package ###
########################################################################

# using the rcmip6 package after a conversation with @mdsumner
# https://twitter.com/mdsumner/status/1529025976025649153

# the rcmip6 package allows you to query the cmip6 portal
# returns list of urls that can be used to download the data and store locally

# need to think carefully about the environmental covariates used in the analysis
# some are easily available at monthly temporal resolution others not

# base initial exploration on the Freer 2019 Diversity and Distributions paper
# https://doi.org/10.1111/ddi.12934

#############################
### Variables of interest ###
#############################

# tos - temperature at surface
# sos - salinity at surface
# zos - sea surface height
# siconc - sea ice area fraction (concentration)
# mlotst - mixed layer depth

####################
### Set up query ###
####################

#remotes::install_github("eliocamp/rcmip6", force = T)
require(rcmip6)
require(tidyverse)

# initially focus on Last Glacial Maximum
query <- list(
  type          = "Dataset",
  experiment_id = "lgm",
  variable_id   = c("tos", "sos", "zos", "siconc", "mlotst"),
  replica       = "false",
  latest        = "true",
  project       = "CMIP6",
  frequency     = "mon",                          
  table_id      = c("Omon", "SImon")
)

# check for availability on portal
results <- cmip_search(query)

# output summary table
df <- cmip_simplify(results) %>% tibble() %>% dplyr::select(source_id, experiment_id, member_id, variable_id, grid_label, nominal_resolution)
df <- df %>% pivot_wider(names_from = variable_id, values_from = variable_id, values_fill = list(variable_id = NA))
df %>% arrange(source_id)

# download files required
cmip_root_set("/Volumes/Nifty 128/ANTSIE data/") # specify root directory
options(timeout = 1200)            # downloading from server so specify how long to keep connection open (will fail if not defined)
files <- cmip_download(results)    # trigger download

# files will be in various scripts
# reprocess to regular 1x1 degree grid using cdo bash script

# ends