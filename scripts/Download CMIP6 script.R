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
# O2os - surface dissolved oxygen concentration
# intpp - Primary Organic Carbon Production by All Types of Phytoplankton
# siconc - sea ice area fraction (concentration)

# o2 - mean dissolved oxygen concentration through water column - then pull 200m
# mlotst - could use Ocean Mixed Layer Thickness Defined by Sigma T instead of trying to extract 200m contour for o2 and temperature?

####################
### Set up query ###
####################

# initially focus on MIROC-ES2L
# this has several products for lgm and holocene - many other modelling groups have not run this
# five variables from one model and three scenarios
query <- list(
  type          = "Dataset",
  source_id     = "MIROC-ES2L",
  grid_label    = "gr1",
  experiment_id = c("lgm", "holocene", "historical"),
  variable_id   = c("tos", "sos", "o2os", "intpp", "siconc"),
  replica       = "false",
  latest        = "true",
  project       = "CMIP6",
  frequency     = "mon",                          
  table_id      = c("Omon", "SImon"),
  member_id     = "r1i1p1f2"
)

# check for availability on portal
results <- cmip_search(query)
as.data.frame(results)[,1:14]

# download products
cmip_root_set("data/")             # specify root directory
options(timeout = 1200)            # downloading from server so specify how long to keep connection open (will fail if not defined)
files <- cmip_download(results)    # trigger download

# check products stored locally
foo <- cmip_available(root = cmip_root_get()) %>% as_tibble

# ends