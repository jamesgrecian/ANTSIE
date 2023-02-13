###################################################
### Process species presence data for modelling ###
###################################################

# updated 2023-02-13

# load libraries
require(tidyverse)
require(sf)
require(lubridate)

# load raw data
myctophids <- read_csv("data/species presence data/myctophids.csv")
krill <- read_csv("data/species presence data/krill.csv")
cephalopods <- read_csv("data/species presence data/cephalopods.csv")

# only consider austral spring/ summer... October to March
myctophids <- myctophids %>% filter(month %in% c(10, 11, 12, 1, 2, 3)) # Jen already filtered the data...
krill <- krill %>% mutate(DATE = dmy(DATE))
krill <- krill %>% filter(month(DATE) %in% c(10, 11, 12, 1, 2, 3)) # loose about ~500 records (~3.5%)
krill <- krill %>% filter(year(DATE) > 1970) # start at 1975 - looses around 3500 records
# date data not available for cephalopods

# format to combined data frame to pass to analysis
# group, species, date, lon, lat, 

myctophids <- myctophids %>% dplyr::select(NAME, eventDate, LON, LAT)

krill <- krill %>% dplyr::select(DATE, LONGITUDE, LATITUDE, STANDARDISED_KRILL_UNDER_1M2, NUMBER_OF_SALPS_UNDER_1M2)
names(krill)[4:5] <- c("Euphausia superba", "Salpidae")
krill <- krill %>% pivot_longer(cols = 4:5, names_to = "species", values_to = "count_under_1m3")
krill <- krill %>% dplyr::select(species, DATE, LONGITUDE, LATITUDE, count_under_1m3)
krill <- krill %>% drop_na(count_under_1m3) # drop cases where no krill/salps counted # 5.5k empty rows dropped

cephalopods <- cephalopods %>% dplyr::select(ScientificName_accepted, Longitude, Latitude)

# unify column names
myctophids$group <- "myctophids"
names(myctophids) <- c("species", "date", "lon", "lat", "group")
myctophids <- myctophids %>% dplyr::select("group", "species", "date", "lon", "lat")

krill$group <- krill$species
krill <- krill %>% dplyr::select("group", "species", "DATE", "LONGITUDE", "LATITUDE")
names(krill) <- c("group", "species", "date", "lon", "lat")

names(cephalopods) <- c("species", "lon", "lat")
cephalopods$group <- "cephalopods"
cephalopods$date <- NA
cephalopods <- cephalopods %>% dplyr::select("group", "species", "date", "lon", "lat")

# combine...
dat <- rbind(myctophids, krill, cephalopods)
saveRDS(dat, "data/combined_presence_data.rds")
