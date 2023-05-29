#########################################################
### Species-specific covariate selection with weights ###
#########################################################

# revise the species-specific covariate selection code to include weights
# https://www.tidyverse.org/blog/2022/05/case-weights/#getting-feedback
# write custom `manual_gam_selection` function
# uncertain that tidymodels is applying the weights correctly across folds...
# https://github.com/tidymodels/hardhat/issues/240

# Find the optimal number of environmental covariates to describe the habitat of each species
# Create candidate covariate sets for all possible combinations of covariates between 3 and 8
# only consider covariate combinations that aren't colinear

# Fit a GAM to each using 10-fold cross validation on the spatialsample folds
# calculate the weights on the fly dependent on the ratio of presence to absence points in each fold

######################################
### Libraries and custom functions ###
######################################

library(tidymodels)
library(sf)
library(spatialsample)
library(mgcv)

source("R/all_combs.R")
source("R/check_cor.R")
source("R/filter_vars.R")
source("R/tidy_covariate_selection.R")

############
### Data ###
############

# load original data
dat <- readRDS("data/presence_absence_data_10k_with_covariates_2023-05-23.rds")

# load pre-folded data
# NB this is an sf object outputed by spatialsample
folds <- readRDS("data/folds_weights.rds")

#######################################
### Set up candidate covariate sets ###
#######################################

# use a GAM to choose which is the best set of environmental covariates
form <- PresAbs ~ s(sst, bs = "ts", k = 5) + s(sst_grad, bs = "ts", k = 5) + s(sal, bs = "ts", k = 5) + s(ssh, bs = "ts", k = 5) + s(ssh_grad, bs = "ts", k = 5) + s(mld, bs = "ts", k = 5) + s(bat, bs = "ts", k = 5) + s(sic, bs = "ts", k = 5)
new_f <- all_combs(form, min.terms = 3) # calculate all possible combinations between 3 and 8 (with no interactions)

# initialise empty list to store results in
tidy_results <- list()

################################
###  1. Electrona antarctica ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[1]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
rev_f <- filter_vars(rev_f, "sal,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[1]] <- manual_gam_selection(rev_f, folds[[1]])

################################
###  2. Electrona carlsbergi ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[2]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[2]] <- manual_gam_selection(rev_f, folds[[2]])


tidy_results[[2]] |> arrange(desc(mean))

rev_f[[79]]

################################
###  3. Gymnoscopelus bolini ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[3]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[3]] <- manual_gam_selection(rev_f, folds[[3]])

#################################
###  4. Gymnoscopelus braueri ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[4]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[4]] <- manual_gam_selection(rev_f, folds[[4]])

#################################
###  5. Gymnoscopelus fraseri ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[5]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[5]] <- manual_gam_selection(rev_f, folds[[5]])

##################################
###  6. Gymnoscopelus nicholsi ###
##################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[6]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[6]] <- manual_gam_selection(rev_f, folds[[6]])

#######################################
###  7. Gymnoscopelus opisthopterus ###
#######################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[7]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[7]] <- manual_gam_selection(rev_f, folds[[7]])

####################################
###  8. Krefftichthys anderssoni ###
####################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[8]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[8]] <- manual_gam_selection(rev_f, folds[[8]])

#################################
###  9. Protomyctophum bolini ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[9]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[9]] <- manual_gam_selection(rev_f, folds[[9]])

###################################
### 10. Protomyctophum tenisoni ###
###################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[10]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[10]] <- manual_gam_selection(rev_f, folds[[10]])

#############################
### 11. Euphausia superba ###
#############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[11]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[11]] <- manual_gam_selection(rev_f, folds[[11]])

####################
### 12. Salpidae ###
####################

# check for correlations greater than .7
dat |> filter(species == unique(species)[12]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[12]] <- manual_gam_selection(rev_f, folds[[12]])

#####################################
### 13. Alluroteuthis antarcticus ###
#####################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[13]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[13]] <- manual_gam_selection(rev_f, folds[[13]])

##################################
### 14. Bathyteuthis abyssicol ###        
##################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[14]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[14]] <- manual_gam_selection(rev_f, folds[[14]])

###########################################
### 15. Slosarczykovia circumantarctica ###
###########################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[15]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[15]] <- manual_gam_selection(rev_f, folds[[15]])

#################################
### 16. Galiteuthis glacialis ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[16]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[16]] <- manual_gam_selection(rev_f, folds[[16]])

###############################
### 17. Gonatus antarcticus ###
###############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[17]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[17]] <- manual_gam_selection(rev_f, folds[[17]])

###################################
### 18. Histioteuthis atlantica ###
###################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[18]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[18]] <- manual_gam_selection(rev_f, folds[[18]])

###################################
### 19. Histioteuthis eltaninae ###
###################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[19]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[19]] <- manual_gam_selection(rev_f, folds[[19]])

################################
### 20. Kondakovia longimana ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[20]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[20]] <- manual_gam_selection(rev_f, folds[[20]])

#######################
### 21. Loligo gahi ###
#######################

#tidy_results[[21]] <- NULL

#############################
### 22. Martialia hyadesi ###
#############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[22]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[22]] <- manual_gam_selection(rev_f, folds[[22]])

######################################
### 23. Mesonychoteuthis hamiltoni ###
######################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[23]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[23]] <- manual_gam_selection(rev_f, folds[[23]])

####################################
### 24. Psychroteuthis glacialis ###
####################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[24]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[24]] <- manual_gam_selection(rev_f, folds[[24]])

################################
### 25. Todarodes filippovae ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[25]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[25]] <- manual_gam_selection(rev_f, folds[[25]])

###############################
### 26. Moroteuthis robsoni ###
###############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[26]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[26]] <- manual_gam_selection(rev_f, folds[[26]])

##############################
### 27. Moroteuthis ingens ###
##############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[27]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
rev_f <- filter_vars(rev_f, "sal,", "ssh,")  
length(rev_f)

# run manual selection avoiding using `fit_resamples`
tidy_results[[27]] <- manual_gam_selection(rev_f, folds[[27]])

####################
### Save results ###
####################

saveRDS(tidy_results, "data/tidy_results_weights.rds")

