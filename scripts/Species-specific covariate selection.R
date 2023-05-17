############################################
### Species-specific covariate selection ###
############################################

# 2023-04-11
# updated 2023-05-20

# Find the optimal number of environmental covariates to describe the habitat of each species
# Create candidate covariate sets for all possible combinations of covariates between 3 and 8
# Fit a GAM to each using 10-fold cross validation on the spatialsample folds
# Calculate AUC, AIC and TSS for each fold
# Output a tidy dataframe that can be ranked by chosen criteria#

# load libraries
library(tidyverse)
library(sf)
library(tidymodels)
library(mgcv)
library(purrr)
library(furrr)
library(spatialsample)
library(future)

# use future package to parallelise
plan(multisession) # use all 8 cores

# load spatial cross-validation functions
source("R/check_cor.R")
source("R/filter_vars.R")
source("R/tidy_models.R")
source("R/tidy_helpers.R")
source("R/all_combs.R")

# load original data
dat <- readRDS("data/presence_absence_data_10k_with_covariates_2023-04-05.rds")

# only consider cephalopod species included in the Xavier paper
drop_id <- c("Illex argentinus", "Lycoteuthis lorigera", "Teuthowenia pellucida",
             "Mastigoteuthis psychrophila", "Abraliopsis gilchristi", "Semirossia patagonica",
             "Batoteuthis skolops", "Chiroteuthis veranyi", "Galiteuthis suhmi",
             "Grimpoteuthis megaptera", "Histioteuthis miranda", "Moroteuthis knipovitchi",
             "Parateuthis tunicata", "Taningia danae")
dat <- dat |> filter(!species %in% drop_id)

# load pre-folded data
# NB this is an sf object outputed by spatialsample
folds <- readRDS("data/folds.rds")

# define maximum formula
form <- PresAbs ~ s(sst, bs = "cr", k = 5) + s(sst_grad, bs = "cr", k = 5) + s(sal, bs = "cr", k = 5) + s(ssh, bs = "cr", k = 5) + s(ssh_grad, bs = "cr", k = 5) + s(mld, bs = "cr", k = 5) + s(bat, bs = "cr", k = 5) + s(sic, bs = "cr", k = 5)
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

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[1]] <- tidy_models(spatial_blocks = folds[[1]], formula_list = rev_f)

################################
###  2. Electrona carlsbergi ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[2]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[2]] <- tidy_models(spatial_blocks = folds[[2]], formula_list = rev_f)

################################
###  3. Gymnoscopelus bolini ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[3]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[3]] <- tidy_models(spatial_blocks = folds[[3]], formula_list = rev_f)

#################################
###  4. Gymnoscopelus braueri ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[4]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[4]] <- tidy_models(spatial_blocks = folds[[4]], formula_list = rev_f)

#################################
###  5. Gymnoscopelus fraseri ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[5]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[5]] <- tidy_models(spatial_blocks = folds[[5]], formula_list = rev_f)

##################################
###  6. Gymnoscopelus nicholsi ###
##################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[6]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[6]] <- tidy_models(spatial_blocks = folds[[6]], formula_list = rev_f)

#######################################
###  7. Gymnoscopelus opisthopterus ###
#######################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[7]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[7]] <- tidy_models(spatial_blocks = folds[[7]], formula_list = rev_f)

####################################
###  8. Krefftichthys anderssoni ###
####################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[8]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[8]] <- tidy_models(spatial_blocks = folds[[8]], formula_list = rev_f)

#################################
###  9. Protomyctophum bolini ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[9]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[9]] <- tidy_models(spatial_blocks = folds[[9]], formula_list = rev_f)

###################################
### 10. Protomyctophum tenisoni ###
###################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[10]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[10]] <- tidy_models(spatial_blocks = folds[[10]], formula_list = rev_f)

#############################
### 11. Euphausia superba ###
#############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[11]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[11]] <- tidy_models(spatial_blocks = folds[[11]], formula_list = rev_f)

####################
### 12. Salpidae ###
####################

# check for correlations greater than .7
dat |> filter(species == unique(species)[12]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[12]] <- tidy_models(spatial_blocks = folds[[12]], formula_list = rev_f)

#####################################
### 13. Alluroteuthis antarcticus ###
#####################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[13]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[13]] <- tidy_models(spatial_blocks = folds[[13]], formula_list = rev_f)

##################################
### 14. Bathyteuthis abyssicol ###        
##################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[14]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[14]] <- tidy_models(spatial_blocks = folds[[14]], formula_list = rev_f)

###########################################
### 15. Slosarczykovia circumantarctica ###
###########################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[15]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[15]] <- tidy_models(spatial_blocks = folds[[15]], formula_list = rev_f)

#################################
### 16. Galiteuthis glacialis ###
#################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[16]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[16]] <- tidy_models(spatial_blocks = folds[[16]], formula_list = rev_f)

###############################
### 17. Gonatus antarcticus ###
###############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[17]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[17]] <- tidy_models(spatial_blocks = folds[[17]], formula_list = rev_f)

###################################
### 18. Histioteuthis atlantica ###
###################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[18]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[18]] <- tidy_models(spatial_blocks = folds[[18]], formula_list = rev_f)

###################################
### 19. Histioteuthis eltaninae ###
###################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[19]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[19]] <- tidy_models(spatial_blocks = folds[[19]], formula_list = rev_f)

################################
### 20. Kondakovia longimana ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[20]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[20]] <- tidy_models(spatial_blocks = folds[[20]], formula_list = rev_f)

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

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[22]] <- tidy_models(spatial_blocks = folds[[22]], formula_list = rev_f)

######################################
### 23. Mesonychoteuthis hamiltoni ###
######################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[23]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[23]] <- tidy_models(spatial_blocks = folds[[23]], formula_list = rev_f)

####################################
### 24. Psychroteuthis glacialis ###
####################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[24]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[24]] <- tidy_models(spatial_blocks = folds[[24]], formula_list = rev_f)

################################
### 25. Todarodes filippovae ###
################################

# check for correlations greater than .7
dat |> filter(species == unique(species)[25]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[25]] <- tidy_models(spatial_blocks = folds[[25]], formula_list = rev_f)

###############################
### 26. Moroteuthis robsoni ###
###############################

# check for correlations greater than .7
dat |> filter(species == unique(species)[26]) |> dplyr::select(9:16) |> check_cor(threshold = .7)

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
length(rev_f)

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[26]] <- tidy_models(spatial_blocks = folds[[26]], formula_list = rev_f)

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

# fit all models in rev_f to k folds and store results as nested tibble
tidy_results[[27]] <- tidy_models(spatial_blocks = folds[[27]], formula_list = rev_f)

saveRDS(tidy_results, "data/tidy_results.rds")

#####################
### Check results ###
#####################

# tidy results is now a list containing all the summary tibbles for each species and each model
# should be possible to do something like this code below
# but for each element of the list in one go...

#tidy_results[[1]] <- E_ant_results
#tidy_results[[2]] <- E_car_results
#tidy_results[[3]] <- G_bol_results
#tidy_results[[4]] <- G_bra_results
#tidy_results[[5]] <- G_fra_results

tidy_summary <- tidy_results |> 
  bind_rows(.id = "id") |> # convert from list to tibble but add an id column for each species
  unnest(cols = results) |> # pull out the results
  group_by(id, model) |> # group by species and model
  summarise(formula = unique(formula), # summarise to produce table with formula
            AUC_mean = mean(AUC), # mean AUC, TSS (and AIC for comparison
            TSS_mean = mean(TSS),
            AIC_mean = mean(AIC)) |>
  arrange(desc(AUC_mean), desc(TSS_mean)) |> # sort by AUC
  slice(1:round((.25*n()))) # slice to keep the top 25% of models based on the number of input formula

# check the number of rows differs per group
# i.e that n() is working as a way to take quantile of formulas
tidy_summary |> summarise(n = n()) |> print(n = Inf)

summary_formulas <- tidy_summary |> ungroup() |> pull(formula)

table(summary_formulas) |> as_tibble() %>% arrange(desc(n))




tidy_summary <- tidy_results |> 
  bind_rows(.id = "id") |> # convert from list to tibble but add an id column for each species
  unnest(cols = results) |> # pull out the results
  group_by(id, model) |> # group by species and model
  summarise(formula = unique(formula), # summarise to produce table with formula
            AUC_mean = mean(AUC), # mean AUC, TSS (and AIC for comparison
            TSS_mean = mean(TSS),
            AIC_mean = mean(AIC)) |>
  arrange(desc(AUC_mean), desc(TSS_mean)) |> # sort by AUC
  slice(1:round((.25*n()))) # slice to keep the top 25% of models based on the number of input formula


# arrange by AIC -> sst + sst_grad + ssh_grad + mld + bat + sic -> n = 25

# arrange by AUC -> sal + ssh + mld + bat + sic -> n = 20


# compare against previous analysis that included myctophids, krill, salps and three cephalopod species
foo <- tidy_results[c(1:12, 14, 16, 27)] 
# format summary output
foo <- foo |> 
  dplyr::bind_rows(.id = "id") |>
  unnest(cols = results) |>
  group_by(id, model) |>
  summarise(formula = unique(formula),
            AUC_mean = mean(AUC),
            TSS_mean = mean(TSS),
            AIC_mean = mean(AIC)) |>
  arrange(desc(AUC_mean), desc(TSS_mean)) |>
  slice(1:round((.25*length(rev_f))))

foo_formulas <- foo |> ungroup() |> pull(formula)
table(foo_formulas) |> as_tibble() %>% arrange(desc(n))

# arrange by AUC -> sst + sst_grad + mld + bat + sic -> n = 13




tidy_results_new[[2]] |> 
  unnest(cols = results) |>
  group_by(model) |>
  summarise(formula = unique(formula),
            AUC_mean = mean(AUC),
            TSS_mean = mean(TSS),
            AIC_mean = mean(AIC)) |>
  arrange(desc(AUC_mean), desc(TSS_mean)) |>
  slice(1:round((.25*n()))) # slice to keep the top 25% of models based on the number of input formula



### interesting that the AUC and TSS scores are the same whether predictions are response scale or model scale











# can you convert this to a function?
# something like...
# tidy_models <- function(spatial_blocks, formula_list){
# check that the inputs are right
# double check the names being passed between functions
# why did it go wrong before?

# set up receiving list
results_summary <- list()

# set up progress bar
pb1 <- txtProgressBar(min = 1, max = length(rev_f), style = 3)

for(i in 1:length(rev_f)){
  setTxtProgressBar(pb1, i) # update progress bar
  f <- rev_f[[i]]

# This is the key paragraph
# take the folded dataset
# fit a gam to each fold using the Analysis set
# pull out and store the coefficient info for later
# extract the predicted values using the Assessment set
test_models <- test |>
  mutate(model = future_map(splits, fit_gam, .options = furrr_options(packages = "sf", seed = T)),
         coef_info = future_map(model, tidy, .options = furrr_options(packages = "broom", seed = T)),
         results = future_map(model, glance, .options = furrr_options(packages = "broom", seed = T)),
         preds = future_map(splits, pred_gam, .options = furrr_options(seed = T)))

# Calculate Area Under the Curve
out <- test_models |>
  unnest(preds) |>
  group_by(id) |>
  roc_auc(factor(PresAbs), preds, event_level = "second") |>
  # roc_curve(factor(PresAbs), .pred, event_level = "second") |> autoplot() will plot the roc curve coloured by fold
  dplyr::select(id, roc_auc = .estimate)

# Calculate the True Skill Statistic
TSS <- test_models |>
  unnest(preds) |>
  group_by(id) |>
  group_map(~tss(.x)) |>
  unlist()

# Combine AUC and TSS
out <- out |>
  mutate(tss = TSS)

# pull other stats from nested glance object
stats <- test_models |> 
  unnest(results) |>
  dplyr::select(id, logLik, AIC, BIC)

# combine all stats
out <- out |>
  left_join(stats, by = join_by(id))

# Append to the test model object
test_models <- test_models |> 
  left_join(out, by = join_by(id))

results_summary[[i]] <- test_models |>
  dplyr::select(AUC = roc_auc, TSS = tss, logLik, AIC, BIC) |>
  mutate(formula = paste(all.vars(attr(terms(f), "variables"))[-1], collapse = " + "),
         model = i) |>
  relocate(model, formula) %>%
  nest(results = c(AUC, TSS, logLik, AIC, BIC))
}

results_summary <- results_summary |> bind_rows()
proc.time() - ptm


# format summary output
results_summary |> 
  unnest(cols = results) |>
  group_by(model) |>
  summarise(AUC_mean = mean(AUC),
            TSS_mean = mean(TSS),
            AIC_mean = mean(AIC)) |>
  arrange(desc(AUC_mean), desc(TSS_mean)) |>
  slice(1:round((.25*length(rev_f))))
  

training(foo)
foo <- folds[[1]]
f <- rev_f[[127]]

m <- gam(f, data = foo$splits[[1]] %>% analysis, family = binomial)
preds1 = predict(m, foo$splits[[1]] %>% assessment, type = "response")
preds2 = predict(m, foo$splits[[1]] %>% assessment)

preds2 = predict(m, foo$splits[[1]] %>% assessment, type = "terms") %>% data.frame()
#preds2$intercept <- as.numeric(m$coefficients[1])
preds2 <- rowSums(preds2)




# what is the correct predictions from a species distribution GAM...?
# type = response?
# some strange prediction without the intercept?


foo_models2 <- foo |>
  mutate(model = future_map(splits, fit_gam, .options = furrr_options(packages = "sf", seed = T)),
         coef_info = future_map(model, tidy, .options = furrr_options(packages = "broom", seed = T)),
         results = future_map(model, glance, .options = furrr_options(packages = "broom", seed = T)),
         preds = future_map(splits, pred_gam, .options = furrr_options(seed = T)))

# Calculate Area Under the Curve
foo_out2 <- foo_models2 |>
  unnest(preds) |>
  group_by(id) |>
  roc_auc(factor(PresAbs), preds, event_level = "second") |>
  # roc_curve(factor(PresAbs), .pred, event_level = "second") |> autoplot() will plot the roc curve coloured by fold
  dplyr::select(id, roc_auc = .estimate)




# Calculate the True Skill Statistic
TSS <- foo_models2 |>
  unnest(preds) |>
  group_by(id) |>
  group_map(~tss(.x)) |>
  unlist()

# Combine AUC and TSS
foo_out2 <- foo_out2 |>
  mutate(tss = TSS)

# pull other stats from nested glance object
stats <- spatial_blocks_models |> 
  unnest(results) |>
  dplyr::select(id, logLik, AIC, BIC)



# remember AUC sst c(0.709, 0.824, 0.755, 0.784, 0.687, 0.805, 0.784, 0.766, 0.795, 0.771)

x <- foo_models |> unnest(preds) |> filter(id == "Fold01")
x <- foo_models2 |> unnest(preds) |> filter(id == "Fold01")


# predicted values for known presences
# predicted values for known absences
p <- as.numeric(unlist(x[x$PresAbs == 1, 7]))
a <- as.numeric(unlist(x[x$PresAbs == 0, 7]))
e <- dismo::evaluate(p, a)
max(e@TPR + e@TNR-1)

foo |>
  unnest(preds) |>
  group_by(id) |>
  roc_auc(factor(PresAbs), preds, event_level = "second") |>
  # roc_curve(factor(PresAbs), .pred, event_level = "second") |> autoplot() will plot the roc curve coloured by fold
  dplyr::select(id, roc_auc = .estimate)

# Calculate the True Skill Statistic
TSS <- spatial_blocks_models |>
  unnest(preds) |>
  group_by(id) |>
  group_map(~tss(.x)) |>
  unlist()




# function to fit GAM and return predictions of model from testing dataset
# essentially doing fit_gam and pred_gam in one function
# quicker? but can't figure out how to return everything I need into nested tibble
# pred_gam <- function(splits){
#   # identify the assessment set
#   m <- gam(f, data = analysis(splits), family = "binomial")
#   holdout <- assessment(splits)
#   tibble::tibble(
#     PresAbs = holdout$PresAbs, 
#     preds = predict(m, holdout)
#   )
# }
#
# might be possible to improve speed by parallelising each iteration of formula to a gam fit
# rather than paralellising each k fold within the iteration of the loop
# would need the foreach package and doFuture package
# results_summary <- foreach(i = 1:2, .combine = c) %dofuture% {
# however lots of issues with how chunks are split and shared
# invested time not worth speed up efficiency...?

  
  