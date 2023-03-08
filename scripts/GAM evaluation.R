#################################################
### Evaluate GAM models of habitat preference ###
#################################################

# 2023-02-13
# Amended from the GAM evaluation script written at Tor 2022-12-24
# Contemporary covariates are now processed and appended in a separate script.

# Load in data frame containing combined species presence and absence data with covariates included
# Run GAM models and test validity using the spatial k-folds cross validation procedure

# load libraries
require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(raster)
require(mgcv)
require(gratia)
require(blockCV)

# load spatial cross-validation functions
source("R/spatial_kfolds.r")
source("R/kfoldCV.r")
source("R/rank_models.r")
source("R/check_cor.R")
source("R/filter_vars.R")

# load data
dat <- readRDS("data/presence_absence_data_with_covariates.rds")
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# think carefully about what's next....

# how to run models - are GAMs best
# what spatial CV actually does
# go back through this....
# remember the Titley paper...

# go with a variation on the Titley paper
# for each species check the total number of combinations of covariates
# remove those that are colinear
# then fit GAMs for each of the covariate combinations
# use spatial CV algorithm to rank the models - select the best one for prediction?



myct <- dat %>% filter(group == "myctophids")
#myct <- myct %>% filter(bat > -7000)
#myct <- myct %>% filter(sal > 32.5)

# define maximum formula
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sst_grad,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(ssh_grad,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5) + s(sic,  bs = "cr", k = 5)
# calculate all possible combinations (with no interactions)
new_f <- all_combs(f)


############################
### Electrona antarctica ###
############################

E_ant <- myct %>% filter(species == "Electrona antarctica")
check_cor(E_ant[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
rev_f <- filter_vars(rev_f, "sal,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = E_ant, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AIC, AUC, TSS)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)



############################
### Electrona carlsbergi ###
############################

E_car <- myct %>% filter(species == "Electrona carlsbergi")
check_cor(E_car[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
rev_f <- filter_vars(rev_f, "sal,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = E_car, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

############################
### Gymnoscopelus bolini ###
############################

G_bol <- myct %>% filter(species == "Gymnoscopelus bolini")
check_cor(G_bol[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = G_bol, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

#############################
### Gymnoscopelus braueri ###
#############################

G_bra <- myct %>% filter(species == "Gymnoscopelus braueri")
check_cor(G_bra[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = G_bra, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

#############################
### Gymnoscopelus fraseri ###
#############################

G_fra <- myct %>% filter(species == "Gymnoscopelus fraseri")
check_cor(G_fra[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = G_fra, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

##############################
### Gymnoscopelus nicholsi ###
##############################

G_nic <- myct %>% filter(species == "Gymnoscopelus nicholsi")
check_cor(G_nic[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = G_nic, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

###################################
### Gymnoscopelus opisthopterus ###
###################################

G_opi <- myct %>% filter(species == "Gymnoscopelus opisthopterus")
check_cor(G_opi[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  
rev_f <- filter_vars(rev_f, "sal,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = G_opi, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

################################
### Krefftichthys anderssoni ###
################################

K_and <- myct %>% filter(species == "Krefftichthys anderssoni")
check_cor(K_and[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = K_and, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

#############################
### Protomyctophum bolini ###
#############################

P_bol <- myct %>% filter(species == "Protomyctophum bolini")
check_cor(P_bol[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "sal,")  
rev_f <- filter_vars(rev_f, "sst,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = P_bol, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)

###############################
### Protomyctophum tenisoni ###
###############################

### Fails - only 143 observations ###

P_ten <- myct %>% filter(species == "Protomyctophum tenisoni")
check_cor(P_ten[,9:16], threshold = .7) # check for correlations greater than .7

# filter formula list to remove pairs of covariates highlighted as being colinear
rev_f <- filter_vars(new_f, "sst,", "ssh,")  

# fit GAM for each of the remaining formulas with k = 6 k-folds spatial cross validation
ptm <- proc.time()
output <- rank_combs(data = P_ten, f = rev_f, k = 6, type = "blockCV")
proc.time() - ptm

# rank models based on AUC and TSS
trial <- output %>% dplyr::select(model, formula, AUC, TSS, AIC)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))
trial %>% head(10)



###
###
###





""        
[11] "Euphausia superba"               "Salpidae"                        


"Abraliopsis gilchristi"          "Alluroteuthis antarcticus"       "Bathyteuthis abyssicola"        
[16] "Batoteuthis skolops"             "Slosarczykovia circumantarctica" "Galiteuthis glacialis"           "Galiteuthis suhmi"               "Gonatus antarcticus"            
[21] "Grimpoteuthis megaptera"         "Histioteuthis atlantica"         "Histioteuthis eltaninae"         "Histioteuthis miranda"           "Illex argentinus"               
[26] "Kondakovia longimana"            "Loligo gahi"                     "Lycoteuthis lorigera"            "Martialia hyadesi"               "Mastigoteuthis psychrophila"    
[31] "Mesonychoteuthis hamiltoni"      "Parateuthis tunicata"            "Psychroteuthis glacialis"        "Semirossia patagonica"           "Taningia danae"                 
[36] "Teuthowenia pellucida"           "Todarodes filippovae"            "Moroteuthis robsoni"             "Moroteuthis ingens"              "Chiroteuthis veranyi"           
[41] "Moroteuthis knipovitchi"  









# so come up with a way to generate all combinations and then filter them to exclude sst and ssh combos

# or eddie kenetic energy
# or geostrophic current anomaly

# quick demonstration of analysis approach
myct <- dat %>% filter(group == "myctophids")
#myct <- myct %>% filter(bat > -7000)
#myct <- myct %>% filter(sal > 32.5)

# fit maximum GAM to each species
# then do k-folds cross validation to select covariates
E_ant <- myct %>% filter(species == "Electrona antarctica")

# first run full model
m <- gam(new_f[[1]], E_ant, family = binomial, method = "REML", select = T)
foo <- kfoldCV(E_ant, new_f[[1]], k)
  
# create output tibble
  out <- tibble_row(model = 1,
                    formula = paste(substr(attr(terms(new_f[[1]]), "term.labels"), 3, 5), collapse = " + "),
                    AIC = AIC(m),
                    AUC = list(foo[[1]]$AUC),
                    TSS = list(foo[[1]]$TSS),
                    spatial_cv = list(foo[[2]]))
  
  # then loop through combination of terms
  for(i in 2:length(new_f)){
    m <- gam(new_f[[i]], data, family = binomial, method = "REML", select = T)
    foo <- kfoldCV(data, new_f[[i]], k)
    
    # append results to output dataframe
    out <- out %>%
      bind_rows(
        tibble_row(model = i+1,
                   formula = paste(substr(attr(terms(new_f[[i]]), "term.labels"), 3, 5), collapse = " + "),
                   AIC = AIC(m),
                   AUC = list(foo[[1]]$AUC),
                   TSS = list(foo[[1]]$TSS),
                   spatial_cv = list(foo[[2]])))
  }


  tri <- out %>% dplyr::select(model, formula, AIC, AUC, TSS)
  tri <- tri %>%
    unnest(c(AUC, TSS)) %>%
    group_by(model) %>%
    mutate(AUC = mean(AUC),
           TSS = mean(TSS)) %>%
    slice(1) %>%
    ungroup() %>% 
    arrange(desc(AUC), desc(TSS))
  
  









                
for (i in 1:nrow(test)){
  print(cor(dat[c(test$V1[i])], dat[c(test$V2[i])]))
}
  
  
cor(dat[test[[3]]])

combn(variables, 8, simplify = F)
combn(variables, 7, simplify = F)
combn(variables, 6, simplify = F)
combn(variables, 5, simplify = F)
combn(variables, 4, simplify = F)
combn(variables, 3, simplify = F)
combn(variables, 2, simplify = F)
combn(variables, 1, simplify = F)



# quick demonstration of analysis approach
myct <- dat %>% filter(group == "myctophids")
#myct <- myct %>% filter(bat > -7000)
#myct <- myct %>% filter(sal > 32.5)

# fit maximum GAM to each species
# then do k-folds cross validation to select covariates

ptm <- proc.time()
E_ant <- myct %>% filter(species == "Electrona antarctica")
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5) + s(sic,  bs = "cr", k = 5)
output <- rank_combs(E_ant, f, k = 6)
proc.time() - ptm

trial <- output %>% dplyr::select(model, formula, AIC, AUC, TSS)
trial <- trial %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))



# want to minimise AIC, maximise the AUC and TSS
trial4 <- output %>% dplyr::select(model, formula, AIC, AUC, TSS)
trial4 <- trial4 %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))


# generate plot for visualising the spatial cv
# define the bounding box
test <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = -20)) %>% 
  sf::st_as_sfc() %>%
  st_set_crs(4326)
# split the bounding box into k segments
test <- test %>% st_make_grid(n = c(k, 1))

# project the box into polar projection
prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
test <- test %>% st_transform(prj)
test <- test %>% st_as_sf

ggplot() + 
  labs(title = "Coefficient of Variation",
       subtitle = "For each spatial k-fold") +
  theme_bw() +
  geom_sf(aes(), fill = NA, data = test) +
  geom_point(aes(x = x, y = y, colour = cv),
             data = output %>% unnest(spatial_cv)) +
  coord_sf(xlim = c(min(output %>% unnest(spatial_cv) %>% pull(x)),
                    max(output %>% unnest(spatial_cv) %>% pull(x))),
           ylim = c(min(output %>% unnest(spatial_cv) %>% pull(y)),
                    max(output %>% unnest(spatial_cv) %>% pull(y))),
           crs = prj) +
  scale_colour_viridis_c() +
  xlab(NULL) + ylab(NULL) +
  facet_wrap(~ model)

###
### Second species...
###

E_car <- myct %>% filter(species == "Electrona carlsbergi")
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
output <- rank_models(E_car, f, k = 6)

E_car_output <- rank_combs(E_car, f, k = 6)

# want to minimise AIC, maximise the AUC and TSS
trial_2 <- E_car_output %>% dplyr::select(model, formula, AIC, AUC, TSS)
trial_2 <- trial_2 %>%
  unnest(c(AUC, TSS)) %>%
  group_by(model) %>%
  mutate(AUC = mean(AUC),
         TSS = mean(TSS)) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(desc(AUC), desc(TSS))



# all terms
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop bat
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop sid
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop mld
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop ssh
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(sal,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop sal
f <- PresAbs ~ s(sst,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)

# drop sst
f <- PresAbs ~ s(sal,  bs = "cr", k = 5) + s(ssh,  bs = "cr", k = 5) + s(mld,  bs = "cr", k = 5) + s(sid,  bs = "cr", k = 5) + s(bat,  bs = "cr", k = 5)
foo <- kfoldCV(E_car, f, k = 6)
m <- gam(f, data = E_car, family = binomial, method = "REML", select = T)
AIC(m)
appraise(m)





