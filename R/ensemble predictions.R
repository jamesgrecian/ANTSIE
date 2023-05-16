########################################
### Attempt at ensemble modelling... ###
########################################

# 2023-04-19

# re-run tidy modelling but checking that min terms is set to 3
# re-run tidy modelling to double check no change with predict type = reponse

# previous iterations :->
# arrange by AIC -> sst + sst_grad + ssh_grad + mld + bat + sic -> n = 25
# arrange by AUC -> sal + ssh + mld + bat + sic -> n = 20

# so need a routine that takes all the covariate data
# stacks the rasters
# generates spatial predictions for the whole surface

# think about how to do this with all the data in folds...
# weighting by AUC?
# how to propogate uncertainty?

# how to combine the different models???

#tidy_models fits all the models but doesn't return them
#object would become too large
# so have to run again?!

# also create seperate script that combines raster stack of environmental data
# and save that in data folder
# this can then be called when needed without having to re-run script...

# set up an RF model using randomForest
# mtry -> vary the number of parameters randomly tested and sampled at each branch between 1 and 3
# fix the minimum number of data points allowed in a node (prevent branch splitting further) to 5
# allow the number of trees to vary between 500 and 2500 in increments of 500?
# rank all these models by AUC and select the best
# both Reisinger and Tickly use randomForest package

# Here's a good example
# https://towardsdatascience.com/dials-tune-and-parsnip-tidymodels-way-to-create-and-tune-model-parameters-c97ba31d6173

require(tidymodels)
require(future)
require(sf)
require(spatialsample)
require(DALEXtra)

plan(multisession)

folds <- readRDS("data/folds.rds")

# set up random forest specification
# should this be a logistic regression?
rf_spec <- 
  rand_forest(mtry = tune(),         # vary number of predictors to randomly sample at each split
              min_n = 10,            # min number of data points in a node for node to be split further
              trees = tune(),        # vary number of trees to contain in the ensemble
              mode = "classification") |>  
  set_engine("randomForest",
             num.threads = 8)

# create workflow and specify model formula
rf_wf <- 
  workflow() |>
  add_model(rf_spec) |>
  add_formula(factor(PresAbs) ~ sal + ssh + mld + bat + sic)

# tune the randomForest model manually
# allow mtry and number of trees to vary...
manual_tune <- rf_wf |>
  tune_grid(resamples = folds[[1]],
            grid = expand.grid(
              mtry = c(2:4),                    # allow the number of covariates to randomly sample at each split to vary between 2 and n-1
              trees = seq(1000, 2500, by = 500) # allow the number of trees to vary between 1000 and 2500
            ))

# rank models by AUC
manual_tune |> 
  collect_metrics() |>
  filter(.metric == "roc_auc") |>
  arrange(desc(mean))

# plot to check
manual_tune |> autoplot(type = "marginals")

# select best model...
best_model <- manual_tune |> select_best(metric = "roc_auc")

# finalise workflow ready to for final fit
final_model <- rf_wf |> finalize_workflow(best_model)

# how to check model?
# how about a plot for each covariate
# each plot containing a line per block
# should also be able to generate spatial predictions for each...

# fit the final model to each of the 10 folds using fit_resamples
# specify control argument to extract model fit for each fold
#final_model_fits <- final_model |>
#  fit_resamples(folds[[1]],
#                control = control_resamples(extract = function (x) extract_fit_parsnip(x)))

# fit model to each fold - not sure fit_resamples with the extract argument is quite right
final_model_fits_list <- lapply(folds[[1]]$splits, FUN = function(x) fit(final_model, analysis(x)))

# partial dependence plots depend on using the DALEX::explain_tidymodels functions
# there may be a more straightforward method to this..?
# can make N small - in which case the line is the average for a small sample of the data
# make N big then less stochasticity BUT run time larger...
explain_wrapper <- function(input_model, input_data){
  pred <- input_data |> analysis() |> st_drop_geometry() |> dplyr::select(sal, ssh, mld, bat, sic)
  resp <- input_data |> analysis() |> st_drop_geometry() |> pull(PresAbs)
  pred_out <- explain_tidymodels(input_model, data = pred, y = resp) |> model_profile(N = 1000, type = "partial")
  pred_out <- pred_out[[2]] |> as_tibble()
  return(pred_out)
}

# map the explain_wrapper function to each model and data set
# out01 <- explain_wrapper(final_model_fits_list[[1]], final_model_fits$splits[[1]])
pdp_preds <- map2(final_model_fits_list, folds[[1]]$splits, explain_wrapper)
pdp_preds <- pdp_preds |> bind_rows(.id = "id")
names(pdp_preds) <- c("id", "var_name", "label", "x", "yhat", "ids")

p1 <- ggplot() + 
  theme_bw() +
  geom_line(aes(x = x, y = yhat, group = id), alpha = 0.5, 
            data = pdp_preds) +
  facet_wrap(~ `var_name`, scales = "free_x") +
  ylim(0, 1) +
  xlab(NULL) + ylab("Average prediction") +
  ggtitle("Partial Dependence Profiles")
p1
  

### generate predictions
### https://github.com/tidymodels/planning/issues/26

### now need to generate 10 spatial predictions
### could simply pass the covariate data to the predict function
### or should this be some other weird explainer thing?
explainer_generator <- function(input_model, input_data){
  pred <- input_data |> analysis() |> st_drop_geometry() |> dplyr::select(sal, ssh, mld, bat, sic)
  resp <- input_data |> analysis() |> st_drop_geometry() |> pull(PresAbs)
  explainer <- explain_tidymodels(input_model, data = pred, y = resp)
  return(explainer)
}
explainers <- map2(final_model_fits_list, folds[[1]]$splits, explainer_generator)


covs <- readRDS("data/covariate_stack.rds")
test <- rasterToPoints(covs) |> as_tibble()
test <- test |> dplyr::select(-x, -y, -sst, -sst_grad, -ssh_grad)
test <- test |> 
  dplyr::mutate(sal = replace_na(sal, mean(sal, na.rm = T)),
                ssh = replace_na(ssh, mean(ssh, na.rm = T)),
                mld = replace_na(mld, mean(mld, na.rm = T)),
                sic = replace_na(sic, mean(sic, na.rm = T)),
                bat = replace_na(bat, mean(bat, na.rm = T)))
test <- rep(list(test), 10)

test_preds <- map2(explainers, test, predict)


terra::predict(covs, explainers[[1]]) |> plot()



test_preds <- lapply(test_preds, tibble) |> bind_rows(.id = "id")
test_preds <- test_preds |> group_by(id)

trial <- rasterToPoints(covs) |> 
  as_tibble() |>
  dplyr::select(x, y, sal, ssh, mld, sic, bat)
trial <- trial[rep(1:nrow(trial),10),]
test_preds <- test_preds |> cbind(trial)
names(test_preds)[2] <- "preds"

test_preds$preds[is.na(test_preds$sal)] <- NA
test_preds$preds[is.na(test_preds$ssh)] <- NA
test_preds$preds[is.na(test_preds$mld)] <- NA
test_preds$preds[is.na(test_preds$sic)] <- NA
test_preds$preds[is.na(test_preds$bat)] <- NA


test_preds <- test_preds |> dplyr::select(id, x, y, preds)
test_preds <- test_preds |> ungroup()
preds_stack <- stack(
  test_preds |> filter(id == 1) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 2) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 3) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 4) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 5) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 6) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 7) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 8) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 9) |> dplyr::select(-id) |> rasterFromXYZ(),
  test_preds |> filter(id == 10) |> dplyr::select(-id) |> rasterFromXYZ()
)

plot(preds_stack)

out_auc <- final_model |> fit_resamples(folds[[1]])
auc_weights <- out_auc |> 
  unnest(.metrics) |>
  filter(.metric == "roc_auc") |>
  pull(.estimate)

preds_rf_raster <- weighted.mean(preds_stack, w = auc_weights)

plot(preds_rf_raster)

# ideally want to output the predictions for each fold as a dataframe
# along with the model? and AUC
# some nested tibble?
# do something with the out_auc tibble?











# when trying to use the DALEX package it gets very complicated
# each time you run the explain_tidymodels function the predictions vary slightly...



# fit the final model to each of the 10 folds using fit_resamples
# specify control argument to extract model fit for each fold
final_fit_models <- final_model |>
  fit_resamples(folds[[1]],
                control = control_resamples(extract = function (x) extract_fit_parsnip(x)))

# pull model fits from the final fit object
ffm_extracts <- final_fit_models |> unnest(.extracts) |> pull(.extracts)

explain_wrapper <- function(input_model, input_data){
  pred <- input_data |> analysis() |> st_drop_geometry() |> dplyr::select(sal, ssh, mld, bat, sic)
  resp <- input_data |> analysis() |> st_drop_geometry() |> pull(PresAbs)
  pred_out <- explain_tidymodels(input_model, data = pred, y = resp) |>
    model_profile(N = 100,
                  type = "accumulated")
  pred_out <- pred_out[[2]] |>
    as_tibble()
  return(pred_out)
}

out01 <- explain_wrapper(ffm_extracts[[1]], final_fit_models$splits[[1]])
ggplot() + 
  geom_line(aes(x = `_x_`, y = `_yhat_`, group = `_vname_`), data = out01, alpha = 0.5) +
  geom_line(aes(x = `_x_`, y = `_yhat_`, group = `_vname_`), data = preds_01, alpha = 0.5) +
  facet_wrap(~ `_vname_`, scales = "free_x") +
  ylim(0, 1) + xlab(NULL) + ylab(NULL) +
  theme_bw()

p2 <- ggplot() + 
  geom_line(aes(x = `_x_`, y = `_yhat_`, group = `_vname_`), data = preds_01, alpha = 0.5) +
  facet_wrap(~ `_vname_`, scales = "free_x") +
  ylim(0, 1) + xlab(NULL) + ylab(NULL) +
  theme_bw()

require(patchwork)
p1 + p2

out02 <- model_profile(out01)

preds_01$agr_profiles %>% dplyr::select(agr_profiles)

preds_01$agr_profiles |> as_tibble() |> mutate(Name = "Fold01")

foo01


# write wrapper function to output `explainer` predictions object
explain_wrapper <- function(splits, .extracts){
  pred <- analysis(splits) |> st_drop_geometry() |> dplyr::select(sal, ssh, mld, bat, sic)
  resp <- analysis(splits) |> st_drop_geometry() |> pull(PresAbs)
  input_model <- .extracts |> pull(.extracts)
  explain_tidymodels(input_model, data = pred, y = resp)
}

explainer_preds <- list()
for(i in 1:10){
  explainer_preds[[i]] <- explain_wrapper(final_fit_models[1,]$splits[[1]], final_fit_models[1,]$.extracts[[1]])
}



explainer_rf_01 <- explain_tidymodels(foo01,
                                      data = folds[[1]]$splits[[1]] |> analysis() |> dplyr::select(-PresAbs),
                                      y = folds[[1]]$splits[[1]] |> analysis() |> dplyr::select(PresAbs))
preds_01 <- model_profile(explainer_rf_01,
                          N = 100,
                          variables = c("sal", "ssh", "mld", "bat", "sic"),
                          type = "accumulated")

model_profile(explainer_preds[[1]],
                          N = 100,
                          variables = c("sal", "ssh", "mld", "bat", "sic"),
                          type = "accumulated")



final_fit_models[1,]$.extracts[[1]] |> pull(.extracts)

lala <- final_fit_models |> mutate(explain_preds = map(splits, .extracts, explain_wrapper))





### boosted regression trees
# boosted regression treees were generated using the gbm r package
# a similar cross-validation approach was used to parameterise the BRT model
# shrinkage parameter set to 0.001
# number of trees set to 5000
# tree complexity allowed to vary between 1 and 4
# tree complexity that minimused summed error across the testing data blocks was used to fit final set of 10 models


#gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 5), 
#n.trees = (1:10)*1000, 
#shrinkage = c(0.1, 0.5, 0.01),
#n.minobsinnode = 20)

# Ryan sets min obs in node at 20
# Tickley sets min obs in node at 10
# mtry (number of predictors at each node) 


# randomForests
# trees - start with 10 x number of covariates as the number of trees
# mtry - controls the split-variable randomization
# tree complexity/ node size
# node size = 1 for classification, node size = 5 for 
# min_n (minimum number of data points required to keep splitting) 

plan(multisession)
brt_spec <-
  boost_tree(mtry = tune(),               # number of covariates to randomly sample at each split
             trees = tune(),              # number of trees was set at 5000
             min_n = 10,                   # min number of data points in a node for it to be split further
             tree_depth = tune(),              # is this the level of interactions?
             learn_rate = .001,           # learning rate is also known as the shrinkage parameter
             sample_size = .8) |>         # proportion of data to sample at each step
  set_engine("xgboost",
             num.threads = 8,
             counts = T,
             importance = "impurity") |>  # gini impurity
  set_mode("classification")

brt_wf <- 
  workflow() |>
  add_model(brt_spec) |>
  add_formula(factor(PresAbs) ~ sal + ssh + mld + bat + sic) # formula selected based on GAM model selection

dials_random <- grid_random(
  mtry(c(2, 4)),         # keep mtry the same as before
  trees(c(1000, 10000)), # tickly 5000 Ryan 10000
  tree_depth(c(1, 4)),   # number of branches allowed in tree
  size = 10)

# manually tune the brt model
tune_brt <- brt_wf |>
  tune_grid(
    resamples = folds[[1]],
    grid = dials_random)

# rank all these models by AUC
tune_brt |> 
  collect_metrics() |>
  filter(.metric == "roc_auc") |>
  arrange(desc(mean))

tune_brt |> autoplot(type = "marginals")

best_model_brt <- tune_brt |> select_best(metric = "roc_auc")

# finalise workflow ready to for final fit
final_model_brt <- brt_wf |> finalize_workflow(best_model_brt)

# fit model to each fold - not sure fit_resamples with the extract argument is quite right
brt_final_model_fits_list <- lapply(folds[[1]]$splits, FUN = function(x) fit(final_model_brt, analysis(x)))
#[14:01:48] WARNING: src/learner.cc:767: 
#Parameters: { "importance", "num_threads" } are not used.



# map the explain_wrapper function to each model and data set
# out01 <- explain_wrapper(final_model_fits_list[[1]], final_model_fits$splits[[1]])
pdp_preds_brt <- map2(brt_final_model_fits_list, folds[[1]]$splits, explain_wrapper)
pdp_preds_brt <- pdp_preds_brt |> bind_rows(.id = "id")
names(pdp_preds_brt) <- c("id", "var_name", "label", "x", "yhat", "ids")

p2 <- ggplot() + 
  theme_bw() +
  geom_line(aes(x = x, y = yhat, group = id), alpha = 0.5, 
            data = pdp_preds_brt) +
  facet_wrap(~ `var_name`, scales = "free_x") +
  ylim(0, 1) +
  xlab(NULL) + ylab("Average prediction") +
  ggtitle("Partial Dependence Profiles")
p2


explainer_generator <- function(input_model, input_data){
  pred <- input_data |> analysis() |> st_drop_geometry() |> dplyr::select(sal, ssh, mld, bat, sic)
  resp <- input_data |> analysis() |> st_drop_geometry() |> pull(PresAbs)
  explainer <- explain_tidymodels(input_model, data = pred, y = resp)
  return(explainer)
}

brt_explainers <- map2(brt_final_model_fits_list, folds[[1]]$splits, explainer_generator)
brt_test_preds <- map2(brt_explainers, test, predict)

brt_test_preds <- lapply(brt_test_preds, tibble) |> bind_rows(.id = "id")
brt_test_preds <- brt_test_preds |> group_by(id)

brt_test_preds <- brt_test_preds |> cbind(trial)
names(brt_test_preds)[2] <- "preds"

brt_test_preds$preds[is.na(brt_test_preds$sal)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$ssh)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$mld)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$sic)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$bat)] <- NA


brt_test_preds <- brt_test_preds |> dplyr::select(id, x, y, preds)
brt_test_preds <- brt_test_preds |> ungroup()
brt_preds_stack <- stack(
  brt_test_preds |> filter(id == 1) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 2) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 3) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 4) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 5) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 6) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 7) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 8) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 9) |> dplyr::select(-id) |> rasterFromXYZ(),
  brt_test_preds |> filter(id == 10) |> dplyr::select(-id) |> rasterFromXYZ()
)

plot(brt_preds_stack)

out_auc_brt <- final_model_brt |> fit_resamples(folds[[1]])
auc_weights_brt <- out_auc_brt |> 
  unnest(.metrics) |>
  filter(.metric == "roc_auc") |>
  pull(.estimate)

preds_brt_raster <- weighted.mean(brt_preds_stack, w = auc_weights_brt)

plot(preds_rf_raster)
plot(preds_brt_raster)
plot(final_map)
final_map <- weighted.mean(stack(preds_rf_raster, preds_brt_raster),
                           w = c(mean(auc_weights), mean(auc_weights_brt)))

# load in shapefile for background mask, clip and project
world_shp <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = 4326) %>% sf::st_as_sfc()
world_shp <- world_shp %>% sf::st_crop(CP)
world_shp <- world_shp %>% sf::st_transform(projection(folds[[1]]$splits[[1]] |> assessment()))
world_shp <- world_shp %>% st_union()

p3 <- ggplot() + 
  theme_bw() +
  geom_raster(aes(x = x, y = y, fill = layer),
              data = rasterToPoints(preds_rf_raster) |> as_tibble()) +
  geom_sf(aes(), data = world_shp) +
  coord_sf(xlim = c(-5300000, 5300000), ylim = c(-5300000, 5300000)) +
#  scale_fill_viridis_c("", limits = c(0, 1)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Random Forest")
  
p4 <- ggplot() + 
  theme_bw() +
  geom_raster(aes(x = x, y = y, fill = layer),
              data = rasterToPoints(preds_brt_raster) |> as_tibble()) +
  geom_sf(aes(), data = world_shp) +
  coord_sf(xlim = c(-5300000, 5300000), ylim = c(-5300000, 5300000)) +
  scale_fill_viridis_c("", limits = c(0, 1)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Boosted Regression Tree")

p5 <- ggplot() +
  theme_bw() +
  geom_raster(aes(x = x, y = y, fill = layer),
              data = rasterToPoints(final_map) |> as_tibble()) +
  geom_sf(aes(), data = world_shp) +
  coord_sf(xlim = c(-5300000, 5300000), ylim = c(-5300000, 5300000)) +
  scale_fill_viridis_c("", limits = c(0, 1)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Average weighted by AUC")


p3 + p4 + p5 + plot_layout(guides = 'collect')

quartz(width = 11, height = 5)
print(p3 + p4 + p5 + plot_layout(guides = 'collect'))
quartz.save(file = "first ensemble attempt.jpg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()







# ideally want to output the predictions for each fold as a dataframe
# along with the model? and AUC
# some nested tibble?
# do something with the out_auc tibble?

# nicer plots...
# could try this https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/

source("R/discrete_gradient.R")
source("R/polar_mask.R")
polar_buffer <- polar_mask(radius_size = 5750000)

p3 <- ggplot() + 
  ggtitle("Random Forest") +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  xlab(NULL) + ylab(NULL) +
  geom_raster(aes(x = x, y = y, fill = layer),
              data = rasterToPoints(preds_rf_raster) |> as_tibble()) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(data = polar_buffer$mask, fill = "white", color = NA) +
  geom_sf(data = polar_buffer$mask, fill = NA, color = "grey40", size = 0.5/.pt) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE, crs = crs_polar, ndiscr = 1000) +
#  geom_text(aes(x = c(0, 6100000, 0, -6100000),
#                y = c(6000000, 0, -6000000, 0)),
#            label = c("0°", "90° E", "180°", "90° W")) +
  scale_fill_discrete_gradient("Habitat Suitability",
                               colours = viridis::viridis(10),
                               bins = 10,
                               limits = c(0, 1),
                               breaks = seq(0, 1, 0.2),
                               labels = seq(0, 1, 0.2),
                               guide = guide_colourbar(nbin = 500,
                                                       raster = T,
                                                       frame.colour = "grey40",
                                                       ticks.colour = "grey40",
                                                       frame.linewidth = .1,
                                                       barwidth = 20,
                                                       barheight = .5,
                                                       direction = "horizontal",
                                                       title.position = "top", #or "right"
                                                       title.theme = element_text(hjust = 0.5,
                                                                                  size = 10))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "vertical")

p4 <- ggplot() + 
  ggtitle("Boosted Regression Tree") +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  xlab(NULL) + ylab(NULL) +
  geom_raster(aes(x = x, y = y, fill = layer),
              data = rasterToPoints(preds_brt_raster) |> as_tibble()) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(data = polar_buffer$mask, fill = "white", color = NA) +
  geom_sf(data = polar_buffer$mask, fill = NA, color = "grey40", size = 0.5/.pt) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE, crs = crs_polar, ndiscr = 1000) +
  #  geom_text(aes(x = c(0, 6100000, 0, -6100000),
  #                y = c(6000000, 0, -6000000, 0)),
  #            label = c("0°", "90° E", "180°", "90° W")) +
  scale_fill_discrete_gradient("Habitat Suitability",
                               colours = viridis::viridis(10),
                               bins = 10,
                               limits = c(0, 1),
                               breaks = seq(0, 1, 0.2),
                               labels = seq(0, 1, 0.2),
                               guide = guide_colourbar(nbin = 500,
                                                       raster = T,
                                                       frame.colour = "grey40",
                                                       ticks.colour = "grey40",
                                                       frame.linewidth = .1,
                                                       barwidth = 20,
                                                       barheight = .5,
                                                       direction = "horizontal",
                                                       title.position = "top", #or "right"
                                                       title.theme = element_text(hjust = 0.5,
                                                                                  size = 10))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "vertical")

p5 <- ggplot() + 
  ggtitle("Average weighted by AUC") +
  theme_void(base_size = 10,
             base_family = "Helvetica Neue") +
  xlab(NULL) + ylab(NULL) +
  geom_raster(aes(x = x, y = y, fill = layer),
              data = rasterToPoints(final_map) |> as_tibble()) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(data = polar_buffer$mask, fill = "white", color = NA) +
  geom_sf(data = polar_buffer$mask, fill = NA, color = "grey40", size = 0.5/.pt) +
  coord_sf(xlim = c(-6400000, 6400000),
           ylim = c(-6400000, 6400000),
           expand = FALSE, crs = crs_polar, ndiscr = 1000) +
  #  geom_text(aes(x = c(0, 6100000, 0, -6100000),
  #                y = c(6000000, 0, -6000000, 0)),
  #            label = c("0°", "90° E", "180°", "90° W")) +
  scale_fill_discrete_gradient("Habitat Suitability",
                               colours = viridis::viridis(10),
                               bins = 10,
                               limits = c(0, 1),
                               breaks = seq(0, 1, 0.2),
                               labels = seq(0, 1, 0.2),
                               guide = guide_colourbar(nbin = 500,
                                                       raster = T,
                                                       frame.colour = "grey40",
                                                       ticks.colour = "grey40",
                                                       frame.linewidth = .1,
                                                       barwidth = 20,
                                                       barheight = .5,
                                                       direction = "horizontal",
                                                       title.position = "top", #or "right"
                                                       title.theme = element_text(hjust = 0.5,
                                                                                  size = 10))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.box = "vertical")

quartz(width = 11, height = 5)
print(p3 + p4 + p5 + plot_layout(guides = 'collect') & theme(legend.position = "bottom"))
quartz.save(file = "second ensemble attempt.jpg",
            type = "jpeg",
            dev = dev.cur(),
            dpi = 500)
dev.off()





#When fitting GLMs, we optimized the combinationof polynomial model terms to maximize model performance in terms of AUC
# for each species, as follows. GLMs were used to fit up to and including third-order polynomials for the
# five predictor variables, resulting in 243 candidatemodel formulations.
# Models were fitted to nine blocks of data, with theremaining block used as a testing dataset to evaluate AUC.
# This was thenrepeated for each of the 10 data blocks.
# The combination of polynomialterms that maximized AUC across the 10 model
# fittings was used to fit a finalset of 10 models.

source("R/glm_poly.R")
poly_f <- glm_poly(cov_terms = c("sal", "ssh", "mld", "sic", "bat"), poly_degree = 3)

glm_model <- logistic_reg() |> 
  set_engine("glm") |> 
  set_mode("classification")

# use workflow_set to fit glm to all the models in the poly_f list
# fit using resamples to give AUC across the folds
glm_workflows <- 
  workflow_set(preproc = poly_f, 
               models = list(logistic = glm_model)) |>
  workflow_map("fit_resamples", resamples = folds[[1]])

# rank these models by AUC
glm_workflows |> 
  collect_metrics() |> 
  filter(.metric == "roc_auc") |> 
  arrange(desc(mean))

# plot auc and accuaracy for all 243 models
glm_workflows |> autoplot(type = "marginals")

#pull out workflow
final_glm_workflow <- glm_workflows |> extract_workflow(id = "formula_238_logistic")

# fit model to each fold - not sure fit_resamples with the extract argument is quite right
glm_final_model_fits_list <- lapply(folds[[1]]$splits, FUN = function(x) fit(final_glm_workflow, analysis(x)))

# map the explain_wrapper function to each model and data set
# out01 <- explain_wrapper(final_model_fits_list[[1]], final_model_fits$splits[[1]])
pdp_preds_glm <- map2(glm_final_model_fits_list, folds[[1]]$splits, explain_wrapper)
pdp_preds_glm <- pdp_preds_glm |> bind_rows(.id = "id")
names(pdp_preds_glm) <- c("id", "var_name", "label", "x", "yhat", "ids")

ggplot() + 
  theme_bw() +
  geom_line(aes(x = x, y = yhat, group = id), alpha = 0.5, 
            data = pdp_preds_glm) +
  facet_wrap(~ `var_name`, scales = "free_x") +
  ylim(0, 1) +
  xlab(NULL) + ylab("Average prediction") +
  ggtitle("Partial Dependence Profiles")
p2


covs2 <- subset(covs, c(3, 4, 6, 7, 8))
glm_test_preds <- lapply(glm_explainers, function(x) terra::predict(covs2, x))

par(mfrow = c(2, 5))
plot(glm_test_preds[[1]])
plot(glm_test_preds[[2]])
plot(glm_test_preds[[3]])
plot(glm_test_preds[[4]])
plot(glm_test_preds[[5]])
plot(glm_test_preds[[6]])
plot(glm_test_preds[[7]])
plot(glm_test_preds[[8]])
plot(glm_test_preds[[9]])
plot(glm_test_preds[[10]])




terra::predict(covs2, glm_explainers[[1]])

|> plot()


glm_explainers <- map2(glm_final_model_fits_list, folds[[1]]$splits, explainer_generator)

glm_test_preds <- map2(glm_explainers, covs2, terra::predict)






brt_test_preds <- lapply(brt_test_preds, tibble) |> bind_rows(.id = "id")
brt_test_preds <- brt_test_preds |> group_by(id)

brt_test_preds <- brt_test_preds |> cbind(trial)
names(brt_test_preds)[2] <- "preds"

brt_test_preds$preds[is.na(brt_test_preds$sal)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$ssh)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$mld)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$sic)] <- NA
brt_test_preds$preds[is.na(brt_test_preds$bat)] <- NA





# Update your workflow and fit: 
final_glm <- 
  final_glm_workflow %>% 
  finalize_workflow(best_results) %>% 
  fit(data = Chicago)



# If there are tuning parameters, get the best results
best_results <- 
  chi_features_res %>% 
  extract_workflow_set_result(id = "plus_pca_lm") %>% 
  select_best(metric = "rmse")

# Update your workflow and fit: 
fitted_workflow_object <- 
  workflow_object %>% 
  finalize_workflow(best_results) %>% 
  fit(data = Chicago)

fitted_workflow_object




# finalise workflow ready to for final fit
final_model_glm <- glm_workflows |> finalize_workflow(best_model_glm)

best_results <- 
  chi_features_res %>% 
  extract_workflow_set_result(id = "plus_pca_lm") %>% 
  select_best(metric = "rmse")

# Update your workflow and fit: 
fitted_workflow_object <- 
  workflow_object %>% 
  finalize_workflow(best_results) %>% 
  fit(data = Chicago)





# manually tune the brt model
tune_brt <- brt_wf |>
  tune_grid(
    resamples = folds[[1]],
    grid = dials_random)

# rank all these models by AUC
tune_brt |> 
  collect_metrics() |>
  filter(.metric == "roc_auc") |>
  arrange(desc(mean))

tune_brt |> autoplot(type = "marginals")

best_model_brt <- tune_brt |> select_best(metric = "roc_auc")

# finalise workflow ready to for final fit
final_model_brt <- brt_wf |> finalize_workflow(best_model_brt)

# fit model to each fold - not sure fit_resamples with the extract argument is quite right
brt_final_model_fits_list <- lapply(folds[[1]]$splits, FUN = function(x) fit(final_model_brt, analysis(x)))
#[14:01:48] WARNING: src/learner.cc:767: 
#Parameters: { "importance", "num_threads" } are not used.








glm_tune |> pull(1)



glm_workflows |> unnest(info) |> pull(preproc)



###
### Old scripts...???
###


# manually tune the brt model
manual_tune_brt <- brt_wf |>
  tune_grid(resamples = folds[[1]],
            grid = expand.grid(
              mtry = c(2:4),                    # allow the number of predictors at each split to vary between 1 and 3
              trees = seq(2500, 10000, by = 2500), # allow the number of trees to vary between 1000 and 2500
              tree_depth = tree_depth(range = c(1, 4))),
            control = control_grid(save_pred = T))
#            control = control_grid(save_workflow = T)) # save workflow to pass to fit_best

manual_tune_brt %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:trees) %>%
  pivot_longer(mtry:trees,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")


# rank all these models by AUC
manual_tune_brt |> 
  collect_metrics() |>
  filter(.metric == "roc_auc") |>
  arrange(desc(mean))

manual_tune_brt |> autoplot(type = "marginals")


# how to actually plot the partial plots from this model to see what it's doing?


manual_tune_brt |> select_best(metric = "roc_auc")

    extract_fit_parsnip() %>% 
  vip(num_features = 5)

last_rf_mod <- 
  rand_forest(mtry = 8, min_n = 7, trees = 1000) %>% 
  set_engine("ranger", num.threads = cores, importance = "impurity") %>% 
  set_mode("classification")

# the last workflow
last_rf_workflow <- 
  rf_workflow %>% 
  update_model(last_rf_mod)

# the last fit
set.seed(345)
last_rf_fit <- 
  last_rf_workflow %>% 
  last_fit(splits)




best_model <- select_best(manual_tune_brt, metric = "roc_auc")
final_wf <- finalize_workflow(brt_wf, best_model)
final_fit <- fit_resamples(final_wf,
                           folds[[1]])





xgb_spec <- boost_tree(
  trees = 1000, 
  tree_depth = tune(),
  min_n = tune(), 
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), 
  mtry = tune(),         ## randomness
  learn_rate = tune(),                         ## step size
) %>% 
  set_engine("xgboost") %>% 
  set_mode("classification")


# go wild with tuning
wild_brt_spec <-
  boost_tree(mtry = tune(),              
             trees = tune(),          
             min_n = tune(),           
             tree_depth = tune(),       
             learn_rate = tune(),
             loss_reduction = tune(),                  
             sample_size = tune()) |>    
  set_engine("xgboost",
             num.threads = 8,
             counts = T,
             importance = "impurity") |> # gini impurity
  set_mode("classification")

wild_xgb_grid <- grid_latin_hypercube(
  mtry = mtry(range = c(1, 5)),
  trees = trees(range = c(1000, 10000)),
  min_n(),
  tree_depth(),
  learn_rate(), 
  loss_reduction(),
  sample_size = sample_prop(),
  size = 100
)

wild_brt_wf <- 
  workflow() |>
  add_model(wild_brt_spec) |>
  add_formula(factor(PresAbs) ~ sal + ssh + mld + bat + sic) # formula selected based on GAM model selection

# manually tune the brt model
wild_tune_brt <- wild_brt_wf |>
  tune_grid(
    resamples = folds[[1]],
    grid = wild_xgb_grid,
    control = control_grid(save_pred = T)) # save predictions for later plotting

wild_tune_brt |>
  collect_metrics() %>% print(n = Inf)


manual_tune_brt %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")


show_best(wild_tune_brt)






















# fit workflow to blocked samples
brt_fit_rs <- 
  brt_wf |> 
  tune::fit_resamples(folds[[1]])

# check summary
brt_fit_rs |> 
  collect_metrics()










impurity

### need to carefully think about parameters to pass/ test with rf and brt
### can use the `tune` parameter to decide on what might be best...
### look at literature -
### Ryan's paper PEIfuture
### Tickly paper and R scripts....

# ryan random forests

#rfGrid <-  data.frame(mtry = c(3, 4, 5),
#                      splitrule = "gini",
#                      min.node.size = c(1))



# randomForests
# trees - start with 10 x number of covariates as the number of trees
# mtry - controls the split-variable randomization
# mtry (number of predictors at each node) 
# tree complexity/ node size
# node size = 1 for classification, node size = 5 for 
# min_n (minimum number of data points required to keep splitting) 



rfw <- workflow() %>%
  add_model(rf_mod) %>%
  add_recipe(wt_rec) 

rf_fit <- fit(rfw,  d_train)

rf_fit %>%
  extract_fit_parsnip() %>%
  vip()





rf_reg_spec <- 
  rand_forest(trees = 200, min_n = 5) |> 
  set_engine("randomForest") |>
  set_mode("classification") %>%   # This model can be used for classification or regression, so set mode
  
  
  rf_reg_spec

rf_reg_fit <- rf_reg_spec %>% fit(PresAbs ~ ., data = gbm_dat_1)

rf_reg_fit

rf_mod <- 
  rand_forest(trees = 200) %>% 
  set_engine("randomForest") %>% 
  set_mode("classification")

rf_wf <- 
  workflow() %>%
  add_model(rf_mod) %>%
  add_formula(PresAbs ~ sal + ssh + mld + bat + sic)

rf_fit_rs <- 
  rf_wf %>% 

  
  
  
  
  
  
    tune::fit_resamples(folds[[1]])






















source("R/tidy_helpers.R")
source("R/tidy_models.R")


# GLM
# GAMS
# RF
# BRTs

# given a formula
# given a dataset
# given a prediction set of new data
# fit model to dataset
# generate predictions for new data
# do this for each fold

foo <- folds[[1]] #|> filter(id == "Fold01")
for (i in 1:10){
  m <- gam(f, data = foo$splits[[i]] %>% analysis, family = "binomial")
  r <- terra::predict(object = subset(covs, c(3, 4, 6, 8, 7)),
                      model = m,
                      fun = predict,
                      type = "response")
  par(ask = T)
  plot(r)
}

# need to tweak this to include AUC calculation
# then work out best way to store output

# then work out how the hell to run GLM, BRT, RF etc on the same data...


library(tidymodels)
library(spatialsample)

workflow() |> 
  add_model(linear_reg()) |> 
  add_formula(canopy_area_2019 ~ land_area * mean_temp) |> 
  fit_resamples(spatial_block_cv(spatialsample::boston_canopy)) |> 
  collect_metrics()

workflow() |> 
  add_model(gen_additive_mod("regression")) |> 
  add_formula(canopy_area_2019 ~ land_area * mean_temp) |> 
  tune::fit_resamples(spatial_block_cv(spatialsample::boston_canopy)) |> 
  collect_metrics()




gam_spec <-
  gen_additive_mod() |>
  set_engine("mgcv") |>
  set_mode("regression")

base_wf <- workflow() |> 
  add_model(gam_spec)

base_wf |> 
  fit(canopy_area_2019 ~ s(land_area) + s(mean_temp), data = spatialsample::boston_canopy)

  
)
model_spec <- linear_reg()

base_wf |>
  add_model(gam_spec) |>
  fit(spatialsample::boston_canopy)
  
  add_model(gen_additive_mod(mode = "regression",
                             engine = "mgcv")) |>
  fit(data = spatialsample::boston_canopy)
  



  add_model() |> 
  add_formula(canopy_area_2019 ~ land_area * mean_temp) |> 
  fit_resamples(vfold_cv(spatialsample::boston_canopy)) |> 
  collect_metrics()

  
gam_wf <- workflow() |>
  add_variables(outcomes = c(PresAbs), predictors = c(sal, ssh, mld, bat, sic)) |>
  add_model(spec = gam_spec,
            formula = factor(PresAbs) ~ s(sal) + s(ssh) + s(mld) + s(bat) + s(sic))

test_gam <- fit_resamples(gam_wf, folds[[1]])


gen_additive_mod(mode = "regression",
                 engine = "mgcv")|>
  fit(canopy_area_2019 ~ land_area * mean_temp, data = spatialsample::boston_canopy))





install.packages("randomForest")
require(randomForest)


rf_reg_spec <- 
  rand_forest(trees = 200, min_n = 5) |> 
  set_engine("randomForest") |>
  set_mode("classification") %>%   # This model can be used for classification or regression, so set mode
  
  
  rf_reg_spec

rf_reg_fit <- rf_reg_spec %>% fit(PresAbs ~ ., data = gbm_dat_1)

rf_reg_fit

rf_mod <- 
  rand_forest(trees = 200) %>% 
  set_engine("randomForest") %>% 
  set_mode("classification")

rf_wf <- 
  workflow() %>%
  add_model(rf_mod) %>%
  add_formula(factor(PresAbs) ~ sal + ssh + mld + bat + sic)

rf_fit_rs <- 
  rf_wf %>% 
  tune::fit_resamples(folds[[1]])


workflow() |> 
  add_model(linear_reg()) |> 
  add_formula(PresAbs ~ sal + ssh + mld + bat + sic) |> 
  fit_resamples(folds[[1]]) |> 
  collect_metrics()
#> # A tibble: 2 × 6
#>   .metric .estimator       mean     n    std_err .config             
#>   <chr>   <chr>           <dbl> <int>      <dbl> <chr>               
#> 1 rmse    standard   379263.       10 13950.     Preprocessor1_Model1
#> 2 rsq     standard        0.354    10     0.0237 Preprocessor1_Model1
#> 
#> 

foo <- folds[[1]]



library(tidymodels)
library(spatialsample)

spatial_block_cv(spatialsample::boston_canopy)

workflow() |> 
  add_model(linear_reg()) |> 
  add_formula(canopy_area_2019 ~ land_area * mean_temp) |> 
  fit_resamples(spatial_block_cv(spatialsample::boston_canopy)) |> 
  collect_metrics()

workflow() |> 
  add_model(linear_reg()) |> 
  add_formula(canopy_area_2019 ~ land_area * mean_temp) |> 
  fit_resamples(vfold_cv(spatialsample::boston_canopy)) |> 
  collect_metrics()



folds <- readRDS("/Users/home/ANTSIE/data/folds.rds")

library(tidymodels)
library(spatialsample)

gam_wf <- workflow() |>
  add_variables(outcomes = c(PresAbs), predictors = c(sal, ssh, mld, bat, sic)) |>
  add_model(spec = gam_spec,
            formula = factor(PresAbs) ~ s(sal) + s(ssh) + s(mld) + s(bat) + s(sic))
test_gam <- fit_resamples(gam_wf, folds[[1]])



folds[[2]] %>% filter(model == 148) %>% unnest(results)


f <- PresAbs ~ s(sal, bs = "cr", k = 5) + s(ssh, bs = "cr", k = 5) + s(mld, bs = "cr", k = 5) + s(bat, bs = "cr", k = 5) + s(sic, bs = "cr", k = 5)

raster::predict


m <- gam(f, data = analysis(splits), family = "binomial")
preds = predict(m, covs, type = "response"))






require(gbm)
mySummary  <- function(data, lev = NULL, model = NULL){
  a1 <- defaultSummary(data, lev, model)
  b1 <- twoClassSummary(data, lev, model)
  c1 <- prSummary(data, lev, model)
  out <- c(a1, b1, c1)
  out}

require(caret)
library(foreach)
library(parallel)
library(doParallel)

tc <- trainControl(method = "cv",
                   number = 10,
                   search = "grid",
                   classProbs = TRUE,
                   allowParallel = TRUE,
                   summaryFunction = mySummary,
                   sampling = "down")
gbmGrid <-  expand.grid(interaction.depth = 1, 
                        n.trees = (1:10)*1000, 
                        shrinkage = c(0.1, 0.5, 0.01),
                        n.minobsinnode = 20)

clust <- makeCluster(detectCores() - 1) # leave 1 core for OS
registerDoParallel(clust)

gbm_dat_1 <- foo$splits[[i]] %>% analysis %>% st_drop_geometry() %>% dplyr::select(PresAbs, sal, ssh, mld, bat, sic)
gbm_dat_1 <- gbm_dat_1 %>% mutate(PresAbs = as.factor(PresAbs))
MOD.gbm <- train(x = gbm_dat_1 %>% dplyr::select(sal, ssh, mld, bat, sic) %>% data.frame(),
                 y = gbm_dat_1 %>% pull(PresAbs),
                 method = "gbm",
                 metric = "ROC",
                 trControl = tc,
                 tuneGrid = gbmGrid)



gbm_m <- gbm(formula = sal + ssh + mld + bat + sic,
    distribution="binomial",
    data=dat.blockNo,
    n.trees=5000,
    n.minobsinnode = 10,
    interaction.depth=tc,
    shrinkage=0.001,
    bag.fraction=0.5,
    verbose=F,
    keep.data=F,
    train.fraction=tfrac)

GBM_eco <- function(data.model, outDir, plotPath, species, eval=TRUE,
                    tree.comp = c(1:3), learn.rate = c(0.01,0.001),
                    PA, blocks, clim.var, distribution="bernoulli")
  




m, 

tidy_predictions <- function(spatial_blocks, formula, raster_stack){
  
  # set up receiving list
  results_summary <- list()
  
  # set up progress bar
  pb1 <- txtProgressBar(min = 1, max = length(formula_list), style = 3)
  
  # loop through formulas in list
  for(i in 1:length(formula_list)){
    

    # function to fit gam to each fold using analysis data
    # return spatial predictions
    return_spatial_predictions <- function(splits){
      m <- gam(f, data = analysis(splits), family = "binomial")
      preds = predict(m, covs, type = "response"))
      )
    }
    
    setTxtProgressBar(pb1, i) # update progress bar
    f <- formula_list[[i]] # iterate through formulas
    
    # This is the key paragraph
    # take the folded dataset
    # fit a gam to each fold using the Analysis set
    # pull out and store the coefficient info for later
    # extract the predicted values using the Assessment set
    spatial_blocks_models <- spatial_blocks |>
      mutate(model = future_map(splits, fit_gam, .options = furrr_options(packages = "sf", seed = T)),
             coef_info = future_map(model, tidy, .options = furrr_options(packages = "broom", seed = T)),
             results = future_map(model, glance, .options = furrr_options(packages = "broom", seed = T)),
             preds = future_map(splits, pred_gam, .options = furrr_options(seed = T)))
    
    # Calculate Area Under the Curve
    out <- spatial_blocks_models |>
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
    
    
    
#    require(spatialsample)
#    require(sf)
#    require(parsnip)
#    require(workflows)
#    
#    gam_spec <- 
#      gen_additive_mod() |>
#      set_engine("mgcv", family = "binomial") |>
#      set_mode("regression")
#    
#    gam_wf <- workflow() |>
#      add_variables(outcomes = c(PresAbs), predictors = c(sal, ssh, mld, bat, sic)) |>
#      add_model(spec = gam_spec,
#                formula = f)
#    
#    test_gam <- fit_resamples(gam_wf, folds[[1]])
#    