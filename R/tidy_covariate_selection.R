#############################################
### Functions for GAM covariate selection ###
#############################################

# 2023-04-26

# two options for species-specific covariate selection
# the first uses weights implemented through hardhat::importance_weights
# theses are then passed to mgcv through fit_resamples
# there are potential issues with the AUC values returned by this...

# the second breaks apart the folds object and fits a GAM to each
# the weights are calculated independently on each fold and passed to mgcv
# I believe these AUC values over the first set
# there are small difference between the two -
# but these lead to different covariate selection

tidy_gam_selection <- function(formula_list, dat, folds){
  
  # set up receiving list
  results_summary <- list()
  
  # set up progress bar
  pb1 <- txtProgressBar(min = 1, max = length(formula_list), style = 3)
  
  # generate formula list copy with covariates and no smooth info
  f_no_s <- lapply(seq_along(formula_list), function(j) {
    reformulate(paste(all.vars(formula_list[[j]])[-1], collapse = " + "), response = all.vars(formula_list[[j]])[1])
  })
  
  # define gam model
  gam_spec <- 
    gen_additive_mod() |> 
    set_engine("mgcv", method = "REML") |> 
    set_mode("classification")
  
  # loop through formulas, fit GAMs and save AUC
  for (i in 1:length(formula_list)){
    # Set up case weights as a recipe step
    weights_recipe <- recipes::recipe(f_no_s[[i]], 
                                      data = sf::st_drop_geometry(dat)) |> 
      recipes::step_mutate(cwts = hardhat::importance_weights(ifelse(PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0))),
                           role = "case_weights")     # Need to set the "case_weights" role explicitly:
    
    setTxtProgressBar(pb1, i) # update progress bar
    
    gam_wf <- 
      workflow(preprocessor = weights_recipe) |> 
      add_model(gam_spec,
                formula = formula_list[[i]]) |> 
      add_case_weights(cwts) |>
      fit_resamples(folds)
    
    results_summary[[i]] <- gam_wf |> collect_metrics() |> mutate(model = i)
  }
  results_summary <- results_summary |> bind_rows() # convert from list to tibble
  return(results_summary)
}



manual_gam_selection <- function(formula_list, folds){
  
  # set up receiving list
  results_summary <- list()
  
  # set up progress bar
  pb1 <- txtProgressBar(min = 1, max = length(formula_list), style = 3)
  
  for(i in 1:length(formula_list)){
    
    setTxtProgressBar(pb1, i) # update progress bar
    
    # now fit to folds...
    auc_values <- vector("numeric", length = 10)
    
    for(j in 1:10) {
      check_gam <- folds$splits[[j]] |> 
        analysis() |> 
        mutate(cwts = ifelse(PresAbs == 1, 1, sum(PresAbs == 1) / sum(PresAbs == 0)))
      
      untidy_gam <- gam(formula_list[[i]],
                        data = check_gam,
                        family = "binomial",
                        weights = cwts)
      
      auc <- yardstick::roc_auc_vec(
        (folds$splits[[j]] |> assessment())$PresAbs,
        predict(
          untidy_gam, 
          folds$splits[[j]] |> assessment(), 
          type = "response"),
        event_level = "second")
      
      auc_values[j] <- auc
    }
    
    results_summary[[i]] <- tibble(model = i,
                                   terms = paste0(all.vars(formula_list[[i]])[-1], collapse = " + "),
                                   metric = "roc_auc",
                                   mean = mean(auc_values),
                                   n = 10,
                                   std_err = sd(auc_values)/sqrt(10))
  }
  results_summary <- results_summary |> bind_rows() # convert from list to tibble
  return(results_summary)
}




