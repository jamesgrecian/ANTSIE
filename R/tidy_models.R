#########################################################################
### Function to fit list of GAM formulas to spatially blocked dataset ###
#########################################################################

# 2023-04-14

# This function aides the ranking of candidate covariates during model selection

# The output from the spatialsample package is an rsample object
# This function takes that object as input alongside a list of formula
# Function then iterates along all the formulas provided in the list (typically > 100)
# During each iteration, focal model will be fit to the K folds provided in the spatial_blocks object
# The predictions from each fold are outputed on the response scale
# These are then used to calculate AUC, TSS and AIC for each fold
# Function then outputs nested tibble containing:
# splits - from rsample that were inputted
# id - of each fold
# model - model object
# coef_inf - a tidy tibble of model coefficients
# results - a tibble of model performance metrics
# preds - a tibble containing inputted PresAbs and model predictions - to pass to AUC and TSS statistics        

tidy_models <- function(spatial_blocks, formula_list){
  
  # set up receiving list
  results_summary <- list()
  
  # set up progress bar
  pb1 <- txtProgressBar(min = 1, max = length(formula_list), style = 3)
  
  # loop through formulas in list
  for(i in 1:length(formula_list)){
    
    # function to fit gam to each fold
    # see https://www.tidymodels.org/learn/statistics/bootstrap/#bootstrapping-models
    fit_gam <- function(splits) {
      m <- gam(f, data = analysis(splits), family = "binomial")
    }
    
    # function to fit gam to each fold using analysis data
    # returns predictions using assessment data from same fold
    # see https://spatialsample.tidymodels.org/articles/spatialsample.html
    pred_gam <- function(splits){
      # identify the assessment set
      m <- gam(f, data = analysis(splits), family = "binomial")
      holdout <- assessment(splits)
      tibble::tibble(
        PresAbs = holdout$PresAbs, 
        preds = as.numeric(predict(m, holdout, type = "response"))
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
    
    # Combine AUC and TSS
    out <- out |>
      mutate(tss = TSS)
    
    # pull other stats from nested glance object
    stats <- spatial_blocks_models |> 
      unnest(results) |>
      dplyr::select(id, logLik, AIC, BIC)
    
    # combine all stats
    out <- out |>
      left_join(stats, by = join_by(id))
    
    # Append to the test model object
    spatial_blocks_models <- spatial_blocks_models |> 
      left_join(out, by = join_by(id))
    
    # format output for summary table
    # nested tibble containing
    # model id, formula, AUC, TSS, logLik, AIC and BIC
    results_summary[[i]] <- spatial_blocks_models |>
      dplyr::select(AUC = roc_auc, TSS = tss, logLik, AIC, BIC) |>
      mutate(formula = paste(all.vars(attr(terms(f), "variables"))[-1], collapse = " + "),
             model = i) |>
      relocate(model, formula) %>%
      nest(results = c(AUC, TSS, logLik, AIC, BIC))
  }
  results_summary <- results_summary |> bind_rows() # convert from list to tibble
  return(results_summary)
}
