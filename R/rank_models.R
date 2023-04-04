##########################
### Spatial CV ranking ###
##########################

# 2022-12-24

# This is a wrapper function bringing together kfoldCV and spatial_kfolds
# Input is a focal dataset, model formula and number of spatial k-folds

# Function will fit full model and subsequent models with one term removed
# Function will return nested tibble containing all metrics for evaluation

rank_models <- function(data, f, k, type){
  
  # how many terms does the formula have?
  n.terms <- length(attr(terms(f), "term.labels"))
  
  # first run full model
  m <- gam(f,
           data,
           family = binomial,
           select = T,
           method = "REML")
  foo <- kfoldCV(data, f, k, type)
  
  # create output tibble
  out <- tibble_row(model = 1,
                    formula = paste(substr(attr(terms(f), "term.labels"), 3, 5), collapse = " + "),
                    AIC = AIC(m),
                    AUC = list(foo[[1]]$AUC),
                    TSS = list(foo[[1]]$TSS),
                    spatial_cv = list(foo[[2]]))
  
  # then loop through combination of terms, only remove one each time
  for(i in 1:n.terms){
    new_f <- formula(terms.formula(drop.terms(terms(f), dropx = i, keep.response = T)))
    m <- gam(new_f,
             data,
             family = binomial,
             select = T,
             method = "REML")
    foo <- kfoldCV(data, new_f, k, type)
    
    # append results to output dataframe
    out <- out %>%
      bind_rows(
        tibble_row(model = i+1,
                   formula = paste(substr(attr(terms(new_f), "term.labels"), 3, 5), collapse = " + "),
                   AIC = AIC(m),
                   AUC = list(foo[[1]]$AUC),
                   TSS = list(foo[[1]]$TSS),
                   spatial_cv = list(foo[[2]])))
  }
  return(out)
}


# given a 
rank_combs <- function(data, f, k, type){

  # first run full model
  m <- gam(f[[1]],
           data,
           family = binomial,
           select = T,
           method = "REML")
  foo <- kfoldCV(data, f[[1]], k, type)
  
  # create output tibble
  out <- tibble_row(model = 1,
                    formula = paste(all.vars(attr(terms(f[[1]]), "variables"))[-1], collapse = " + "),
                    AIC = AIC(m),
                    AUC = list(foo[[1]]$AUC),
                    TSS = list(foo[[1]]$TSS),
                    spatial_cv = list(foo[[2]]))
  
  # then loop through combination of terms
  for(j in 2:length(f)){
    m <- gam(f[[j]],
             data,
             family = binomial,
             select = T,
             method = "REML")
    foo <- kfoldCV(data, f[[j]], k, type)
    
    # append results to output dataframe
    out <- out %>%
      bind_rows(
        tibble_row(model = j,
                   formula = paste(all.vars(attr(terms(f[[j]]), "variables"))[-1], collapse = " + "),
                   AIC = AIC(m),
                   AUC = list(foo[[1]]$AUC),
                   TSS = list(foo[[1]]$TSS),
                   spatial_cv = list(foo[[2]])))
  }
  return(out)
}

