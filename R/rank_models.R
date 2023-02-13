##########################
### Spatial CV ranking ###
##########################

# 2022-12-24

# This is a wrapper function bringing together kfoldCV and spatial_kfolds
# Input is a focal dataset, model formula and number of spatial k-folds

# Function will fit full model and subsequent models with one term removed
# Function will return nested tibble containing all metrics for evaluation

rank_models <- function(data, f, k){
  
  # how many terms does the formula have?
  n.terms <- length(attr(terms(f), "term.labels"))
  
  # first run full model
  m <- gam(f, data, family = binomial, method = "REML", select = T)
  foo <- kfoldCV(data, f, k)
  
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
    m <- gam(new_f, data, family = binomial, method = "REML", select = T)
    foo <- kfoldCV(data, new_f, k)
    
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


# we might need all combinations of the covariates
# this function breaks down the formula and returns all possible combinations
# does not work on interactions
all_combs <- function(f){
  
  new_f <- f
  
  # how many terms does the formula have?
  n.terms <- length(attr(terms(f), "term.labels"))
  
  for(i in 1:n.terms){
    t <- combn(attr(terms(f), "term.labels"), m = i, simplify = F)
    
    for(j in 1:length(t)){
      t2 <- t[[j]]
      temp_f <- reformulate(paste(t2, collapse = " + "), response = all.vars(attr(terms(f), "variables"))[1])
      new_f <- c(new_f, temp_f)
    }
  }
  return(new_f)
}

# given a 
rank_combs <- function(data, f, k){
  
  # how many terms does the formula have?
  n.terms <- length(attr(terms(f), "term.labels"))
  
  # generate list of all combinations
  new_f <- all_combs(f)
  
  # first run full model
  m <- gam(new_f[[1]], data, family = binomial, method = "REML", select = T)
  foo <- kfoldCV(data, new_f[[1]], k)
  
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
  return(out)
}

