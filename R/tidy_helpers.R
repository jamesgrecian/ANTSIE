####################
### Tidy helpers ###
####################

# 2023-04-14

# Three helper functions to use with the tidy model function


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

# function to calculate true skill statistic
# based on dismo package
tss <- function(x){
  p <- as.numeric(x$preds[x$PresAbs == 1]) # predicted values for known presences 
  a <- as.numeric(x$preds[x$PresAbs == 0]) # predicted values for known absences
  e <- dismo::evaluate(p, a)
  max(e@TPR + e@TNR-1)
}
