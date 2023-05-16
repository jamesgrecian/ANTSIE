################################################################################
### Helper function to filter list of formulas based on two input covariates ###
################################################################################

# input is a list of formulas - the output from 'all_combs'
# and two variables that have been found to be colinear
# output will be all combinations that DO NOT contain those two variables

filter_vars <- function(f, V1, V2){
  
  id1 <- grep(V1, f) # find indices of V1
  id2 <- grep(V2, f) # find indices of V2
  
  ids <- id1[id1 %in% id2] # find indices that match both
  idx <- 1:length(f) # dummy index for length of function
  kp <- setdiff(idx, ids) # take difference to give indices of those we keep
  out <- f[kp] # subset formula by indices of those we want to keep
  cat("dropped", length(f) - length(out), "variable sets", "\n")
  return(out)
}
