##########################################################
### Function to generate all combinations of covarites ###
##########################################################

# WJ Grecian
# 2023-03-07

# we might need all combinations of the covariates
# this function breaks down the formula and returns all possible combinations
# does not work on interactions
all_combs <- function(f, min.terms = 1){
  
  new_f <- f
  
  # how many terms does the formula have?
  n.terms <- length(attr(terms(f), "term.labels"))
  
  for(i in min.terms:n.terms){
    t <- combn(attr(terms(f), "term.labels"), m = i, simplify = F)
    
    for(j in 1:length(t)){
      t2 <- t[[j]]
      temp_f <- reformulate(paste(t2, collapse = " + "), response = all.vars(attr(terms(f), "variables"))[1])
      new_f <- c(new_f, temp_f)
    }
  }
  new_f <- new_f[-1]
  return(new_f)
}
