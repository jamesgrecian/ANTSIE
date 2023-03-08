###############################################################
### Helper function to check correlation between covariates ###
###############################################################

# 2023-02-15

# function will return correlation betwen all inputted columns
# only pass the environmental data columns to the function
# output is then filtered by threshold to show those covariates with high correlation

check_cor <- function(dat, threshold = .7){
  
  # pull names for variables
  variables <- names(dat)
  
  # calculate all possible pairs
  test <- combn(variables, 2, simplify = T)
  test <- t(test) %>% as_tibble()
  
  # calculate all pairwise correlations
  test <- test %>% rowwise() %>% mutate(correlation = cor(dat[c(V1)], dat[c(V2)]))
  
  # return only those that exceed threshold
  out <- test %>% filter(correlation > threshold)
  return(out)
}

# ends