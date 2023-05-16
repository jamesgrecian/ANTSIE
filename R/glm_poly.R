########################################################################################
### Function to generate all possible combinations of n polynomials and n covariates ###
########################################################################################

# 2023-05-03
# code based on discussions with Christine Howard

# function to output all possible polynomial combinations of formulas
# input is a vector of covariate terms that will match the data and the degree of polynomial fits
glm_poly <- function(cov_terms, poly_degree) {
  combs <- as.data.frame(matrix(data = 1:poly_degree, nrow = poly_degree, ncol = length(cov_terms))) # n terms x n poly dataframe
  colnames(combs) <- cov_terms
  combs <- as.matrix(expand.grid(combs))   # calculate all possible combinations of n covariates and n polynomial terms
  
  # use combs to generate a list of formulas
  # output will be same length as nrow of input combs
  poly_form <- lapply(1:nrow(combs), function(i) {
    formula(paste("PresAbs ~ ", paste0("poly(", cov_terms, ", degree = ", combs[i, ], ")", collapse = " + ")))
  })
  return(poly_form)
}


