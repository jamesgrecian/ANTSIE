#############################################################
### Function to perform k-folds cross-validation on a GAM ###
#############################################################

# 2022-12-24

# This function is based on Abrahms
# https://github.com/briana-abrahms/DynamicEnsembleSDM/blob/master/model_evaluation.R
# modified for a spatial k folds procedure - see spatial_kfolds R function

# Given an input dataframe, a formula, and k number of folds
# Function will split the data into k spatial clusters
# Fit the model on k-1 clusters and test on remaining
# Extract AUC and TSS statistics
# Repeat this for all k folds
# Return dataframe containing test statistics

# depends on spatial_kfolds, dismo and mgcv

kfoldCV <- function(data, f, k){
  
  input <- data
  input$K <- spatial_kfolds(input, k) #randomly allocate k groups
  
  # output dataframe
  evaluations <- tibble(K = rep(0, k),
                        AUC = rep(0, k),
                        TSS = rep(0, k))
  
  # store spatial coefficient of variation for each fold (standard error/ predicted value for each point)
  # use this to check there isn't a specific region with issues
  spatial_cv <- tibble(x = input$x,
                       y = input$y,
                       K = input$K,
                       cv = rep(NA, nrow(input)))
  
  # set up progress bar
  pb1 <- txtProgressBar(min = 1, max = k, style = 3)
  
  # loop through k folds
  for (i in 1:k){
    setTxtProgressBar(pb1, i) # update progress bar
    # divide
    train <- input[input$K != i,]
    test <- input[input$K == i,]
    # run model
    model_k <- mgcv::gam(formula = f,
                         data = train,
                         family = binomial,
                         select = T,
                         method = "REML")
    # generate predictions
    preds <- mgcv::predict.gam(model_k, test, se = T, type = "response") %>% as_tibble
    
    # spatial coefficient of variation
    cv <- mgcv::predict.gam(model_k, test, se = T, type = "response") %>% as_tibble
    cv <- cv$se.fit/cv$fit
    
    # combine predictions with observations from testing data set
    d <- preds %>% mutate(PresAbs = test$PresAbs)

    # pull out predicted presence and absences for model evaluation
    pres <- d %>% filter(PresAbs == 1) %>% pull(fit) %>% as.numeric()
    abs <- d %>% filter(PresAbs == 0) %>% pull(fit) %>% as.numeric()
    e <- dismo::evaluate(p = pres, a = abs)
    
    # store summary statistics
    evaluations[i, 1] <- k
    evaluations[i, 2] <- e@auc
    evaluations[i, 3] <- max(e@TPR + e@TNR-1)
    
    spatial_cv$cv[spatial_cv$K == i] <- cv
  }
  
#  # generate plot for visualising the spatial cv
#  # define the bounding box
#  test <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = -20)) %>% 
#    sf::st_as_sfc() %>%
#    st_set_crs(4326)
#  # split the bounding box into k segments
#  test <- test %>% st_make_grid(n = c(k, 1))
  
#  # project the box into polar projection
#  prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#  test <- test %>% st_transform(prj)
#  test <- test %>% st_as_sf
#  
#  p1 <- ggplot() +
#    labs(title = "Coefficient of Variation",
#            subtitle = "For each spatial k-fold") +
#    theme_bw() +
#    geom_sf(aes(), fill = NA, data = test) +
#    geom_point(aes(x = x, y = y, colour = cv), data = spatial_cv) +
#    coord_sf(xlim = c(min(spatial_cv$x), max(spatial_cv$x)),
#             ylim = c(min(spatial_cv$y), max(spatial_cv$y)),
#             crs = prj) +
#    scale_colour_viridis_c() +
#    xlab(NULL) + ylab(NULL)
#  print(p1)
#  cat("\nmean AUC =", round(mean(foo[[1]]$AUC), 3))
#  cat("\nmean TSS =", round(mean(foo[[1]]$TSS), 3), "\n")

  return(list(evaluations, spatial_cv))
  
}
