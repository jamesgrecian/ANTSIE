##################
### cmip_stack ###
##################

# 2023-02-09

# helper function to load cmip6 data
# input path to folder
# variable of interest to search for i.e. tos, sos, mlotst 
# time period i.e. historical or lgm

cmip_stack <- function(path, var, period){
  fn <- list.files(path, full.names = T, pattern = "grd")
  fn <- fn[grepl("stere", fn)]
  fn <- fn[grepl(var, fn)]
  fn <- fn[grepl(period, fn)]
  nms <- sapply(strsplit(fn, "_"), "[[", 2)
  stk <- stack(fn)
  names(stk) <- nms
  return(stk)
}

#ends