# Compute Bounds for Treatment Effect Distribution
# This module handles bound computation for ECDFs and Makarov bounds

#' Compute bounds for F_1(x) (treatment group ECDF)
#'
#' @param x Point(s) to evaluate
#' @param bounds Pre-computed bounds list
#' @param obs_Y1 Sorted observed outcomes in treatment group
#' @return Matrix with lower and upper bounds
ecdf_Y1_bound <- function(x, bounds, obs_Y1) {
  u_bound_1 <- bounds[[1]]
  l_bound_1 <- bounds[[2]]
  
  # Find closest index
  idx <- sapply(x, \(z) which.min(abs(obs_Y1 - z)))
  return(cbind(l_bound_1[idx], u_bound_1[idx]))
}

#' Compute bounds for F_0(x) (control group ECDF)
#'
#' @param x Point(s) to evaluate
#' @param bounds Pre-computed bounds list
#' @param obs_Y0 Sorted observed outcomes in control group
#' @return Matrix with lower and upper bounds
ecdf_Y0_bound <- function(x, bounds, obs_Y0) {
  u_bound_0 <- bounds[[3]]
  l_bound_0 <- bounds[[4]]
  
  # Find closest index
  idx <- sapply(x, \(z) which.min(abs(obs_Y0 - z)))
  return(cbind(l_bound_0[idx], u_bound_0[idx]))
}

#' Compute Makarov bounds for individual treatment effect at quantile t
#'
#' @param t Quantile value
#' @param bounds Pre-computed ECDF bounds
#' @param obs_Y1 Sorted observed outcomes in treatment group
#' @param obs_Y0 Sorted observed outcomes in control group
#' @return Vector with (lower_bound, upper_bound, t)
marakov_bounds_t <- function(t, bounds, obs_Y1, obs_Y0) {
  x <- seq(-50, 50)
  Y1_bound <- ecdf_Y1_bound(x, bounds, obs_Y1)
  Y0_bound <- ecdf_Y0_bound(x - t, bounds, obs_Y0)

  argument <- cbind(pmax(Y1_bound[,1] - Y0_bound[,1], 0),
                    pmin(Y1_bound[,2] - Y0_bound[,2], 1))
  
  l_bound <- max(argument)
  u_bound <- min(argument + 1)
  
  c(l_bound, u_bound, t)
}

#' Get Makarov bounds for all quantiles
#'
#' @param bounds Pre-computed ECDF bounds
#' @param obs_Y1 Sorted observed outcomes in treatment group
#' @param obs_Y0 Sorted observed outcomes in control group
#' @param data True data with tau values
#' @return Matrix with bounds for all t values
get_marakov_bounds_all_t <- function(bounds, obs_Y1, obs_Y0, data) {
  t <- seq(-20, 80, length.out = 1000)
  
  res_bounds <- matrix(ncol = 3)
  for(i in t) {
    res_bounds <- rbind(res_bounds, marakov_bounds_t(i, bounds, obs_Y1, obs_Y0))
  }
  
  res_bounds <- res_bounds[2:nrow(res_bounds),]

  plot.ecdf(data$tau, main = "true ECDF of individual tau")
  points(res_bounds[,3], res_bounds[,2], col = "red")
  points(res_bounds[,3], res_bounds[,1], col = "blue")
  
  return(res_bounds)
}

