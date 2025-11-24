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
  
  n <- length(obs_Y1)
  # OPTIMIZATION: Use findInterval instead of sapply + which
  idx <- findInterval(x, obs_Y1)
  below_min <- (idx == 0)
  above_max <- (idx == n & x > obs_Y1[n])
  
  idx[below_min] <- n + 1
  idx[above_max] <- n + 2

  
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
  n <- length(obs_Y0)
  # OPTIMIZATION: Use findInterval instead of sapply + which
  idx <- findInterval(x, obs_Y0)
  below_min <- (idx == 0)
  above_max <- (idx == n & x > obs_Y0[n])
  
  idx[below_min] <- n + 1
  idx[above_max] <- n + 2
  return(cbind(l_bound_0[idx], u_bound_0[idx]))
}

#' Compute Makarov bounds for individual treatment effect at quantile t
#'
#' @param t Quantile value
#' @param bounds Pre-computed ECDF bounds
#' @param obs_Y1 Sorted observed outcomes in treatment group
#' @param obs_Y0 Sorted observed outcomes in control group
#' @param minY Optional: pre-computed minimum of obs_Y1 (for efficiency)
#' @param maxY Optional: pre-computed maximum of obs_Y1 (for efficiency)
#' @return Vector with (lower_bound, upper_bound, t)
marakov_bounds_t <- function(t, bounds, obs_Y1, obs_Y0, minY = NULL, maxY = NULL) {
  # OPTIMIZATION: Use pre-computed min/max if provided, else compute once
  if (is.null(minY)) minY <- min(obs_Y1)
  if (is.null(maxY)) maxY <- max(obs_Y1)
  
  L <- 15
  min_grid <- minY - L - abs(t)
  max_grid <- maxY + L + abs(t)
  
  x <- seq(min_grid, max_grid, by = 0.005)

  Y1_bound <- ecdf_Y1_bound(x, bounds, obs_Y1)
  Y0_bound <- ecdf_Y0_bound(x - t, bounds, obs_Y0)
  lower_arg <- pmax(Y1_bound[,1] - Y0_bound[,2],0)
  upper_arg <- pmin(Y1_bound[,2] - Y0_bound[,1] + 1, 1)
  # OPTIMIZATION: Use vector operations more efficiently
  l_bound <- max(lower_arg)  # Explicitly use first column for clarity
  u_bound <- min(upper_arg)  # Explicitly use second column
  
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
  bounds <- lapply(bounds, \(z) c(z, 0,1))
  L <- 10
  t <- seq(-2*L, 2*L, length.out = 1000) # L is assumed bound on PO

  # CRITICAL OPTIMIZATION: Pre-allocate result matrix instead of using rbind in loop
  # rbind in loop is O(n²) - this is O(n)
  n_t <- length(t)
  res_bounds <- matrix(nrow = n_t, ncol = 3)
  
  # OPTIMIZATION: Compute min/max once outside the loop
  minY <- min(obs_Y1)
  maxY <- max(obs_Y1)
  
  for(i in seq_along(t)) {
    # Debug condition (keep for now)
    #if(abs(t[i] + 5) < 0.1){
    # if(length(obs_Y1) + length(obs_Y0) == 1e4){browser()}
    #}
    # OPTIMIZATION: Fill pre-allocated matrix instead of rbind
    res_bounds[i, ] <- marakov_bounds_t(t[i], bounds, obs_Y1, obs_Y0, minY, maxY)
  }
  
  # OPTIMIZATION: No need to remove first row since we pre-allocated correctly
  return(res_bounds)
}