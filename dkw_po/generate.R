# Generate Data and Compute Initial Bounds
# This module handles data generation and initial bound computation
# install.packages("DeclareDesign) good for potential outcome simualtion
library(DeclareDesign)
library(MASS)
#' Special ECDF function for 1/(pn). Returns a modified ecdf
#' 
#' @parm data vector
#' @parm p of BRE
ecdf_deterministic <- function(data, n, p = 0.5) {

  force(n); force(p)  # ensure they're captured
  
  res <- function(xs) {
    sapply(xs, function(x) {
      less_than_x <- sum(data < x)
      min((1 / (n * p)) * less_than_x,1)
    })
  }
  
  return(res)
}

#' Generate true sample data with varying levels of heterogeneity and
#' outcome distribution modelling
#'
#' @param m Sample size
#' @return Data frame with potential outcomes
generate_true_sample <- function(m, p, tau, rho) {

  design <-
    declare_model(
      N = m,
      draw_multivariate(
        c(Y_Z_0, Y_Z_1) ~ mvrnorm(n = m, mu = c(0,tau), Sigma = matrix(c(1,rho,rho,1),2))
      )
    ) +
    declare_assignment(Z = complete_ra(N = m, m = m/2)) + 
    declare_measurement(Y = reveal_outcomes(Y ~ Z))
  data_true <- draw_data(design)
  data_true[["tau"]] <- data_true$Y_Z_1 - data_true$Y_Z_0
  return(data_true)
}


#' Get initial bounds for ECDFs
#'
#' @param m Sample size
#' @param alpha Significance level
#' @return List containing bounds, data, and observed outcomes
get_bounds <- function(m = 300, p = 0.6, alpha = 0.1, 
                       rho = 0.4,tau = 1) {
  
  # Generate true sample data
  data_truth <- generate_true_sample(m, p = p, tau = tau, rho = rho)

  
  # Observed sample
  obs_Y0 <- data_truth$Y[data_truth$Z == 0]
  obs_Y1 <- data_truth$Y[data_truth$Z == 1]

  # ECDF 0 and evaluation

  #ecdf_fn0 <- ecdf_deterministic(obs_Y0, n = m,p = 1-p)
  ecdf_fn0 <- ecdf(obs_Y0)
  obs_Y0_sorted <- sort(obs_Y0)
  ecdf_Y0_sorted <- ecdf_fn0(obs_Y0_sorted)
  
  # ECDF 1 and evaluation
  #ecdf_fn1 <- ecdf_deterministic(obs_Y1,n = m, p = p)
  ecdf_fn1 <- ecdf(obs_Y1)
  obs_Y1_sorted <- sort(obs_Y1)
  ecdf_Y1_sorted <- ecdf_fn1(obs_Y1_sorted)
  
  # Compute DKW bounds
  # regular eps (tight)
  eps1 <- sqrt(-log((alpha/2)/2) * (1/(2 * m/2)))
  eps0 <- sqrt(-log((alpha/2)/2) * (1/(2 * m/2)))
  eps <- eps1 + eps0


  
  u_bound_1 <- pmin(ecdf_Y1_sorted +eps, 1)
  l_bound_1 <- pmax(ecdf_Y1_sorted -eps, 0)
  u_bound_0 <- pmin(ecdf_Y0_sorted +eps, 1)
  l_bound_0 <- pmax(ecdf_Y0_sorted -eps, 0)
  
  bounds <- list(u_bound_1,
                 l_bound_1,
                 u_bound_0,
                 l_bound_0)
  
  return(list(bounds = bounds, 
              data = data_truth, 
              obs_Y1 = obs_Y1,
              obs_Y0 = obs_Y0))
}


