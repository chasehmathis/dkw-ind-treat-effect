# Generate Data and Compute Initial Bounds
# This module handles data generation and initial bound computation
# install.packages("DeclareDesign) good for potential outcome simualtion
library(DeclareDesign)

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

#' Generate true sample data
#'
#' @param m Sample size
#' @return Data frame with potential outcomes
generate_true_sample <- function(m, p = 0.5) {
  design <-
    declare_model(
      N = m,
      # baseline covariates
      X1 = rnorm(N, 0, 1),
      X2 = sample(letters[1:3], N, replace = TRUE),
      # potential outcomes depending on covariates + treatment
      potential_outcomes(
        Y ~ rlnorm(N, meanlog = 0.3 + 0.5*X1 + 0.2*(X2 == "b") + 1.3*Z,
                   sdlog = 0.5)
      )
    ) +
    declare_assignment(Z = rbinom(n = N, size = 1, prob = p)) +
    declare_measurement(Y = reveal_outcomes(Y ~ Z))
  
  data_true <- draw_data(design)
  data_true[["tau"]] <- data_true$Y_Z_1 - data_true$Y_Z_0
  return(data_true)
}

#' Generate permutation data
#'
#' @param nperm Number of permutations
#' @param data_truth True data
#' @param j Treatment assignment (0 or 1)
#' @return Matrix of permuted outcomes
generate_perm <- function(nperm, data_truth, j = 0) {
  Ys <- matrix(nrow = nperm, ncol = nrow(data_truth)/2)
  for(i in 1:nperm) {
    new_Z <- sample(data_truth$Z)
    new_Y <- data_truth$Y_Z_1 * new_Z + (1-new_Z) * data_truth$Y_Z_0
    Ys[i,] <- new_Y[new_Z == j]
  }
  return(Ys)
}

#' Get initial bounds for ECDFs
#'
#' @param m Sample size
#' @param alpha Significance level
#' @param nsim Number of simulations (currently unused, kept for compatibility)
#' @return List containing bounds, data, and observed outcomes
get_bounds <- function(m = 3000, p = 0.6, alpha = 0.1, nsim = 1e3) {
  
  # Generate true sample data
  data_truth <- generate_true_sample(m, p = p)
  
  # Observed sample
  obs_Y0 <- data_truth$Y[data_truth$Z == 0]
  obs_Y1 <- data_truth$Y[data_truth$Z == 1]

  # ECDF 0 and evaluation
  ecdf_fn0 <- ecdf_deterministic(obs_Y0, n = m,p = 1-p)
  obs_Y0_sorted <- sort(obs_Y0)
  ecdf_Y0_sorted <- ecdf_fn0(obs_Y0_sorted)
  
  # ECDF 1 and evaluation
  ecdf_fn1 <- ecdf_deterministic(obs_Y1,n = m, p = p)
  obs_Y1_sorted <- sort(obs_Y1)
  ecdf_Y1_sorted <- ecdf_fn1(obs_Y1_sorted)
  
  # Compute DKW bounds
  eps1 <- sqrt(-log((alpha/2)/4)*8/(m*p))
  eps0 <- sqrt(-log((alpha/2)/4)*8/(m*(1-p)))
  eps <- eps1 + eps0
  browser()
  # regular eps (tight)
  # eps1 <- sqrt(-log((alpha/2)/2) * (1/(2 * m * p)))
  # eps0 <- sqrt(-log((alpha/2)/2) * (1/(2 * m * (1-p))))
  # eps <- eps1 + eps0
  
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

