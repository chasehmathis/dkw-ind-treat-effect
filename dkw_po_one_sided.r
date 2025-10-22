#dkw_po one sided
# Two-Sided
library(DeclareDesign)

get_bounds <- function(m = 200, alpha = 0.1, nsim = 1e3) {
  
  # Helper to draw sample from specified distribution
  generate_true_sample <- function(m) {
    
    design <-
      declare_model(
        N = m, 
        potential_outcomes(Y ~ rnorm(N, mean = 20*Z, sd = 1))
        #potential_outcomes(Y ~ rlnorm(N, meanlog = 1*Z, sdlog = 1.2))
      ) +
      declare_assignment(Z = complete_ra(N, prob = 1/2)) +
      declare_measurement(Y = reveal_outcomes(Y ~ Z))
    
    data_true <- draw_data(design)
    data_true[["tau"]] <- data_true$Y_Z_1 - data_true$Y_Z_0
    return(data_true)
  }
  
  generate_perm <- function(nperm, data_truth, j = 0){
    Ys <- matrix(nrow = nperm, ncol = m/2)
    for(i in 1:nperm){
      new_Z <- sample(data_truth$Z)
      new_Y <- data_truth$Y_Z_1 * new_Z + (1-new_Z) * data_truth$Y_Z_0
      Ys[i,] <- new_Y[new_Z == j]
    }
    return(Ys)
  }
  data_truth <- generate_true_sample(m)
  
  # Y0s <- generate_perm(nperm = 1e3, data_truth = data_truth)
  # Y1s <- generate_perm(nperm = 1e3, data_truth = data_truth, 1)
  # true sample
  obs_Y0 <- data_truth$Y[data_truth$Z == 0]
  obs_Y1 <- data_truth$Y[data_truth$Z == 1]
  
  # ECDF 0 and evaluation
  ecdf_fn0 <- ecdf(obs_Y0)
  obs_Y0_sorted<- sort(obs_Y0)
  ecdf_Y0_sorted <- ecdf_fn0(obs_Y0_sorted)
  
  # ECDF 1 and evaluation
  ecdf_fn1 <- ecdf(obs_Y1)
  obs_Y1_sorted<- sort(obs_Y1)
  ecdf_Y1_sorted <- ecdf_fn1(obs_Y1_sorted)
  
  # ECDF tau
  obs_tau <- obs_Y1- obs_Y0
  
  eps <- sqrt(-log(alpha/4)*2/m) # change for union bound
  
  u_bound_1 <- pmin(ecdf_Y1_sorted + eps, 1)
  l_bound_1 <- pmax(ecdf_Y1_sorted - eps, 0)
  
  u_bound_0 <- pmin(ecdf_Y0_sorted + eps, 1)
  l_bound_0 <- pmax(ecdf_Y0_sorted - eps, 0)
  
  bounds <- list(u_bound_1,
                 l_bound_1,
                 u_bound_0,
                 l_bound_0)
  
  return(list(bounds = bounds, data = data_truth, obs_Y1 = obs_Y1,
              obs_Y0 = obs_Y0))
}

# assume have marginals 

result <- get_bounds()
bounds <- result$bounds; data <- result$data
obs_Y1 <- sort(result$obs_Y1); obs_Y0 <- sort(result$obs_Y0)

ecdf_Y1_bound <- function(x){
  
  u_bound_1 <- bounds[[1]]; l_bound_1 <- bounds[[2]]
  # closest_idx
  idx <- sapply(x, \(z) which.min(abs(obs_Y1 - z)))
  return(cbind(l_bound_1[idx], u_bound_1[idx]))
}

ecdf_Y0_bound <- function(x){
  
  u_bound_0 <- bounds[[3]]; l_bound_0 <- bounds[[4]]
  # closest_idx
  idx <- sapply(x, \(z) which.min(abs(obs_Y0 - z)))
  return(cbind(l_bound_0[idx], u_bound_0[idx]))
}

marakov_bounds_t <- function( t){
  x <- seq(-50, 50)
  Y1_bound <- ecdf_Y1_bound(x)
  Y0_bound <- ecdf_Y0_bound(x-t)
  
  argument <- cbind(pmax(Y1_bound[,1] - Y0_bound[,1], 0),
                    pmin(Y1_bound[,2] - Y0_bound[,2], 1))
  
  l_bound <- max(argument)
  u_bound <- min(argument + 1)
  c(l_bound, u_bound, t)
  
  # F1 <- ecdf(data$Y_Z_1); F0 <- ecdf(data$Y_Z_0)
  # diffs <- sapply(x, function(z) F1(z) - F0(z - t))
  
  # GML <- max(0, max(diffs))          # lower bound
  # GMU <- min(1, min(diffs + 1, 1))   # upper bound
  # c(lower = GML, upper = GMU, t)
}

marakov_bounds_t(6)
get_marakov_bounds_all_t <- function() {
  t <- seq(-20, 80, length.out = 1000)
  
  res_bounds <- matrix( ncol = 3)
  for(i in t){
    res_bounds <- rbind(res_bounds, marakov_bounds_t(i))
  }
  
  res_bounds <- res_bounds[2:nrow(res_bounds),]
  
  plot.ecdf(data$tau)
  points(res_bounds[,3], res_bounds[,2],col = "red")
  points(res_bounds[,3], res_bounds[,1], col = "blue")
  return(res_bounds)
}


diff_bounds <- get_marakov_bounds_all_t()


library(RIQITE) # xinran et al method

ci6 = ci_quantile( Z=data$Z, Y=data$Y,
                   alternative="greater", 
                   method.list=list(name="Stephenson", s=10),
                   nperm=1e5, alpha=0.1 )
plot_quantile_CIs(ci6, main="Xinran Method: Stephenson Rank=10")

browser()