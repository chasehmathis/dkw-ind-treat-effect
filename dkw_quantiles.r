library(DeclareDesign)

plot_ecdf_with_sims <- function(m = 200, alpha = 0.05, nsim = 1e3) {
  theta_m <- sin(seq(m))
  
  # Helper to draw sample from specified distribution
  generate_true_sample <- function() {

    design <-
      declare_model(
        N = m, 
        potential_outcomes(Y ~ rlnorm(N, meanlog = 1*Z, sdlog = 1.2))
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
  data_truth <- generate_true_sample()
  
  Y0s <- generate_perm(nperm = 1e3, data_truth = data_truth)
  Y1s <- generate_perm(nperm = 1e3, data_truth = data_truth, 1)
  # true sample
  true_Y0 <- data_truth$Y[data_truth$Z == 0]
  true_Y1 <- data_truth$Y[data_truth$Z == 1]
  
  # ECDF 0 and evaluation
  ecdf_fn0 <- ecdf(true_Y0)
  true_Y0_sorted<- sort(true_Y0)
  ecdf_Y0_sorted <- ecdf_fn0(true_Y0_sorted)
  
  # ECDF 1 and evaluation
  ecdf_fn1 <- ecdf(true_Y1)
  true_Y1_sorted<- sort(true_Y1)
  ecdf_Y1_sorted <- ecdf_fn1(true_Y1_sorted)
  
  #ECDF tau and evaluation
  ecdf_fntau <- ecdf(true_Y1 - true_Y0)
  true_tau_sorted <- sort(true_Y1 - true_Y0)
  ecdf_tau_sorted <- ecdf_fntau(true_tau_sorted)
  
  eps <- sqrt(-log(alpha/4)*2/m) # change for union bound
  
  u_bound <- pmin(ecdf_tau_sorted + eps, 1)
  l_bound <- pmax(ecdf_tau_sorted - eps, 0)
  
  # Plot all the simulated ECDFs (step lines), then the main ECDF with its bounds
  plot(true_tau_sorted, ecdf_tau_sorted, type="n", ylim=c(0,1), 
       main="ECDF with Confidence Bounds and Simulations",
       ylab="ECDF", xlab="x")
  
  # Overlay all simulated ECDFs in light gray
  for(i in 1:nsim){
    tau_sim <- Y1s[i,] - Y0s[i,]
    ecdf_sim <- ecdf(tau_sim)
    lines(true_tau_sorted, ecdf_sim(true_tau_sorted), type="s", col=rgb(0.7,0.7,0.7,0.15))
  }
  
  # Plot actual ECDF of the new sample
  lines(true_tau_sorted, ecdf_tau_sorted, type="s", lwd=2, col="black")
  
  # Plot bounds
  lines(true_tau_sorted, u_bound, col="red", lty=2, lwd=2)
  lines(true_tau_sorted, l_bound, col="blue", lty=2, lwd=2)
  
  return(list(true_tau_sorted, l_bound, u_bound, data_truth))
}

res <- plot_ecdf_with_sims(m = 200)
tau_sorted <- res[[1]]
l_bound <- res[[2]]; u_bound <- res[[3]]
data_truth <- res[[4]]
# observe median is
tau_median <- tau_sorted[50]

# give confidence interval for median
qtls <- ecdf(tau_sorted)(tau_sorted)
print(paste(
  tau_median - tau_sorted[round(l_bound[50] * (m/2))],
  tau_median + tau_sorted[round(u_bound[50] * (m/2))]
)
)

# give confidence interval for max
tau_max <- tau_sorted[100]
print( paste ( 
  tau_max - tau_sorted[round(l_bound[100] * m/2)],
  tau_max - tau_sorted[round(u_bound[100] * m/2)]
  )
)

library(RIQITE) # xinrans method

ci6 = ci_quantile( Z=data_truth$Z, Y=data_truth$Y,
                   alternative="greater", method.list=list(name="Stephenson", s=10),
                   nperm=1e5, alpha=0.05 )
plot_quantile_CIs(ci6, main="s=6", k_start = n-37)

