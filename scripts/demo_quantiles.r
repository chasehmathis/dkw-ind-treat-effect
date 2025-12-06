# Demo: Quantile Confidence Intervals with DKW Bounds

library(dkw)
library(DeclareDesign)
library(RIQITE)

set.seed(123)

m <- 200
alpha <- 0.05
nsim <- 1000

# Generate data with lognormal outcomes
design <- declare_model(
  N = m,
  potential_outcomes(Y ~ rlnorm(N, meanlog = 1 * Z, sdlog = 1.2))
) +
  declare_assignment(Z = complete_ra(N, prob = 1/2)) +
  declare_measurement(Y = reveal_outcomes(Y ~ Z))

data_truth <- draw_data(design)
data_truth$tau <- data_truth$Y_Z_1 - data_truth$Y_Z_0

# Generate permutations for simulation
generate_perm <- function(nperm, data_truth, j = 0) {
  Ys <- matrix(nrow = nperm, ncol = m / 2)
  for (i in 1:nperm) {
    new_Z <- sample(data_truth$Z)
    new_Y <- data_truth$Y_Z_1 * new_Z + (1 - new_Z) * data_truth$Y_Z_0
    Ys[i, ] <- new_Y[new_Z == j]
  }
  return(Ys)
}

Y0s <- generate_perm(nperm = nsim, data_truth = data_truth)
Y1s <- generate_perm(nperm = nsim, data_truth = data_truth, j = 1)

# Observed data
true_Y0 <- data_truth$Y[data_truth$Z == 0]
true_Y1 <- data_truth$Y[data_truth$Z == 1]
true_tau_sorted <- sort(true_Y1 - true_Y0)
ecdf_tau <- ecdf(true_tau_sorted)(true_tau_sorted)

# DKW bounds
eps <- sqrt(-log(alpha / 4) * 2 / m)
u_bound <- pmin(ecdf_tau + eps, 1)
l_bound <- pmax(ecdf_tau - eps, 0)

# Plot ECDF with simulations
plot(true_tau_sorted, ecdf_tau, type = "n", ylim = c(0, 1),
     main = "ECDF with DKW Bounds and Permutation Simulations",
     ylab = "ECDF", xlab = "Treatment Effect")

# Overlay simulated ECDFs
for (i in 1:nsim) {
  tau_sim <- Y1s[i, ] - Y0s[i, ]
  ecdf_sim <- ecdf(tau_sim)
  lines(true_tau_sorted, ecdf_sim(true_tau_sorted),
        type = "s", col = rgb(0.7, 0.7, 0.7, 0.15))
}

# Plot actual ECDF and bounds
lines(true_tau_sorted, ecdf_tau, type = "s", lwd = 2, col = "black")
lines(true_tau_sorted, u_bound, col = "red", lty = 2, lwd = 2)
lines(true_tau_sorted, l_bound, col = "blue", lty = 2, lwd = 2)
legend("bottomright", c("ECDF", "Upper Bound", "Lower Bound"),
       col = c("black", "red", "blue"), lty = c(1, 2, 2), lwd = 2)

# Median confidence interval
cat("\nMedian treatment effect:", true_tau_sorted[m / 4], "\n")

