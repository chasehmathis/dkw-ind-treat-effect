# Demo: One-Sided Makarov Bounds with Large Treatment Effect

library(dkw)
library(RIQITE)

set.seed(123)

# Parameters - large treatment effect for clear visualization
m <- 200
alpha <- 0.1
tau_true <- 20  # Large effect

# Generate data
cat("Generating data with large treatment effect (tau =", tau_true, ")...\n")
result <- get_bounds(m = m, p = 0.5, alpha = alpha, tau = tau_true, rho = 0.5)

bounds <- result$bounds
data <- result$data
obs_Y1 <- sort(result$obs_Y1)
obs_Y0 <- sort(result$obs_Y0)

# Extend bounds for edge cases
bounds_ext <- lapply(bounds, function(z) c(z, 0, 1))

# Compute Makarov bounds
cat("Computing Makarov bounds...\n")
diff_bounds <- get_makarov_bounds(bounds_ext, obs_Y1, obs_Y0,
                                   t_range = c(-10, 50), n_points = 500)

# Plot results
cat("Plotting...\n")


# Plot 1: Makarov bounds vs true CDF
plot(ecdf(data$tau), main = "Makarov Bounds vs True CDF",
     xlab = "Treatment Effect", ylab = "CDF")
points(diff_bounds[, 3], diff_bounds[, 2], col = "red", pch = 16, cex = 0.5)
points(diff_bounds[, 3], diff_bounds[, 1], col = "blue", pch = 16, cex = 0.5)
legend("bottomright", c("True CDF", "Upper Bound", "Lower Bound"),
       col = c("black", "red", "blue"), pch = c(NA, 16, 16), lty = c(1, NA, NA))

# Plot 2: RIQITE comparison
cat("Computing RIQITE bounds...\n")
ci_riqite <- ci_quantile(
  Z = data$Z,
  Y = data$Y,
  alternative = "greater",
  method.list = list(name = "Stephenson", s = 10),
  nperm = 1e4,
  alpha = alpha,
  k.vec = floor(seq(0.1, 1, by = 0.1)*m)
)
ci <- extract_quantile_ci(diff_bounds)
ci$k <- ci$k * m

plot_quantile_ci(ci_riqite, main = "RIQITE Method (Stephenson s=10)",
                ci, k_start = 0, xlim = c(15,22))

# Summary statistics
cat("\n--- Summary ---\n")
cat("True ATE:", mean(data$tau), "\n")
cat("Estimated ATE:", mean(obs_Y1) - mean(obs_Y0), "\n")
cat("True median tau:", median(data$tau), "\n")

# Extract quantile CIs
cat("\nMakarov Quantile CIs:\n")
print(ci)

cat("\nDone!\n")
