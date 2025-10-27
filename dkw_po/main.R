# Main Execution Script
# This script orchestrates the DKW bound computation and comparison with RIQITE

source("generate.R")
source("bounds.R")

# Generate data and compute initial bounds
cat("Generating data and computing initial bounds...\n")
result <- get_bounds()

# Extract results
bounds <- result$bounds
data <- result$data
obs_Y1 <- sort(result$obs_Y1)
obs_Y0 <- sort(result$obs_Y0)

# Compute Makarov bounds for all t values
cat("Computing Makarov bounds for all quantiles...\n")

diff_bounds <- get_marakov_bounds_all_t(bounds, obs_Y1, obs_Y0, data)

# Compare with RIQITE method
cat("Computing bounds using RIQITE method...\n")
library(RIQITE)

ci6 <- ci_quantile(Z = data$Z, 
                   Y = data$Y,
                   alternative = "greater", 
                   method.list = list(name = "Stephenson", s = 10),
                   nperm = 1e5, 
                   alpha = 0.1)

# Prepare plotting
sorted_percentiles <- c(1, sort(unique(diff_bounds[,1]), decreasing = TRUE))

taus <- sapply(sorted_percentiles[-1], function(perc) {
  min(diff_bounds[diff_bounds[,1] == perc, 3, drop = FALSE])
})

xmax <- max(taus) + 10

# Create comparison plot
par(mfrow = c(1, 2))
plot_quantile_CIs(ci6, 
                  main = "Xinran Method: Stephenson Rank=10",
                  k_start = 0)

plot(NA, 
     ylab = "k", 
     xlab = expression("lower" ~ "confidence" ~ "limit" ~ "for" ~ tau[(k)]),
     ylim = c(0, nrow(data)), 
     xlim = c(-20, xmax))

for (idx in seq_along(sorted_percentiles[-1])) {
  this_percentile <- sorted_percentiles[idx]
  next_percentile <- sorted_percentiles[idx+1]
  go_through <- unique(round(seq(next_percentile, this_percentile, 
                                   length.out = 100)*nrow(data)))
  
  tau_i <- taus[idx]
  for (k in go_through) {
    lines(c(tau_i, xmax + 1), c(k, k), col = "grey")
  }
}

