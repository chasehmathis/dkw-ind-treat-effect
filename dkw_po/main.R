# Main Execution Script
# This script orchestrates the DKW bound computation and comparison with RIQITE
setwd("~/Research/dkw/dkw_po")
source("generate.R")
source("bounds.R")

# Generate data and compute initial bounds
cat("Generating data and computing initial bounds...\n")
# high level params for stratifying results
# Define parameters as a list in R
params <- list(
  m = 500,     # Sample size
  p = 0.6,      # Probability parameter
  alpha = 0.1,   # Significance level
  heterogeneity = 1,
  outcome_model <- rlnorm
  
)
result <- do.call(get_bounds, params)

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
sorted_percentiles_lower <- c( sort(unique(diff_bounds[,1]), decreasing = TRUE))
sorted_percentiles_upper <- c(sort(unique(diff_bounds[,2]), decreasing = TRUE))

grid <- expand.grid(sorted_percentiles_lower, sorted_percentiles_upper)
X <- apply(grid, 1, \(x) which(x[1] == diff_bounds[,1] & x[2] == diff_bounds[,2]))
# get min and max
min_max <- function(data){
  if(nrow(data) == 0){return(NA)}
  min <- apply(data, 2, min)[3]
  max <- apply(data, 2, max)[3]
  c(data[1,1:2], min,max)
}

# prepare for ranks 
diff_bounds[,1:2] <- params$m*diff_bounds[,1:2]
plot_bounds <- sapply(X, \(z) min_max(diff_bounds[z, , drop = F]))
plot_bounds <-plot_bounds[!is.na(plot_bounds)]
mintaus <- min(diff_bounds[,3]);
maxtaus <- max(diff_bounds[,3])

# Create comparison plot
par(mfrow = c(1, 2))
plot_quantile_CIs(ci6, 
                  main = "Recent Method (Caughey et al.)",
                  k_start = 0)

plot(NA, 
     xlab = "k", 
     ylab = expression("confidence"~"limit"~"for" ~ tau[(k)]),
     xlim = c(0, nrow(data)), 
     ylim = c(mintaus, maxtaus),
     main = "Our Closed Form Method")

for (idx in seq_along(plot_bounds)) {

  this_bound <- plot_bounds[[idx]]
  go_through <- seq(this_bound[3], this_bound[4], length.out = 1e2)
  
  for (k in go_through) {
    if(idx != length(plot_bounds)){
      points(this_bound[1], k, pch = 16, col = "black", cex = 0.6)
    }
    if(idx != 1){
      points(this_bound[2], k, pch = 16, col = "black", cex = 0.6)
    }
    lines(c(this_bound[1], this_bound[2] + 1), c(k, k), col = "grey")


  }
}





