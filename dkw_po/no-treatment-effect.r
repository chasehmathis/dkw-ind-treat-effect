# no treatment effect

set.seed(123)

source("generate.R")
source("bounds.R")
source("plot_quantile_CI_comparison.r")
library(RIQITE)
library(MASS)
library(ggplot2)
plot.dir <- "../fig/"
size <- 500
params <- list(
  m = size,
  p = 0.5,
  alpha = 0.1,
  tau = 1,
  rho = 0.4
)

result <- do.call(get_bounds, params)
bounds <- result$bounds
data <- result$data
obs_Y1 <- sort(result$obs_Y1)
obs_Y0 <- sort(result$obs_Y0)
# naive bounds
naive_u <- max(obs_Y1) - min(obs_Y0)
naive_l <- min(obs_Y1) - max(obs_Y0)
# Compute Makarov bounds for all t values
diff_bounds <- get_marakov_bounds_all_t(bounds, obs_Y1, obs_Y0, data)

quants_check <- seq(0.1,1,by = 0.1)
quants_true <- quantile(data$tau, quants_check)
conf_int_idx <- lapply(quants_check, \(z) which(diff_bounds[,2] >= z))
# just do lower bound since that is what riqite does well in 
conf_int_us <- lapply(conf_int_idx, \(z)  max( c(min(diff_bounds[z, 3])), naive_l))
conf_int_us$lower <- unlist(conf_int_us)
conf_riqite <- ci_quantile(data$Z, data$Y, quants_check*nrow(data), 
                           nperm = 1e3, alpha = params$alpha)

plot_quantile_CI_comparison(conf_riqite, k_start = 1, result2 = conf_int_us)
