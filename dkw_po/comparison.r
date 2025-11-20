
# Compare our Makarov quantile method and save simulation results

compare_makarov_quantile <- function(params) {
  result <- do.call(get_bounds, params)
  bounds <- result$bounds
  data <- result$data
  obs_Y1 <- sort(result$obs_Y1)
  obs_Y0 <- sort(result$obs_Y0)

  # Compute Makarov bounds for all t values
  diff_bounds <- get_marakov_bounds_all_t(bounds, obs_Y1, obs_Y0, data)

  quant_true <- quantile(data$tau, 0.5)
  conf_int_idx <- which(diff_bounds[,1] < 0.5 & diff_bounds[,2] > 0.5)
  if (length(conf_int_idx) == 0) {
    # If interval not found, set as NA
    conf_int_us <- list(lower = NA, upper = NA)
    width_us <- NA
    cover_us <- NA
  } else {
    conf_int_us <- list(lower = min(diff_bounds[conf_int_idx, 3]), upper = max(diff_bounds[conf_int_idx, 3]))
    width_us <- conf_int_us$upper - conf_int_us$lower
    cover_us <- quant_true > conf_int_us$lower & quant_true < conf_int_us$upper
  }
  return(
    data.frame(
      m = params$m,
      p = params$p,
      alpha = params$alpha,
      heterogeneity = params$heterogeneity,
      tight_bounds = params$tight_bounds,
      quant_true = quant_true,
      ci_lower = conf_int_us$lower,
      ci_upper = conf_int_us$upper,
      width = width_us,
      cover = cover_us
    )
  )
}

set.seed(123) # for reproducibility
NSIM <- 200 # or set to desired number of simulations
m_values <- c(1e2, 5e2, 1e3)
sim_results <- data.frame()

for (size in m_values) {
  params <- list(
    m = size,
    p = 0.5,
    alpha = 0.1,
    heterogeneity = 3,
    outcome_model = rlnorm,
    tight_bounds = FALSE
  )
  for (sim in 1:NSIM) {
    print(sim)
    sim_results <- rbind(sim_results, compare_makarov_quantile(params))
  }
}

# Optionally, save results to a file
write.csv(sim_results, file = "makarov_quantile_sim_results.csv", row.names = FALSE)

