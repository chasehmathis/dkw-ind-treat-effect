
# Compare our Makarov quantile method and save simulation results
set.seed(123)
source("generate.R")
source("bounds.R")
library(RIQITE)
compare_makarov_quantile <- function(params) {
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

  quant_true <- quantile(data$tau, 0.5)
  conf_int_idx <- which(diff_bounds[,1] < 0.5 & diff_bounds[,2] > 0.5)
  conf_int_us <- list(lower = min(diff_bounds[conf_int_idx, 3]), upper = max(diff_bounds[conf_int_idx, 3]))
  conf_riqite <- ci_quantile(data$Z, data$Y, 0.5, nperm = 1e4, alpha = params$alpha)
  if(abs(conf_int_us$lower - naive_l) < 0.25) {
    conf_int_us$lower <- -Inf
  }
  if(abs(conf_int_us$upper - naive_u) < 0.25){
  conf_int_us$upper <- Inf
  }
  width_us <- conf_int_us$upper - conf_int_us$lower
  cover_us <- quant_true > conf_int_us$lower & quant_true < conf_int_us$upper
  width_riqite <- conf_riqite$upper - conf_riqite$lower
  cover_riqite <- quant_true > conf_riqite$lower & quant_true < conf_riqite$upper

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
      cover = cover_us,
      width_riqite = width_riqite,
      cover_riqite = cover_riqite
    )
  )
}

set.seed(123) # for reproducibility
NSIM <- 20 # or set to desired number of simulations
m_values <- c(100,500,1e3,2e3)
sim_results <- data.frame()
for(tight_bounds in c(FALSE, TRUE)){
  for(alpha in c(0.2, 0.3, 0.4)){
    for (size in m_values) {
      params <- list(
        m = size,
        p = 0.5,
        alpha = alpha,
        heterogeneity = 3,
        outcome_model = rlnorm,
        tight_bounds = tight_bounds
      )
      
      for (sim in 1:NSIM) {
        print(sim)
        sim_results <- rbind(sim_results, compare_makarov_quantile(params))
      }
    }
  }
}


sim_results |> 
  dplyr::group_by(tight_bounds,m,alpha) |> 
  dplyr::summarise(
    width_us = mean(width),
    width_riqite = mean(width_riqite == Inf)
  ) |> print(n =200)
# Optionally, save results to a file
write.csv(sim_results, file = "makarov_quantile_sim_results.csv", row.names = FALSE)

