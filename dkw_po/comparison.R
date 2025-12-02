
# Compare our Makarov quantile method and save simulation results
set.seed(123)

source("generate.R")
source("bounds.R")
source("plot_quantile_CI_comparison.r")
library(RIQITE)
library(MASS)
library(ggplot2)
plot.dir <- "../fig/"
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

  quants_check <- seq(0.1,1,by = 0.1)
  quants_true <- quantile(data$tau, quants_check)
  conf_int_idx <- lapply(quants_check, \(z) which(diff_bounds[,2] >= z))
  # just do lower bound since that is what riqite does well in 
  conf_int_us <- lapply(conf_int_idx, \(z)  max( c(min(diff_bounds[z, 3])), naive_l))
  conf_riqite <- ci_quantile(data$Z, data$Y, quants_check*nrow(data), 
                             nperm = 1e3, alpha = params$alpha)

  for(i in 1:length(conf_int_us)){
    if(abs(conf_int_us[[i]] - naive_l) < 0.25) {
      conf_int_us[[i]] <- -Inf
    }
  }

  cover_us <- quants_true > conf_int_us
  cover_riqite <- quants_true > conf_riqite$lower
  width_diff <- unlist(conf_int_us) - conf_riqite$lower

  return(
    data.frame(
      m = params$m,
      p = params$p,
      tau = params$tau,
      rho = params$rho,
      lower_us = unlist(conf_int_us),
      lower_riqite = conf_riqite$lower,
      cover_us = cover_us,
      cover_riqite = cover_riqite,
      quants_check = quants_check
    )
  )
}

set.seed(123) # for reproducibility
NSIM <- 20 # or set to desired number of simulations
m <- c(1e2, 5e2, 1e3)
sim_results <- data.frame()
for(tau in seq(1,9, length = 3)){
for(rho in c(0.1, 0.3, 0.7)){
      for (size in m) {
        params <- list(
          m = size,
          p = 0.5,
          alpha = 0.1,
          tau = tau,
          rho = rho
        )
      
        for (sim in 1:NSIM) {
          print(sim)
          sim_results <- rbind(sim_results, compare_makarov_quantile(params))
        }
  }
}
}






# Optionally, save results to a file
write.csv(sim_results, file = "makarov_quantile_sim_results.csv", row.names = FALSE)

sim_results <- sim_results |> 
  dplyr::mutate(
    percent_improvement = ifelse(
      is.infinite(width_riqite),
      Inf,
      100 * (width_us - width_riqite) / abs(width_riqite)
    )
  )


library(dplyr)
library(purrr)

# Get all unique values for m, rho, and tau from sim_results
unique_m_vals <- sort(unique(sim_results$m))
unique_rho_vals <- sort(unique(sim_results$rho))
unique_tau_vals <- sort(unique(sim_results$tau))

# For each parameter, compute mean lower_us and lower_riqite by quantile
avg_quantile_results_m <- sim_results %>%
  group_by(m, quants_check) %>%
  dplyr::summarise(
    lower_us = mean(lower_us, na.rm = TRUE),
    lower_riqite = mean(lower_riqite, na.rm = TRUE),
    .groups = "drop"
  )

avg_quantile_results_rho <- sim_results %>%
  group_by(rho, quants_check) %>%
  dplyr::summarise(
    lower_us = mean(lower_us, na.rm = TRUE),
    lower_riqite = mean(lower_riqite, na.rm = TRUE),
    .groups = "drop"
  )

avg_quantile_results_tau <- sim_results %>%
  group_by(tau, quants_check) %>%
  dplyr::summarise(
    lower_us = mean(lower_us, na.rm = TRUE),
    lower_riqite = mean(lower_riqite, na.rm = TRUE),
    .groups = "drop"
  )

# Set output directory for plots
if (!dir.exists(plot.dir)) dir.create(plot.dir)

# Plot for each m value and save (High Quality)
plots_by_m <- map(unique_m_vals, function(m_val) {
  m_result <- avg_quantile_results_m %>% filter(m == m_val)
  riqite_list <- list(
    lower = m_result$lower_riqite,
    k = m_result$quants_check * m_val
  )
  us_list <- list(
    lower = m_result$lower_us
  )
  plot_file <- file.path(plot.dir, sprintf("quantile_CI_comparison_m_%s.png", m_val))
  png(
    filename = plot_file,
    width = 2400,           # Higher resolution
    height = 1800,
    units = "px",
    res = 300,              # 300 DPI for publication quality
    pointsize = 14
  )
  plot_quantile_CI_comparison(
    riqite_list,
    k_start = 1,
    result2 = us_list,
    xlim = c(0,6),
    main = paste0("Quantile Confidence Intervals (m = ", m_val, ")")
  )
  dev.off()
})

# Plot for each rho value and save (High Quality)
plots_by_rho <- map(unique_rho_vals, function(rho_val) {
  rho_result <- avg_quantile_results_rho %>% filter(rho == rho_val)
  riqite_list <- list(
    lower = rho_result$lower_riqite,
    k = rho_result$quants_check * 100 # or a typical N, just for display
  )
  us_list <- list(
    lower = rho_result$lower_us
  )
  plot_file <- file.path(plot.dir, sprintf("quantile_CI_comparison_rho_%s.png", gsub("\\.", "_", as.character(rho_val))))
  png(
    filename = plot_file,
    width = 2400,
    height = 1800,
    units = "px",
    res = 300,
    pointsize = 14
  )
  plot_quantile_CI_comparison(
    riqite_list,
    k_start = 1,
    result2 = us_list,
    xlim = c(0,6),
    main = paste0("Quantile Confidence Intervals (rho = ", rho_val, ")")
  )
  dev.off()
})

# Plot for each tau value and save (High Quality)
plots_by_tau <- map(unique_tau_vals, function(tau_val) {
  tau_result <- avg_quantile_results_tau %>% filter(tau == tau_val)
  riqite_list <- list(
    lower = tau_result$lower_riqite,
    k = tau_result$quants_check * 100 # or a typical N, just for display
  )
  us_list <- list(
    lower = tau_result$lower_us
  )
  plot_file <- file.path(plot.dir, sprintf("quantile_CI_comparison_tau_%s.png", gsub("\\.", "_", as.character(tau_val))))
  png(
    filename = plot_file,
    width = 2400,
    height = 1800,
    units = "px",
    res = 300,
    pointsize = 14
  )
  plot_quantile_CI_comparison(
    riqite_list,
    k_start = 1,
    result2 = us_list,
    xlim = c(tau_val-5, tau_val +3),
    main = paste0("Quantile Confidence Intervals (tau = ", tau_val, ")")
  )
  dev.off()
})

# Optionally, print or save plots if running interactively.
# If you want to save each plot, wrap the plotting call in png() or ggsave() as appropriate.
# figure k




