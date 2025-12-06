# Monte Carlo Simulation: Compare Makarov Bounds vs RIQITE

devtools::load_all("./")
library(RIQITE)
library(ggplot2)
library(dplyr)
library(purrr)

set.seed(123)

plot.dir <- "output/plots/"
data.dir <- "output/data/"

# Create directories if needed
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)
if (!dir.exists(data.dir)) dir.create(data.dir, recursive = TRUE)

compare_makarov_quantile <- function(params) {
  result <- get_bounds(
    m = params$m,
    p = params$p,
    alpha = params$alpha,
    tau = params$tau,
    rho = params$rho
  )
  bounds <- result$bounds
  data <- result$data
  obs_Y1 <- sort(result$obs_Y1)
  obs_Y0 <- sort(result$obs_Y0)

  # Naive bounds
  naive_u <- max(obs_Y1) - min(obs_Y0)
  naive_l <- min(obs_Y1) - max(obs_Y0)

  # Extend bounds and compute Makarov bounds
  bounds_ext <- lapply(bounds, function(z) c(z, 0, 1))
  diff_bounds <- get_makarov_bounds(bounds_ext, obs_Y1, obs_Y0)

  quants_check <- seq(0.1, 1, by = 0.1)
  quants_true <- quantile(data$tau, quants_check)
  conf_int_idx <- lapply(quants_check, function(z) which(diff_bounds[, 2] >= z))
  conf_int_us <- lapply(conf_int_idx, function(z) max(c(min(diff_bounds[z, 3])), naive_l))
  conf_riqite <- ci_quantile(data$Z, data$Y, quants_check * nrow(data),
                             nperm = 1e3, alpha = params$alpha)

  for (i in seq_along(conf_int_us)) {
    if (abs(conf_int_us[[i]] - naive_l) < 0.25) {
      conf_int_us[[i]] <- -Inf
    }
  }

  cover_us <- quants_true > unlist(conf_int_us)
  cover_riqite <- quants_true > conf_riqite$lower

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
}

# Simulation parameters
NSIM <- 20
m_vals <- c(100, 500, 1000)
tau_vals <- seq(1, 9, length = 3)
rho_vals <- c(0.1, 0.3, 0.7)

# Run simulations
cat("Running Monte Carlo simulations...\n")
sim_results <- data.frame()

for (tau in tau_vals) {
  for (rho in rho_vals) {
    for (size in m_vals) {
      params <- list(m = size, p = 0.5, alpha = 0.1, tau = tau, rho = rho)
      for (sim in 1:NSIM) {
        cat(sprintf("tau=%.1f, rho=%.1f, m=%d, sim=%d\n", tau, rho, size, sim))
        sim_results <- rbind(sim_results, compare_makarov_quantile(params))
      }
    }
  }
}

# Save results
write.csv(sim_results, file = file.path(data.dir, "makarov_quantile_sim_results.csv"),
          row.names = FALSE)

# Aggregate and plot
unique_m_vals <- sort(unique(sim_results$m))
unique_rho_vals <- sort(unique(sim_results$rho))
unique_tau_vals <- sort(unique(sim_results$tau))

avg_m <- sim_results %>%
  group_by(m, quants_check) %>%
  summarise(lower_us = mean(lower_us, na.rm = TRUE),
            lower_riqite = mean(lower_riqite, na.rm = TRUE), .groups = "drop")

avg_rho <- sim_results %>%
  group_by(rho, quants_check) %>%
  summarise(lower_us = mean(lower_us, na.rm = TRUE),
            lower_riqite = mean(lower_riqite, na.rm = TRUE), .groups = "drop")

avg_tau <- sim_results %>%
  group_by(tau, quants_check) %>%
  summarise(lower_us = mean(lower_us, na.rm = TRUE),
            lower_riqite = mean(lower_riqite, na.rm = TRUE), .groups = "drop")

# Generate plots
walk(unique_m_vals, function(m_val) {
  m_result <- filter(avg_m, m == m_val)
  riqite_list <- list(lower = m_result$lower_riqite, k = m_result$quants_check * m_val)
  us_list <- list(lower = m_result$lower_us)
  plot_file <- file.path(plot.dir, sprintf("quantile_CI_comparison_m_%s.png", m_val))
  png(plot_file, width = 2400, height = 1800, units = "px", res = 300, pointsize = 14)
  plot_quantile_ci(riqite_list, k_start = 1, result2 = us_list, xlim = c(0, 6),
                   main = paste0("Quantile CIs (m = ", m_val, ")"))
  dev.off()
  cat("Saved:", plot_file, "\n")
})

walk(unique_rho_vals, function(rho_val) {
  rho_result <- filter(avg_rho, rho == rho_val)
  riqite_list <- list(lower = rho_result$lower_riqite, k = rho_result$quants_check * 100)
  us_list <- list(lower = rho_result$lower_us)
  plot_file <- file.path(plot.dir, sprintf("quantile_CI_comparison_rho_%s.png",
                                           gsub("\\.", "_", as.character(rho_val))))
  png(plot_file, width = 2400, height = 1800, units = "px", res = 300, pointsize = 14)
  plot_quantile_ci(riqite_list, k_start = 1, result2 = us_list, xlim = c(0, 6),
                   main = paste0("Quantile CIs (rho = ", rho_val, ")"))
  dev.off()
  cat("Saved:", plot_file, "\n")
})

walk(unique_tau_vals, function(tau_val) {
  tau_result <- filter(avg_tau, tau == tau_val)
  riqite_list <- list(lower = tau_result$lower_riqite, k = tau_result$quants_check * 100)
  us_list <- list(lower = tau_result$lower_us)
  plot_file <- file.path(plot.dir, sprintf("quantile_CI_comparison_tau_%s.png",
                                           gsub("\\.", "_", as.character(tau_val))))
  png(plot_file, width = 2400, height = 1800, units = "px", res = 300, pointsize = 14)
  plot_quantile_ci(riqite_list, k_start = 1, result2 = us_list,
                   xlim = c(tau_val - 5, tau_val + 3),
                   main = paste0("Quantile CIs (tau = ", tau_val, ")"))
  dev.off()
  cat("Saved:", plot_file, "\n")
})

cat("Done!\n")
