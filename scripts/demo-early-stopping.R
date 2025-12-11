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

early_stopping_ite <- function(params) {
  # Generate data
  
  data <- generate_potential_outcomes(
    m = params$m,
    p = params$p,
    tau = params$tau,
    rho = params$rho
  )

  for(t in seq(20, m, by = 20)){
    sub_data <- data[(1:t), ]
    # Extract observed outcomes by treatment status
    obs_Y1 <- sub_data$Y[sub_data$Z == 1]
    obs_Y0 <- sub_data$Y[sub_data$Z == 0]
    
    # Compute quantile CIs using the new convenience function
    quants_check <- seq(0.1, 1, by = 0.1)

    conf_int_us <- quantile_ci(obs_Y1, obs_Y0,
                               quantiles = quants_check,
                               alpha = params$alpha, time_uniform = TRUE)
    
    
    # True quantiles
    quants_true <- quantile(sub_data$tau, quants_check)
    
    # Coverage
    cover_us <- quants_true > conf_int_us$lower
    print(conf_int_us[5,]$lower)
    print(conf_int_us[8,]$lower)
  }
  
  
  
  
  data.frame(
    m = params$m,
    p = params$p,
    tau = params$tau,
    rho = params$rho,
    lower_us = conf_int_us$lower,
    lower_riqite = conf_riqite$lower,
    cover_us = cover_us,
    cover_riqite = cover_riqite,
    quants_check = quants_check
  )
}

# Simulation parameters
NSIM <- 20
m_vals <- c(500, 1000)
tau_vals <- c(2); #seq(1, 9, length = 3)
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
        sim_results <- rbind(sim_results, early_stopping_ite(params))
      }
    }
  }
}