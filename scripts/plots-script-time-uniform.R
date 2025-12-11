# Generate CDF plots using time-uniform bounds across parameter grid

library(dkw)
library(ggplot2)

set.seed(123)

plot.dir <- "output/plots/"
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

cat("Generating CDF plots with time-uniform bounds...\n")

plot_all_time_uniform <- function(params) {
  result <- get_time_uniform_bounds(
    m = params$m,
    p = params$p,
    alpha = params$alpha,
    tau = params$tau,
    rho = params$rho,
    delta = params$delta
  )

  bounds <- result$bounds
  data <- result$data
  obs_Y1 <- sort(result$obs_Y1)
  obs_Y0 <- sort(result$obs_Y0)

  # Extend bounds and compute Makarov bounds
  bounds_ext <- lapply(bounds, function(z) c(z, 0, 1))
  cat("Computing Makarov bounds with time-uniform CDF bands...\n")
  diff_bounds <- get_makarov_bounds(bounds_ext, obs_Y1, obs_Y0)

  fname <- paste0("cdf-tu-m-", params$m, "-rho-", params$rho,
                  "-p-", params$p, "-tau-", params$tau, ".png")

  plot_title <- bquote("Time-Uniform Bounds: ECDF of individual " * tau ~
    "(" * m == .(params$m) * "," ~ p == .(params$p) * "," ~
    rho == .(params$rho) * "," ~ tau == .(params$tau) * ")")

  # Prepare data for ggplot
  ecdf_vals <- ecdf(data$tau)
  ecdf_data <- data.frame(
    tau = sort(unique(data$tau)),
    ecdf = ecdf_vals(sort(unique(data$tau)))
  )

  bounds_df <- data.frame(
    tau = diff_bounds[, 3],
    lower = diff_bounds[, 1],
    upper = diff_bounds[, 2]
  )

  gg_cdf <- ggplot() +
    geom_step(data = ecdf_data, aes(x = tau, y = ecdf),
              color = "#2C3E50", linewidth = 1.1) +
    geom_point(data = bounds_df, aes(x = tau, y = lower),
               color = "#2980B9", size = 1, alpha = 0.85, shape = 16) +
    geom_point(data = bounds_df, aes(x = tau, y = upper),
               color = "#E74C3C", size = 1, alpha = 0.85, shape = 16) +
    labs(
      title = plot_title,
      x = expression("Treatment Effect (" * tau * ")"),
      y = expression("Empirical Cumulative Distribution Function " * (widehat(F)(tau))),
      caption = "Blue: lower Makarov bound; Red: upper Makarov bound"
    ) +
    theme_grey() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14),
      axis.title = element_text(size = 12)
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(-5, 6))

  ggsave(filename = file.path(plot.dir, fname), plot = gg_cdf,
         width = 7, height = 5, dpi = 300)
  cat("Saved:", fname, "\n")
}

# Generate plots for all parameter combinations
m_vals <- c(100, 500, 1000)
rho_vals <- c(0.1, 0.3, 0.7)

for (rho in rho_vals) {
  for (m in m_vals) {
    params <- list(m = m, p = 0.5, alpha = 0.1, tau = 1, rho = rho, delta = 0.5)
    plot_all_time_uniform(params)
  }
}


cat("Done!\n")
