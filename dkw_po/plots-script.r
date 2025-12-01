#plots quad box
set.seed(123)
setwd("~/Research/dkw/dkw_po")
source("generate.R")
source("bounds.R")
plot.dir <- "../fig/"
# Create fig directory if it doesn't exist
if (!dir.exists(plot.dir)) {
  dir.create(plot.dir, recursive = TRUE)
}
# Generate data and compute initial bounds
cat("Generating data and computing initial bounds...\n")
# high level params for stratifying results
# Define parameters as a list in R

plot_all <- function(params){

  result <- do.call(get_bounds, params)
 
  # Extract results
  bounds <- result$bounds
  data <- result$data
  obs_Y1 <- sort(result$obs_Y1)
  obs_Y0 <- sort(result$obs_Y0)

  ATE <- mean(obs_Y1) - mean(obs_Y0)
  # Compute Makarov bounds for all t values
  cat("Computing Makarov bounds for all quantiles...\n")
  
  diff_bounds <- get_marakov_bounds_all_t(bounds, obs_Y1, obs_Y0, data)

  library(ggplot2)

  fname <- paste0("cdf-m-", params$m, "-rho-", 
                  params$rho, "-p-", params$p, "-tau-", params$tau, ".png")
  # Dynamic title using current parameters
  plot_title <- bquote("True ECDF of individual " * tau ~ 
      "(" * m == .(params$m) * "," ~ p == .(params$p) * "," 
      ~ rho == .(params$rho) * "," ~ tau == .(params$tau) * ")")

  # Prepare ECDF data for ggplot
  ecdf_vals <- ecdf(data$tau)
  ecdf_data <- data.frame(
    tau = sort(unique(data$tau)),
    ecdf = ecdf_vals(sort(unique(data$tau)))
  )

  # Prepare bounds for ggplot
  # diff_bounds columns: lower, upper, tau_grid (assumed from code)
  bounds_df <- data.frame(
    tau = diff_bounds[,3],
    lower = diff_bounds[,1],
    upper = diff_bounds[,2]
  )

  gg_cdf <- ggplot() +
    geom_step(data = ecdf_data, aes(x = tau, y = ecdf), color = "#2C3E50", size = 1.1) +
    geom_point(data = bounds_df, aes(x = tau, y = lower), 
               color = "#2980B9", size = 1, alpha = 0.85, shape = 16) +
    geom_point(data = bounds_df, aes(x = tau, y = upper), 
               color = "#E74C3C", size = 1, alpha = 0.85, shape = 16) +
    labs(
      title = plot_title,
      x = expression("Treatment Effect (" * tau * ")"),
      y = expression("Empirical Cumulative Distribution Function"~(widehat(F)(tau))),
      caption = "Blue: lower Makarov bound; Red: upper Makarov bound"
    ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 16, family = "serif", face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          plot.caption = element_text(size = 10, color = "gray50")
    ) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) + 
    scale_x_continuous(limits = c(min(ecdf_data$tau) - 1, max(ecdf_data$tau) + 1))

  ggsave(
    filename = paste0(plot.dir, fname),
    plot = gg_cdf,
    width = 7,
    height = 5,
    dpi = 300
  )

  cat("Saved plot to:", paste0(plot.dir, fname), "\n")

  
}



# FIGURE 1
m <- c(1e2, 5e2, 1e3)
for(tau in seq(1,9, length = 3)){
    for(rho in c(0.1, 0.3, 0.7)){
      for(size in m){
        
        params <- list(
          m = size,
          p = p,
          alpha = 0.1,
          tau =tau,
          rho = rho
        )
        plot_all(params)
        
      }
  }
}






