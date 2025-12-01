# Illustrating conservativeness of finite population inference with worst-case bounds:

library(DeclareDesign)
library(MASS)
library(ggplot2)

set.seed(123)
m <- 50
tau <- 1
rho <- -0.4
p <- 0.5
NSIM <- 500

# --- Design for super-population draw ---
design <- 
  declare_model(
    N = m,
    Z = rbinom(n = N, size = 1, prob = p),
    draw_multivariate(
      c(Y_Z_0, Y_Z_1) ~ mvrnorm(
        n = m, 
        mu = c(0, tau), 
        Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
      )
    )
  ) + 
  declare_measurement(Y = reveal_outcomes(Y ~ Z))

pop_dat <- draw_data(design)

# --- Simulation: collect ECDFs for both finite- and super-pop draws ---
finite_ecdfs <- vector("list", NSIM)
pop_draw_ecdfs <- vector("list", NSIM)
all_control_vals <- numeric()
all_popdraw_vals <- numeric()

for (s in seq_len(NSIM)) {
  # (A) Finite population randomization inference: re-randomize Z
  data_sim <- pop_dat
  data_sim$Z <- rbinom(m, 1, p)
  obs_control <- data_sim$Y_Z_0[data_sim$Z == 0]
  finite_ecdfs[[s]] <- ecdf(obs_control)
  all_control_vals <- c(all_control_vals, obs_control)
  
  # (B) Super-population draw + randomization
  draw_dat <- draw_data(design)
  obs_control_popdraw <- draw_dat$Y[draw_dat$Z == 0]
  pop_draw_ecdfs[[s]] <- ecdf(obs_control_popdraw)
  all_popdraw_vals <- c(all_popdraw_vals, obs_control_popdraw)
}

# --- Evaluation grid for plotting CDFs ---
all_vals <- c(all_control_vals, all_popdraw_vals)
min_val <- min(all_vals, na.rm = TRUE)
max_val <- max(all_vals, na.rm = TRUE)
eval_grid <- sort(unique(all_vals))
if (length(eval_grid) < 25) {
  eval_grid <- seq(min_val, max_val, length.out = m)
}

# --- Create ECDF Matrices (faster than per-simulation data frames) ---
finite_cdfs_matrix <- sapply(finite_ecdfs, function(f) f(eval_grid))
popdraw_cdfs_matrix <- sapply(pop_draw_ecdfs, function(f) f(eval_grid))

# --- Prepare tidy data frames for ggplot (vectorized, more efficient) ---
to_long_df <- function(cdf_matrix, sim_type) {
  data.frame(
    value = rep(eval_grid, NSIM),
    cdf = as.vector(cdf_matrix),
    sim = rep(seq_len(NSIM), each = length(eval_grid)),
    type = sim_type
  )
}
plot_df <- to_long_df(finite_cdfs_matrix, "Finite-pop assignment")
plot_df_popdraw <- to_long_df(popdraw_cdfs_matrix, "Super-population")
plot_df_combined <- rbind(plot_df, plot_df_popdraw)

# --- Calculate worst-case (envelope) bounds for each method ---
worst_case_bounds_df <- data.frame(
  value = eval_grid,
  finite_lower = apply(finite_cdfs_matrix, 1, min),
  finite_upper = apply(finite_cdfs_matrix, 1, max),
  popdraw_lower = apply(popdraw_cdfs_matrix, 1, min),
  popdraw_upper = apply(popdraw_cdfs_matrix, 1, max)
)

# --- Plot ECDFs and envelopes ---
p <- ggplot() +
  # Finite population worst-case bounds: dotted black
  geom_line(
    data = worst_case_bounds_df, aes(x = value, y = finite_lower),
    color = "#111111", linetype = "dotted", linewidth = 1.1
  ) +
  geom_line(
    data = worst_case_bounds_df, aes(x = value, y = finite_upper),
    color = "#111111", linetype = "dotted", linewidth = 1.1
  ) +
  # Super-population worst-case bounds: solid blue
  geom_line(
    data = worst_case_bounds_df, aes(x = value, y = popdraw_lower),
    color = "#1362b1", linetype = "solid", linewidth = 1.1
  ) +
  geom_line(
    data = worst_case_bounds_df, aes(x = value, y = popdraw_upper),
    color = "#1362b1", linetype = "solid", linewidth = 1.1
  ) +
  labs(
    title = "Empirical CDFs with Worst-case Bounds: Finite vs. Super-population",
    subtitle = "Finite-pop assignment bands: dotted black; Super-population bands: solid blue",
    x = "Value", y = "Cumulative Probability"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, family = "serif"),
    plot.subtitle = element_text(size = 12, family = "serif"),
  panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.caption = element_text(size = 10, color = "gray50")
)

plot.dir <- "../fig/"
if (!dir.exists(plot.dir)) {
  dir.create(plot.dir, recursive = TRUE)
}
ggsave(
  filename = file.path(plot.dir, "finite-vs-superpop-worstcase-bands.png"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 320
)
