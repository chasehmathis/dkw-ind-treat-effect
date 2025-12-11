# Quantile-uniform confidence sequence II (Theorem S1)

# Helper functions
logit <- function(p) log(p / (1 - p))
logit_inv <- function(l) exp(l) / (1 + exp(l))

# r_{p,t} from (S8)
r_pt <- function(p, t) {
  if (p >= 0.5) {
    return(p)
  } else {
    return(pmin(0.5, logit_inv(logit(p) + sqrt(2.1 / t))))
  }
}
r_pt <- Vectorize(r_pt)

# ell(p, t) from (S9)
ell_pt <- function(p, t, alpha = 0.05) {
  1.4 * log(log(2.1 * t)) + 1.4 * log(sqrt(t) * abs(logit(p)) + 1) + log(72 / alpha)
}

# g_tilde_t(p) from (S10)
g_tilde <- function(p, t, delta = 0.5, alpha = 0.05) {
  r <- r_pt(p, t)
  ell <- ell_pt(p, t, alpha)
  delta * sqrt(2.1 * t * r * (1 - r)) + 1.5 * sqrt(r * (1 - r) * t * ell) + 0.81 * ell
}

# Confidence band half-widths (asymmetric)
# Lower: g_tilde(1-p) / t
# Upper: g_tilde(p) / t
band_lower <- function(p, t, delta = 0.5, alpha = 0.05) {
  g_tilde(1 - p, t, delta, alpha) / t
}

band_upper <- function(p, t, delta = 0.5, alpha = 0.05) {
  g_tilde(p, t, delta, alpha) / t
}

# Plot with example data
set.seed(42)
n <- 500
X_stream <- rnorm(n, mean = 0, sd = 1)

par(mfrow = c(2, 2))


for (t in c(50, 100, 250, 500)) {
  X_t <- X_stream[1:t]
  F_hat <- ecdf(X_t)
  
  # Evaluate on a fine grid
  x_grid <- seq(min(X_t) - 0.5, max(X_t) + 0.5, length.out = 500)
  p_hat <- F_hat(x_grid)
  true_p <- pnorm(x_grid)
  
  # Bands (clip p_hat away from 0 and 1 for numerical stability)
  p_clipped <- pmax(0.001, pmin(0.999, p_hat))
  lower <- pmax(0, p_hat - band_lower(p_clipped, t))
  upper <- pmin(1, p_hat + band_upper(p_clipped, t))
  
  # Plot
  plot(x_grid, p_hat, type = "l", lwd = 2,
       xlab = "x", ylab = "F(x)",
       main = paste0("Theorem S1: t = ", t),
       ylim = c(0, 1))
  
  # Band as lines
  lines(x_grid, lower, col = "red", lwd = 1.5, lty = 2)
  lines(x_grid, upper, col = "red", lwd = 1.5, lty = 2)
  
  # True CDF
  lines(x_grid, true_p, col = "blue", lwd = 2)
  
  # Check coverage
  inside <- all(true_p >= lower & true_p <= upper)
  legend("bottomright",
         c("Empirical", "Band", "True", paste("Covered:", inside)),
         col = c("black", "red", "blue", NA), lty = c(1, 2, 1, NA), lwd = 2,
         cex = 0.7)
}