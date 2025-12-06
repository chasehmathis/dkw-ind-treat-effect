#' Plot Quantile Confidence Interval Comparison
#'
#' Creates a plot comparing lower confidence limits for treatment effect
#' quantiles, optionally overlaying two methods for comparison.
#'
#' @param result A list with at least a \code{lower} component containing
#'   the lower confidence limits. Optionally includes a \code{k} component
#'   for custom y-axis labels.
#' @param k_start Integer or NULL. Starting rank for display. If NULL,
#'   automatically determined from non-missing values. Default is NULL.
#' @param main Character or NULL. Plot title. Default is NULL.
#' @param xlim Numeric vector of length 2 or NULL. X-axis limits.
#'   Default is NULL (auto-determined).
#' @param result2 A second result list for comparison, or NULL. Default is NULL.
#'
#' @return Invisibly returns the input \code{result}.
#'
#' @details
#' This function visualizes confidence intervals for ordered treatment effects
#' (rank statistics). The x-axis shows the lower confidence limit for
#' \eqn{\tau_{(k)}}, and the y-axis shows the rank k.
#'
#' When \code{result2} is provided, both sets of confidence limits are
#' plotted, with the first in black and the second in red.
#'
#' @examples
#' # Create mock results
#' result1 <- list(lower = c(-Inf, -2, -1, 0, 0.5, 1, 1.5, 2, 2.5, 3),
#'                 k = 1:10)
#' result2 <- list(lower = c(-Inf, -1.5, -0.5, 0.2, 0.8, 1.2, 1.8, 2.2, 2.7, 3.2))
#'
#' plot_quantile_ci(result1, main = "Method Comparison", result2 = result2)
#'
#' @export
plot_quantile_ci <- function(result, k_start = NULL, main = NULL,
                             xlim = NULL, result2 = NULL) {
  ci_limit <- result$lower
  n <- length(ci_limit)

  # Determine starting point
  if (is.null(k_start)) {
    k_start <- n - sum(!is.nan(ci_limit) & !is.infinite(ci_limit))
  }

  ylim <- c(k_start, n + 1)

  # Determine x-axis limits
  if (is.null(xlim)) {
    finite_vals <- ci_limit[is.finite(ci_limit)]
    if (length(finite_vals) > 0) {
      xlim <- range(finite_vals) * 1.1
    } else {
      xlim <- c(-1, 1)
    }
  }

  # Y-axis labels
  y_axis_labels <- if (!is.null(result$k)) result$k else seq_along(ci_limit)

  # Create plot
  graphics::plot(NA,
       ylab = "k",
       xlab = expression("lower" ~ "confidence" ~ "limit" ~ "for" ~ tau[(k)]),
       ylim = ylim,
       xlim = xlim,
       main = main,
       yaxt = "n")

  graphics::axis(2, at = seq_along(y_axis_labels), labels = y_axis_labels, las = 1)

  # Draw horizontal lines
  min_finite <- if (any(ci_limit > -Inf)) min(ci_limit[ci_limit > -Inf]) - 100 else -100
  max_ci <- if (any(is.finite(ci_limit))) max(ci_limit[is.finite(ci_limit)]) + 10 else 10

  for (k in seq_along(ci_limit)) {
    graphics::lines(c(max(ci_limit[k], min_finite), max_ci),
          rep(k, 2), col = "grey")
  }

  # Plot points for first result
  graphics::points(ci_limit, seq_along(ci_limit), pch = 20, col = 1)

  # Plot second result if provided
  # TOOD What if different k
  if (!is.null(result2)) {
    ci_limit2 <- result2$lower
    graphics::points(ci_limit2, seq_along(ci_limit2), pch = 20, col = 2)
    graphics::legend("topleft", legend = c("RIQITE", "Makarov"),
           col = c(1, 2), pch = 20, bty = "n")
  }

  # Add reference line at zero
  graphics::abline(v = 0, lty = 2)

  invisible(result)
}


#' Plot Makarov Bounds on Treatment Effect CDF
#'
#' Creates a plot showing the Makarov confidence bounds on the cumulative
#' distribution function of individual treatment effects.
#'
#' @param bounds_matrix Matrix from \code{\link{get_makarov_bounds}} or
#'   \code{\link{makarov_bounds}}.
#' @param true_tau Numeric vector or NULL. True individual treatment effects
#'   for plotting the empirical CDF. Default is NULL.
#' @param main Character. Plot title. Default is "Makarov Bounds on Treatment Effect CDF".
#' @param xlab Character. X-axis label. Default is "Treatment Effect (tau)".
#' @param ylab Character. Y-axis label. Default is "Cumulative Probability".
#' @param lower_col Color for lower bounds. Default is "blue".
#' @param upper_col Color for upper bounds. Default is "red".
#' @param ... Additional arguments passed to \code{plot}.
#'
#' @return Invisibly returns the bounds matrix.
#'
#' @examples
#' set.seed(123)
#' result <- get_bounds(m = 200, p = 0.5, alpha = 0.1, rho = 0.3, tau = 2)
#' bounds_ext <- lapply(result$bounds, function(z) c(z, 0, 1))
#' makarov <- get_makarov_bounds(bounds_ext, sort(result$obs_Y1), sort(result$obs_Y0))
#'
#' plot_makarov_bounds(makarov, true_tau = result$data$tau)
#'
#' @export
plot_makarov_bounds <- function(bounds_matrix, true_tau = NULL,
                                main = "Makarov Bounds on Treatment Effect CDF",
                                xlab = "Treatment Effect (tau)",
                                ylab = "Cumulative Probability",
                                lower_col = "blue", upper_col = "red", ...) {
  # Handle list input
  if (is.list(bounds_matrix) && "bounds_matrix" %in% names(bounds_matrix)) {
    bounds_matrix <- bounds_matrix$bounds_matrix
  }

  # Plot true ECDF if available
  if (!is.null(true_tau)) {
    graphics::plot(stats::ecdf(true_tau), main = main, xlab = xlab, ylab = ylab, ...)
  } else {
    t_range <- range(bounds_matrix[, "t"])
    graphics::plot(NA, xlim = t_range, ylim = c(0, 1),
         main = main, xlab = xlab, ylab = ylab, ...)
  }

  # Add bounds
  graphics::points(bounds_matrix[, "t"], bounds_matrix[, "lower"],
         col = lower_col, pch = 16, cex = 0.3)
  graphics::points(bounds_matrix[, "t"], bounds_matrix[, "upper"],
         col = upper_col, pch = 16, cex = 0.3)

  # Add legend
  graphics::legend("bottomright",
         legend = c("Lower bound", "Upper bound"),
         col = c(lower_col, upper_col),
         pch = 16,
         bty = "n")

  invisible(bounds_matrix)
}


#' Plot DKW Confidence Bands
#'
#' Creates a plot showing DKW confidence bands around an empirical CDF.
#'
#' @param x Numeric vector. Sorted sample values.
#' @param Fhat Numeric vector. ECDF values at x.
#' @param lower Numeric vector. Lower confidence band values.
#' @param upper Numeric vector. Upper confidence band values.
#' @param main Character. Plot title. Default is "DKW Confidence Bands".
#' @param xlab Character. X-axis label. Default is "x".
#' @param ylab Character. Y-axis label. Default is "F(x)".
#' @param band_col Color for confidence bands. Default is "red".
#' @param ... Additional arguments passed to \code{plot}.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' set.seed(123)
#' x <- sort(rnorm(100))
#' Fhat <- ecdf(x)(x)
#' bands <- dkw_band(Fhat, alpha = 0.05, n = 100)
#'
#' plot_dkw_bands(x, Fhat, bands$lower, bands$upper)
#'
#' @export
plot_dkw_bands <- function(x, Fhat, lower, upper,
                           main = "DKW Confidence Bands",
                           xlab = "x", ylab = "F(x)",
                           band_col = "red", ...) {
  graphics::plot(x, Fhat, type = "s", lwd = 2, ylim = c(0, 1),
       main = main, xlab = xlab, ylab = ylab, ...)
  graphics::lines(x, lower, col = band_col, lwd = 2, lty = 2)
  graphics::lines(x, upper, col = band_col, lwd = 2, lty = 2)

  graphics::legend("bottomright",
         legend = c("Empirical CDF", "DKW bands"),
         col = c("black", band_col),
         lty = c(1, 2),
         lwd = 2,
         bty = "n")

  invisible(NULL)
}
