plot_quantile_CI_comparison <- function (result, k_start = NULL, main = NULL, 
                                         xlim=NULL, result2=NULL) 
{
  ci.limit = result$lower
  n = length(ci.limit)
  if (is.null(k_start)) {
    k_start = n - sum(!is.nan(ci.limit) & !is.infinite(ci.limit))
  }
  ylim = c(k_start, n + 1)
  if (is.null(xlim)){
    xlim = range(ci.limit[ci.limit > -Inf]) * 1.1
  }
  y_axis_labels <- if (!is.null(result$k)) result$k else seq_along(ci.limit)
  plot(NA, ylab = "k", xlab = expression("lower" ~ "confidence" ~ 
                                           "limit" ~ "for" ~ tau[(k)]), ylim = ylim, xlim = xlim, 
       main = main, yaxt = "n")
  axis(2, at = seq_along(y_axis_labels), labels = y_axis_labels, las = 1)
  for (k in 1:length(ci.limit)) {
    lines(c(max(ci.limit[k], min(ci.limit[ci.limit > -Inf]) - 
                  100), max(ci.limit) + 10), rep(k, 2), col = "grey")
  }
  points(ci.limit, seq_along(ci.limit), pch = 20, col = 1)
  if (!is.null(result2)) {
    ci.limit2 = result2$lower
    points(ci.limit2, seq_along(ci.limit2), pch = 20, col = 2)
    legend("topleft", legend = c("riqite", "us"),
           col = c(1,2), pch = 20, bty = "n")
  }
  abline(v = 0, lty = 2)
  invisible(result)
}
