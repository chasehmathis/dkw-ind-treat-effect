#' dkw: DKW and Makarov Bounds for Individual Treatment Effect Distributions
#'
#' The dkw package provides tools for non-parametric inference on individual
#' treatment effects in randomized experiments. It implements Dvoretzky-Kiefer-
#' Wolfowitz (DKW) confidence bands for empirical cumulative distribution
#' functions and Makarov bounds for treatment effect distributions.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{makarov_bounds}}}{Compute Makarov bounds on the treatment
#'     effect distribution from raw treatment and control outcomes.}
#'   \item{\code{\link{get_makarov_bounds}}}{Compute Makarov bounds from
#'     pre-computed DKW bounds.}
#'   \item{\code{\link{extract_quantile_ci}}}{Extract quantile confidence
#'     intervals from Makarov bounds.}
#'   \item{\code{\link{dkw_band}}}{Compute standard DKW confidence bands.}
#'   \item{\code{\link{hybrid_band}}}{Compute hybrid Simes/Dempster/DKW bands.}
#' }
#'
#' @section Data Generation:
#' \describe{
#'   \item{\code{\link{generate_potential_outcomes}}}{Simulate treatment/control
#'     data with bivariate normal potential outcomes.}
#'   \item{\code{\link{get_bounds}}}{Generate data and compute DKW bounds
#'     in one step.}
#' }
#'
#' @section Visualization:
#' \describe{
#'   \item{\code{\link{plot_makarov_bounds}}}{Plot Makarov bounds on
#'     treatment effect CDF.}
#'   \item{\code{\link{plot_quantile_ci}}}{Plot quantile confidence interval
#'     comparison.}
#'   \item{\code{\link{plot_dkw_bands}}}{Plot DKW confidence bands.}
#' }
#'
#' @section Background:
#' The DKW inequality provides non-parametric confidence bands for the
#' empirical CDF that are valid uniformly over all points. The Makarov
#' bounds combine these marginal CDF bounds to provide valid confidence
#' sets for the distribution of treatment effects \eqn{\tau_i = Y_i(1) - Y_i(0)}.
#'
#' The key insight is that while individual treatment effects are never
#' observed (due to the fundamental problem of causal inference), bounds
#' on their distribution can be derived from the marginal distributions
#' of observed outcomes.
#'
#' @references
#' Dvoretzky, A., Kiefer, J., & Wolfowitz, J. (1956). Asymptotic minimax
#' character of the sample distribution function and of the classical
#' multinomial estimator. The Annals of Mathematical Statistics, 27(3), 642-669.
#'
#' Makarov, G. D. (1982). Estimates for the distribution function of the sum
#' of two random variables with given marginal distributions. Theory of
#' Probability & Its Applications, 26(4), 803-806.
#'
#' Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz
#' inequality. The Annals of Probability, 18(3), 1269-1283.
#'
#' @docType package
#' @name dkw-package
#' @aliases dkw
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats ecdf quantile uniroot
#' @importFrom graphics abline axis legend lines plot points
## usethis namespace: end
NULL
