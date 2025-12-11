#' Compute DKW Confidence Bands for Empirical CDF
#'
#' Computes Dvoretzky-Kiefer-Wolfowitz (DKW) confidence bands for the
#' empirical cumulative distribution function. These bands provide
#' non-parametric confidence sets for the true CDF.
#'
#' @param Fhat Numeric vector. Empirical CDF values at which to compute bounds.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#' @param n Integer. Sample size used to compute the ECDF.
#'
#' @return A list with components:
#'   \describe{
#'     \item{lower}{Numeric vector. Lower confidence band values.}
#'     \item{upper}{Numeric vector. Upper confidence band values.}
#'     \item{epsilon}{Numeric. The DKW bandwidth parameter.}
#'   }
#'
#' @details
#' The DKW inequality (with Massart's tight constant) provides:
#' \deqn{P\left(\sup_x |F_n(x) - F(x)| > \epsilon\right) \leq 2\exp(-2n\epsilon^2)}
#'
#' Inverting gives the bandwidth:
#' \deqn{\epsilon = \sqrt{\frac{\log(2/\alpha)}{2n}}}
#'
#' The confidence band is then:
#' \deqn{[\max(0, \hat{F}(x) - \epsilon), \min(1, \hat{F}(x) + \epsilon)]}
#'
#' @references
#' Dvoretzky, A., Kiefer, J., & Wolfowitz, J. (1956). Asymptotic minimax
#' character of the sample distribution function and of the classical
#' multinomial estimator. The Annals of Mathematical Statistics, 27(3), 642-669.
#'
#' Massart, P. (1990). The tight constant in the Dvoretzky-Kiefer-Wolfowitz
#' inequality. The Annals of Probability, 18(3), 1269-1283.
#'
#' @examples
#' # Generate sample and compute ECDF
#' set.seed(123)
#' x <- rnorm(100)
#' Fhat <- ecdf(x)(sort(x))
#'
#' # Compute 95% DKW bands
#' bands <- dkw_band(Fhat, alpha = 0.05, n = 100)
#'
#' # Plot
#' plot(sort(x), Fhat, type = "s", ylim = c(0, 1))
#' lines(sort(x), bands$lower, col = "red", lty = 2)
#' lines(sort(x), bands$upper, col = "red", lty = 2)
#'
#' @export
dkw_band <- function(Fhat, N, finite_pop = FALSE, alpha = 0.05) {
  if (alpha <= 0 || alpha >= 1) {
    stop("Significance level alpha must be in (0, 1)")
  }
  n <- length(Fhat)
  epsilon <- sqrt(log(2 / alpha) / (2 * n))
  if(finite_pop){
    epsilon <- epsilon * sqrt(1 - (n - 1) / N)
  }

  list(
    lower = pmax(Fhat - epsilon, 0),
    upper = pmin(Fhat + epsilon, 1),
    epsilon = epsilon
  )
}


#' Compute DKW Bounds for Treatment and Control Groups
#'
#' Computes DKW confidence bands separately for treatment and control groups,
#' with optional finite population correction.
#'
#' @param obs_Y1 Numeric vector. Observed outcomes for treatment group.
#' @param obs_Y0 Numeric vector. Observed outcomes for control group.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.1.
#' @param finite_pop Logical. If TRUE, apply finite population correction
#'   to the DKW bandwidth. Default is TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{bounds}{List of four vectors: upper and lower bounds for
#'       treatment group (indices 1-2) and control group (indices 3-4).}
#'     \item{obs_Y1_sorted}{Sorted treatment group outcomes.}
#'     \item{obs_Y0_sorted}{Sorted control group outcomes.}
#'     \item{ecdf_Y1}{ECDF values for treatment group.}
#'     \item{ecdf_Y0}{ECDF values for control group.}
#'     \item{epsilon1}{DKW bandwidth for treatment group.}
#'     \item{epsilon0}{DKW bandwidth for control group.}
#'   }
#'
#' @details
#' When \code{finite_pop = TRUE}, the DKW bandwidth is multiplied by a
#' finite population correction factor:
#' \deqn{\epsilon_{adj} = \epsilon \cdot \sqrt{1 - \frac{n_g - 1}{m}}}
#'
#' where \eqn{n_g} is the group sample size and \eqn{m} is the total sample size.
#'
#' @examples
#' # Simulate treatment/control data
#' set.seed(123)
#' Y1 <- rnorm(50, mean = 1)
#' Y0 <- rnorm(50, mean = 0)
#'
#' # Compute DKW bounds
#' result <- compute_dkw_bounds(Y1, Y0, alpha = 0.1)
#'
#' @export
compute_dkw_bounds <- function(obs_Y1, obs_Y0, alpha = 0.1, finite_pop = FALSE) {
  n1 <- length(obs_Y1)
  n0 <- length(obs_Y0)
  m <- n1 + n0

  # Sort observations
  obs_Y1_sorted <- sort(obs_Y1)
  obs_Y0_sorted <- sort(obs_Y0)

  # Compute ECDFs
  ecdf_fn1 <- stats::ecdf(obs_Y1)
  ecdf_fn0 <- stats::ecdf(obs_Y0)

  ecdf_Y1 <- ecdf_fn1(obs_Y1_sorted)
  ecdf_Y0 <- ecdf_fn0(obs_Y0_sorted)

  
  # Compute DKW bandwidths
  # bound_1 <- dkw_band(ecdf_Y1, N = m, finite_pop = finite_pop, alpha = alpha)
  # bound_0 <- dkw_band(ecdf_Y0, N = m, finite_pop = finite_pop, alpha = alpha)
  # compute stronger hybrid bands
  bound_1 <- hybrid_band(ecdf_Y1, m = m, finite_pop = finite_pop, alpha = alpha)
  bound_0 <- hybrid_band(ecdf_Y0, m = m, finite_pop = finite_pop, alpha = alpha)
  # Compute bounds
  u_bound_1 <- bound_1$upper
  l_bound_1 <- bound_1$lower
  u_bound_0 <- bound_0$upper
  l_bound_0 <- bound_0$lower

  list(
    bounds = list(u_bound_1, l_bound_1, u_bound_0, l_bound_0),
    obs_Y1_sorted = obs_Y1_sorted,
    obs_Y0_sorted = obs_Y0_sorted,
    ecdf_Y1 = ecdf_Y1,
    ecdf_Y0 = ecdf_Y0
  )
}


#' Generate Data and Compute DKW Bounds
#'
#' Convenience function that generates simulated potential outcomes data
#' and computes DKW bounds for treatment and control group ECDFs.
#'
#' @param m Integer. Sample size. Default is 300.
#' @param p Numeric in (0, 1). Treatment proportion. Default is 0.5.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.1.
#' @param rho Numeric in (-1, 1). Correlation between potential outcomes.
#'   Default is 0.4.
#' @param tau Numeric. True average treatment effect. Default is 1.
#' @param finite_pop Logical. Apply finite population correction. Default is TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{bounds}{DKW bounds list (see \code{\link{compute_dkw_bounds}}).}
#'     \item{data}{Generated data frame with potential outcomes.}
#'     \item{obs_Y1}{Treatment group observed outcomes.}
#'     \item{obs_Y0}{Control group observed outcomes.}
#'   }
#'
#' @examples
#' # Generate data and bounds
#' result <- get_bounds(m = 200, p = 0.5, alpha = 0.1, rho = 0.3, tau = 2)
#'
#' # Access bounds
#' str(result$bounds)
#'
#' # True treatment effects
#' summary(result$data$tau)
#'
#' @export
get_bounds <- function(m = 300, p = 0.5, alpha = 0.1,
                       rho = 0.4, tau = 1, finite_pop = FALSE) {
  # Generate data
data <- generate_potential_outcomes(m = m, p = p, tau = tau, rho = rho)

  # Extract observed outcomes by treatment status
  obs_Y0 <- data$Y[data$Z == 0]
  obs_Y1 <- data$Y[data$Z == 1]

  # Compute DKW bounds
  dkw_result <- compute_dkw_bounds(obs_Y1, obs_Y0, alpha = alpha, finite_pop = finite_pop)

  list(
    bounds = dkw_result$bounds,
    data = data,
    obs_Y1 = obs_Y1,
    obs_Y0 = obs_Y0
  )
}
