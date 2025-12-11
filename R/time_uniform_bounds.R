#' Compute r_{p,t} for Time-Uniform DKW Bounds
#'
#' Computes the adjusted probability value r_{p,t} used in the time-uniform
#' confidence sequence construction (Equation S8 from the supplement).
#'
#' @param p Numeric. Probability value in (0, 1).
#' @param t Integer. Current sample size (time point).
#'
#' @return Numeric. The adjusted probability r_{p,t}.
#'
#' @details
#' For p >= 0.5, r_{p,t} = p. For p < 0.5:
#' \deqn{r_{p,t} = \min\left(0.5, \text{logit}^{-1}\left(\text{logit}(p) + \sqrt{\frac{2.1}{t}}\right)\right)}
#'
#' This adjustment ensures the confidence sequence remains valid uniformly
#' over time while maintaining good coverage properties.
#'
#' Uses R's built-in \code{qlogis} (logit) and \code{plogis} (inverse logit) functions.
#'
#' @examples
#' # For p >= 0.5, returns p unchanged
#' r_pt(0.7, 100)
#'
#' # For p < 0.5, returns adjusted value
#' r_pt(0.3, 100)
#'
#' @seealso \code{\link{time_uniform_band}} for the main confidence band function
#' @keywords internal
#' @export
r_pt <- function(p, t) {
  if (p >= 0.5) {
    return(p)
  } else {
    return(pmin(0.5, plogis(qlogis(p) + sqrt(2.1 / t))))
  }
}

# Vectorize r_pt for use with vectors
r_pt <- Vectorize(r_pt)


#' Compute ell(p, t) for Time-Uniform DKW Bounds
#'
#' Computes the logarithmic term ell(p, t) used in the time-uniform confidence
#' sequence construction (Equation S9 from the supplement).
#'
#' @param p Numeric. Probability value in (0, 1).
#' @param t Integer. Current sample size (time point).
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#'
#' @return Numeric. The ell(p, t) value.
#'
#' @details
#' The ell function captures the time and probability-dependent penalty:
#' \deqn{\ell(p, t) = 1.4 \log(\log(2.1 \cdot t)) + 1.4 \log\left(\sqrt{t} \cdot |\text{logit}(p)| + 1\right) + \log\left(\frac{72}{\alpha}\right)}
#'
#' Uses R's built-in \code{qlogis} function for the logit transformation.
#'
#' @examples
#' ell_pt(0.5, 100, alpha = 0.05)
#' ell_pt(0.25, 100, alpha = 0.05)
#'
#' @seealso \code{\link{time_uniform_band}} for the main confidence band function
#' @keywords internal
#' @export
ell_pt <- function(p, t, alpha = 0.05) {
  1.4 * log(log(2.1 * t)) + 1.4 * log(sqrt(t) * abs(qlogis(p)) + 1) + log(72 / alpha)
}


#' Compute g_tilde(p, t) for Time-Uniform DKW Bounds
#'
#' Computes the g_tilde function used in constructing time-uniform confidence
#' sequences for the empirical CDF (Equation S10 from the supplement).
#'
#' @param p Numeric. Probability value in (0, 1).
#' @param t Integer. Current sample size (time point).
#' @param delta Numeric. Tuning parameter controlling the width-time tradeoff.
#'   Default is 0.5.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#'
#' @return Numeric. The g_tilde(p, t) value.
#'
#' @details
#' The g_tilde function combines three terms to achieve time-uniform validity:
#' \deqn{\tilde{g}_t(p) = \delta\sqrt{2.1 \cdot t \cdot r(1-r)} + 1.5\sqrt{r(1-r) \cdot t \cdot \ell} + 0.81 \cdot \ell}
#'
#' where r = r_{p,t} and \eqn{\ell} = ell(p, t).
#'
#' @examples
#' g_tilde(0.5, 100, delta = 0.5, alpha = 0.05)
#'
#' @seealso \code{\link{time_uniform_band}} for the main confidence band function
#' @keywords internal
#' @export
g_tilde <- function(p, t, delta = 0.5, alpha = 0.05) {
  r <- r_pt(p, t)
  ell <- ell_pt(p, t, alpha)
  delta * sqrt(2.1 * t * r * (1 - r)) + 1.5 * sqrt(r * (1 - r) * t * ell) + 0.81 * ell
}


#' Compute Lower Band Half-Width for Time-Uniform DKW
#'
#' Computes the lower confidence band half-width for the time-uniform
#' confidence sequence.
#'
#' @param p Numeric vector. Probability values (typically from ECDF).
#' @param t Integer. Current sample size (time point).
#' @param delta Numeric. Tuning parameter. Default is 0.5.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#'
#' @return Numeric vector. Lower band half-widths (to subtract from ECDF).
#'
#' @details
#' The lower band half-width is computed as:
#' \deqn{\text{lower width} = \frac{\tilde{g}_t(1-p)}{t}}
#'
#' The asymmetry arises from the asymmetric treatment of probabilities
#' near 0 and 1 in the r_{p,t} function.
#'
#' @examples
#' band_lower(0.5, 100)
#' band_lower(c(0.25, 0.5, 0.75), 100)
#'
#' @seealso \code{\link{band_upper}}, \code{\link{time_uniform_band}}
#' @keywords internal
#' @export
band_lower <- function(p, t, delta = 0.5, alpha = 0.05) {
  g_tilde(1 - p, t, delta, alpha) / t
}


#' Compute Upper Band Half-Width for Time-Uniform DKW
#'
#' Computes the upper confidence band half-width for the time-uniform
#' confidence sequence.
#'
#' @param p Numeric vector. Probability values (typically from ECDF).
#' @param t Integer. Current sample size (time point).
#' @param delta Numeric. Tuning parameter. Default is 0.5.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#'
#' @return Numeric vector. Upper band half-widths (to add to ECDF).
#'
#' @details
#' The upper band half-width is computed as:
#' \deqn{\text{upper width} = \frac{\tilde{g}_t(p)}{t}}
#'
#' @examples
#' band_upper(0.5, 100)
#' band_upper(c(0.25, 0.5, 0.75), 100)
#'
#' @seealso \code{\link{band_lower}}, \code{\link{time_uniform_band}}
#' @keywords internal
#' @export
band_upper <- function(p, t, delta = 0.5, alpha = 0.05) {
  g_tilde(p, t, delta, alpha) / t
}


#' Time-Uniform Confidence Bands for Empirical CDF
#'
#' Computes time-uniform (anytime-valid) confidence bands for the empirical
#' cumulative distribution function. These bands remain valid uniformly over
#' all sample sizes, enabling sequential monitoring without alpha-spending.
#'
#' @param Fhat Numeric vector. Empirical CDF values at which to compute bounds.
#' @param t Integer. Current sample size (time point). If NULL, uses length(Fhat).
#' @param delta Numeric. Tuning parameter controlling early vs. late stopping
#'   power tradeoff. Smaller values give tighter bands at large t but wider
#'   bands at small t. Default is 0.5.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#'
#' @return A list with components:
#'   \describe{
#'     \item{lower}{Numeric vector. Lower confidence band values, clipped to [0,1].}
#'     \item{upper}{Numeric vector. Upper confidence band values, clipped to [0,1].}
#'     \item{lower_width}{Numeric vector. Raw lower half-widths before clipping.}
#'     \item{upper_width}{Numeric vector. Raw upper half-widths before clipping.}
#'     \item{t}{Integer. Sample size used.}
#'     \item{delta}{Numeric. Delta parameter used.}
#'     \item{alpha}{Numeric. Significance level used.}
#'   }
#'
#' @details
#' This implements the quantile-uniform confidence sequence from Theorem S1
#' of the supplement. Unlike standard DKW bands which are only valid at a
#' fixed sample size, these bands maintain (1-alpha) coverage probability
#' uniformly over all sample sizes t = 1, 2, 3, ...
#'
#' The bands are asymmetric: the lower half-width uses g_tilde(1-p)/t while
#' the upper half-width uses g_tilde(p)/t. This asymmetry provides better
#' coverage properties near the boundaries (p near 0 or 1).
#'
#' Key properties:
#' \itemize{
#'   \item \strong{Time-uniform validity}: For any stopping rule,
#'     \eqn{P(\exists t: F(x) \notin [L_t(x), U_t(x)]) \leq \alpha}
#'   \item \strong{Asymptotic tightness}: Width shrinks at rate O(1/sqrt(t))
#'   \item \strong{Continuous monitoring}: Can peek at data at any time
#' }
#'
#' @references
#' Howard, S. R., Ramdas, A., McAuliffe, J., & Sekhon, J. (2021).
#' Time-uniform, nonparametric, nonasymptotic confidence sequences.
#' The Annals of Statistics, 49(2), 1055-1080.
#'
#' @examples
#' # Generate sample and compute ECDF
#' set.seed(123)
#' x <- rnorm(100)
#' Fhat <- ecdf(x)(sort(x))
#'
#' # Compute time-uniform bands
#' bands <- time_uniform_band(Fhat, t = 100, alpha = 0.05)
#'
#' # Plot
#' plot(sort(x), Fhat, type = "s", ylim = c(0, 1),
#'      main = "Time-Uniform Confidence Band")
#' lines(sort(x), bands$lower, col = "red", lty = 2)
#' lines(sort(x), bands$upper, col = "red", lty = 2)
#' lines(sort(x), pnorm(sort(x)), col = "blue")
#' legend("bottomright", c("ECDF", "Band", "True CDF"),
#'        col = c("black", "red", "blue"), lty = c(1, 2, 1))
#'
#' @seealso \code{\link{dkw_band}} for fixed-sample DKW bands,
#'   \code{\link{plot_time_uniform_bands}} for visualization
#' @export
time_uniform_band <- function(Fhat, t = NULL, delta = 0.5, alpha = 0.05) {
  if (is.null(t)) {
    t <- length(Fhat)
  }

  if (alpha <= 0 || alpha >= 1) {
    stop("Significance level alpha must be in (0, 1)")
  }

  if (t < 1) {
    stop("Sample size t must be at least 1")
  }

  # Clip Fhat away from 0 and 1 for numerical stability
  p_clipped <- pmax(0.001, pmin(0.999, Fhat))

  # Compute half-widths
  lower_width <- band_lower(p_clipped, t, delta, alpha)
  upper_width <- band_upper(p_clipped, t, delta, alpha)

  list(
    lower = pmax(0, Fhat - lower_width),
    upper = pmin(1, Fhat + upper_width),
    lower_width = lower_width,
    upper_width = upper_width,
    t = t,
    delta = delta,
    alpha = alpha
  )
}


#' Compute Time-Uniform DKW Bounds for Treatment and Control Groups
#'
#' Computes time-uniform confidence bands separately for treatment and control
#' groups, providing anytime-valid inference for sequential experiments.
#'
#' @param obs_Y1 Numeric vector. Observed outcomes for treatment group.
#' @param obs_Y0 Numeric vector. Observed outcomes for control group.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.1.
#' @param delta Numeric. Tuning parameter for time-uniform bands. Default is 0.5.
#'
#' @return A list with components:
#'   \describe{
#'     \item{bounds}{List of four vectors: upper and lower bounds for
#'       treatment group (indices 1-2) and control group (indices 3-4).}
#'     \item{obs_Y1_sorted}{Sorted treatment group outcomes.}
#'     \item{obs_Y0_sorted}{Sorted control group outcomes.}
#'     \item{ecdf_Y1}{ECDF values for treatment group.}
#'     \item{ecdf_Y0}{ECDF values for control group.}
#'   }
#'
#' @details
#' Unlike standard DKW bounds which are only valid at a fixed sample size,
#' time-uniform bounds maintain (1-alpha) coverage probability uniformly
#' over all sample sizes. This enables valid sequential monitoring without
#' alpha-spending adjustments.
#'
#' The bounds use the construction from Howard et al. (2021) based on
#' sub-exponential confidence sequences.
#'
#' @examples
#' # Simulate treatment/control data
#' set.seed(123)
#' Y1 <- rnorm(50, mean = 1)
#' Y0 <- rnorm(50, mean = 0)
#'
#' # Compute time-uniform bounds
#' result <- compute_time_uniform_bounds(Y1, Y0, alpha = 0.1)
#'
#' @seealso \code{\link{compute_dkw_bounds}} for fixed-sample bounds,
#'   \code{\link{time_uniform_band}} for the underlying band computation
#' @export
compute_time_uniform_bounds <- function(obs_Y1, obs_Y0, alpha = 0.1, delta = 0.5) {
  n1 <- length(obs_Y1)
  n0 <- length(obs_Y0)

  # Sort observations
  obs_Y1_sorted <- sort(obs_Y1)
  obs_Y0_sorted <- sort(obs_Y0)

  # Compute ECDFs
  ecdf_fn1 <- stats::ecdf(obs_Y1)
  ecdf_fn0 <- stats::ecdf(obs_Y0)

  ecdf_Y1 <- ecdf_fn1(obs_Y1_sorted)
  ecdf_Y0 <- ecdf_fn0(obs_Y0_sorted)

  # Compute time-uniform bands
  bound_1 <- time_uniform_band(ecdf_Y1, t = n1, delta = delta, alpha = alpha)
  bound_0 <- time_uniform_band(ecdf_Y0, t = n0, delta = delta, alpha = alpha)

  # Extract bounds in same format as compute_dkw_bounds
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


#' Generate Data and Compute Time-Uniform Bounds
#'
#' Convenience function that generates simulated potential outcomes data
#' and computes time-uniform bounds for treatment and control group ECDFs.
#'
#' @param m Integer. Sample size. Default is 300.
#' @param p Numeric in (0, 1). Treatment proportion. Default is 0.5.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.1.
#' @param rho Numeric in (-1, 1). Correlation between potential outcomes.
#'   Default is 0.4.
#' @param tau Numeric. True average treatment effect. Default is 1.
#' @param delta Numeric. Tuning parameter for time-uniform bands. Default is 0.5.
#'
#' @return A list with components:
#'   \describe{
#'     \item{bounds}{Time-uniform bounds list (see \code{\link{compute_time_uniform_bounds}}).}
#'     \item{data}{Generated data frame with potential outcomes.}
#'     \item{obs_Y1}{Treatment group observed outcomes.}
#'     \item{obs_Y0}{Control group observed outcomes.}
#'   }
#'
#' @examples
#' # Generate data and bounds
#' result <- get_time_uniform_bounds(m = 200, p = 0.5, alpha = 0.1, rho = 0.3, tau = 2)
#'
#' # Access bounds
#' str(result$bounds)
#'
#' @seealso \code{\link{get_bounds}} for fixed-sample DKW bounds
#' @export
get_time_uniform_bounds <- function(m = 300, p = 0.5, alpha = 0.1,
                                     rho = 0.4, tau = 1, delta = 0.5) {
  # Generate data
  data <- generate_potential_outcomes(m = m, p = p, tau = tau, rho = rho)

  # Extract observed outcomes by treatment status
  obs_Y0 <- data$Y[data$Z == 0]
  obs_Y1 <- data$Y[data$Z == 1]

  # Compute time-uniform bounds
  tu_result <- compute_time_uniform_bounds(obs_Y1, obs_Y0, alpha = alpha, delta = delta)

  list(
    bounds = tu_result$bounds,
    data = data,
    obs_Y1 = obs_Y1,
    obs_Y0 = obs_Y0
  )
}


#' Compute Makarov Bounds Using Time-Uniform CDF Bands
#'
#' Computes Makarov bounds on the treatment effect distribution using
#' time-uniform confidence bands instead of standard DKW bands.
#'
#' @param Y1 Numeric vector. Observed outcomes for treatment group.
#' @param Y0 Numeric vector. Observed outcomes for control group.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.1.
#' @param delta Numeric. Tuning parameter for time-uniform bands. Default is 0.5.
#' @param t_range Numeric vector of length 2. Range of treatment effect values.
#'   Default is NULL (auto-determined from data).
#' @param n_points Integer. Number of evaluation points. Default is 1000.
#'
#' @return A list with components:
#'   \describe{
#'     \item{bounds_matrix}{Matrix of Makarov bounds with columns (lower, upper, t).}
#'     \item{tu_bounds}{Time-uniform bounds for marginal ECDFs.}
#'     \item{alpha}{Significance level used.}
#'     \item{delta}{Delta parameter used.}
#'   }
#'
#' @details
#' This function combines time-uniform confidence bands with Makarov's
#' inequality to provide anytime-valid bounds on the treatment effect
#' distribution. Unlike standard Makarov bounds based on DKW, these bounds
#' remain valid under continuous monitoring.
#'
#' @examples
#' set.seed(42)
#' Y1 <- rnorm(100, mean = 1.5, sd = 1)
#' Y0 <- rnorm(100, mean = 0, sd = 1)
#'
#' result <- makarov_bounds_time_uniform(Y1, Y0, alpha = 0.1)
#'
#' plot(result$bounds_matrix[, "t"], result$bounds_matrix[, "lower"],
#'      type = "l", col = "blue", ylim = c(0, 1),
#'      xlab = "Treatment Effect", ylab = "CDF")
#' lines(result$bounds_matrix[, "t"], result$bounds_matrix[, "upper"],
#'       col = "red")
#'
#' @seealso \code{\link{makarov_bounds}} for fixed-sample version
#' @export
makarov_bounds_time_uniform <- function(Y1, Y0, alpha = 0.1, delta = 0.5,
                                         t_range = NULL, n_points = 1000) {
  # Compute time-uniform bounds
  tu_result <- compute_time_uniform_bounds(Y1, Y0, alpha = alpha, delta = delta)

  # Auto-determine t_range if not provided
  if (is.null(t_range)) {
    y_range <- range(c(Y1, Y0))
    spread <- diff(y_range)
    t_range <- c(-spread - 5, spread + 5)
  }

  # Extend bounds for edge cases
  bounds_ext <- lapply(tu_result$bounds, function(z) c(z, 0, 1))

  # Compute Makarov bounds using existing function
  bounds_matrix <- get_makarov_bounds(
    bounds = bounds_ext,
    obs_Y1 = tu_result$obs_Y1_sorted,
    obs_Y0 = tu_result$obs_Y0_sorted,
    t_range = t_range,
    n_points = n_points
  )

  list(
    bounds_matrix = bounds_matrix,
    tu_bounds = tu_result,
    alpha = alpha,
    delta = delta
  )
}
