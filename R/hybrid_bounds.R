#' Simes Upper Band Piece
#'
#' Computes the Simes upper confidence band piece for empirical CDFs,
#' based on Proposition C.3.
#'
#' @param Fhat Numeric vector. Empirical CDF values in [0, 1].
#' @param delta Numeric in (0, 1). Significance level for this component.
#' @param m Integer. Sample size.
#' @param k Integer. Order parameter. Default is \code{floor(m/2)}.
#'
#' @return Numeric vector of upper bound values.
#'
#' @details
#' Uses log-space arithmetic for numerical stability when computing
#' binomial coefficients.
#'
#' @keywords internal
h_simes <- function(Fhat, delta, m, k = floor(m / 2)) {
  Fhat <- pmin(pmax(Fhat, 0), 1)

  b <- numeric(m + 2L)
  b[1L] <- 0
  b[m + 2L] <- 1

  logC_mk <- lchoose(m, k)
  delta_root <- delta^(1 / k)

  for (i in 1:m) {
    if (i < k) {
      ratio_1k <- 0
    } else {
      logCi <- lchoose(i, k)
      log_ratio <- (logCi - logC_mk) / k
      ratio_1k <- exp(log_ratio)
    }
    b[m + 2L - i] <- 1 - delta_root * ratio_1k
  }

  j <- ceiling(m * Fhat)
  j <- pmin(pmax(j, 0L), m + 1L)
  b[j + 1L]
}


#' Dempster Line-Crossing Formula
#'
#' Computes the Dempster line-crossing probability (Proposition C.4).
#'
#' @param a Numeric in (0, 1). Intercept parameter.
#' @param b Numeric in (0, 1). Slope parameter.
#' @param m Integer. Sample size.
#'
#' @return Numeric. Crossing probability.
#'
#' @keywords internal
Delta_Dempster <- function(a, b, m) {
  jmax <- floor(m * (1 - b))
  j <- 0:jmax

  comb <- choose(m, j)
  term1 <- a + (1 - a) * j / ((1 - b) * m)
  term2 <- 1 - a - (1 - a) * j / ((1 - b) * m)

  pow1 <- numeric(length(j))
  pow1[j == 0] <- 1 / term1[j == 0]
  pow1[j > 0] <- term1[j > 0]^(j[j > 0] - 1)

  pow2 <- term2^(m - j)

  a * sum(comb * pow1 * pow2)
}


#' Solve for Dempster Parameter
#'
#' Solves for parameter a given b and delta using the Dempster formula.
#'
#' @param b Numeric. Slope parameter.
#' @param delta Numeric. Target probability.
#' @param m Integer. Sample size.
#'
#' @return Numeric. Solved value of a.
#'
#' @keywords internal
a_dempster <- function(b, delta, m) {
  f <- function(a) Delta_Dempster(a, b, m) - delta
  stats::uniroot(f, c(1e-6, 1 - 1e-6))$root
}


#' Dempster Lower Bound
#'
#' Computes the pure Dempster lower confidence bound for empirical CDFs.
#'
#' @param Fhat Numeric vector. Empirical CDF values.
#' @param delta Numeric. Significance level.
#' @param m Integer. Sample size.
#' @param b Numeric. Slope parameter. Default is \code{5/m}.
#'
#' @return Numeric vector of lower bound values.
#'
#' @keywords internal
dempster_lower <- function(Fhat, delta, m, b = 5 / m) {
  a <- a_dempster(b, delta, m)
  coef <- (1 - a) / (1 - b)
  pmax(0, coef * pmax(0, Fhat - b))
}


#' Hybrid Upper Confidence Band
#'
#' Computes a hybrid upper confidence band combining Simes and DKW methods.
#' This band is sharper than either method alone in different regions.
#'
#' @param Fhat Numeric vector. Empirical CDF values.
#' @param delta Numeric in (0, 1). Significance level.
#' @param m Integer. Sample size.
#' @param k Integer. Order parameter for Simes. Default is \code{floor(m/2)}.
#'
#' @return Numeric vector of upper bound values.
#'
#' @details
#' The hybrid band takes the minimum of the Simes and DKW upper bounds:
#' \deqn{\bar{F}_{hybrid}(x) = \min(\bar{F}_{Simes}(x), \bar{F}_{DKW}(x))}
#'
#' This exploits the fact that Simes provides tighter bounds near the tails
#' while DKW is tighter in the middle of the distribution.
#'
#' @examples
#' set.seed(123)
#' x <- sort(runif(100))
#' Fhat <- ecdf(x)(x)
#'
#' upper <- hybrid_upper(Fhat, delta = 0.05, m = 100)
#' plot(x, Fhat, type = "s")
#' lines(x, upper, col = "red")
#'
#' @export
hybrid_upper <- function(Fhat, delta, m, k = floor(m / 2)) {
  up_simes <- h_simes(Fhat, delta / 2, m, k)
  up_dkw <- Fhat + sqrt(log(2 / delta) / (2 * m))
  pmin(up_simes, up_dkw, 1)
}


#' Hybrid Lower Confidence Band
#'
#' Computes a hybrid lower confidence band combining Dempster and DKW methods.
#'
#' @param Fhat Numeric vector. Empirical CDF values.
#' @param delta Numeric in (0, 1). Significance level.
#' @param m Integer. Sample size.
#' @param b Numeric. Slope parameter for Dempster. Default is \code{5/m}.
#'
#' @return Numeric vector of lower bound values.
#'
#' @details
#' The hybrid band takes the maximum of the Dempster and DKW lower bounds:
#' \deqn{\underline{F}_{hybrid}(x) = \max(\underline{F}_{Dempster}(x), \underline{F}_{DKW}(x))}
#'
#' @examples
#' set.seed(123)
#' x <- sort(runif(100))
#' Fhat <- ecdf(x)(x)
#'
#' lower <- hybrid_lower(Fhat, delta = 0.05, m = 100)
#' plot(x, Fhat, type = "s")
#' lines(x, lower, col = "blue")
#'
#' @export
hybrid_lower <- function(Fhat, delta, m, b = 5 / m) {
  low_demp <- dempster_lower(Fhat, delta, m, b)
  eps <- sqrt(log(2 / delta) / (2 * m))
  low_dkw <- Fhat - eps
  pmax(0, pmax(low_demp, low_dkw))
}


#' Compute Hybrid Confidence Bands
#'
#' Computes hybrid confidence bands for empirical CDFs combining
#' Simes (upper) and Dempster (lower) methods with DKW.
#'
#' @param Fhat Numeric vector. Empirical CDF values.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#' @param m Integer. Sample size.
#' @param k Integer. Order parameter for Simes. Default is \code{floor(m/2)}.
#' @param b Numeric. Slope parameter for Dempster. Default is \code{5/m}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{lower}{Numeric vector of lower bounds.}
#'     \item{upper}{Numeric vector of upper bounds.}
#'     \item{lower_dkw}{Standard DKW lower bounds (for comparison).}
#'     \item{upper_dkw}{Standard DKW upper bounds (for comparison).}
#'   }
#'
#' @details
#' This function computes confidence bands that are provably valid and
#' often tighter than standard DKW bands, particularly in the tails of
#' the distribution.
#'
#' The bands are computed at significance level \code{alpha/8} for each
#' component to account for the combination of methods.
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x <- sort(runif(n))
#' Fhat <- ecdf(x)(x)
#'
#' bands <- hybrid_band(Fhat, alpha = 0.05, m = n)
#'
#' plot(x, Fhat, type = "s", ylim = c(0, 1))
#' lines(x, bands$lower, col = "red", lwd = 2)
#' lines(x, bands$upper, col = "red", lwd = 2)
#' lines(x, bands$lower_dkw, col = "gray", lty = 2)
#' lines(x, bands$upper_dkw, col = "gray", lty = 2)
#' legend("bottomright",
#'        c("ECDF", "Hybrid", "DKW"),
#'        col = c("black", "red", "gray"),
#'        lty = c(1, 1, 2))
#'
#' @export
hybrid_band <- function(Fhat,m, finite_pop = FALSE, alpha = 0.05, k = floor(m / 2), b = 5 / m) {
  delta <- alpha / 8
  up_simes <- h_simes(Fhat, delta, m, k)
  lower_simes <- 1-h_simes(1-Fhat, delta, m, k)
  
  upper_d <- 1-dempster_lower(1-Fhat, delta, m, b)
  lower_d <- dempster_lower(Fhat, delta, m, b)
  
  dkw <- dkw_band(Fhat, N = m, alpha = alpha, finite_pop = finite_pop)

  upper_h <- pmin(up_simes, upper_d, dkw$upper)
  lower_h <- pmax(lower_simes, lower_d, dkw$lower)
  
  list(
    lower = lower_h,
    upper = upper_h,
    lower_dkw = dkw$lower,
    upper_dkw = dkw$upper
  )
}
