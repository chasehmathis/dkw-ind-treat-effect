#' Generate Simulated Treatment/Control Data with Potential Outcomes
#'
#' Generates a dataset with bivariate normal potential outcomes for use in
#' treatment effect estimation. Uses the DeclareDesign framework for
#' principled simulation of randomized experiments.
#'
#' @param m Integer. Sample size (total number of units).
#' @param p Numeric in (0, 1). Proportion of units assigned to treatment.
#'   Default is 0.5.
#' @param tau Numeric. True average treatment effect (difference in means
#'   between treatment and control potential outcomes). Default is 1.
#' @param rho Numeric in (-1, 1). Correlation between potential outcomes
#'   Y(0) and Y(1). Higher values indicate more homogeneous treatment effects.
#'   Default is 0.5.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{ID}{Unit identifier}
#'     \item{Y_Z_0}{Potential outcome under control}
#'     \item{Y_Z_1}{Potential outcome under treatment}
#'     \item{Z}{Treatment assignment indicator (0 = control, 1 = treatment)}
#'     \item{Y}{Observed outcome (revealed based on Z)}
#'     \item{tau}{True individual treatment effect (Y_Z_1 - Y_Z_0)}
#'   }
#'
#' @details
#' The potential outcomes are drawn from a bivariate normal distribution:
#' \deqn{(Y(0), Y(1))^T \sim N\left(\begin{pmatrix} 0 \\ \tau \end{pmatrix},
#' \begin{pmatrix} 1 & \rho \\ \rho & 1 \end{pmatrix}\right)}
#'
#' Treatment is assigned via complete random assignment with exactly \code{m * p}

#' units receiving treatment.
#'
#' @examples
#' # Generate data with moderate correlation
#' data <- generate_potential_outcomes(m = 100, p = 0.5, tau = 2, rho = 0.3)
#' head(data)
#'
#' # Check that observed outcomes match potential outcomes
#' all(data$Y == ifelse(data$Z == 1, data$Y_Z_1, data$Y_Z_0))
#'
#' @export
generate_potential_outcomes <- function(m, p = 0.5, tau = 1, rho = 0.5) {
  if (!requireNamespace("DeclareDesign", quietly = TRUE)) {
    stop("Package 'DeclareDesign' is required. Install with: install.packages('DeclareDesign')")
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install with: install.packages('MASS')")
  }

  # Validate inputs
  if (m < 2) stop("Sample size m must be at least 2")
  if (p <= 0 || p >= 1) stop("Treatment proportion p must be in (0, 1)")
  if (abs(rho) >= 1) stop("Correlation rho must be in (-1, 1)")

  design <- DeclareDesign::declare_model(
    N = m,
    fabricatr::draw_multivariate(
      c(Y_Z_0, Y_Z_1) ~ MASS::mvrnorm(
        n = m,
        mu = c(0, tau),
        Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
      )
    )
  ) +
    DeclareDesign::declare_assignment(Z = randomizr::complete_ra(N = m, m = round(m * p))) +
    DeclareDesign::declare_measurement(Y = fabricatr::reveal_outcomes(Y ~ Z))

  data <- DeclareDesign::draw_data(design)
  data$tau <- data$Y_Z_1 - data$Y_Z_0

  return(data)
}


#' Create Modified ECDF with Custom Scaling
#'
#' Creates an empirical cumulative distribution function with scaling
#' factor 1/(n*p), useful for certain finite population corrections.
#'
#' @param data Numeric vector. The sample data.
#' @param n Integer. Population size for scaling.
#' @param p Numeric in (0, 1). Treatment proportion for scaling. Default is 0.5.
#'
#' @return A function that evaluates the modified ECDF at given points.
#'
#' @details
#' For each point x, the modified ECDF returns:
#' \deqn{\hat{F}_{mod}(x) = \min\left(\frac{1}{np} \sum_{i=1}^{n} \mathbf{1}(X_i < x), 1\right)}
#'
#' This is primarily used internally for finite population inference.
#'
#' @keywords internal
ecdf_scaled <- function(data, n, p = 0.5) {
  force(n)
  force(p)

  function(xs) {
    vapply(xs, function(x) {
      less_than_x <- sum(data < x)
      min((1 / (n * p)) * less_than_x, 1)
    }, numeric(1))
  }
}
