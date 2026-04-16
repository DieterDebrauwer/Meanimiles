
#'
#' @title Estimate a univariate meanimile
#'
#' @description
#' Calculates the empirical meanimile for a given (univariate) data sample, weight function,
#' and loss derivative.
#'
#' @param d_tau A function representing the density function of the distributional weight function D_tau.
#' @param tau A numeric risk level/index (e.g., 0.95). Argument of d_tau function.
#' @param lossDerivative A function defining the derivative of the loss function.
#' @param delta Optional numeric parameter for the loss function (default: NULL).
#' @param sample A numeric vector of univariate data.
#' @param grid_size An integer defining the resolution of the fallback grid (default: 1000).
#'
#' @return A numeric scalar representing the estimated meanimile.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' d_tau_expectile <- function(u, tau){1}
#' lossDerivativeExpectiles = function(x, c, delta) { -2 * ifelse(x > c, delta, 1 - delta) * (x - c) }
#' MeanimileEstimator(d_tau_expectile, tau=0, lossDerivativeExpectiles,delta=0.95, sample = x)
#'
#' @references
#' D. Debrauwer, I. Gijbels, and K. Herrmann. (2026)
#' On a general class of functionals: Statistical inference and application to risk measures
#' Electronic Journal of Statistics
#' doi:https://doi.org/10.1214/25-EJS2391

MeanimileEstimator <- function(d_tau, tau, lossDerivative, delta = NULL, sample, grid_size = 1000) {

  if (!is.numeric(sample)) stop("'sample' must be a numeric vector.")

  x_sorted <- sort(sample[!is.na(sample)])
  n <- length(x_sorted)

  if (n == 0) stop("No non-missing data provided in 'sample'.")
  if (!is.function(d_tau)) stop("'d_tau' must be a function.")
  if (!is.function(lossDerivative)) stop("'lossDerivative' must be a function.")

  ranks <- seq_len(n) / (n + 1)
  d_tauWeights <- d_tau(ranks, tau)

  if (is.null(delta)) {
    lambdaHat <- function(c) sum(d_tauWeights * lossDerivative(x_sorted, c))
  } else {
    lambdaHat <- function(c) sum(d_tauWeights * lossDerivative(x_sorted, c, delta))
  }

  res <- tryCatch({
    root_res <- stats::uniroot(lambdaHat, interval = c(min(x_sorted), max(x_sorted)), extendInt = "yes")
    root_res$root
  }, error = function(e) {
    grid <- seq(min(x_sorted), max(x_sorted), length.out = grid_size)
    grid_evals <- vapply(grid, lambdaHat, FUN.VALUE = numeric(1))
    grid[which.min(abs(grid_evals))]
  })

  return(res)
}

