
#' @title  Estimate portfolio meanimile
#'
#' @description
#' Estimates the risk of an aggregated portfolio using the meanimile framework.
#' Supports parametric and nonparametric modeling
#' of marginal distributions and copula.
#'
#' @param data A numeric matrix or data.frame (n observations x d assets).
#' @param weights A numeric vector of portfolio weights summing to 1.
#' @param d_tau A function representing the density function of the distributional weight function D_tau.
#' @param tau A numeric risk level/index (e.g., 0.95). Argument of d_tau function.
#' @param lossDerivative A function defining the derivative of the loss function.
#' @param delta Optional numeric parameter for the loss function (default: NULL).
#' @param copula_obj An unfitted 'copula' object (e.g., \code{copula::gumbelCopula(dim=3)}).
#'   If \code{NULL}, an Empirical Beta Copula is used.
#' @param margin_dists A character vector of base R distribution root names for \code{fitdistrplus::fitdist}
#'   (e.g., \code{c("norm", "exp", "weibull")}). If \code{NULL}, empirical margins are used.
#' @param N Integer. Number of Monte Carlo simulations to draw (default: 100 000).
#' @param grid_size An integer defining the resolution of the fallback grid (default: 1000).
#'
#' @return  A numeric scalar representing the estimated meanimile.
#' @export
#' @examples
#'
#' gumbel.cop <- copula::gumbelCopula(10/9, dim=3)
#' myMvdX <- copula::mvdc(copula=gumbel.cop, margins=c("norm", "norm","exp"),
#'                paramMargins=list(list(mean=0, sd=1),
#'                                  list(mean=0, sd=1),
#'                                  list(rate=1)))
#' X <- copula::rMvdc(500,myMvdX)
#'
#' d_tau_expectile <- function(u, tau){1}
#' lossDerivativeExpectiles = function(x, c, delta) { -2 * ifelse(x > c, delta, 1 - delta) * (x - c) }
#' gumbelCopula=copula::gumbelCopula(dim=3)
#' margin_dists=c("norm","norm","exp")
#' weights=c(1/3,1/3,1/3)
#' PortfolioMeanimile(data=X, weights, d_tau_expectile, tau=0, lossDerivativeExpectiles, delta = 0.95,
#'   copula_obj = gumbelCopula, margin_dists = margin_dists, N = 100000, grid_size = 1000)
#'
#'@references
#'D. Debrauwer and I. Gijbels (2026)
#'Copula-based estimation of meanimiles of aggregated risks.
#'Metrika
#'doi:https://doi.org/10.1007/s00184-026-01022-9
PortfolioMeanimile <- function(data, weights, d_tau, tau, lossDerivative, delta = NULL,
                               copula_obj = NULL, margin_dists = NULL, N = 100000, grid_size = 1000) {

  # Input Validation
  data <- as.matrix(data)
  if (!is.numeric(data)) stop("'data' must be a numeric matrix or data.frame.")
  if (any(is.na(data))) stop("'data' contains missing values. Please handle NAs before fitting.")

  n <- nrow(data)
  d <- ncol(data)

  if (d<=1) stop("Dimension must be at least two.")
  if (!is.numeric(weights)) stop("'weights' must be numeric.")
  if (length(weights) != d) stop(sprintf("Length of weights (%d) must match number of assets (%d).", length(weights), d))
  if (abs(sum(weights) - 1) > 1e-5) stop("Portfolio weights must sum to 1.")
  if (N <= 0 || grid_size <= 0) stop("'N' and 'grid_size' must be positive integers.")

  is_param_copula <- !is.null(copula_obj)
  is_param_margins <- !is.null(margin_dists)

  if (is_param_margins && length(margin_dists) != d) {
    stop("If 'margin_dists' is provided, it must have a distribution name for every column in 'data'.")
  }

  # Estimate Margins & Compute Pseudo-Observations
  U_data <- matrix(NA, nrow = n, ncol = d)
  fitted_margin_params <- list()

  for (j in 1:d) {
    if (is_param_margins) {

      # User passes the actual root (e.g., "norm", "weibull", "exp")
      r_root <- margin_dists[j]

      # Fit safely using fitdistrplus
      fit <- tryCatch({
        fitdistrplus::fitdist(data[, j], distr = r_root)
      }, error = function(e) {
        stop(sprintf("Failed to fit '%s' on column %d. Check your data or distribution choice. Error: %s", r_root, j, e$message))
      })

      fitted_margin_params[[j]] <- as.list(fit$estimate)

      # Transform to Uniform
      p_func_name <- paste0("p", r_root)
      if (!exists(p_func_name, mode = "function")) {
        stop(sprintf("Probability function %s not found.", p_func_name))
      }
      p_func <- get(p_func_name, mode = "function")

      U_data[, j] <- do.call(p_func, c(list(q = data[, j]), fitted_margin_params[[j]]))

    } else {
      # Nonparametric Margins
      U_data[, j] <- rank(data[, j]) / (n + 1)
    }
  }

  # Estimate Copula & Draw N Samples
  if (is_param_copula) {

    if (!inherits(copula_obj, "copula")) stop("'copula_obj' must be an object of class 'copula'.")

    fitted_copula <- tryCatch({
      copula::fitCopula(copula_obj, data = U_data, method = "mpl",estimate.variance=F)@copula
    }, error = function(e) {
      stop(paste("Failed to fit parametric copula. Error:", e$message))
    })

    U_sim <- copula::rCopula(N, fitted_copula)

  } else {
    # Nonparametric Copula (Empirical Beta Copula with smoothing)
    ecop <- copula::empCopula(U_data, smoothing = "beta")
    U_sim <- copula::rCopula(N, copula = ecop)
  }

  # Transform U_sim back to X_sim via Inverse Margins
  X_sim <- matrix(NA, nrow = N, ncol = d)

  for (j in 1:d) {
    if (is_param_margins) {
      q_func_name <- paste0("q", margin_dists[j])
      if (!exists(q_func_name, mode = "function")) stop(sprintf("Quantile function %s not found.", q_func_name))
      q_func <- get(q_func_name, mode = "function")

      X_sim[, j] <- do.call(q_func, c(list(p = U_sim[, j]), as.list(fitted_margin_params[[j]])))
    } else {
      X_sim[, j] <- stats::quantile(data[, j], probs = U_sim[, j], type = 1, names = FALSE)
    }
  }

  # Aggregate Portfolio and Estimate
  Z_sim <- as.vector(X_sim %*% weights)

  result <- MeanimileEstimator(
    d_tau = d_tau,
    tau = tau,
    lossDerivative = lossDerivative,
    delta = delta,
    sample = Z_sim,
    grid_size = grid_size
  )

  return(result)
}




