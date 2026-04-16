#' @title Estimate nonparametric conditional meanimile
#'
#' @description
#'    Estimates a conditional meanimile in nonparametric way.
#'    This function acts as a unified wrapper, allowing the user to seamlessly switch
#'    between the square loss formulation and the general loss formulation
#'    as described in Debrauwer et al. (2026).
#'
#' @param X A numeric matrix of covariates (n x d).
#' @param Y A numeric vector of the response variable.
#' @param X0 A numeric matrix of evaluation points (m x d) where the curve should be estimated.
#' @param d_tau A function representing the density function of the
#'    distributional weight function D_tau.
#' @param tau A numeric risk level/index. Argument of d_tau function.
#' @param h_x_loc A numeric vector of local spatial bandwidths for the covariates X.
#' @param h_x_cdf A numeric vector of bandwidths for the covariates X used in the CDF estimation.
#' @param h_y_cdf A numeric scalar bandwidth for the response variable Y used in the CDF estimation.
#' @param loss_type Character string specifying the loss formulation (\code{"square"} or \code{"general"}).
#' @param lossDerivative A function defining the derivative of the loss function. Required for \code{"general"}.
#' @param delta Optional numeric parameter for the loss function. Only for \code{"general"}.
#'
#' @return A numeric vector of length m, containing the estimated meanimile predictions at each point in X0.
#' @export
#' @examples
#' set.seed(123)
#' n <- 300
#'
#' # 1. Simulate 1D Non-linear Data (Sine Wave)
#' X <- matrix(runif(n, -2, 2), ncol = 1)
#' Y <- sin(X[, 1] * 1.5) + rnorm(n, sd = 0.4)
#'
#' # 2. Define Evaluation Grid
#' X0 <- matrix(seq(-2, 2, length.out = 50), ncol = 1)
#'
#' # 3. Define Functions
#' d_tau_minvar <- function(u, tau) {
#'   (tau + 1) * u^tau
#' }
#' loss_deriv_expectile <- function(x, c, delta) {
#'   -2 * ifelse(x > c, delta, 1 - delta) * (x - c)
#' }
#'
#' # 4. Fit Square Loss Nonparametric
#' preds_sq <- MeanimileNonparametricRegression(
#'   X = X, Y = Y, X0 = X0, d_tau = d_tau_minvar, tau = 2,
#'   h_x_loc = 0.4, h_x_cdf = 0.4, h_y_cdf = 0.4, loss_type = "square"
#' )
#'
#' # 5. Fit General Loss Nonparametric (Expectile)
#' preds_gen <- MeanimileNonparametricRegression(
#'   X = X, Y = Y, X0 = X0, d_tau = d_tau_minvar, tau = 2,
#'   h_x_loc = 0.4, h_x_cdf = 0.4, h_y_cdf = 0.4, loss_type = "general",
#'   lossDerivative = loss_deriv_expectile, delta = 0.9
#' )
MeanimileNonparametricRegression <- function(X, Y, X0, d_tau, tau,
                                             h_x_loc, h_x_cdf, h_y_cdf,
                                             loss_type = c("square", "general"),
                                             lossDerivative = NULL, delta = NULL) {
  loss_type <- match.arg(loss_type)
  X <- as.matrix(X)
  X0 <- as.matrix(X0)

  if (ncol(X) != ncol(X0)) stop("X and X0 must have the same number of columns.")

  if (length(h_x_loc) == 1) h_x_loc <- rep(h_x_loc, ncol(X))
  if (length(h_x_cdf) == 1) h_x_cdf <- rep(h_x_cdf, ncol(X))

  if (loss_type == "square") {
    return(fit_nonpar_square(X, Y, X0, d_tau, tau, h_x_loc, h_x_cdf, h_y_cdf))
  } else {
    if (is.null(lossDerivative)) stop("lossDerivative required for general loss.")
    return(fit_nonpar_general(X, Y, X0, d_tau, tau, h_x_loc, h_x_cdf, h_y_cdf, lossDerivative, delta))
  }
}

# -------------------------------------------------------------------------
# Helper: Spatial Kernel Weights
# -------------------------------------------------------------------------
compute_spatial_weights <- function(X, x0, h_x_loc) {
  W_spatial <- rep(1, nrow(X))
  for (j in seq_len(ncol(X))) {
    W_spatial <- W_spatial * stats::dnorm((X[, j] - x0[j]) / h_x_loc[j])
  }
  return(W_spatial)
}

# -------------------------------------------------------------------------
# Internal Engine: Square Loss
# -------------------------------------------------------------------------
fit_nonpar_square <- function(X, Y, X0, d_tau, tau, h_x_loc, h_x_cdf, h_y_cdf) {
  n <- nrow(X)
  m <- nrow(X0)
  preds <- numeric(m)

  for (i in seq_len(m)) {
    x0 <- X0[i, ]

    x0_mat <- matrix(x0, nrow = n, ncol = length(x0), byrow = TRUE)
    F_hat <- ConditionalCDF(
      Y_eval = Y, X_eval = x0_mat,
      Y_train = Y, X_train = X,
      h_y = h_y_cdf, h_x = h_x_cdf
    )

    W_cdf <- d_tau(F_hat, tau)
    valid <- !is.na(W_cdf)

    # Spatial weights
    W_spatial <- compute_spatial_weights(X, x0, h_x_loc)
    W_total <- W_cdf * W_spatial

    # Local Linear Math
    X_loc <- -sweep(X, 2, x0, FUN = "-")
    X_loc <- cbind(1, X_loc)

    fit <- stats::lm.wfit(
      x = X_loc[valid, , drop = FALSE],
      y = Y[valid], w = W_total[valid]
    )
    preds[i] <- fit$coefficients[1]
  }

  return(preds)
}

# -------------------------------------------------------------------------
# Internal Engine: General Loss
# -------------------------------------------------------------------------
fit_nonpar_general <- function(X, Y, X0, d_tau, tau, h_x_loc, h_x_cdf, h_y_cdf, lossDerivative, delta) {
  n <- nrow(X)
  m <- nrow(X0)
  preds <- numeric(m)

  solve_local_root <- function(Y_val, X_loc_val, W_total_val) {
    obj_fun <- function(beta) {
      preds <- as.vector(X_loc_val %*% beta)
      if (is.null(delta)) {
        derivs <- lossDerivative(Y_val, preds)
      } else {
        derivs <- lossDerivative(Y_val, preds, delta)
      }
      eq_sums <- colSums(X_loc_val * (W_total_val * derivs))
      return(sum(eq_sums^2))
    }

    start_beta <- stats::lm.wfit(x = X_loc_val, y = Y_val, w = W_total_val)$coefficients
    start_beta[is.na(start_beta)] <- 0


    opt_res <- stats::optim(par = start_beta, fn = obj_fun, method = "BFGS")
    return(opt_res$par[1])
  }

  for (i in seq_len(m)) {
    x0 <- X0[i, ]

    x0_mat <- matrix(x0, nrow = n, ncol = length(x0), byrow = TRUE)
    F_hat <- ConditionalCDF(
      Y_eval = Y, X_eval = x0_mat,
      Y_train = Y, X_train = X,
      h_y = h_y_cdf, h_x = h_x_cdf
    )

    W_cdf <- d_tau(F_hat, tau)
    valid <- !is.na(W_cdf)

    # Spatial weights
    W_spatial <- compute_spatial_weights(X, x0, h_x_loc)
    W_total <- W_cdf * W_spatial

    X_loc <- -sweep(X, 2, x0, FUN = "-")
    X_loc <- cbind(1, X_loc)[valid, , drop = FALSE]

    preds[i] <- solve_local_root(Y[valid], X_loc, W_total[valid])
  }

  return(preds)
}
