#' @title Estimate linear conditional meanimile
#'
#' @description
#'    Estimates the parameters of a linear conditional meanimile model using a two-fold
#'    sample splitting procedure to mitigate bias from the conditional CDF estimation.
#'    This function acts as a unified wrapper, allowing the user to seamlessly switch
#'    between the square loss formulation and the general loss formulation.
#'
#' @param X A numeric matrix of covariates (n x d).
#' @param Y A numeric vector of the response variable.
#' @param d_tau A function representing the density function of the
#'    distributional weight function D_tau.
#' @param tau A numeric risk level/index. Argument of d_tau function.
#' @param h_x A numeric vector of bandwidths for the covariates X.
#' @param h_y A numeric scalar bandwidth for the response variable Y.
#' @param loss_type A character string specifying the loss formulation. Must be either \code{"square"} (default) or \code{"general"}.
#' @param lossDerivative A function defining the derivative of the loss function. Required only if \code{loss_type = "general"}.
#' @param delta Optional numeric parameter for the loss function. Used only if \code{loss_type = "general"}.
#'
#' @return A numeric vector of estimated regression coefficients (including the intercept).
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 400
#' X <- matrix(rnorm(n * 2), ncol = 2)
#' Y <- 0.9 * X[, 1] + 0.5 * X[, 2] + rnorm(n)
#' d_tau_minvar <- function(u, tau) {
#'   (tau + 1) * u^tau
#' }
#' lossDerivativeExpectiles <- function(x, c, delta) {
#'   -2 * ifelse(x > c, delta, 1 - delta) * (x - c)
#' }
#'
#' # 1. Square Loss
#' MeanimileLinearRegression(X, Y,
#'   d_tau = d_tau_minvar, tau = 2,
#'   h_x = c(0.4, 0.4), h_y = 0.4, loss_type = "square"
#' )
#'
#' # 2. General Loss
#' MeanimileLinearRegression(X, Y,
#'   d_tau = d_tau_minvar, tau = 2,
#'   h_x = c(0.4, 0.4), h_y = 0.4, loss_type = "general",
#'   lossDerivative = lossDerivativeExpectiles, delta = 0.9
#' )
MeanimileLinearRegression <- function(X, Y, d_tau, tau, h_x, h_y,
                                      loss_type = c("square", "general"),
                                      lossDerivative = NULL, delta = NULL) {
  # Ensure the user picks a valid option, defaults to "square"
  loss_type <- match.arg(loss_type)

  # Route to the correct internal engine
  if (loss_type == "square") {
    return(fit_linear_square(X, Y, d_tau, tau, h_x, h_y))
  } else if (loss_type == "general") {
    if (is.null(lossDerivative)) {
      stop("You must provide 'lossDerivative' when using loss_type = 'general'.")
    }
    return(fit_linear_general(X, Y, d_tau, tau, h_x, h_y, lossDerivative, delta))
  }
}

# -------------------------------------------------------------------------
# Internal Engine: Square Loss
# -------------------------------------------------------------------------
fit_linear_square <- function(X, Y, d_tau, tau, h_x, h_y) {
  X <- as.matrix(X)
  n <- nrow(X)
  X_G <- cbind(1, X)

  # SAFE COLUMN NAMING
  if (is.null(colnames(X))) {
    colnames(X_G) <- c("(Intercept)", paste0("X", seq_len(ncol(X))))
  } else {
    colnames(X_G) <- c("(Intercept)", colnames(X))
  }

  indices <- sample(seq_len(n))
  split_point <- floor(n / 2)
  idx_S1 <- indices[1:split_point]
  idx_S2 <- indices[(split_point + 1):n]

  # Fold 1
  F_hat_S2 <- ConditionalCDF(Y[idx_S2], X[idx_S2, , drop = FALSE], Y[idx_S1], X[idx_S1, , drop = FALSE], h_y, h_x)
  W_S2 <- d_tau(F_hat_S2, tau)
  valid_S2 <- !is.na(W_S2)
  fit_S2 <- stats::lm.wfit(
    x = X_G[idx_S2, , drop = FALSE][valid_S2, , drop = FALSE],
    y = Y[idx_S2][valid_S2], w = W_S2[valid_S2]
  )
  beta_2 <- fit_S2$coefficients

  # Fold 2
  F_hat_S1 <- ConditionalCDF(Y[idx_S1], X[idx_S1, , drop = FALSE], Y[idx_S2], X[idx_S2, , drop = FALSE], h_y, h_x)
  W_S1 <- d_tau(F_hat_S1, tau)
  valid_S1 <- !is.na(W_S1)
  fit_S1 <- stats::lm.wfit(
    x = X_G[idx_S1, , drop = FALSE][valid_S1, , drop = FALSE],
    y = Y[idx_S1][valid_S1], w = W_S1[valid_S1]
  )
  beta_1 <- fit_S1$coefficients

  return((beta_1 + beta_2) / 2)
}


fit_linear_general <- function(X, Y, d_tau, tau, h_x, h_y, lossDerivative, delta) {
  X <- as.matrix(X)
  n <- nrow(X)
  X_G <- cbind(1, X)

  #  SAFE COLUMN NAMING
  if (is.null(colnames(X))) {
    colnames(X_G) <- c("(Intercept)", paste0("X", seq_len(ncol(X))))
  } else {
    colnames(X_G) <- c("(Intercept)", colnames(X))
  }


  indices <- sample(seq_len(n))
  split_point <- floor(n / 2)
  idx_S1 <- indices[1:split_point]
  idx_S2 <- indices[(split_point + 1):n]


  # Helper function to solve the estimating equations for a given fold
  solve_fold <- function(idx_eval, idx_train) {
    F_hat <- ConditionalCDF(
      Y[idx_eval], X[idx_eval, , drop = FALSE],
      Y[idx_train], X[idx_train, , drop = FALSE], h_y, h_x
    )
    W <- d_tau(F_hat, tau)
    valid <- !is.na(W)

    Y_val <- Y[idx_eval][valid]
    X_val <- X_G[idx_eval, , drop = FALSE][valid, , drop = FALSE]
    W_val <- W[valid]

    # We want to find beta such that the sum of the estimating equations is zero.
    obj_fun <- function(beta) {
      preds <- as.vector(X_val %*% beta)
      if (is.null(delta)) {
        derivs <- lossDerivative(Y_val, preds)
      } else {
        derivs <- lossDerivative(Y_val, preds, delta)
      }
      eq_sums <- colSums(X_val * (W_val * derivs))
      return(sum(eq_sums^2))
    }

    start_beta <- stats::lm.fit(x = X_val, y = Y_val)$coefficients
    start_beta[is.na(start_beta)] <- 0

    opt_res <- stats::optim(par = start_beta, fn = obj_fun, method = "BFGS")
    return(opt_res$par)
  }

  beta_2 <- solve_fold(idx_eval = idx_S2, idx_train = idx_S1)
  beta_1 <- solve_fold(idx_eval = idx_S1, idx_train = idx_S2)

  return((beta_1 + beta_2) / 2)
}
