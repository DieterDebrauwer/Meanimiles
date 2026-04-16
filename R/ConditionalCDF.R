#' @title Estimate conditional cumulative distribution function
#'
#' @description
#'    Computes a doubly smoothed Nadaraya-Watson estimator for the conditional
#'    cumulative distribution function.
#'
#' @param Y_eval A numeric vector of Y values where the CDF should be evaluated.
#' @param X_eval A numeric matrix of X values where the CDF should be evaluated.
#' @param Y_train A numeric vector of training Y values.
#' @param X_train A numeric matrix of training X values.
#' @param h_y A numeric scalar representing the bandwidth for the
#'    response variable Y.
#' @param h_x A numeric vector representing the bandwidths for the covariates X.
#'
#' @return A numeric vector of estimated CDF probabilities.
#' @export
ConditionalCDF <- function(Y_eval, X_eval, Y_train, X_train, h_y, h_x) {
  X_eval <- as.matrix(X_eval)
  X_train <- as.matrix(X_train)

  n_eval <- nrow(X_eval)
  n_train <- nrow(X_train)
  d <- ncol(X_train)

  if (length(h_x) == 1) h_x <- rep(h_x, d)

  # Initialize the weight matrix with 1s (size: n_eval x n_train)
  W <- matrix(1, nrow = n_eval, ncol = n_train)

  # Vectorized computation for all covariate dimensions
  for (k in seq_len(d)) {
    # outer() creates an n_eval x n_train matrix of all pairwise differences
    K_k <- stats::dnorm(outer(X_eval[, k], X_train[, k], "-") / h_x[k])
    W <- W * K_k # Element-wise multiplication
  }

  # Vectorized computation for the Y response
  # We use outer() here too to ensure we get the full n_eval x n_train matrix
  W_y <- stats::pnorm(outer(Y_eval, Y_train, "-") / h_y)

  # Nadaraya-Watson combination using rowSums
  numerator <- rowSums(W * W_y)
  denominator <- rowSums(W)

  # Calculate final CDF and handle division by zero
  F_hat <- numerator / denominator
  F_hat[denominator < 1e-12] <- 0

  return(F_hat)
}
