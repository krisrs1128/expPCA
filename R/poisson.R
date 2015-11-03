
################################################################################
# Exponential Family PCA for Poisson Data
################################################################################

# optim-funs -------------------------------------------------------------------
#' @title Version of x * log(x) continuous at 0
z_log_z <- function(z) {
  res <- z * log(z)
  res[z == 0] <- 0
  return (res)
}

#' @title Bregman Loss for Poisson data
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the bregman loss of x given the current parameters.
poisson_obj <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- rep(1, length(x))
  }
  theta <- a * v + s
  sum(exp(theta) - x * theta + z_log_z(x) - x) +
    lambda * sum(exp(theta) - mu0 * theta + z_log_z(mu0) - mu0)
}

#' @title Gradient of Bregman Loss for Poisson data with respect v
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the gradient of the bregman loss with respect to the
#' parameter v.
poisson_grad_v <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- rep(1, length(x))
  }
  theta <- a * v + s
  sum(a * exp(theta) - x * a * theta) +
    lambda * sum(a * exp(theta) - mu0 * a * theta)
}

#' @title Gradient of Bregman Loss for Poisson data with respect a
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the gradient of the bregman loss with respect to the
#' parameter v.
poisson_grad_a <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- rep(1, length(x))
  }
  theta <- a * v + s
  sum(v * exp(theta) - x * v * theta) +
    lambda * sum(v * exp(theta) - mu0 * v * theta)
}

# poisson-pca ------------------------------------------------------------------
#' @title Poisson PCA
#' @param X The n x p Poisson matrix with samples along rows which we want to
#' decompose using exponential family PCA.
#' @param n_comp How many principal components should we return?
#' @param n_cycle How many times should the iterative optimization pass through
#' across each component?
#' @param n_iter The maximum number of iterations to run the PCA.
#' @param eps The convergence criterion. If the mean change in the scores is
#' less than eps, we return.
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return A list containing the converged values of A and V.
#' @references Collins, Michael, Sanjoy Dasgupta, and Robert E. Schapire. "A generalization of principal components analysis to the exponential family." Advances in neural information processing systems. 2001.
#' @export
poisson_exp_pca <- function(X, n_comp = 2, n_cycle = 30, n_iter = 30, eps = 1e-4,
                            lambda = 0.001, mu0 = NULL) {

  args <- list(obj_fun = poisson_obj,
               obj_grad_a = poisson_grad_a,
               obj_grad_v = poisson_grad_v,
               X = X,
               n_comp = n_comp,
               n_cycle = n_cycle,
               n_iter = n_iter,
               eps = eps,
               lambda = lambda,
               mu0 = mu0)

  do.call(exp_pca, args)
}
