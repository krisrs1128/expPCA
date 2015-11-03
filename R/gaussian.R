
################################################################################
# Exponential Family PCA for Gaussian Data
################################################################################

# optim-funs -------------------------------------------------------------------
#' @title Bregman Loss in Gaussian case
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the bregman loss of x given the current parameters.
gaussian_obj <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- rep(0, length(x))
  }
  theta <- a * v + s
  .5 * (sum((x - theta) ^ 2) + lambda * sum(mu0 - theta) ^ 2)
}

#' @title Gradient of Bregman Loss in Gaussian Case with respect to v
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the gradient of the bregman loss with respect to the
#' parameter v.
gaussian_grad_v <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- rep(0, length(x))
  }
  theta <- a * v + s
  - sum(a * ((x - theta) + lambda * (mu0 - theta)))
}

#' @title Gradient of Bregman Loss in Gaussian Case with respect to a
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the gradient of the bregman loss with respect to the
#' parameter a.
gaussian_grad_a <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- rep(0, length(x))
  }
  theta <- a * v + s
  - sum(v * ((x - theta) + lambda * (mu0 - theta)))
}

# gaussian-pca ----------------------------------------------------------------
#' @title Gaussian PCA
#' @param X The n x p binary matrix with samples along rows which we want to
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
gaussian_exp_pca <- function(X, n_comp = 2, n_cycle = 30, n_iter = 30, eps = 1e-4,
                             lambda = 0.01, mu0 = NULL) {
  args <- list(obj_fun = gaussian_obj,
               obj_grad_a = gaussian_grad_a,
               obj_grad_v = gaussian_grad_v,
               X = X,
               n_comp = n_comp,
               n_cycle = n_cycle,
               n_iter = n_iter,
               eps = eps,
               lambda = lambda,
               mu0 = mu0)
  do.call(exp_pca, args)
}
