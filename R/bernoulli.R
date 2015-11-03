
################################################################################
# Exponential Family PCA for Bernoulli Data
################################################################################

# optim-funs -------------------------------------------------------------------
#' @title Logistic Function
#' @param x Input vector or scale for logistic fun.
#' @return v Value of logistic function at input.
r <- function(x) {
  exp(-x) / (1 + exp(-x))
}

#' @title Bregman Loss in Bernoulli case
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the bregman loss of x given the current parameters.
bern_obj <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  z <- 2 * x - 1
  sum(log(1 + exp(- z * (a * v + s)))) + lambda * sum(log(1 + exp(- mu0 * (a * v + s))))
}

#' @title Gradient of Bregman Loss in Bernoulli Case with respect to v
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the gradient of the bregman loss with respect to the
#' parameter v.
bern_grad_v <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  z <- 2 * x - 1
  sum(- z * a * r(z * (a * v + s))) +  lambda * sum(- mu0 * v * r(mu0 * (a * v + s)))
}

#' @title Gradient of Bregman Loss in Bernoulli Case with respect to a
#' @param x A single row of column of x, for which we want to minimize distance
#' to the mean.
#' @param a Either a vector or scale representing the latent scores.
#' @param v Either a vector or scale representing the latent loadings.
#' @param s Either the s^th row or column of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return The value of the gradient of the bregman loss with respect to the
#' parameter a.
bern_grad_a <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  z <- 2 * x - 1
  sum(- z * v * r(z * (a * v + s))) +  lambda * sum(- mu0 * a * r(mu0 * (a * v + s)))
}

# bernoulli-pca ----------------------------------------------------------------
#' @title Bernoulli PCA
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
bern_exp_pca <- function(X, n_comp = 2, n_cycle = 30, n_iter = 30, eps = 1e-4,
                         lambda = 0.1, mu0 = NULL) {
  args <- list(obj_fun = bern_obj,
               obj_grad_a = bern_grad_a,
               obj_grad_v = bern_grad_v,
               X = X,
               n_comp = n_comp,
               n_cycle = n_cycle,
               n_iter = n_iter,
               eps = eps,
               lambda = lambda,
               mu0 = mu0)
  do.call(exp_pca, args)
}
