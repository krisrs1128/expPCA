
################################################################################
# Perform exponential family PCA using a bernoulli model
################################################################################

# optim-funs -------------------------------------------------------------------
#' @title Logistic Function
#' @param x Input vector or scale for logistic fun.
#' @return v Value of logistic function at input.
r <- function(x) {
  exp(-x) / (1 + exp(-x))
}

#' @title Bregman Loss in Bernoulli case
#' @export
bern_obj <- function(x, a, v, s, lambda, mu0) {
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  z <- 2 * x - 1
  sum(log(1 + exp(- z * (a * v + s)))) + lambda * sum(log(1 + exp(- mu0 * (a * v + s))))
}

#' @title Update one element of scores A
#' @param obj_fun Version of loss given by bregman divergence between data and
#' link to use in current updates.
#' @param x_i The i^th sample of X.
#' @param a_ic The initial value of a_ic to start the search from.
#' @param v_c The c^th component of V.
#' @param s_i The s^th row of of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimizatyion.
#' @param mu0 The value to regularize towards.
#' @return a_ic The updated version of a_ic.
#' @export
update_a_ic <- function(obj_fun, x_i, a_ic, v_c, s_i, lambda, mu0) {
  f_min_a <- function(a) { obj_fun(x_i, a, v_c, s_i, lambda, mu0) }
  optim_a <- optim(a_ic, f_min_a, method = "BFGS")
  if(optim_a$convergence == 1) warning("failed to converge")
  return (optim_a$par)
}

#' @title Update one element of loadings V
#' @param obj_fun Version of loss given by bregman divergence between data and
#' link to use in current updates.
#' @param x_j The j^th dimension of X.
#' @param v_jc The initial value of v_jc to start the search from.
#' @param a_c The c^th component of A.
#' @param s_j The s^th dimension of the matrix AV^T - a_{c}v_{c}^{T}.
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return v_jc The updated version of v_jc.
#' @export
update_v_jc <- function(obj_fun, x_j, v_jc, a_c, s_j, lambda, mu0) {
  f_min_v <- function(v) { obj_fun(x_j, a_c, v, s_j, lambda, mu0) }
  optim_v <- optim(v_jc, f_min_v, method = "BFGS")
  if(optim_v$convergence == 1) warning("failed to converge")
  return (optim_v$par)
}

#' @title Bernoulli PCA with one component
#' @param X The n x p binary matrix with samples along rows which we want to
#' decompose using exponential family PCA.
#' @param iter_max The maximum number of iterations to run the PCA.
#' @param eps The convergence criterion. If the mean change in the scores is
#' less than eps, we return.
#' @return A list containing the converged values of a and v, along with traces
#' of their values across the optimization [A and V, the columns are values
#' across iterations].
#' @references Collins, Michael, Sanjoy Dasgupta, and Robert E. Schapire. "A generalization of principal components analysis to the exponential family." Advances in neural information processing systems. 2001.
#' @export
bern_exp_pca <- function(X, n_comp = 2, n_cycle = 30, n_iter = 30, eps = 1e-4,
                         lambda = 0.1, mu0 = NULL) {
  # initialize score and loadings to 0
  A <- matrix(0, nrow(X), n_comp)
  V <- matrix(0, ncol(X), n_comp)

  # cycle through the n_comp components n_cycle times
  for(cur_cycle in seq_len(n_cycle)) {
    for(cur_comp in seq_len(n_comp)) {

      # Prepare optimization over c^th component, see section 4 in reference
      V[, cur_comp] <- rnorm(nrow(V))
      S <- A %*% t(V) - A[, cur_comp] %*% t(V[, cur_comp])

      # alternating minimization
      for(cur_iter in seq_len(n_iter)) {
        cat(sprintf("cycle %g \t component %g \t iteration %g \n", cur_cycle, cur_comp, cur_iter))
        for(i in seq_len(nrow(X))) {
          A[i, cur_comp] <- update_a_ic(bern_obj, X[i, ], A[i, cur_comp],
                                        V[, cur_comp], S[i, ], lambda, mu0)
        }
        for(j in seq_len(ncol(X))) {
          V[j, cur_comp] <- update_v_jc(bern_obj, X[, j], V[j, cur_comp],
                                        A[, cur_comp], S[, j], lambda, mu0)
        }
      }
    }
  }
  list(A = A, V = V)
}
