
################################################################################
# Updates used in all exponential family PCA methodsb
################################################################################

# alternating-updates ----------------------------------------------------------
#' @title Update one element of scores A
#' @param obj_fun Version of loss given by bregman divergence between data and
#' link to use in current updates.
#' @param obj_grad_a Gradient of the provided loss, with respect to a.
#' @param x_i The i^th sample of X.
#' @param a_ic The initial value of a_ic to start the search from.
#' @param v_c The c^th component of V.
#' @param s_i The s^th row of the matrix AV^T - a_{c}v_{c}^{T}
#' @param lambda The regularization parameter in the optimizatyion.
#' @param mu0 The value to regularize towards.
#' @return a_ic The updated version of a_ic.
#' @export
update_a_ic <- function(obj_fun, obj_grad_a, x_i, a_ic, v_c, s_i, lambda, mu0) {
  f_min_a <- function(a) { obj_fun(x_i, a, v_c, s_i, lambda, mu0) }
  f_grad_a <- function(a) { obj_grad_a(x_i, a, v_c, s_i, lambda, mu0) }
  optim_a <- optim(a_ic, f_min_a, method = "BFGS")
  if(optim_a$convergence == 1) warning("failed to converge")
  return (optim_a$par)
}

#' @title Update one element of loadings V
#' @param obj_fun Version of loss given by bregman divergence between data and
#' link to use in current updates.
#' @param obj_grad_v Gradient of the provided loss, with respect to v.
#' @param x_j The j^th dimension of X.
#' @param v_jc The initial value of v_jc to start the search from.
#' @param a_c The c^th component of A.
#' @param s_j The s^th dimension of the matrix AV^T - a_{c}v_{c}^{T}.
#' @param lambda The regularization parameter in the optimization.
#' @param mu0 The value to regularize towards.
#' @return v_jc The updated version of v_jc.
#' @export
update_v_jc <- function(obj_fun, obj_grad_v, x_j, v_jc, a_c, s_j, lambda, mu0) {
  f_min_v <- function(v) { obj_fun(x_j, a_c, v, s_j, lambda, mu0) }
  f_grad_v <- function(v) { obj_grad_v(x_j, a_c, v, s_j, lambda, mu0) }
  optim_v <- optim(v_jc, f_min_v, method = "BFGS")
  if(optim_v$convergence == 1) warning("failed to converge")
  return (optim_v$par)
}

# alternating-minimization -----------------------------------------------------
#' @title Exponential Family PCA via Alternating Minimization
#' @description General exponential-family PCA algorithm, within which different
#' family-specific bregman losses can be input.
#' @param obj_fun The main bregman loss function to attempt to minimize. See
#' bern_breg() for an example of the required input / output format.
#' @param obj_grad_a The derivative of the bregman loss with respect to the
#' scores. See bern_brad_grad_a() for an example of the required input / output
#' format.
#' @param obj_grad_v The derivative of the bregman loss with respect to the
#' loadings. See bern_brad_grad_a() for an example of the required input /
#' output format.
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
exp_pca <- function(obj_fun, obj_grad_a, obj_grad_v, X, n_comp, n_cycle, n_iter,
                    eps, lambda, mu0) {
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
          A[i, cur_comp] <- update_a_ic(obj_fun, obj_grad_a, X[i, ],
                                        A[i, cur_comp], V[, cur_comp], S[i, ],
                                        lambda, mu0)
        }
        for(j in seq_len(ncol(X))) {
          V[j, cur_comp] <- update_v_jc(obj_fun, obj_grad_v, X[, j],
                                        V[j, cur_comp], A[, cur_comp], S[, j],
                                        lambda, mu0)
        }
      }
    }
  }
  list(A = A, V = V)
}
