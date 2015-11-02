
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

#' @title Bernoulli Bregman Divergence Objective as a Function of parameter a
#' @description To fit Exponential Family PCA, we need to alternate
#' n + d optimizations, n with respect to the scores and d with respect to the
#' loadings. This function helps return the function to optimize in each of the
#' n optimizations across scores.
#' @param x A vector representing a single sample from the matrix we are
#' decomposing with PCA.
#' @param v The current value of the loadings to keep fixed during the
#' optimization over scores.
#' @return The function to optimize over a, corresponding to the current values
#' of the loadings and the i^th sample of X.
#' @export
bern_breg_obj_a <- function(x, v, mu0 = NULL) {
  z <- 2 * x - 1
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  function(a, lambda = 0.01) {
    sum(log(1 + exp(- z * a * v))) + lambda * sum(log(1 + exp(- mu0 * a * v)))
  }
}

#' @title Gradient of Bregman Divergence with respect to parameter a
#' @param x A vector representing a single sample from the matrix we are
#' decomposing with PCA.
#' @param v The current value of the loadings to keep fixed during the
#' optimization over scores.
#' @return Gradient of the function to optimize over a, corresponding to the
#' current values of the loadings and the i^th sample of X.
#' @export
bern_breg_grad_a <- function(x, v, mu0 = NULL) {
  z <- 2 * x - 1
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  function(a, lambda = 0.1) {
    sum(- z * v * r(z * a * v)) +  lambda * sum(- mu0 * v * r(z * a * v))
  }
}

#' @title Bernoulli Bregman Divergence Objective as a Function of parameter v
#' @description To fit Exponential Family PCA, we need to alternate
#' n + d optimizations, n with respect to the scores and d with respect to the
#' loadings. This function helps return the function to optimize in each of the
#' d optimizations across loadings coordinates.
#' @param x A vector representing a single dimension from the matrix we are
#' decomposing with PCA.
#' @param a The current value of the scores to keep fixed during the
#' optimization over loadings.
#' @return The function to optimize over v, corresponding to the current values
#' of the scores and the i^th dimension of X.
bern_breg_obj_v <- function(x, a, mu0 = NULL) {
  z <- 2 * x - 1
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  function(v, lambda = 0.1) {
    sum(log(1 + exp(- z * a * v))) + lambda * sum(log(1 + exp(- mu0 * a * v)))
  }
}

#' @title Gradient of Bregman Divergence with respect to parameter v
#' @param x A vector representing a single dimension from the matrix we are
#' decomposing with PCA.
#' @param a The current value of the scores to keep fixed during the
#' optimization over loadings.
#' @return The gradient of the function to optimize over v, corresponding to the
#' current values of the scores and the i^th dimension of X.
#' @export
bern_breg_grad_v <- function(x, a, mu0 = NULL) {
  z <- 2 * x - 1
  if(is.null(mu0)) {
    mu0 <- 2 * rep(0, length(x)) - 1
  }
  function(v, lambda = 0.1) {
    sum(- z * a * r(z * a * v)) + lambda * sum(- mu0 * a * r(z * a * v))
  }
}

# alternating-minimization -----------------------------------------------------
#' @title Minimization over scores in Bernoulli PCA
#' @param X The n x p binary matrix with samples along rows which we want to
#' decompose using exponential family PCA.
#' @param a_cur The current values of the latent scores. Used to initialize the
#' search for the next iteration's scores.
#' @param v_cur The current values of the loadings.
#' @return The scores at the next iteration.
#' @export
a_step <- function(X, a_cur, v_cur, lambda = 0.1) {
  a_next <- a_cur
  for(i in seq_len(nrow(X))) {
    f_min_a <-  function(a) { bern_breg_obj_a(X[i, ], v_cur)(a, lambda) }
    f_grad_a <- function(a) { bern_breg_grad_a(X[i, ], v_cur)(a, lambda) }
    optim_a <- optim(a_cur[i], f_min_a, f_grad_a, method = "BFGS")
    if(optim_a$convergence == 1) warning("failed to converge")
    a_next[i] <- optim_a$par
  }
  return (a_next)
}

#' @title Minimization over loadings in Bernoulli PCA
#' @param X The n x p binary matrix with samples along rows which we want to
#' decompose using exponential family PCA.
#' @param a_cur The current values of the latent scores.
#' @param v_cur The current values of the loadings. Used to initialize the
#' search for the next iteration's loadings.
#' @return The loadings at the next iteration.
#' @export
v_step <- function(X, a_cur, v_cur, lambda = 0.1) {
  v_next <- v_cur
  for(j in seq_len(ncol(X))) {
    f_min_v <- function(v) { bern_breg_obj_v(X[, j], a_cur)(v, lambda) }
    f_grad_v <- function(x) { bern_breg_grad_v(X[, j], a_cur)(v, lambda) }
    optim_v <- optim(v_cur[j], f_min_v, f_grad_v, method = "BFGS")
    if(optim_v$convergence == 1) warning("failed to converge")
    v_next[j] <- optim_v$par
  }
  return (v_next)
}

#' @title Bernoulli PCA with one component
#' @param X The n x p binary matrix with samples along rows which we want to
#' decompose using exponential family PCA.
#' @param a0 Initial values for the latent scores. If not supplied, initialized
#' with standard normals.
#' @param v0 Initial values for loading vector. If not supplied, initialized
#' with standard normals.
#' @param iter_max The maximum number of iterations to run the PCA.
#' @param eps The convergence criterion. If the mean change in the scores is
#' less than eps, we return.
#' @return A list containing the converged values of a and v, along with traces
#' of their values across the optimization [A and V, the columns are values
#' across iterations].
#' @references Collins, Michael, Sanjoy Dasgupta, and Robert E. Schapire. "A generalization of principal components analysis to the exponential family." Advances in neural information processing systems. 2001.
bern_exp_pca <- function(X, a0 = NULL, v0 = NULL, iter_max = 100, eps = 1e-4,
                         lambda = 0.1) {
  # initialize score and loadings, if not supplied
  if(is.null(a0)) {
    a0 <- rnorm(nrow(X))
  }
  if(is.null(v0)) {
    v0 <- rnorm(ncol(X))
  }

  # initialize trace of scores / loadings
  A <- cbind(a0, matrix(0, length(a0), iter_max))
  V <- cbind(v0, matrix(0, length(v0), iter_max))

  # alternating minimization
  for(iter in seq_len(iter_max)) {
    cat(sprintf("iteration %g \n", iter))
    A[, iter + 1] <- a_step(X, A[, iter], V[, iter], lambda)
    V[, iter + 1] <- v_step(X, A[, iter + 1], V[, iter], lambda)

    tol <- mean(abs(A[, iter + 1] - A[, iter]))
    if(tol < eps) break
  }

  list(a = A[, iter + 1], v = V[, iter + 1],
       A = A[, seq_len(iter + 1)], V = V[, seq_len(iter + 1)])
}
