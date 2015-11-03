
################################################################################
# Simulation utiltiies for testing exponential family PCA.
################################################################################

#' @title Simulate test data for bernoulli PCA
#' @param k The number of latent binary classes.
#' @param n The number of samples to generate.
#' @param d The dimension of the samples.
#' @param p The probability of a 1 in each coordinate of each latent class.
#' @param eps The probability of flipping the sign in the observed data.
#' @export
generate_bern_data <- function(k = 4, n = 100, d = 10, p = 0.5, eps = 0.05) {
  Z <- matrix(rbinom(d * k, 1, p), k, d)
  copies <- sample(1:k, n, replace = TRUE)
  X <- (Z[copies, ] + matrix(rbinom(n * d, 1, eps), n, d)) %% 2
  return (list(X = X, copies = copies))
}

#' @title Simulate test data for usual Gaussian PCA
#' @param k The latent subspace dimension.
#' @param n The number of samples to generate.
#' @param d The dimension of the samples.
#' @param eps The variance of the noise around the latent subspace.
#' @return A list with the following elements, \cr
#'   $X: The simulated gaussian data.
#'   $V: The latent subspace.
#'   $A: The latent scores.
generate_gaussian_data <- function(k = 4, n = 100, d = 10, eps = 0.05) {
  V <- matrix(rnorm(k * d), d, k)
  V <- qr.Q(qr(V))
  A <- matrix(rnorm(n * k), n, k)

  X <- A %*% t(V) + matrix(rnorm(n * d, 0, eps), n, d)
  return (list(X = X, V = V, A = A))
}

#' @title Simulate test data for Poisson PCA
#' @param k The latent subspace dimension.
#' @param n The number of samples to generate.
#' @param d The dimension of the samples.
#' @param eps The variance of the noise around the latent subspace.
#' @return A list with the following elements, \cr
#'   $X: The simulated poisson data.
#'   $V: The latent subspace in the natural parameter space.
#'   $A: The latent scores in the natural parameter space.
generate_poisson_data <- function(k = 4, n = 100, d = 10, eps = 0.05) {
  V <- matrix(rnorm(k * d), d, k)
  V <- qr.Q(qr(V))
  A <- matrix(rnorm(n * k), n, k)

  X <- matrix(0, n, d)

  Lambda <- exp(A %*% t(V))
  for(i in 1:nrow(X)) {
    for(j in 1:ncol(X)) {
      X[i, j] <- rpois(1, Lambda[i, j])
    }
  }

  return (list(X = X, V = V, A = A))
}
