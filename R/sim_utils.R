
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
generate_data <- function(k = 4, n = 100, d = 10, p = 0.5, eps = 0.05) {
  Z <- matrix(rbinom(d * k, 1, p), k, d)
  copies <- sample(1:k, n, replace = TRUE)
  X <- (Z[copies, ] + matrix(rbinom(n * d, 1, eps), n, d)) %% 2
  return (list(X = X, copies = copies))
}
