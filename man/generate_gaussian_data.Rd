% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/sim_utils.R
\name{generate_gaussian_data}
\alias{generate_gaussian_data}
\title{Simulate test data for usual Gaussian PCA}
\usage{
generate_gaussian_data(k = 4, n = 100, d = 10, eps = 0.05)
}
\arguments{
\item{k}{The latent subspace dimension.}

\item{n}{The number of samples to generate.}

\item{d}{The dimension of the samples.}

\item{eps}{The variance of the noise around the latent subspace.}
}
\value{
A list with the following elements, \cr
  $X: The simulated gaussian data.
  $V: The latent subspace.
  $A: The latent scores.
}
\description{
Simulate test data for usual Gaussian PCA
}

