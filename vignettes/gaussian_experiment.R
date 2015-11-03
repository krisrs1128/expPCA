
################################################################################
# Experiment using exponential family PCA for binary data data
################################################################################

## ---- libraries ----
library("ggplot2")
library("reshape2")
library("plyr")
library("dplyr")
library("expPCA")
theme_set(theme_bw())

## ---- setup ----
n <- 1000
k <- 2
d <- 15
X_data <- generate_gaussian_data(k, n, d, .5, 0.1)
heatmap(X_data$X)

## ---- pca ----
pca_res <- gaussian_exp_pca(X_data$X, 2, 5, 10)

## ---- vis-result ----
A_est <- data.frame(pca_res$A, ix = 1:nrow(pca_res$A))
A_true <- data.frame(X_data$A, ix = 1:nrow(X_data$A))
ggplot(A_est) +
  geom_text(aes(x = X1, y = X2, label = ix)) +
  ggtitle("Estimated Scores")
ggplot(A_true) +
  geom_text(aes(x = X1, y = X2, label = ix)) +
  ggtitle("True Scores")
