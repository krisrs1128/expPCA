
################################################################################
#
################################################################################

## ---- libraries ----
library("ggplot2")
library("reshape2")
library("plyr")
library("dplyr")
library("expPCA")
theme_set(theme_bw())

## ---- setup ----
n <- 100
k <- 2
d <- 5
X_data <- generate_poisson_data(k, n, d)
heatmap(X_data$X)

## ---- pca ----
pca_res <- poisson_exp_pca(X_data$X, 2, 2, 5)

## ---- vis-result ----
graphics.off()
A_est <- data.frame(pca_res$A, ix = 1:nrow(pca_res$A))
A_true <- data.frame(X_data$A, ix = 1:nrow(X_data$A))
ggplot(A_est) +
  geom_text(aes(x = X1, y = X2, label = ix)) +
  ggtitle("Estimated Scores")
dev.new()
ggplot(A_true) +
  geom_text(aes(x = X1, y = X2, label = ix)) +
  ggtitle("True Scores")

