
################################################################################
# Experiment using exponential family PCA for binary data data
################################################################################

## ---- libraries ----
library("ggplot2")
library("reshape2")
library("plyr")
library("dplyr")

## ---- setup ----
n <- 1000
k <- 5
d <- 50
iter_max <- 20
X_data <- generate_data(k, n, d, .5, 0.1)
heatmap(X_data$X)

## ---- pca ----
pca_res <- bern_exp_pca(X_data$X, iter_max = iter_max)
a_df <- data.frame(a = pca_res$a, label = X_data$copies)

## ---- study-pca ----
ggplot(a_df) +
  geom_jitter(aes(x = a, y = 0, col = as.factor(label))) +
  ylim(c(-5, 5)) +
  ggtitle("Exponential Family PCA Scores (1-d)")

ggplot(a_df) +
  geom_jitter(aes(x = a, y = 0, col = as.factor(label))) +
  ylim(c(-5, 5)) +
  facet_wrap(~label) +
  ggtitle("Exponential Family PCA Scores (1-d)")

Z <- data.frame(X_data$X, label = X_data$copies)
Z$ix <- 1:nrow(Z)
mZ <- melt(Z, id.vars = c("ix", "label"))
mZ$ix <- as.factor(mZ$ix)

ggplot(mZ) +
  geom_tile(aes(x = variable, y = ix, fill = as.factor(value))) +
  theme(axis.text = element_blank())

ggplot(mZ) +
  geom_tile(aes(x = variable, y = ix, fill = as.factor(value))) +
  facet_wrap(~label, scales = "free_y") +
  theme(axis.text = element_blank())
