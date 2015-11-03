
################################################################################
# Experiment using exponential family PCA for binary data
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
k <- 5
d <- 15
X_data <- generate_bern_data(k, n, d, .5, 0.1)
heatmap(X_data$X)

## ---- pca ----
pca_res <- bern_exp_pca(X_data$X, 2, 5, 10)
a_df <- data.frame(a = pca_res$A, label = X_data$copies)
head(a_df)

## ---- study-pca ----
ggplot(a_df) +
  geom_jitter(aes(x = a.1, y = a.2, col = as.factor(label)), alpha = 0.5) +
  ylim(c(-5, 5)) +
  ggtitle("Exponential Family PCA Scores (2-d)")

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
